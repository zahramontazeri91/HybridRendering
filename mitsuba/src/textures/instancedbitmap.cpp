/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/mipmap.h>
#include <mitsuba/hw/renderer.h>
#include <mitsuba/hw/gputexture.h>
#include <mitsuba/hw/gpuprogram.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

class InstanceBitmapTexture : public Texture2D {
public:
	/* Store texture data using half precision, but perform computations in
	   single/double precision based on compilation flags. The following
	   generates efficient implementations for both luminance and RGB data */
	typedef TSpectrum<Float, 1> Color1;
	typedef TSpectrum<Float, 3> Color3;
	typedef TSpectrum<half, 1>  Color1h;
	typedef TSpectrum<half, 3>  Color3h;
	typedef TMIPMap<Color1, Color1h> MIPMap1;
	typedef TMIPMap<Color3, Color3h> MIPMap3;

	InstanceBitmapTexture(const Properties &props) : Texture2D(props) {
		uint64_t timestamp = 0;
		bool tryReuseCache = false;
		fs::path cacheFile;
		ref<Bitmap> bitmap;

		m_channel = boost::to_lower_copy(props.getString("channel", ""));

		if (props.hasProperty("bitmap")) {
			/* Support initialization via raw data passed from another plugin */
			bitmap = reinterpret_cast<Bitmap *>(props.getData("bitmap").ptr);
		} else {
			m_filename = Thread::getThread()->getFileResolver()->resolve(
				props.getString("filename"));

			Log(EInfo, "Loading texture \"%s\"", m_filename.filename().string().c_str());
			if (!fs::exists(m_filename))
				Log(EError, "Texture file \"%s\" could not be found!", m_filename.string().c_str());

			boost::system::error_code ec;
			timestamp = (uint64_t) fs::last_write_time(m_filename, ec);
			if (ec.value())
				Log(EError, "Could not determine modification time of \"%s\"!", m_filename.string().c_str());

			cacheFile = m_filename;

			if (m_channel.empty())
				cacheFile.replace_extension(".mip");
			else
				cacheFile.replace_extension(formatString(".%s.mip", m_channel.c_str()));

			tryReuseCache = fs::exists(cacheFile) && props.getBoolean("cache", true);
		}

		std::string filterType = boost::to_lower_copy(props.getString("filterType", "ewa"));
		std::string wrapMode = props.getString("wrapMode", "repeat");
		m_wrapModeU = parseWrapMode(props.getString("wrapModeU", wrapMode));
		m_wrapModeV = parseWrapMode(props.getString("wrapModeV", wrapMode));

		m_gamma = props.getFloat("gamma", 0);

		if (filterType == "ewa")
			m_filterType = EEWA;
		else if (filterType == "bilinear")
			m_filterType = EBilinear;
		else if (filterType == "trilinear")
			m_filterType = ETrilinear;
		else if (filterType == "nearest")
			m_filterType = ENearest;
		else
			Log(EError, "Unknown filter type '%s' -- must be "
				"'ewa', 'trilinear', or 'nearest'!", filterType.c_str());

		m_maxAnisotropy = props.getFloat("maxAnisotropy", 20);

		if (m_filterType != EEWA)
			m_maxAnisotropy = 1.0f;

		if (tryReuseCache && MIPMap3::validateCacheFile(cacheFile, timestamp,
				Bitmap::ERGB, m_wrapModeU, m_wrapModeV, m_filterType, m_gamma)) {
			/* Reuse an existing MIP map cache file */
			m_mipmap3 = new MIPMap3(cacheFile, m_maxAnisotropy);
		} else if (tryReuseCache && MIPMap1::validateCacheFile(cacheFile, timestamp,
				Bitmap::ELuminance, m_wrapModeU, m_wrapModeV, m_filterType, m_gamma)) {
			/* Reuse an existing MIP map cache file */
			m_mipmap1 = new MIPMap1(cacheFile, m_maxAnisotropy);
		} else {
			if (bitmap == NULL) {
				/* Load the input image if necessary */
				ref<Timer> timer = new Timer();
				ref<FileStream> fs = new FileStream(m_filename, FileStream::EReadOnly);
				bitmap = new Bitmap(Bitmap::EAuto, fs);
				if (m_gamma != 0)
					bitmap->setGamma(m_gamma);
				Log(EDebug, "Loaded \"%s\" in %i ms", m_filename.filename().string().c_str(),
					timer->getMilliseconds());
			}

			Bitmap::EPixelFormat pixelFormat;
			if (!m_channel.empty()) {
				/* Create a texture from a certain channel of an image */
				pixelFormat = Bitmap::ELuminance;
				bitmap = bitmap->extractChannel(findChannel(bitmap, m_channel));
				if (m_channel == "a")
					bitmap->setGamma(1.0f);
			} else {
				switch (bitmap->getPixelFormat()) {
					case Bitmap::ELuminance:
					case Bitmap::ELuminanceAlpha:
						pixelFormat = Bitmap::ELuminance;
						break;
					case Bitmap::ERGB:
					case Bitmap::ERGBA:
						pixelFormat = Bitmap::ERGB;
						break;
					default:
						Log(EError, "The input image has an unsupported pixel format!");
						return;
				}
			}

			/* (Re)generate the MIP map hierarchy; downsample using a
			    2-lobed Lanczos reconstruction filter */
			Properties rfilterProps("lanczos");
			rfilterProps.setInteger("lobes", 2);
			ref<ReconstructionFilter> rfilter = static_cast<ReconstructionFilter *> (
				PluginManager::getInstance()->createObject(
				MTS_CLASS(ReconstructionFilter), rfilterProps));
			rfilter->configure();

			/* Potentially create a new MIP map cache file */
			bool createCache = !cacheFile.empty() && props.getBoolean("cache",
				bitmap->getSize().x * bitmap->getSize().y > 1024*1024);
			
			if (pixelFormat == Bitmap::ELuminance)
				m_mipmap1 = new MIPMap1(bitmap, pixelFormat, Bitmap::EFloat,
					rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
					createCache ? cacheFile : fs::path(), timestamp);
			else 
				m_mipmap3 = new MIPMap3(bitmap, pixelFormat, Bitmap::EFloat,
					rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
					createCache ? cacheFile : fs::path(), timestamp);
		}
		
		///tiling purpose
		m_divideReso.x = props.getInteger("divideX");
        m_divideReso.y = props.getInteger("divideY");
        m_reso.x = props.getInteger("tileX");
        m_reso.y = props.getInteger("tileY"); 
		m_blockFile = props.getString("blockFile"); 
		ref<FileStream> fs = new FileStream(m_blockFile, FileStream::EReadWrite);
		bool constructBlockID = 1;
		
#if 0
		std::cout<<"@@@@@@@ constructing the dat file ... @@@@@@@ " <<std::endl;		
		fs->flush();
  		fs->writeInt(m_reso.x*m_reso.y);
  		for (int j=0; j<m_reso.y; j++)
			for (int i=0; i<m_reso.x; i++) {
				int temp = ( (j%m_divideReso.y)*m_divideReso.x + i%m_divideReso.x );
				fs->writeInt(temp);
			}  
#else		
 		m_sz = fs->readInt();
		m_blockID.assign(m_reso.x*m_reso.y,0);
		for (int i=0; i<m_reso.y; i++) {
			for (int j=0; j<m_reso.x; j++) {
				m_blockID.at(i*m_reso.x + j) = fs->readInt();
			}
			if (i != m_reso.y-1)
				for (int k=m_reso.x; k<400; k++) 
					fs->readInt();
		} 
#endif		
	}

	static int findChannel(const Bitmap *bitmap, const std::string channel) {
		int found = -1;
		std::string channelNames;
		for (int i=0; i<bitmap->getChannelCount(); ++i) {
			std::string name = boost::to_lower_copy(bitmap->getChannelName(i));
			if (name == channel)
				found = i;
			channelNames += name;
			if (i + 1 < bitmap->getChannelCount())
				channelNames += std::string(", ");
		}

		if (found == -1) {
			Log(EError, "Channel \"%s\" not found! Must be one of: [%s]",
				channel.c_str(), channelNames.c_str());
		}

		return found;
	}

	inline ReconstructionFilter::EBoundaryCondition parseWrapMode(const std::string &wrapMode) {
		if (wrapMode == "repeat")
			return ReconstructionFilter::ERepeat;
		else if (wrapMode == "clamp")
			return ReconstructionFilter::EClamp;
		else if (wrapMode == "mirror")
			return ReconstructionFilter::EMirror;
		else if (wrapMode == "zero" || wrapMode == "black")
			return ReconstructionFilter::EZero;
		else if (wrapMode == "one" || wrapMode == "white")
			return ReconstructionFilter::EOne;
		else
			Log(EError, "Unknown wrap mode '%s' -- must be "
				"'repeat', 'clamp', 'black', or 'white'!", wrapMode.c_str());
		return ReconstructionFilter::EZero; // make gcc happy
	}

	InstanceBitmapTexture(Stream *stream, InstanceManager *manager)
	 : Texture2D(stream, manager) {
		m_filename = stream->readString();
		Log(EDebug, "Unserializing texture \"%s\"", m_filename.filename().string().c_str());
		m_filterType = (EMIPFilterType) stream->readUInt();
		m_wrapModeU = (ReconstructionFilter::EBoundaryCondition) stream->readUInt();
		m_wrapModeV = (ReconstructionFilter::EBoundaryCondition) stream->readUInt();
		m_gamma = stream->readFloat();
		m_maxAnisotropy = stream->readFloat();
		m_channel = stream->readString();

		size_t size = stream->readSize();
		ref<MemoryStream> mStream = new MemoryStream(size);
		stream->copyTo(mStream, size);
		mStream->seek(0);
		ref<Bitmap> bitmap = new Bitmap(Bitmap::EAuto, mStream);
		if (m_gamma != 0)
			bitmap->setGamma(m_gamma);

		/* Downsample using a 2-lobed Lanczos reconstruction filter */
		Properties rfilterProps("lanczos");
		rfilterProps.setInteger("lobes", 2);
		ref<ReconstructionFilter> rfilter = static_cast<ReconstructionFilter *> (
			PluginManager::getInstance()->createObject(
			MTS_CLASS(ReconstructionFilter), rfilterProps));
		rfilter->configure();

		Bitmap::EPixelFormat pixelFormat;
		if (!m_channel.empty()) {
			/* Create a texture from a certain channel of an image */
			pixelFormat = Bitmap::ELuminance;
			bitmap = bitmap->extractChannel(findChannel(bitmap, m_channel));
			if (m_channel == "a")
				bitmap->setGamma(1.0f);
		} else {
			switch (bitmap->getPixelFormat()) {
				case Bitmap::ELuminance:
				case Bitmap::ELuminanceAlpha:
					pixelFormat = Bitmap::ELuminance;
					break;
				case Bitmap::ERGB:
				case Bitmap::ERGBA:
					pixelFormat = Bitmap::ERGB;
					break;
				default:
					Log(EError, "The input image has an unsupported pixel format!");
					return;
			}
		}

		if (pixelFormat == Bitmap::ELuminance)
			m_mipmap1 = new MIPMap1(bitmap, pixelFormat, Bitmap::EFloat,
				rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
				fs::path(), 0);
		else
			m_mipmap3 = new MIPMap3(bitmap, pixelFormat, Bitmap::EFloat,
				rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
				fs::path(), 0);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture2D::serialize(stream, manager);
		stream->writeString(m_filename.string());
		stream->writeUInt(m_filterType);
		stream->writeUInt(m_wrapModeU);
		stream->writeUInt(m_wrapModeV);
		stream->writeFloat(m_gamma);
		stream->writeFloat(m_maxAnisotropy);

		if (!m_filename.empty() && fs::exists(m_filename)) {
			/* We still have access to the original image -- use that, since
			   it is probably much smaller than the in-memory representation */
			ref<Stream> is = new FileStream(m_filename, FileStream::EReadOnly);
			stream->writeString(m_channel);
			stream->writeSize(is->getSize());
			is->copyTo(stream);
		} else {
			/* No access to the original image anymore. Create an EXR image
			   from the top MIP map level and serialize that */
			ref<MemoryStream> mStream = new MemoryStream();
			ref<Bitmap> bitmap = m_mipmap1.get() ?
				m_mipmap1->toBitmap() : m_mipmap3->toBitmap();
			bitmap->write(Bitmap::EOpenEXR, mStream);

			stream->writeString("");
			stream->writeSize(mStream->getSize());
			stream->write(mStream->getData(), mStream->getSize());
		}
	}

	Spectrum eval(const Point2 &uv) const {
		/* There are no ray differentials to do any kind of
		   prefiltering. Evaluate the full-resolution texture */

				///tiling purpose		
		bool read_from_file = false;

#if 0
		std::cout<< "@@@ not reading from input@@@" <<std::endl;
		int tileX = 4;
		int tileY = 6;
		int divideX = 6;
		int divideY = 6;

		Point2 uv_tiled(uv);
		std::vector<int> m_blockID(tileX*tileY);
		m_blockID.assign(tileX*tileY,0);
  		for (int j=0; j<tileY; j++)
			for (int i=0; i<tileX; i++) 
				m_blockID.at(j*tileX + i) = ( (j%divideY)*divideX + i%divideX );
	
		
		float ix = std::floor(uv.x*tileX);
		float iy = std::floor(uv.y*tileY);
		int idx = m_blockID[int(iy*tileX + ix)]; 

		int bx = idx % divideX;
		int by = std::floor(idx/divideX);
		
		uv_tiled.x = float(bx)/float(divideX) + (float(uv.x*tileX)-ix)/float(divideX);
		uv_tiled.y = float(by)/float(divideY) + (float(uv.y*tileY)-iy)/float(divideY); 		
 			
#else		
		float ix = std::floor(uv.x*m_reso.x);
		float iy = std::floor(uv.y*m_reso.y);
		int idx = m_blockID[int(iy*m_reso.x + ix)]; 

		int bx = idx % m_divideReso.x;
		int by = std::floor(idx/m_divideReso.x);
		
		Point2 uv_tiled(uv);
		uv_tiled.x = float(bx)/float(m_divideReso.x) + (float(uv.x*m_reso.x)-ix)/float(m_divideReso.x);
		uv_tiled.y = float(by)/float(m_divideReso.y) + (float(uv.y*m_reso.y)-iy)/float(m_divideReso.y); 
#endif	

		Spectrum result;
		if (m_mipmap3.get()) {
			Color3 value;
			if (m_mipmap3->getFilterType() != ENearest)
				value = m_mipmap3->evalBilinear(0, uv_tiled);
			else
				value = m_mipmap3->evalBox(0, uv_tiled);
			result.fromLinearRGB(value[0], value[1], value[2]);
		} else {
			Color1 value;
			if (m_mipmap1->getFilterType() != ENearest)
				value = m_mipmap1->evalBilinear(0, uv_tiled);
			else
				value = m_mipmap1->evalBox(0, uv_tiled);
			result = Spectrum(value[0]);
		}
		stats::filteredLookups.incrementBase();

		return result;
	}

	void evalGradient(const Point2 &uv, Spectrum *gradient) const {
		/* There are no ray differentials to do any kind of
		   prefiltering. Evaluate the full-resolution texture */
		if (m_mipmap3.get()) {
			Color3 result[2];
			if (m_mipmap3->getFilterType() != ENearest) {
				m_mipmap3->evalGradientBilinear(0, uv, result);
				gradient[0].fromLinearRGB(result[0][0], result[0][1], result[0][2]);
				gradient[1].fromLinearRGB(result[1][0], result[1][1], result[1][2]);
			} else {
				gradient[0] = gradient[1] = Spectrum(0.0f);
			}
		} else {
			Color1 result[2];
			if (m_mipmap1->getFilterType() != ENearest) {
				m_mipmap1->evalGradientBilinear(0, uv, result);
				gradient[0] = Spectrum(result[0][0]);
				gradient[1] = Spectrum(result[1][0]);
			} else {
				gradient[0] = gradient[1] = Spectrum(0.0f);
			}
		}
		stats::filteredLookups.incrementBase();
	}

	ref<Bitmap> getBitmap(const Vector2i &/* unused */) const {
	
		ref<Bitmap> untiled_bitmap =  m_mipmap1.get() ? m_mipmap1->toBitmap() : m_mipmap3->toBitmap();
		int width_untiled = untiled_bitmap->getWidth();
		int height_untiled = untiled_bitmap->getHeight();
	int tiled_reso_y = (m_reso.y  * width_untiled ) / m_divideReso.y ;
		int tiled_reso_x = (m_reso.x  * height_untiled ) / m_divideReso.x;
		ref<Bitmap> tiled_bitmap = new Bitmap(Bitmap::ELuminance, Bitmap::EUInt16, Vector2i(tiled_reso_x,tiled_reso_y), 1, 0); 
		//ref<Bitmap> tiled_bitmap = new Bitmap(Bitmap::ELuminance, Bitmap::EUInt16, Vector2i((457*4)/6,(671*6)/6), 1, 0);
		
		int width = tiled_bitmap->getWidth();
		int height = tiled_bitmap->getHeight();
		Point2 uv;
		for (int i=0; i< height-1; i++)
			for (int j=0; j < width-1; j++) {
				uv.x = (j+0.5f)/width;
				uv.y = (i+0.5f)/height;
				tiled_bitmap->setPixel(Point2i(j,i), eval(uv, Vector2(0,0), Vector2(0,0) ) );
			}
			
		return tiled_bitmap;		
	}

	Spectrum eval(const Point2 &uv, const Vector2 &d0, const Vector2 &d1) const {
		stats::filteredLookups.incrementBase();
		++stats::filteredLookups;

		///tiling purpose		
		bool read_from_file = false;

#if 0
		std::cout<< "@@@ not reading from input@@@" <<std::endl;
		int tileX = 4;
		int tileY = 6;
		int divideX = 6;
		int divideY = 6;

		Point2 uv_tiled(uv);
		std::vector<int> m_blockID(tileX*tileY);
		m_blockID.assign(tileX*tileY,0);
  		for (int j=0; j<tileY; j++)
			for (int i=0; i<tileX; i++) 
				m_blockID.at(j*tileX + i) = ( (j%divideY)*divideX + i%divideX );
	
		
		float ix = std::floor(uv.x*tileX);
		float iy = std::floor(uv.y*tileY);
		int idx = m_blockID[int(iy*tileX + ix)]; 

		int bx = idx % divideX;
		int by = std::floor(idx/divideX);
		
		uv_tiled.x = float(bx)/float(divideX) + (float(uv.x*tileX)-ix)/float(divideX);
		uv_tiled.y = float(by)/float(divideY) + (float(uv.y*tileY)-iy)/float(divideY); 		
 			
#else		
		float ix = std::floor(uv.x*m_reso.x);
		float iy = std::floor(uv.y*m_reso.y);
		int idx = m_blockID[int(iy*m_reso.x + ix)]; 

		int bx = idx % m_divideReso.x;
		int by = std::floor(idx/m_divideReso.x);
		
		Point2 uv_tiled(uv);
		uv_tiled.x = float(bx)/float(m_divideReso.x) + (float(uv.x*m_reso.x)-ix)/float(m_divideReso.x);
		uv_tiled.y = float(by)/float(m_divideReso.y) + (float(uv.y*m_reso.y)-iy)/float(m_divideReso.y); 
#endif		

		Spectrum result;
		if (m_mipmap3.get()) {
			Color3 value = m_mipmap3->eval(uv_tiled, d0, d1);
			result.fromLinearRGB(value[0], value[1], value[2]);
		} else {
			Color1 value = m_mipmap1->eval(uv_tiled, d0, d1);
			result = Spectrum(value[0]);
		}
		return result;
	}

	Spectrum getAverage() const {
		Spectrum result;
		if (m_mipmap3.get()) {
			Color3 value = m_mipmap3->getAverage();
			result.fromLinearRGB(value[0], value[1], value[2]);
		} else {
			Color1 value = m_mipmap1->getAverage();
			result = Spectrum(value[0]);
		}
		return result;
	}

	Spectrum getMaximum() const {
		Spectrum result;
		if (m_mipmap3.get()) {
			Color3 value = m_mipmap3->getMaximum();
			result.fromLinearRGB(value[0], value[1], value[2]);
		} else {
			Color1 value = m_mipmap1->getMaximum();
			result = Spectrum(value[0]);
		}
		return result;
	}

	Spectrum getMinimum() const {
		Spectrum result;
		if (m_mipmap3.get()) {
			Color3 value = m_mipmap3->getMinimum();
			result.fromLinearRGB(value[0], value[1], value[2]);
		} else {
			Color1 value = m_mipmap1->getMinimum();
			result = Spectrum(value[0]);
		}
		return result;
	}

	bool isConstant() const {
		return false;
	}

	bool usesRayDifferentials() const {
		return true;
	}

	bool isMonochromatic() const {
		return m_mipmap1.get() != NULL;
	}

	Vector3i getResolution() const {
		if (m_mipmap3.get()) {
			return Vector3i(
				m_mipmap3->getWidth(),
				m_mipmap3->getHeight(),
				1
			);
		} else {
			return Vector3i(
				m_mipmap1->getWidth(),
				m_mipmap1->getHeight(),
				1
			);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "InstanceBitmapTexture[" << endl
			<< "  filename = \"" << m_filename.string() << "\"," << endl;

		if (m_mipmap3.get())
			oss << "  mipmap = " << indent(m_mipmap3.toString()) << endl;
		else
			oss << "  mipmap = " << indent(m_mipmap1.toString()) << endl;

		oss << "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	ref<MIPMap1> m_mipmap1;
	ref<MIPMap3> m_mipmap3;
	EMIPFilterType m_filterType;
	ReconstructionFilter::EBoundaryCondition m_wrapModeU;
	ReconstructionFilter::EBoundaryCondition m_wrapModeV;
	Float m_gamma, m_maxAnisotropy;
	std::string m_channel;
	fs::path m_filename;
	
	///tiling purpose
    Vector2i m_divideReso;
    Vector2i m_reso;
	std::vector<int> m_blockID;
	std::string m_blockFile;	
	int m_sz;
};

// ================ Hardware shader implementation ================
class InstanceBitmapTextureShader : public Shader {
public:
	InstanceBitmapTextureShader(Renderer *renderer, const std::string &filename,
			const InstanceBitmapTexture::MIPMap1* mipmap1,
			const InstanceBitmapTexture::MIPMap3* mipmap3,
			const Point2 &uvOffset, const Vector2 &uvScale,
			ReconstructionFilter::EBoundaryCondition wrapModeU,
			ReconstructionFilter::EBoundaryCondition wrapModeV,
			Float maxAnisotropy)
		: Shader(renderer, ETextureShader), m_uvOffset(uvOffset), m_uvScale(uvScale) {

		ref<Bitmap> bitmap = mipmap1 ? mipmap1->toBitmap() : mipmap3->toBitmap();
		m_gpuTexture = renderer->createGPUTexture(filename, bitmap);

		switch (wrapModeU) {
			case ReconstructionFilter::EClamp:
				m_gpuTexture->setWrapType(GPUTexture::EClampToEdge);
				break;
			case ReconstructionFilter::EMirror:
				m_gpuTexture->setWrapType(GPUTexture::EMirror);
				break;
			case ReconstructionFilter::ERepeat:
				m_gpuTexture->setWrapType(GPUTexture::ERepeat);
				break;
			case ReconstructionFilter::EZero:
				m_gpuTexture->setWrapType(GPUTexture::EClampToBorder);
				m_gpuTexture->setBorderColor(Color3(0.0f));
				break;
			case ReconstructionFilter::EOne:
				m_gpuTexture->setWrapType(GPUTexture::EClampToBorder);
				m_gpuTexture->setBorderColor(Color3(1.0f));
				break;
			default:
				Log(EError, "Unknown wrap mode!");
		}

		switch (mipmap1 ? mipmap1->getFilterType() : mipmap3->getFilterType()) {
			case ENearest:
				m_gpuTexture->setFilterType(GPUTexture::ENearest);
				break;
			case EBilinear:
				m_gpuTexture->setFilterType(GPUTexture::ELinear);
				m_gpuTexture->setMipMapped(false);
				break;
			default:
				m_gpuTexture->setFilterType(GPUTexture::EMipMapLinear);
				break;
		}

		m_gpuTexture->setMaxAnisotropy(maxAnisotropy);
		m_gpuTexture->setMaxAnisotropy(maxAnisotropy);
		m_gpuTexture->initAndRelease();
	}

	void cleanup(Renderer *renderer) {
		m_gpuTexture->cleanup();
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform sampler2D " << evalName << "_texture;" << endl
			<< "uniform vec2 " << evalName << "_uvOffset;" << endl
			<< "uniform vec2 " << evalName << "_uvScale;" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv) {" << endl
			<< "    return texture2D(" << evalName << "_texture, vec2(" << endl
			<< "          uv.x * " << evalName << "_uvScale.x + " << evalName << "_uvOffset.x," << endl
			<< "          uv.y * " << evalName << "_uvScale.y + " << evalName << "_uvOffset.y)).rgb;" << endl
			<< "}" << endl;
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_texture", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_uvOffset", false));
		parameterIDs.push_back(program->getParameterID(evalName + "_uvScale", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs,
		int &textureUnitOffset) const {
		m_gpuTexture->bind(textureUnitOffset++);
		program->setParameter(parameterIDs[0], m_gpuTexture.get());
		program->setParameter(parameterIDs[1], m_uvOffset);
		program->setParameter(parameterIDs[2], m_uvScale);
	}

	void unbind() const {
		m_gpuTexture->unbind();
	}

	MTS_DECLARE_CLASS()
private:
	ref<GPUTexture> m_gpuTexture;
	Point2 m_uvOffset;
	Vector2 m_uvScale;
	
};

Shader *InstanceBitmapTexture::createShader(Renderer *renderer) const {
	return new InstanceBitmapTextureShader(renderer, m_filename.filename().string(),
			m_mipmap1.get(), m_mipmap3.get(), m_uvOffset, m_uvScale,
			m_wrapModeU, m_wrapModeV, m_maxAnisotropy);
}

MTS_IMPLEMENT_CLASS_S(InstanceBitmapTexture, false, Texture2D)
MTS_IMPLEMENT_CLASS(InstanceBitmapTextureShader, false, Shader)
MTS_EXPORT_PLUGIN(InstanceBitmapTexture, "Bitmap texture");
MTS_NAMESPACE_END
