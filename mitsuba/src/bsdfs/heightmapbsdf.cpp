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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{diffuse}{Smooth diffuse material}
 * \order{1}
 * \icon{bsdf_diffuse}
 * \parameters{
 *     \parameter{reflectance}{\Spectrum\Or\Texture}{
 *       Specifies the diffuse albedo of the
 *       material \default{0.5}
 *     }
 * }
 *
 * \renderings{
 *     \rendering{Homogeneous reflectance, see \lstref{diffuse-uniform}}
 *         {bsdf_diffuse_plain}
 *     \rendering{Textured reflectance, see \lstref{diffuse-textured}}
 *         {bsdf_diffuse_textured}
 * }
 *
 * The smooth diffuse material (also referred to as ``Lambertian'')
 * represents an ideally diffuse material with a user-specified amount of
 * reflectance. Any received illumination is scattered so that the surface
 * looks the same independently of the direction of observation.
 *
 * Apart from a  homogeneous reflectance value, the plugin can also accept
 * a nested or referenced texture map to be used as the source of reflectance
 * information, which is then mapped onto the shape based on its UV
 * parameterization. When no parameters are specified, the model uses the default
 * of 50% reflectance.
 *
 * Note that this material is one-sided---that is, observed from the
 * back side, it will be completely black. If this is undesirable,
 * consider using the \pluginref{twosided} BRDF adapter plugin.
 * \vspace{4mm}
 *
 * \begin{xml}[caption={A diffuse material, whose reflectance is specified
 *     as an sRGB color}, label=lst:diffuse-uniform]
 * <bsdf type="diffuse">
 *     <srgb name="reflectance" value="#6d7185"/>
 * </bsdf>
 * \end{xml}
 *
 * \begin{xml}[caption=A diffuse material with a texture map,
 *     label=lst:diffuse-textured]
 * <bsdf type="diffuse">
 *     <texture type="bitmap" name="reflectance">
 *         <string name="filename" value="wood.jpg"/>
 *     </texture>
 * </bsdf>
 * \end{xml}
 */
class HeightMapBSDF : public BSDF {
public:
	HeightMapBSDF(const Properties &props)
		: BSDF(props) {
		/* For better compatibility with other models, support both
		   'reflectance' and 'diffuseReflectance' as parameter names */
		m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
			props.hasProperty("reflectance") ? "reflectance"
				: "diffuseReflectance", Spectrum(.5f)));
		m_tangentMap = new ConstantSpectrumTexture(
			props.getSpectrum("tangentMap", Spectrum(1.0f)));
		m_Kd = new ConstantSpectrumTexture(
			props.getSpectrum("Kd", Spectrum(0.5f)));
		m_Ks = new ConstantSpectrumTexture(
			props.getSpectrum("Ks", Spectrum(0.2f)));
		m_Beta = new ConstantSpectrumTexture(
			props.getSpectrum("Beta", Spectrum(0.2f)));
	}

	HeightMapBSDF(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_reflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_Kd = static_cast<Texture *>(manager->getInstance(stream));
		m_Ks = static_cast<Texture *>(manager->getInstance(stream));
		m_Beta = static_cast<Texture *>(manager->getInstance(stream));

		configure();
	}

	void configure() {
		/* Verify the input parameter and fix them if necessary */
		m_reflectance = ensureEnergyConservation(m_reflectance, "reflectance", 1.0f);

		m_components.clear();
		if (m_reflectance->getMaximum().max() > 0)
			m_components.push_back(EDiffuseReflection | EFrontSide
				| (m_reflectance->isConstant() ? 0 : ESpatiallyVarying));
			m_usesRayDifferentials = m_reflectance->usesRayDifferentials();

		BSDF::configure();
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		return m_reflectance->eval(its);
	}

	double CalTheta(Vector V, Vector T) const{
		V = normalize(V);
		T = normalize(T);
		double VdotT = dot(V, T);
		double phi = acos(VdotT);
		return phi - 0.5*M_PI;
	}

	double CalPhi(Vector V, Vector T, Vector N) const{
		V = normalize(V);
		T = normalize(T);
		double theta = CalTheta(V, T);
		Vector X = normalize(V - T*sin(theta));
		double XdotN = dot(X, N);
		return acos(XdotN);
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		Spectrum Kd = m_Kd->eval(bRec.its);
		Spectrum Ks = m_Ks->eval(bRec.its);
		Spectrum Beta = m_Beta->eval(bRec.its);

		Spectrum F(0.0f);
		Spectrum tangent = m_tangentMap->eval(bRec.its);

		//Vector t = normalize(Vector(tangent[0], tangent[1], tangent[2]));
		Vector t = normalize(Vector(tangent[0] * 2 - 1, tangent[1] * 2 - 1, tangent[2] * 2 - 1));
#if 0
        Spectrum ret;
        ret.fromLinearRGB(std::abs(t.x), std::abs(t.y), std::abs(t.z));
        return ret/M_PI;
#else

		Vector n = normalize(bRec.its.shFrame.n);

		Vector s = normalize(cross(t, n));
		n = normalize(cross(s, t));

		Frame tFrame;
		tFrame.n = n;
		tFrame.t = t;
		tFrame.s = s;

		Vector t_wi = normalize(bRec.its.shFrame.toWorld(bRec.wi));
		Vector t_wo = normalize(bRec.its.shFrame.toWorld(bRec.wo));

		double theta_i = CalTheta(t_wi, t);
		double theta_o = CalTheta(t_wo, t);
		double phi_i = CalPhi(t_wi, t, n);
		double phi_o = CalPhi(t_wo, t, n);

		double cos_item = dot(t_wo, n) > 0 ? dot(t_wo, n):0;

		F[0] = Kd[0] + (1.0 / pow(cos(theta_o), 2.0))*(Ks[0] / sqrt(2 * Beta[0] * Beta[0] * M_PI))*exp(-0.5*pow((theta_i + theta_o) / Beta[0], 2));
		F[1] = Kd[1] + (1.0 / pow(cos(theta_o), 2.0))*(Ks[1] / sqrt(2 * Beta[1] * Beta[1] * M_PI))*exp(-0.5*pow((theta_i + theta_o) / Beta[1], 2));
		F[2] = Kd[2] + (1.0 / pow(cos(theta_o), 2.0))*(Ks[2] / sqrt(2 * Beta[2] * Beta[2] * M_PI))*exp(-0.5*pow((theta_i + theta_o) / Beta[2], 2));

		//return F *(INV_PI * Frame::cosTheta(bRec.wo));
		return cos_item * 1.0 * F;// *(Frame::cosTheta(bRec.wo));
#endif
		}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;

		return warp::squareToCosineHemispherePdf(bRec.wo);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const {
		if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

		bRec.wo = warp::squareToCosineHemisphere(_sample);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		pdf = warp::squareToCosineHemispherePdf(bRec.wo);
		//return m_reflectance->eval(bRec.its);
		return eval(bRec, ESolidAngle) / pdf;
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
		Float pdf;
		return sample(bRec, pdf, _sample);
	}

	//void addChild(const std::string &name, ConfigurableObject *child) {
	//	if (child->getClass()->derivesFrom(MTS_CLASS(Texture))
	//			&& (name == "reflectance" || name == "diffuseReflectance")) {
	//		m_reflectance = static_cast<Texture *>(child);
	//	} else {
	//		BSDF::addChild(name, child);
	//	}
	//}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "reflectance" || name == "diffuseReflectance")
				m_reflectance = static_cast<Texture *>(child);
			else if (name == "tangentMap")
				m_tangentMap = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		}
		else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_reflectance.get());
	}

	Float getRoughness(const Intersection &its, int component) const {
		return std::numeric_limits<Float>::infinity();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothDiffuse[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  reflectance = " << indent(m_reflectance->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_reflectance;
	ref<Texture> m_tangentMap;
	ref<Texture> m_Kd;
	ref<Texture> m_Ks;
	ref<Texture> m_Beta;
};

// ================ Hardware shader implementation ================

class HeightMapBSDFShader : public Shader {
public:
	HeightMapBSDFShader(Renderer *renderer, const Texture *reflectance)
		: Shader(renderer, EBSDFShader), m_reflectance(reflectance) {
		m_reflectanceShader = renderer->registerShaderForResource(m_reflectance.get());
	}

	bool isComplete() const {
		return m_reflectanceShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_reflectance.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_reflectanceShader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * inv_pi * cosTheta(wo);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_reflectance;
	ref<Shader> m_reflectanceShader;
};

Shader *HeightMapBSDF::createShader(Renderer *renderer) const {
	return new HeightMapBSDFShader(renderer, m_reflectance.get());
}

MTS_IMPLEMENT_CLASS(HeightMapBSDFShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(HeightMapBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(HeightMapBSDF, "Custom BSDF for Height Maps")
MTS_NAMESPACE_END