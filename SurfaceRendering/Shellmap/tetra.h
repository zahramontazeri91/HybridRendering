#ifdef __TETRAHEDRON2_H
#error "badness"
#endif

#ifndef __TETRAHEDRON_H
#define __TETRAHEDRON_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/stream.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/tls.h>
#include <mitsuba/core/triangle.h>

//#define TETRA_VERBOSE

#ifdef DOUBLE_PRECISION
#define TETRA_EPSILON 1e-7
#else
#define TETRA_EPSILON 1e-4
#endif

MTS_NAMESPACE_BEGIN

struct Tetrahedron
{	
	uint32_t idx[4];

	inline AABB getAABB(const Point *positions) const
	{
		AABB result(positions[idx[0]]);
		result.expandBy(positions[idx[1]]);
		result.expandBy(positions[idx[2]]);
		result.expandBy(positions[idx[3]]);
		return result;
	}

	inline bool inside(const Point *positions, const Point &p, Point4 &bb) const
	{
		Vector3 a(positions[idx[0]] - p);
		Vector3 b(positions[idx[1]] - p);
		Vector3 c(positions[idx[2]] - p);
		Vector3 d(positions[idx[3]] - p);

		Float v = std::abs(dot(a - d, cross(b - d, c - d)));
		Float v0 = std::abs(dot(b, cross(c, d)));
		Float v1 = std::abs(dot(a, cross(c, d)));
		Float v2 = std::abs(dot(a, cross(b, d)));
		Float v3 = std::abs(dot(a, cross(b, c)));

		bb.x = v0/v;
		bb.y = v1/v;
		bb.z = v2/v;
		bb.w = v3/v;

		/*
		return bb.x >= 0 && bb.x <= 1 &&
			bb.y >= 0 && bb.y <= 1 &&
			bb.z >= 0 && bb.z <= 1 &&
			bb.w >= 0 && bb.w <= 1;
		*/

		return bb.x + bb.y + bb.z + bb.w < 1.0f + TETRA_EPSILON;
	}

	inline bool rayIntersect(const Ray &ray, Float &mint, Float &maxt, const Point *positions) const 
	{
		mint = -std::numeric_limits<Float>::infinity();
		maxt  = std::numeric_limits<Float>::infinity();
		bool intersect = false;

		//Matrix4x4 A = worldToBarycentricMatrix(positions);

		FOR(vertexIndex, 4) 
		{
			Float u, v, t;
			const Point &a = positions[idx[(vertexIndex+1) % 4]];
			const Point &b = positions[idx[(vertexIndex+2) % 4]];
			const Point &c = positions[idx[(vertexIndex+3) % 4]];
			const Point &d = positions[idx[(vertexIndex+0) % 4]];
			bool hit = Triangle::rayIntersect(a, b, c, ray, u, v, t);
			if (hit)
			{	
				/*			
				SLog(EInfo, "vertexIndex = %d hit, t = %f!", vertexIndex, t);
				Vector4 aa = A*Vector4(a.x, a.y, a.z, 1);
				SLog(EInfo, "aa = %s", aa.toString().c_str());
				Vector4 bb = A*Vector4(b.x, b.y, b.z, 1);
				SLog(EInfo, "bb = %s", bb.toString().c_str());
				Vector4 cc = A*Vector4(c.x, c.y, c.z, 1);
				SLog(EInfo, "cc = %s", cc.toString().c_str());

				Point hitPoint = (1-u-v)*a + u*b + v*c;
				SLog(EInfo, "hitPoint = %s", hitPoint.toString().c_str());
				Vector4 hh = A*Vector4(hitPoint.x, hitPoint.y, hitPoint.z, 1);
				SLog(EInfo, "hh = %s", hh.toString().c_str());

				Point rayPoint = ray(t);
				SLog(EInfo, "rayPoint = %s", rayPoint.toString().c_str());
				*/
				intersect = true;
				Vector n = cross(b-a, c-a);
				int vertexSide = (dot(d-a, n) >= 0) ? 1 : -1;
				int dirSide = (dot(ray.d, n) >= 0) ? 1 : -1;
				bool goingIn = (vertexSide * dirSide) == 1;
				if (goingIn) 
				{
					mint = std::max(mint, t);
				} 
				else 
				{
					maxt = std::min(maxt, t);
				}
			}
		}

		if (!std::isfinite(mint) || !std::isfinite(maxt))
		{
			return false;
		}
		return intersect && mint <= maxt;

		/*
		if (intersect && mint <= maxt + 1e-5f)
		{
			Float mid = (mint + maxt) / 2;
			Float width = std::max(std::abs(maxt - mid), std::abs(mint - mid));
			maxt = mid + width * 1.001f;
			mint = mid - width * 1.001f;
			return true;
		}
		else
		{
			return false;
		}
		*/
	}
	
	inline Transform worldToTexCoordTransform(const Point *positions, const Point *texCoords) const
	{
		Matrix4x4 A = worldToBarycentricMatrix(positions);
		
		Vector3 ta(texCoords[idx[0]]);
		Vector3 tb(texCoords[idx[1]]);
		Vector3 tc(texCoords[idx[2]]);
		Vector3 td(texCoords[idx[3]]);

		Matrix4x4 B;
		B.m[0][0] = ta.x;
		B.m[1][0] = ta.y;
		B.m[2][0] = ta.z;
		B.m[3][0] = 1;

		B.m[0][1] = tb.x;
		B.m[1][1] = tb.y;
		B.m[2][1] = tb.z;
		B.m[3][1] = 1;

		B.m[0][2] = tc.x;
		B.m[1][2] = tc.y;
		B.m[2][2] = tc.z;
		B.m[3][2] = 1;

		B.m[0][3] = td.x;
		B.m[1][3] = td.y;
		B.m[2][3] = td.z;
		B.m[3][3] = 1;
				
		Matrix4x4 M = B*A;		

		return Transform(M);
	}

	inline Matrix4x4 worldToBarycentricMatrix(const Point *positions) const 
	{
		Vector3 a(positions[idx[0]]);
		Vector3 b(positions[idx[1]]);
		Vector3 c(positions[idx[2]]);
		Vector3 d(positions[idx[3]]);		

		Float V = dot(a-d, cross(b-d, c-d));		

		Vector3 r0 = -cross(d-b, c-b) / V;
		Vector3 r1 = -cross(c-a, d-a) / V;
		Vector3 r2 = -cross(d-a, b-a) / V;
		Vector3 r3 = -cross(b-a, c-a) / V;

		Matrix4x4 M;
		M.m[0][0] = r0.x;
		M.m[1][0] = r1.x;
		M.m[2][0] = r2.x;
		M.m[3][0] = r3.x;

		M.m[0][1] = r0.y;
		M.m[1][1] = r1.y;
		M.m[2][1] = r2.y;
		M.m[3][1] = r3.y;

		M.m[0][2] = r0.z;
		M.m[1][2] = r1.z;
		M.m[2][2] = r2.z;
		M.m[3][2] = r3.z;

		M.m[0][3] = -dot(r0, b);
		M.m[1][3] = -dot(r1, a);
		M.m[2][3] = -dot(r2, a);
		M.m[3][3] = -dot(r3, a);

		return M;
	}
};


class TetrahedronMesh
{
public:
	TetrahedronMesh()
	{		
		m_vertexCount = m_tetrahedronCount = 0;
		m_vtxPosition = m_vtxTexcoord = NULL;
		m_vtxNormal = NULL; m_vtxTangent = NULL;
		m_tetra = NULL;

		m_list = NULL; m_tree = NULL;
		m_treeSize = 0;
	}	

	~TetrahedronMesh()
	{
		if ( m_vtxPosition ) delete[] m_vtxPosition;
		if ( m_vtxTexcoord ) delete[] m_vtxTexcoord;
		if ( m_vtxNormal ) delete[] m_vtxNormal;
		if ( m_vtxTangent ) delete[] m_vtxTangent;
		if ( m_tetra ) delete[] m_tetra;

		if ( m_list ) delete[] m_list;
		if ( m_tree ) delete[] m_tree;
	}

	void load(Stream *stream) 
	{
		m_vertexCount = stream->readUInt();
		m_tetrahedronCount = stream->readUInt();

		m_vtxPosition = new Point[m_vertexCount];
		m_vtxTexcoord = new Point[m_vertexCount];
		m_vtxNormal = new Vector[m_vertexCount];
		m_vtxTangent = new TangentSpace[m_vertexCount];
		m_tetra = new Tetrahedron[m_tetrahedronCount];

		stream->readArray((Float *)m_vtxPosition, 3*m_vertexCount);
		stream->readArray((Float *)m_vtxTexcoord, 3*m_vertexCount);
		stream->readArray((Float *)m_vtxNormal, 3*m_vertexCount);
		stream->readArray((Float *)m_vtxTangent, 6*m_vertexCount);

		stream->readArray((uint32_t *)m_tetra, 4*m_tetrahedronCount);

		/*
		uint32_t i;
		for ( i = 0; i < m_vertexCount; ++i )
		{
			m_vtxPosition[i] = Point(stream);
			m_vtxTexcoord[i] = Point(stream);
			m_vtxNormal[i] = Vector(stream);
			m_vtxTangent[i] = TangentSpace(stream);			
		}
		for ( i = 0; i < m_tetrahedronCount; ++i )
		{
			m_tetra[i].idx[0] = stream->readUInt();
			m_tetra[i].idx[1] = stream->readUInt();
			m_tetra[i].idx[2] = stream->readUInt();
			m_tetra[i].idx[3] = stream->readUInt();
		}
		*/

		buildBvh();
	}

	bool load(const char *fname)
	{
		SLog(EInfo, "fname = %s", fname);
		#define VERIFY_VALUE(x, y)	{ if ( (x) != (y) ) return false; }		

		FILE *fin = fopen(fname, "rt");
		if ( fin )
		{
			VERIFY_VALUE(fscanf(fin, "%u %u", &m_vertexCount, &m_tetrahedronCount), 2);

			m_vtxPosition = new Point[m_vertexCount];
			m_vtxTexcoord = new Point[m_vertexCount];
			m_vtxNormal = new Vector[m_vertexCount];
			m_vtxTangent = new TangentSpace[m_vertexCount];
			m_tetra = new Tetrahedron[m_tetrahedronCount];

			uint32_t i;
			for ( i = 0; i < m_vertexCount; ++i )
			{
#ifdef SINGLE_PRECISION
				VERIFY_VALUE(fscanf(fin, "%f %f %f", &m_vtxPosition[i].x, &m_vtxPosition[i].y, &m_vtxPosition[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%f %f %f", &m_vtxTexcoord[i].x, &m_vtxTexcoord[i].y, &m_vtxTexcoord[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%f %f %f", &m_vtxNormal[i].x, &m_vtxNormal[i].y, &m_vtxNormal[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%f %f %f", &m_vtxTangent[i].dpdu.x, &m_vtxTangent[i].dpdu.y, &m_vtxTangent[i].dpdu.z), 3);
				VERIFY_VALUE(fscanf(fin, "%f %f %f", &m_vtxTangent[i].dpdv.x, &m_vtxTangent[i].dpdv.y, &m_vtxTangent[i].dpdv.z), 3);
#else
				VERIFY_VALUE(fscanf(fin, "%lf %lf %lf", &m_vtxPosition[i].x, &m_vtxPosition[i].y, &m_vtxPosition[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%lf %lf %lf", &m_vtxTexcoord[i].x, &m_vtxTexcoord[i].y, &m_vtxTexcoord[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%lf %lf %lf", &m_vtxNormal[i].x, &m_vtxNormal[i].y, &m_vtxNormal[i].z), 3);
				VERIFY_VALUE(fscanf(fin, "%lf %lf %lf", &m_vtxTangent[i].dpdu.x, &m_vtxTangent[i].dpdu.y, &m_vtxTangent[i].dpdu.z), 3);
				VERIFY_VALUE(fscanf(fin, "%lf %lf %lf", &m_vtxTangent[i].dpdv.x, &m_vtxTangent[i].dpdv.y, &m_vtxTangent[i].dpdv.z), 3);
#endif
			}
			for ( i = 0; i < m_tetrahedronCount; ++i )
				VERIFY_VALUE(fscanf(fin, "%u %u %u %u", &m_tetra[i].idx[0], &m_tetra[i].idx[1], &m_tetra[i].idx[2], &m_tetra[i].idx[3]), 4);

			fclose(fin);

			SLog(EInfo, "building BVH");
			buildBvh();
			
			return true;
		}
		else
			return false;

		#undef VERIFY_VALUE
	}

	void serialize(Stream *stream) 
	{
		//SLog(EInfo, "tetMesh vertex count = %d", m_vertexCount);
		//SLog(EInfo, "tetMesh tetrahedron count  = %d", m_tetrahedronCount);
		
		stream->writeUInt(m_vertexCount);
		stream->writeUInt(m_tetrahedronCount);

		stream->writeArray((Float *)&m_vtxPosition[0], 3*m_vertexCount);
		stream->writeArray((Float *)&m_vtxTexcoord[0], 3*m_vertexCount);
		stream->writeArray((Float *)&m_vtxNormal[0], 3*m_vertexCount);
		stream->writeArray((Float *)&m_vtxTangent[0], 6*m_vertexCount);

		stream->writeArray((uint32_t *)&m_tetra[0], 4*m_tetrahedronCount);

		/*
		for ( size_t i = 0; i < m_vertexCount; ++i )		
		{
			m_vtxPosition[i].serialize(stream);
			m_vtxTexcoord[i].serialize(stream);
			m_vtxNormal[i].serialize(stream);
			m_vtxTangent[i].serialize(stream);
			//if (i % 10000 == 0)
			//	SLog(EInfo, "transferred %d vertices", i);
		}
		for ( size_t i = 0; i < m_tetrahedronCount; ++i ) {
			stream->writeUInt(m_tetra[i].idx[0]);
			stream->writeUInt(m_tetra[i].idx[1]);
			stream->writeUInt(m_tetra[i].idx[2]);
			stream->writeUInt(m_tetra[i].idx[3]);
			//if (i % 10000 == 0)
			//	SLog(EInfo, "transferred %d tetrahedra", i);
		}
		*/
	}

	void buildBvh() 
	{
		m_aabb = AABB();
		uint32_t i;
		for ( i = 0; i < m_tetrahedronCount; ++i )
			m_aabb.expandBy(m_tetra[i].getAABB(m_vtxPosition));
		SLog(EInfo, "m_aabb = %s", m_aabb.toString().c_str());

		// build BVH
		m_list = new uint32_t[m_tetrahedronCount];
		m_tree = new _node_type[2*m_tetrahedronCount];

		AABB *abl = new AABB[m_tetrahedronCount], aabb;
		for ( i = 0; i < m_tetrahedronCount; ++i )
		{
			m_list[i] = i;
			abl[i] = m_tetra[i].getAABB(m_vtxPosition);
			aabb.expandBy(abl[i]);
		}
		std::random_shuffle(m_list, m_list + m_tetrahedronCount);

		_node_type &root = m_tree[1];
		root.L = 0; root.R = m_tetrahedronCount;
		root.Lch = root.Rch = -1;
		root.aabb = aabb;
		m_treeSize = 1;
		m_treeDepth = (uint32_t)build(0, 1, abl);
		SAssert( m_treeDepth <= MAX_TREE_DEPTH );

		delete[] abl;
	}

	inline AABB getAABB() const
	{
		return m_aabb;
	}

	bool lookupPoint_BruteForce(const Point &p, Point &tex) const
	{
		Point4 bb;
		for ( uint32_t i = 0; i < m_tetrahedronCount; ++i )
			if ( m_tetra[i].inside(m_vtxPosition, p, bb) )
			{
				tex = m_vtxTexcoord[m_tetra[i].idx[0]]*bb.x
					+ m_vtxTexcoord[m_tetra[i].idx[1]]*bb.y
					+ m_vtxTexcoord[m_tetra[i].idx[2]]*bb.z
					+ m_vtxTexcoord[m_tetra[i].idx[3]]*bb.w;
				return true;
			}
		return false;
	}

	bool lookupPoint(const Point &p, Point &tex) const
	{
		uint32_t id;
		Point4 bb;

		//if ( !lookup(p, 1, id, bb) ) return false;
		bool inside = lookup(p, id, bb);
		if ( !inside ) return false;

		tex = m_vtxTexcoord[m_tetra[id].idx[0]]*bb.x
			+ m_vtxTexcoord[m_tetra[id].idx[1]]*bb.y
			+ m_vtxTexcoord[m_tetra[id].idx[2]]*bb.z
			+ m_vtxTexcoord[m_tetra[id].idx[3]]*bb.w;

		return true;
	}

	bool lookupPoint(const Point &p, Point &tex, Vector &normal, TangentSpace &tangent) const
	{
		uint32_t id;
		Point4 bb;

		//if ( !lookup(p, 1, id, bb) ) return false;
		if ( !lookup(p, id, bb) ) return false;

#ifdef TETRA_VERBOSE
		SLog(EInfo, "id = %d", id);
		SLog(EInfo, "bb = %s", bb.toString().c_str());
		SLog(EInfo, "m_vtxPosition[0] = %s", m_vtxPosition[m_tetra[id].idx[0]].toString().c_str());
		SLog(EInfo, "m_vtxPosition[1] = %s", m_vtxPosition[m_tetra[id].idx[1]].toString().c_str());
		SLog(EInfo, "m_vtxPosition[2] = %s", m_vtxPosition[m_tetra[id].idx[2]].toString().c_str());
		SLog(EInfo, "m_vtxPosition[3] = %s", m_vtxPosition[m_tetra[id].idx[3]].toString().c_str());
		SLog(EInfo, "m_vtxTexCoord[0] = %s", m_vtxTexcoord[m_tetra[id].idx[0]].toString().c_str());
		SLog(EInfo, "m_vtxTexCoord[1] = %s", m_vtxTexcoord[m_tetra[id].idx[1]].toString().c_str());
		SLog(EInfo, "m_vtxTexCoord[2] = %s", m_vtxTexcoord[m_tetra[id].idx[2]].toString().c_str());
		SLog(EInfo, "m_vtxTexCoord[3] = %s", m_vtxTexcoord[m_tetra[id].idx[3]].toString().c_str());
#endif

		tex = m_vtxTexcoord[m_tetra[id].idx[0]]*bb.x
			+ m_vtxTexcoord[m_tetra[id].idx[1]]*bb.y
			+ m_vtxTexcoord[m_tetra[id].idx[2]]*bb.z
			+ m_vtxTexcoord[m_tetra[id].idx[3]]*bb.w;
		normal = m_vtxNormal[m_tetra[id].idx[0]]*bb.x
			   + m_vtxNormal[m_tetra[id].idx[1]]*bb.y
			   + m_vtxNormal[m_tetra[id].idx[2]]*bb.z
			   + m_vtxNormal[m_tetra[id].idx[3]]*bb.w;
		tangent.dpdu = m_vtxTangent[m_tetra[id].idx[0]].dpdu*bb.x
					 + m_vtxTangent[m_tetra[id].idx[1]].dpdu*bb.y
					 + m_vtxTangent[m_tetra[id].idx[2]].dpdu*bb.z
					 + m_vtxTangent[m_tetra[id].idx[3]].dpdu*bb.w;
		tangent.dpdv = m_vtxTangent[m_tetra[id].idx[0]].dpdv*bb.x
					 + m_vtxTangent[m_tetra[id].idx[1]].dpdv*bb.y
					 + m_vtxTangent[m_tetra[id].idx[2]].dpdv*bb.z
					 + m_vtxTangent[m_tetra[id].idx[3]].dpdv*bb.w;

		return true;
	}

	inline uint32_t getTreeSize() const
	{
		return m_treeSize;
	}

	inline uint32_t getTreeDepth() const
	{
		return m_treeDepth;
	}

    inline uint32_t getTetrahedronCount() const
    {
        return m_tetrahedronCount;
    }

public:
	struct _node_type
	{
		_node_type()
		{
			L = R = 0;
			Lch = Rch = -1;
		}

		_node_type(const _node_type& in)
		{
			L = in.L; R = in.R;
			Lch = in.Lch; Rch = in.Rch;
			aabb = in.aabb;
		}

		uint32_t L, R;
		int Lch, Rch;
		AABB aabb;
	};

	//struct _cache_struct
	//{
	//	inline _cache_struct()
	//	{
	//		data = new uint32_t[MAX_TREE_DEPTH + 1];
	//		data[0] = data[1] = 1;
	//	}

	//	~_cache_struct()
	//	{
	//		delete[] data;
	//	}

	//	inline uint32_t &operator[] (size_t idx)
	//	{
	//		return data[idx];
	//	}

	//	uint32_t *data;
	//};

	static const uint32_t INVALID_INDEX = 0xffffffff;
	static const uint32_t MAX_TREE_DEPTH = 200;

	bool buildSplit(uint32_t r, const AABB *abl, int &L, int &R, int axis) {
		Float pos = 0.5f*(m_tree[r].aabb.min[axis] + m_tree[r].aabb.max[axis]);

		while ( L < R )
		{
			while ( L < R && pos - abl[m_list[L]].getCenter()[axis] > -Epsilon ) ++L;
			while ( L < R && abl[m_list[R - 1]].getCenter()[axis] - pos > -Epsilon ) --R;
			if ( L < R ) std::swap(m_list[L], m_list[R - 1]);
		}

		return !(L == static_cast<int>(m_tree[r].L) || L == static_cast<int>(m_tree[r].R));
	}

	size_t build(uint32_t dep, uint32_t r, const AABB *abl)
	{
        if ( dep == MAX_TREE_DEPTH /*|| m_tree[r].R - m_tree[r].L <= 5*/ )
            return 1;
						
		int L, R;		
		L = m_tree[r].L;
		R = m_tree[r].R;
		bool found = buildSplit(r, abl, L, R, m_tree[r].aabb.getLargestAxis());

		if (!found) {
			FOR(axis, 3) {
				L = m_tree[r].L;
				R = m_tree[r].R;
				if (buildSplit(r, abl, L, R, axis)) {
					found = true;
					break;
				}
			}			
		}		

		if (!found)
			return 1;
		
		
		_node_type nodeL, nodeR;
		int i;

		nodeL.L = m_tree[r].L; nodeL.R = L;
		nodeL.Lch = nodeL.Rch = -1;
		for ( i = m_tree[r].L; i < L; ++i )
			nodeL.aabb.expandBy(abl[m_list[i]]);
		nodeR.L = L; nodeR.R = m_tree[r].R;
		nodeR.Lch = nodeR.Rch = -1;
		for ( i = static_cast<int>(L); i < static_cast<int>(m_tree[r].R); ++i )
			nodeR.aabb.expandBy(abl[m_list[i]]);

		m_tree[++m_treeSize] = nodeL;
		m_tree[++m_treeSize] = nodeR;
		m_tree[r].Lch = static_cast<int>(m_treeSize) - 1;
		m_tree[r].Rch = static_cast<int>(m_treeSize);
		return std::max(build(dep + 1, m_tree[r].Lch, abl), build(dep + 1, m_tree[r].Rch, abl)) + 1;
	}

#if 0
	size_t build(uint32_t dep, uint32_t r, const AABB *abl)
	{
		if ( dep == MAX_TREE_DEPTH /*|| m_tree[r].R - m_tree[r].L <= 5*/ )
			return 1;
		int axis = m_tree[r].aabb.getLargestAxis();
		Float pos = 0.5f*(m_tree[r].aabb.min[axis] + m_tree[r].aabb.max[axis]);
		int L = m_tree[r].L, R = m_tree[r].R;
		while ( L < R )
		{
			while ( L < R && pos - abl[m_list[L]].getCenter()[axis] > -Epsilon ) ++L;
			while ( L < R && abl[m_list[R - 1]].getCenter()[axis] - pos > -Epsilon ) --R;
			if ( L < R ) std::swap(m_list[L], m_list[R - 1]);
		}
		if ( L == static_cast<int>(m_tree[r].L) || L == static_cast<int>(m_tree[r].R) )
			return 1;
		_node_type nodeL, nodeR;
		int i;
		nodeL.L = m_tree[r].L; nodeL.R = L;
		nodeL.Lch = nodeL.Rch = -1;
		for ( i = m_tree[r].L; i < L; ++i )
			nodeL.aabb.expandBy(abl[m_list[i]]);
		nodeR.L = L; nodeR.R = m_tree[r].R;
		nodeR.Lch = nodeR.Rch = -1;
		for ( i = static_cast<int>(L); i < static_cast<int>(m_tree[r].R); ++i )
			nodeR.aabb.expandBy(abl[m_list[i]]);
		m_tree[++m_treeSize] = nodeL;
		m_tree[++m_treeSize] = nodeR;
		m_tree[r].Lch = static_cast<int>(m_treeSize) - 1;
		m_tree[r].Rch = static_cast<int>(m_treeSize);
		return std::max(build(dep + 1, m_tree[r].Lch, abl), build(dep + 1, m_tree[r].Rch, abl)) + 1;
	}
#endif

	uint32_t lookup(const Point &p, uint32_t r, uint32_t &id, Point4 &bb /*, uint32_t *stack = NULL*/) const
	{
		//if ( stack ) *stack = r;
		if ( !m_tree[r].aabb.contains(p) ) return 0;

		uint32_t j;
		if ( m_tree[r].Lch < 0 )
		{
			for ( uint32_t i = m_tree[r].L; i < m_tree[r].R; ++i )
			{
				j = m_list[i];
				if ( m_tetra[j].inside(m_vtxPosition, p, bb) )
				{
					id = j;
#ifndef TETRAHEDRON_MESH_NO_CACHE
                    m_cache_tetra.get() = j;
#endif
					return 1;
				}
			}
		}
		else
		{
			uint32_t ret;

			j = m_tree[r].Lch;
			ret = lookup(p, j, id, bb /*, stack ? stack + 1 : NULL*/);
			if ( ret ) return ret + 1;

			j = m_tree[r].Rch;
			ret = lookup(p, j, id, bb /*, stack ? stack + 1 : NULL*/);
			if ( ret ) return ret + 1;
		}

		return 0;
	}

	inline bool lookup(const Point &p, uint32_t &id, Point4 &bb) const
	{
#ifndef TETRAHEDRON_MESH_NO_CACHE
		uint32_t val = m_cache_tetra.get();
		if ( val <= m_treeSize && m_tetra[val].inside(m_vtxPosition, p, bb) )
		{
			id = val; return true;
		}
#endif

#if 0
		_cache_struct &cache = m_cache_path.get();
		uint32_t i, j;

		i = cache[0];
		if ( j = lookup(p, cache[i], id, bb, &cache[i]) )
		{
			cache[0] = i + j - 1;
			return true;
		}

		for ( i = cache[0] - 1; i > 0; --i )
		{
			j = cache[i];
			cache[i + 1] = m_tree[j].Lch + m_tree[j].Rch - cache[i + 1];
			if ( j = lookup(p, cache[i + 1], id, bb, &cache[i + 1]) )
			{
				cache[0] = i + j;
				return true;
			}
		}

		cache[0] = cache[1] = 1;
		return false;
#endif
		return lookup(p, 1, id, bb);
	}

	AABB m_aabb;

	uint32_t m_vertexCount, m_tetrahedronCount;
	Point *m_vtxPosition, *m_vtxTexcoord;
	Vector *m_vtxNormal;
	TangentSpace *m_vtxTangent;
	Tetrahedron *m_tetra;

	// BVH
	uint32_t *m_list;
	uint32_t m_treeSize, m_treeDepth;
	_node_type *m_tree;

	// per-thread cache
#ifdef TETRAHEDRON_MESH_NO_CACHE
    #pragma message("Tetrahedron mesh compiled without per-thread cache")
#else
	mutable PrimitiveThreadLocal<uint32_t> m_cache_tetra;
    //mutable PrimitiveThreadLocal<uint32_t> m_cache_hit, m_cache_query;
	//mutable PrimitiveThreadLocal<_cache_struct> m_cache_path;
#endif
};


MTS_NAMESPACE_END


#endif
