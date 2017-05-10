#include "StdAfx.h"
#include "MeshHelper.h"
#include <algorithm>
#include <set>

#include "oldobj/GeoObject.h"

static void ImportWavefrontObjFile(LPCTSTR fn, wglh::CwglMesh& mesh)
{
    using namespace num;

	GeoObj::CGeoObject obj;
	if(obj.LoadOBJModel(CString(fn)))
	{
		int i;
		mesh.Destroy();

		mesh.SetMeshRenderMode(GL_TRIANGLES);
		mesh.SetSize_MeshElementIndic(obj.m_nTriangles*3);

		if( !obj.m_bBindTexture )
		{	
			// no texture binded, just copy vertex array
			// Triangle list
			UINT * p = mesh.GetBuffer_MeshElementIndic();
			for(i=0;i<obj.m_nTriangles;i++,p+=3)
			{	memcpy(p,obj.m_pTriangles[i].m_nVertexIndices,sizeof(UINT)*3);
			}
			// vertices array
			mesh.SetSize_Vertex(obj.m_nVertices);
			{
				Vec3f * Pos = mesh.GetBuffer_VertexPosition();
				
				for(i=0;i<obj.m_nVertices;i++,Pos++)
				{
					Pos->CopyFrom((float*)&obj.m_pVertices[i].m_vPosition);
				}
			}
			if(obj.m_bBindNormal)
			{
				Vec3f * Nor = mesh.GetBuffer_VertexNormal();
				for(i=0;i<obj.m_nVertices;i++,Nor++)
				{
					Nor->CopyFrom((float*)&obj.m_pVertices[i].m_vN);
					Nor->Normalize();
				}
			}
		}
		else
		{	// Split shared vertex if texCoord is inconsistent 
			struct VBO_Vertex
			{	UINT	index_to_ObjFile;
				float	TextureCoord[3];
				float	Tangent[3];
				float	Normal[3];
				bool operator <(const VBO_Vertex&x) const
				{	
					if(index_to_ObjFile<x.index_to_ObjFile)
						return true;
					else if(index_to_ObjFile > x.index_to_ObjFile)
						return false;
					else
					{	if( x.TextureCoord[0] - TextureCoord[0] > EPSILON )
							return true;
						else if( x.TextureCoord[0] - TextureCoord[0] < -EPSILON )
							return false;
						if( x.TextureCoord[1] - TextureCoord[1] > EPSILON )
							return true;
						else if( x.TextureCoord[1] - TextureCoord[1] < -EPSILON )
							return false;
						if( x.TextureCoord[2] - TextureCoord[2] > EPSILON )
							return true;
						else if( x.TextureCoord[2] - TextureCoord[2] < -EPSILON )
							return false;
					}
					return false; 
				}
				UINT	push_index;
			};

			std::set<VBO_Vertex>	Vbo_Vertexs;
			UINT submitted_count = 0;
			UINT * pIndex;
			_CheckDump(_T("\tDetecting shared vertices with inconsistent texture coordinate.\n"));
			{	VBO_Vertex			NewAdd;
				UINT * p = mesh.GetBuffer_MeshElementIndic();
				pIndex = p;

				UINT i=0;
				w32::CTimeMeasureDisplay<>	prog(&i,obj.m_nTriangles);
				for(;i<(UINT)obj.m_nTriangles;i++)
				{
					for(UINT j=0;j<3;j++,p++)
					{
						int ind = obj.m_pTriangles[i].m_nVertexIndices[j];
						Vec3f tex = (Vec3f&)obj.m_pTriangles[i].m_vTexture[j];
						Vec3f tang = (Vec3f&)obj.m_pTriangles[i].m_vS[j];
						Vec3f Norm = (Vec3f&)obj.m_pTriangles[i].m_vN[j];
						
						NewAdd.index_to_ObjFile = ind;
						tex.CopyTo(NewAdd.TextureCoord);
						tang.CopyTo(NewAdd.Tangent);
						Norm.CopyTo(NewAdd.Normal);
						NewAdd.push_index = submitted_count;

						std::pair<std::set<VBO_Vertex>::iterator, bool> & ret = Vbo_Vertexs.insert(NewAdd);

						if( ret.second ) //newly added
							*p = submitted_count++;
						else
							*p = ret.first->push_index;
					}
				}
			}
			ASSERT(submitted_count == Vbo_Vertexs.size());
			_CheckDump("\tVertex count increased from "<<obj.m_nVertices<<" to "<<submitted_count<<".\n");
			mesh.SubmitBuffers(TRUE);

			mesh.SetSize_Vertex(submitted_count);
			
			// Position
			{	Vec3f* pPos = mesh.GetBuffer_VertexPosition();
				Vec3f* pNorm = mesh.GetBuffer_VertexNormal();
				Vec3f* pTang = mesh.GetBuffer_VertexTangent();
				Vec3f* pCoord = mesh.GetBuffer_VertexTexCoord();
				ASSERT(pPos && pNorm && pCoord);
				std::set<VBO_Vertex>::iterator	pVbo = Vbo_Vertexs.begin();
				for(;pVbo!=Vbo_Vertexs.end();pVbo++)
				{
					pPos[pVbo->push_index] = (Vec3f&)obj.m_pVertices[pVbo->index_to_ObjFile].m_vPosition;
					pNorm[pVbo->push_index] = (Vec3f&)obj.m_pVertices[pVbo->index_to_ObjFile].m_vN;
					pTang[pVbo->push_index] = (Vec3f&)obj.m_pVertices[pVbo->index_to_ObjFile].m_vS;
					//pNorm[pVbo->push_index].CopyFrom(pVbo->Normal);
					//pTang[pVbo->push_index].CopyFrom(pVbo->Tangent);
					pCoord[pVbo->push_index].CopyFrom(pVbo->TextureCoord);
				}
			}
		}

		//mesh.NormalizePosition();
		mesh.SubmitBuffers(TRUE);
	}
}


CMeshHelper::CMeshHelper(LPCTSTR obj_in, bool format)
{
    bool done = false;

    if ( format )
    {
        ImportWavefrontObjFile(obj_in, mesh);
        done = true;
    }
    else
    {
        done = SUCCEEDED(mesh.Load(obj_in, wglh::CwglMesh::LOADING_FLAG_BUFFER_IN_MEMORY | wglh::CwglMesh::LOADING_FLAG_KEEP_OBJECT_SIZE));
    }

    done = done && mesh.IsNormalBinded() && mesh.IsTextureBinded() && mesh.IsTangentBinded();
    if ( done )
    {
        float area = 0.0f;

        n = mesh.GetVertexCount();
        m = mesh.GetMeshElementIndicCount();
        positions = mesh.GetBuffer_VertexPosition();
        normals = mesh.GetBuffer_VertexNormal();
        tangents = mesh.GetBuffer_VertexTangent();
        texcoords = mesh.GetBuffer_VertexTexCoord();
        faces = mesh.GetBuffer_MeshElementIndic();

        for ( UINT i = 0; i < m; i += 3 )
        {
            USING_EXPRESSION;

            UINT a = faces[i], b = faces[i + 1], c = faces[i + 2];
            num::Vec3f e0, e1, e2;
            e0 = positions[b] - positions[a];
            e1 = positions[c] - positions[a];

            e2.CrossProduct(e0, e1);
            area += 0.5f*sqrt(e2.L2Norm_Sqr());
        }

        printf("Total Surface Area = %.4f, Edge Length = %.4f\n", area, sqrt(area));
    }
    else
        printf("Badness\n");
}


CMeshHelper::~CMeshHelper(void)
{
    mesh.Destroy();
}


void CMeshHelper::GenerateShellMap(float step, LPCTSTR fname) const
{
    UINT i, j, k;

    _MAP map;
    map.clear();
    for ( i = 0; i < m; i += 3 )
    {
        UINT a = faces[i], b = faces[i + 1], c = faces[i + 2];
        map.insert(_MAP::value_type(_PAIR(a, b), i/3));
        map.insert(_MAP::value_type(_PAIR(b, c), i/3));
        map.insert(_MAP::value_type(_PAIR(c, a), i/3));
    }

    _FLAG_MAP flag;
    rt::Buffer<bool> used;

    flag.clear();
    used.SetSize(m/3);

    for ( i = 0; i < m; i += 3 )
    {
        UINT a = faces[i], b = faces[i + 1], c = faces[i + 2];
        char v0, v1, v2;
        _FLAG_MAP::const_iterator it;

        it = flag.find(_PAIR(b, a));
        if ( it != flag.end() )
            v0 = ( it->second ? -1 : 1 );
        else
            v0 = 0;

        it = flag.find(_PAIR(c, b));
        if ( it != flag.end() )
            v1 = ( it->second ? -1 : 1 );
        else
            v1 = 0;

        it = flag.find(_PAIR(a, c));
        if ( it != flag.end() )
            v2 = ( it->second ? -1 : 1 );
        else
            v2 = 0;

        if ( abs(v0 + v1 + v2) == 3 )
        {
            _PAIR p;
            switch ( rand()%3 )
            {
            case 0:
                p.first = b; p.second = a; v0 = -v0; break;
            case 1:
                p.first = c; p.second = b; v1 = -v1; break;
            case 2:
                p.first = a; p.second = c; v2 = -v2; break;
            }

            used.Set(false);
            used[i/3] = true;

            flag.insert(_FLAG_MAP::value_type(_PAIR(a, b), v0 > 0));
            flag.insert(_FLAG_MAP::value_type(_PAIR(b, c), v1 > 0));
            flag.insert(_FLAG_MAP::value_type(_PAIR(c, a), v2 > 0));
            resolveConflict(used, map, flag, map[p], p);
        }
        else
        {
            UINT ord[] = {1, 2, 3, 4, 5, 6};
            std::random_shuffle(ord, ord + 6);
            for ( j = 0; j < 6; ++j )
            {
                k = ord[j];
                if ( (v0 == 0 || (v0 + 1) >> 1 == (k&1)) && (v1 == 0 || (v1 + 1) >> 1 == (k&2) >> 1) && (v2 == 0 || (v2 + 1) >> 1 == (k&4) >> 2) )
                {
                    flag.insert(_FLAG_MAP::value_type(_PAIR(a, b), (k&1) > 0));
                    flag.insert(_FLAG_MAP::value_type(_PAIR(b, c), (k&2) > 0));
                    flag.insert(_FLAG_MAP::value_type(_PAIR(c, a), (k&4) > 0));

                    break;
                }
            }
        }
    }


#ifdef _DEBUG
    for ( i = 0; i < m; i += 3 )
    {
        UINT a = faces[i], b = faces[i + 1], c = faces[i + 2];

        bool v0 = flag[_PAIR(a, b)];
        bool v1 = flag[_PAIR(b, c)];
        bool v2 = flag[_PAIR(c, a)];

        if ( v0 == v1 && v0 == v2 )
            __debugbreak();

        _FLAG_MAP::const_iterator it;

        it = flag.find(_PAIR(b, a));
        if ( it != flag.end() && it->second == v0 )
            __debugbreak();

        it = flag.find(_PAIR(c, b));
        if ( it != flag.end() && it->second == v1 )
            __debugbreak();

        it = flag.find(_PAIR(a, c));
        if ( it != flag.end() && it->second == v2 )
            __debugbreak();
    }
#endif

    /*
    FILE *fout;

    fout = fopen("vertices.dat", "wt");
    for ( i = 0; i < n; ++i )
        fprintf(fout, "%f %f %f\n", positions[i].x, positions[i].y, positions[i].z);
    fclose(fout);

    fout = fopen("edges.dat", "wt");

    _FLAG_MAP::const_iterator it;
    for ( it = flag.begin(); it != flag.end(); ++it )
        fprintf(fout, "%u %u %d\n", it->first.first, it->first.second, it->second);


    fclose(fout);
    */
    rt::BufferEx<num::Vec4u> tetra;
    tetra.SetSize();

    for ( i = 0; i < m; i += 3 )
    {
        UINT a = faces[i], b = faces[i + 1], c = faces[i + 2];
        UINT a0, b0, c0;

        bool v0 = flag[_PAIR(a, b)], v1 = flag[_PAIR(b, c)], v2 = flag[_PAIR(c, a)];
        bool v;

        if ( v0 == v1 )
        {
            v = v0;
            a0 = a; b0 = b; c0 = c;
        }
        else if ( v1 == v2 )
        {
            v = v1;
            a0 = b; b0 = c; c0 = a;
        }
        else
        {
            ASSERT( v2 == v0 );
            v = v2;
            a0 = c; b0 = a; c0 = b;
        }
        ASSERT( flag[_PAIR(a0, b0)] == flag[_PAIR(b0, c0)] );

        if ( v )
        {
            tetra.push_back(num::Vec4u(c0, c0 + n, b0 + n, a0 + n));
            tetra.push_back(num::Vec4u(c0, a0 + n, b0    , b0 + n));
            tetra.push_back(num::Vec4u(c0, a0    , b0    , a0 + n));
        }
        else
        {
            tetra.push_back(num::Vec4u(a0, a0 + n, c0 + n, b0 + n));
            tetra.push_back(num::Vec4u(a0, b0 + n, c0 + n, b0    ));
            tetra.push_back(num::Vec4u(a0, c0 + n, c0    , b0    ));
        }
    }

    /*
    rt::BufferEx<num::Vec3f> tangs;
    //FILE *fin = fopen("tang_plane.dat", "rt");
    FILE *fin = fopen("tang_silk.dat", "rt");
    for ( ; ; )
    {
        num::Vec3f p;
        if ( fscanf(fin, "%f %f %f", &p.x, &p.y, &p.z) != 3 )
            break;
        p.Normalize();
        tangs.push_back(p);
    }
    fclose(fin);
    ASSERT( tangs.GetSize() == static_cast<SIZE_T>(n+n) );
    */
    rt::BufferEx<num::Vec3f> tangs;
    for ( i = 0; i < n ; ++i )
    {
        num::Vec3f p(tangents[i]), q;
        q.CrossProduct(normals[i], tangents[i]);
        tangs.push_back(p);
        tangs.push_back(q);
    }

    FILE *fout = fopen(ATL::CT2A(fname), "wt");
    fprintf(fout, "%u %u\n", n+n, tetra.GetSize());
    for ( i = 0; i < n; ++i )
    {
        fprintf(fout, "%f %f %f\n", positions[i].x, positions[i].y, positions[i].z);
        fprintf(fout, "%f %f 0.0\n", texcoords[i].x, texcoords[i].y);
        fprintf(fout, "%f %f %f\n", normals[i].x, normals[i].y, normals[i].z);
        fprintf(fout, "%f %f %f %f %f %f\n",
            tangs[i + i].x, tangs[i + i].y, tangs[i + i].z,
            tangs[i + i + 1].x, tangs[i + i + 1].y, tangs[i + i + 1].z
        );
    }
    for ( i = 0; i < n; ++i )
    {
        USING_EXPRESSION;
        num::Vec3f p;
        p = positions[i] + step*normals[i];
        fprintf(fout, "%f %f %f\n", p.x, p.y, p.z);
        fprintf(fout, "%f %f 1.0\n", texcoords[i].x, texcoords[i].y);
        fprintf(fout, "%f %f %f\n", normals[i].x, normals[i].y, normals[i].z);
        fprintf(fout, "%f %f %f %f %f %f\n",
            tangs[i + i].x, tangs[i + i].y, tangs[i + i].z,
            tangs[i + i + 1].x, tangs[i + i + 1].y, tangs[i + i + 1].z
        );
    }
    for ( i = 0; i < static_cast<UINT>(tetra.GetSize()); ++i )
        fprintf(fout, "%u %u %u %u\n", tetra[i].x, tetra[i].y, tetra[i].z, tetra[i].w);
    fclose(fout);
}


void CMeshHelper::resolveConflict(rt::Buffer<bool> &used, _MAP &map, _FLAG_MAP &flag, UINT r, const _PAIR &p) const
{
    _FLAG_MAP::iterator it = flag.find(p);
    if ( it == flag.end() ) return;
    if ( it->second != flag[_PAIR(p.second, p.first)] )
        __debugbreak();

    used[r] = true;
    it->second = !it->second;

    UINT a, b, c;
    if ( p.first == faces[3*r] && p.second == faces[3*r + 1] )
    {
        a = faces[3*r + 1]; b = faces[3*r + 2]; c = faces[3*r];
    }
    else if ( p.first == faces[3*r + 1] && p.second == faces[3*r + 2] )
    {
        a = faces[3*r + 2]; b = faces[3*r]; c = faces[3*r + 1];
    }
    else if ( p.first == faces[3*r + 2] && p.second == faces[3*r] )
    {
        a = faces[3*r]; b = faces[3*r + 1]; c = faces[3*r + 2];
    }
    else
        __debugbreak();

    if ( flag[_PAIR(a, b)] != it->second || flag[_PAIR(b, c)] != it->second )
        return;

    _PAIR q0(b, a), q1(c, b), *q;

    if ( map.find(q0) == map.end() || used[map[q0]] )
    {
        if ( map.find(q1) == map.end() || used[map[q1]] )
            __debugbreak();
        q = &q1;
    }
    else if ( map.find(q1) == map.end() || used[map[q1]] )
        q = &q0;
    else if ( rand() % 2 )
        q = &q1;
    else
        q = &q0;

    it = flag.find(_PAIR(q->second, q->first));
    it->second = !it->second;

    resolveConflict(used, map, flag, map[*q], *q);
}


void CMeshHelper::GenerateBoundingMesh(float step, LPCTSTR fname) const
{
    _FLAG_MAP map;
    
    UINT i;
    for ( i = 0; i < m; i += 3 )
    {
        UINT a = faces[i], b = faces[i + 1], c = faces[i + 2];
        map.insert(_FLAG_MAP::value_type(_PAIR(a, b), true));
        map.insert(_FLAG_MAP::value_type(_PAIR(b, c), true));
        map.insert(_FLAG_MAP::value_type(_PAIR(c, a), true));
    }

    FILE *fout = fopen(ATL::CT2A(fname), "wt");

    for ( i = 0; i < n; ++i )
    {
        fprintf(fout, "v %f %f %f\n", positions[i].x, positions[i].y, positions[i].z);
        //fprintf(fout, "vn %f %f %f\n", -normals[i].x, -normals[i].y, -normals[i].z);
    }
    for ( i = 0; i < n; ++i )
    {
        USING_EXPRESSION;
        num::Vec3f p;
        p = positions[i] + step*normals[i];
        fprintf(fout, "v %f %f %f\n", p.x, p.y, p.z);
        //fprintf(fout, "vn %f %f %f\n", normals[i].x, normals[i].y, normals[i].z);
    }

    for ( i = 0; i < m; i += 3 )
    {
        UINT a = faces[i] + 1, b = faces[i + 1] + 1, c = faces[i + 2] + 1;
        fprintf(fout, "f %u %u %u\n", a, c, b);
        fprintf(fout, "f %u %u %u\n", a + n, b + n, c + n);
    }

    for ( i = 0; i < m; i += 3 )
    {
        UINT a = faces[i], b = faces[i + 1], c = faces[i + 2];
        if ( map.find(_PAIR(b, a)) == map.end() )
        {
            fprintf(fout, "f %u %u %u\n", b + 1, a + n + 1, a + 1);
            fprintf(fout, "f %u %u %u\n", b + 1, b + n + 1, a + n + 1);
        }
        if ( map.find(_PAIR(c, b)) == map.end() )
        {
            fprintf(fout, "f %u %u %u\n", c + 1, b + n + 1, b + 1);
            fprintf(fout, "f %u %u %u\n", c + 1, c + n + 1, b + n + 1);
        }
        if ( map.find(_PAIR(a, c)) == map.end() )
        {
            fprintf(fout, "f %u %u %u\n", a + 1, c + n + 1, c + 1);
            fprintf(fout, "f %u %u %u\n", a + 1, a + n + 1, c + n + 1);
        }
    }

    fclose(fout);
}


void CMeshHelper::GenerateSphereMesh(UINT n, LPCTSTR fname)
{
    static const float dx[] = { 0.0f,  0.0f,  0.0f,  0.0f, -1.0f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f};
    static const float dy[] = { 1.0f,  0.0f, -1.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  1.0f,  0.0f, -1.0f};
    static const float dz[] = { 0.0f,  1.0f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f,  0.0f,  0.0f,  0.0f};

    float unit = 1.0f/static_cast<float>(n);
    rt::BufferEx<num::Vec3f> pos, tex;
    rt::BufferEx<num::Vec4u> face;

    UINT i, j, k;
    for ( i = 0; i < 6; ++i )
    {
        num::Vec3f p, q;
        if ( max(dx[i + i], dx[i + i + 1]) > 0.5f )
            p.x = -0.5f;
        else if ( min(dx[i + i], dx[i + i + 1]) < -0.5f )
            p.x = 0.5f;
        else
            p.x = ( i % 2 ? -0.5f : 0.5f );

        if ( max(dy[i + i], dy[i + i + 1]) > 0.5f )
            p.y = -0.5f;
        else if ( min(dy[i + i], dy[i + i + 1]) < -0.5f )
            p.y = 0.5f;
        else
            p.y = ( i % 2 ? -0.5f : 0.5f );

        if ( max(dz[i + i], dz[i + i + 1]) > 0.5f )
            p.z = -0.5f;
        else if ( min(dz[i + i], dz[i + i + 1]) < -0.5f )
            p.z = 0.5f;
        else
            p.z = ( i % 2 ? -0.5f : 0.5f );

        for ( j = 0; j <= n; ++j )
            for ( k = 0; k <= n; ++k )
            {
                USING_EXPRESSION;
                q = p + static_cast<float>(k)*unit*num::Vec3f(dx[i + i], dy[i + i], dz[i + i])
                    + static_cast<float>(j)*unit*num::Vec3f(dx[i + i + 1], dy[i + i + 1], dz[i + i + 1]);

                q.Normalize();
                pos.push_back(q);
                tex.push_back(num::Vec3f(static_cast<float>(k)/n, static_cast<float>(j)/n, 0.0f));
                if ( j < n && k < n )
                {
                    UINT id = static_cast<UINT>(pos.GetSize());
                    face.push_back(num::Vec4u(id, id + 1, id + n + 2, id + n + 1));
                }
            }
    }

    FILE *fout = fopen(ATL::CT2A(fname), "wt");
    for ( i = 0; i < pos.GetSize(); ++i )
    {
        fprintf(fout, "v %.4f %.4f %.4f\n", pos[i].x, pos[i].y, pos[i].z);
        fprintf(fout, "vt %.4f %.4f\n", tex[i].x, tex[i].y);
        fprintf(fout, "vn %.4f %.4f %.4f\n", pos[i].x, pos[i].y, pos[i].z);
    }
    for ( i = 0; i < face.GetSize(); ++i )
        fprintf(fout, "f %u/%u/%u %u/%u/%u %u/%u/%u %u/%u/%u\n",
            face[i].x, face[i].x, face[i].x,
            face[i].y, face[i].y, face[i].y,
            face[i].z, face[i].z, face[i].z,
            face[i].w, face[i].w, face[i].w
        );
    fclose(fout);
}
