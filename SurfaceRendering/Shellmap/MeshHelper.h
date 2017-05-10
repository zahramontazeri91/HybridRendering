#pragma once
#include <1.0/wglh/wglMesh.h>
#include <map>


class CMeshHelper
{
public:
    CMeshHelper(LPCTSTR obj_in, bool format = false);
    ~CMeshHelper(void);

    void GenerateShellMap(float step, LPCTSTR fname) const;
    void GenerateBoundingMesh(float step, LPCTSTR fname) const;
    static void GenerateSphereMesh(UINT n, LPCTSTR fname);

protected:

    typedef std::pair<UINT, UINT> _PAIR;
    typedef std::map<_PAIR, UINT> _MAP;
    typedef std::map<_PAIR, bool> _FLAG_MAP;

    wglh::CwglMesh mesh;
    UINT n, m;
    const num::Vec3f *positions, *normals, *texcoords, *tangents;
    const UINT *faces;

    void resolveConflict(rt::Buffer<bool> &used, _MAP &map, _FLAG_MAP &flag, UINT r, const _PAIR &p) const;
};
