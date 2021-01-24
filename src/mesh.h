#ifndef MESH_H
#define MESH_H

#include <string>
#include <array>
#include <optixpp_namespace.h>

struct material
{
    float impedance, attenuation, mu0, mu1, sigma, specularity, shininess, thickness;
};

struct mesh
{
    const std::string filename;
    const bool is_rigid, is_vascular;
    const std::array<float,3> deltas;
    const bool outside_normals;
    const optix::Material & material_inside;
    const optix::Material & material_outside;
};

#endif // MESH_H
