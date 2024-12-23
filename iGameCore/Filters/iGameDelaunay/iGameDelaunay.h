#ifndef iGame_Delaunay_H
#define iGame_Delaunay_H

#include "iGameFilter.h"
#include "iGameObject.h"
#include "iGamePoints.h"
#include "iGameSurfaceMesh.h"

IGAME_NAMESPACE_BEGIN
class Delaunay : public Filter {
public:
    I_OBJECT(Delaunay);
    static Pointer New() { return new Delaunay; }
    bool Execute() override;
    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }
    Point get_TriFace_Circumcenter(int i);
    double get_TriFace_Area(int i);
    bool is_filp_ok_triangle(IGsize e);

protected:
    SurfaceMesh::Pointer mesh{};
    Delaunay() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
    };
    Delaunay(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }
};
IGAME_NAMESPACE_END
#endif