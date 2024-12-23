#ifndef iGame_Remesh_H
#define iGame_Remesh_H

#include "iGameFilter.h"
#include "iGameObject.h"
#include "iGamePoints.h"
#include "iGameSurfaceMesh.h"

IGAME_NAMESPACE_BEGIN
class Remesh : public Filter {
public:
    I_OBJECT(Remesh);
    static Pointer New() { return new Remesh; }
    bool Execute() override;
    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }
    bool is_collapse_ok(int eid, double low, double high);
    double calculateTargetEdgeLength();
    void split_long_edges(double high);
    void collapse_short_edges(double high, double low);
    void equalize_valences();
    void tangential_relaxation();
    void project_to_surface();
    void get_AABB_tree();

protected:
    SurfaceMesh::Pointer mesh{};
    Remesh() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
    };
    Remesh(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }
};
IGAME_NAMESPACE_END
#endif