#ifndef iGame_ARAPInterpolation_H
#define iGame_ARAPInterpolation_H

#include "iGameFilter.h"
#include "iGameObject.h"
#include "iGamePoints.h"
#include "iGameSurfaceMesh.h"
IGAME_NAMESPACE_BEGIN
class ARAPInterpolation : public Filter {
public:
    I_OBJECT(ARAPInterpolation);
    static Pointer New() { return new ARAPInterpolation; }
    bool Execute() override;
    // void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }

    void SetT(double t) { m_T = t;}

protected:
    SurfaceMesh::Pointer source{}, target{};
    ARAPInterpolation() {
        source = nullptr;
        target = nullptr;
        m_T = 0;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(2);
    };
    ARAPInterpolation(SurfaceMesh::Pointer _mesh, SurfaceMesh::Pointer _mesh2) {
        this->source = _mesh;
        this->target = _mesh2;
    }
    double m_T;
};
IGAME_NAMESPACE_END
#endif