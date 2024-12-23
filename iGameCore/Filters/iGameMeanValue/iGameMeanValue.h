#ifndef iGame_MeanValue_H
#define iGame_MeanValue_H

#include "iGameFilter.h"
#include "iGameObject.h"
#include "iGamePoints.h"
#include "iGameSurfaceMesh.h"

IGAME_NAMESPACE_BEGIN
class MeanValue : public Filter {
public:
    I_OBJECT(MeanValue);
    static Pointer New() { return new MeanValue; }
    bool Execute() override;
    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }

protected:
    SurfaceMesh::Pointer mesh{};
    MeanValue() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
    };
    MeanValue(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }
};
IGAME_NAMESPACE_END
#endif