#ifndef iGame_Arap_H
#define iGame_Arap_H

#include "iGameFilter.h"
#include "iGameObject.h"
#include "iGamePoints.h"
#include "iGameSurfaceMesh.h"

IGAME_NAMESPACE_BEGIN
class Arap : public Filter {
public:
    I_OBJECT(Arap);
    static Pointer New() { return new Arap; }
    bool Execute() override;
    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }

protected:
    SurfaceMesh::Pointer mesh{};
    Arap() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
    };
    Arap(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }
};
IGAME_NAMESPACE_END
#endif