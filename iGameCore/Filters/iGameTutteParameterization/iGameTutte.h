#ifndef iGame_Tutte_H
#define iGame_Tutte_H

#include "iGameFilter.h"
#include "iGameObject.h"
#include "iGamePoints.h"
#include "iGameSurfaceMesh.h"

IGAME_NAMESPACE_BEGIN
class Tutte : public Filter {
public:
    I_OBJECT(Tutte);
    static Pointer New() { return new Tutte; }
    bool Execute() override;
    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }

protected:
    SurfaceMesh::Pointer mesh{};
    Tutte() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
    };
    Tutte(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }

};
IGAME_NAMESPACE_END
#endif