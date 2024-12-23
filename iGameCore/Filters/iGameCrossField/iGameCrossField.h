#ifndef iGame_CrossField_H
#define iGame_CrossField_H

#include "iGameFilter.h"
#include "iGameObject.h"
#include "iGamePoints.h"
#include "iGameSurfaceMesh.h"

IGAME_NAMESPACE_BEGIN
class CrossField : public Filter {
public:
    I_OBJECT(CrossField);
    static Pointer New() { return new CrossField; }
    bool Execute() override;
    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }

    void calculateMeshFaceBase(SurfaceMesh* mesh,std::vector<Point> &f_base); //在每个面上算一组基

protected:
    SurfaceMesh::Pointer mesh{};
    CrossField() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
    };
    CrossField(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }
};
IGAME_NAMESPACE_END
#endif