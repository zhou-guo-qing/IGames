#ifndef iGame_VSA_H
#define iGame_VSA_H

#include "iGameFilter.h"
#include "iGameObject.h"
#include "iGamePoints.h"
#include "iGameSurfaceMesh.h"

IGAME_NAMESPACE_BEGIN
struct Tri {
    int id;
    int tag;
    double weight;

    Tri(int id, int tag, double weight) : tag(tag), id(id), weight(weight){};

    bool operator<(const Tri& a) const { return weight > a.weight; }
};
class VSA : public Filter {
public:
    I_OBJECT(VSA);
    static Pointer New() { return new VSA; }
    bool Execute() override;
    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }

    std::vector<int> getPartition() { return partition; }
    double cal_area(int i);

protected:
    SurfaceMesh::Pointer mesh{};
    VSA() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
        partition.resize(50000);
        for (int i = 0; i < 50000; ++i) {
            partition[i] = 0.0; // 将所有元素初始化为 0.0
        }
    };
    VSA(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }

    std::vector<int> partition;
};
IGAME_NAMESPACE_END
#endif