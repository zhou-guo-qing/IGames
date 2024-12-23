#ifndef iGame_QEM_H
#define iGame_QEM_H

#include "iGameFilter.h"
#include "iGameObject.h"
#include "iGamePoints.h"
#include "iGameSurfaceMesh.h"
#include "iGameVector.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <map>
#include <queue>
#include <time.h>
IGAME_NAMESPACE_BEGIN

struct Edge_priority {
    double cost;
    int eid;
    int state;

    Vector3f NewPoint;
};
struct cmp {
    bool operator()(const Edge_priority& a, const Edge_priority& b) {
        return a.cost > b.cost;
    }
};
class QEM : public Filter {
public:
    I_OBJECT(QEM);
    static Pointer New() { return new QEM; }
    bool Execute() override;
    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }
    
    void cal_Cost(int id);
    void cal_Q(int id);
    void update_Q_and_Cost(int id);
    bool QEM_collapse(const Edge_priority& temp_edge);


protected:
    SurfaceMesh::Pointer mesh{};
    QEM() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
    };
    QEM(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }
    std::priority_queue<Edge_priority, std::vector<Edge_priority>, cmp> Cost;
    std::map<Point, Eigen::Matrix4d> Qv;
    std::map<int, int> State;
    std::map<int, bool> PointDeleted;
};
IGAME_NAMESPACE_END
#endif