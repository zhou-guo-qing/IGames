#include "iGameQEM.h"
IGAME_NAMESPACE_BEGIN

//计算每条边的损失
void QEM::cal_Cost(int id) {

    Edge_priority temp;

    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    auto edge = mesh->GetEdges();
    vcnt = edge->GetCellIds(id, vhs);
    assert(vcnt == 2);

    Point p0 = mesh->GetPoint(vhs[0]);
    Point p1 = mesh->GetPoint(vhs[1]);

    Eigen::Matrix4d Q_plus = Qv[p0] + Qv[p1], Q_solve = Q_plus;
    Q_solve(3, 0) = 0.0, Q_solve(3, 1) = 0.0, Q_solve(3, 2) = 0.0,
               Q_solve(3, 3) = 1.0;
    Point new_point;
    Eigen::Vector4d new_vec;

    if (Q_solve.determinant() == 0) {
        Vector3f temp = (p0 + p1) * 0.5;
        new_point = {temp[0], temp[1], temp[2]};
        new_vec = {new_point[0], new_point[1], new_point[2], 1.0};
    } else {
        Eigen::Vector4d temp = {0.0, 0.0, 0.0, 1.0};
        new_vec = Q_solve.inverse() * temp;
        new_point = {static_cast<float>(new_vec[0]),
                     static_cast<float>(new_vec[1]),
                     static_cast<float>(new_vec[2])};
    }
    temp.cost = new_vec.transpose() * Q_plus * new_vec;
    temp.eid = id;
    temp.state = State[id];
    temp.NewPoint = new_point;
    Cost.push(temp);
}

//计算每个点的Q
void QEM::cal_Q(int id) {
    int id1 = id;
    Point pv = mesh->GetPoint(id1);
    Eigen::Matrix4d Q_temp;
    Q_temp.setZero();

    Vector3f faceNor;
    int vNum = mesh->GetNumberOfPoints();
    auto faces = mesh->GetFaces();

    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];

    double a, b, c, d;
    Eigen::Matrix<double, 4, 1> p;

    if (!mesh->IsBoundaryPoint(id1)) {
        fcnt = mesh->GetPointToNeighborFaces(id1, fhs);
        for (int i = 0; i < fcnt; i++) {

            auto f = mesh->GetFace(fhs[i]);
            faceNor = f->GetNormal();
            a = faceNor[0], b = faceNor[1], c = faceNor[2],
            d = -faceNor.dot(pv);
            p = {a, b, c, d};
            Q_temp += p * p.transpose();
        }
        Qv[pv] = Q_temp;

    } else {
        vcnt = mesh->GetPointToNeighborEdges(id1, vhs);
        for (int j = 0; j < vcnt; j++) {

            fcnt = mesh->GetEdgeToNeighborFaces(vhs[j], fhs);
            if (fcnt == 0) { break; }
            auto f = mesh->GetFace(fhs[0]);
            Vector3f faceNor = f->GetNormal();
            double a = faceNor[0], b = faceNor[1], c = faceNor[2],
                   d = -faceNor.dot(pv);
            Eigen::Matrix<double, 4, 1> p = {a, b, c, d};
            Q_temp += p * p.transpose();
        }
        Qv[pv] = Q_temp;
    }
}


void QEM::update_Q_and_Cost(int id) {
    
    cal_Q(id);
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];
    vcnt = mesh->GetPointToOneRingPoints(id, vhs);
    for (int i = 0; i < vcnt; i++) { cal_Q(vhs[i]); }

    for (int i = 0; i < vcnt; i++) {

        int eid = mesh->GetEdgeIdFormPointIds(id, vhs[i]);
        cal_Cost(eid);
    }
}
int pv[2]; //另外二个点
bool QEM::QEM_collapse(const Edge_priority& temp_edge) {

    int id = temp_edge.eid;
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];

    igIndex fvcnt = 0;
    igIndex fvhs[IGAME_CELL_MAX_SIZE];
    auto edge = mesh->GetEdges();   
    auto faces = mesh->GetFaces();

    vcnt = edge->GetCellIds(id, vhs);

    Point p0 = mesh->GetPoint(vhs[0]);
    Point p1 = mesh->GetPoint(vhs[1]);

    mesh->SetPoint(vhs[1], temp_edge.NewPoint);

    fcnt = mesh->GetPointToNeighborFaces(vhs[0], fhs); //vhs0的1领域面

    for (int i = 0; i < fcnt; i++) {

        fvcnt = faces->GetCellIds(fhs[i], fvhs);

        bool isflag = false;
        assert(fvcnt == 3);
        int num = 0;
        for (int j = 0; j < fvcnt; j++) {

            if (fvhs[j] == vhs[1]) {
                isflag = true;
                break;
            }
            if (fvhs[j] != vhs[0]) {
                pv[num++] = fvhs[j];
                isflag = false;
            }
        }
        if (!isflag) {
            assert(num == 2);
            igIndex face[3]{vhs[1], pv[0], pv[1]};
            mesh->AddFace(face, 3);
        }
    }
    mesh->DeletePoint(vhs[0]);
    PointDeleted[vhs[0]] = true;
    update_Q_and_Cost(vhs[1]);
    return true;
}
bool QEM::Execute() {
    std::cout << "QEM start \n";
    UpdateProgress(0);
    mesh = DynamicCast<SurfaceMesh>(GetInput(0));

    if (mesh == nullptr) {
        std::cout << " mesh is null  \n";
        return false;
    }
    std::cout << "mesh gets \n";
    mesh->RequestEditStatus();

    //初始化Q
    int vNum = mesh->GetNumberOfPoints();
    for (int i = 0; i < vNum; i++) { 
        PointDeleted[i] = false;
        cal_Q(i);
    }
    //初始化cost
    int eNum = mesh->GetNumberOfEdges();
    for (int i = 0; i < eNum; i++) {
        cal_Cost(i);
    }
    // simplification
    int N_V = mesh->GetNumberOfPoints();
    int target_num = std::min((int) (0.5 * N_V), 1000);
    while (N_V > target_num) {
        Edge_priority temp_edge = Cost.top();
        Cost.pop();
        int id = temp_edge.eid;
        int vid[2];
        mesh->GetEdgePointIds(id, vid);
        if (!PointDeleted[vid[0]] && !PointDeleted[vid[1]]) {
            if (id != -1) {
                if (QEM_collapse(temp_edge)) { N_V--; }
            }
        }
    }

    std::cout << "QEM end" << std::endl;
    mesh->GarbageCollection();
    SetOutput(mesh); // 设置输出数据
    return true;
}

IGAME_NAMESPACE_END