#include"iGameTutte.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "iGameVector.h"
#define PI 3.1415926
IGAME_NAMESPACE_BEGIN

bool compareVectors(const Vector3f& v1, const Vector3f& v2) {
    for (int i = 0; i < 3; ++i) {
        if (v1[i] != v2[i]) { return false; }
    }
    return true;
}
bool Tutte::Execute()
{

    std::cout << "Tutte Parameterization start \n";
    UpdateProgress(0);
    mesh = DynamicCast<SurfaceMesh>(GetInput(0));

    if (mesh == nullptr) {
        std::cout << " mesh is null  \n";
        return false;
    }

    std::cout << "mesh gets \n";
    mesh->RequestEditStatus();

    int faceNum = mesh->GetNumberOfFaces();
    int vNum = mesh->GetNumberOfPoints();
    auto faces = mesh->GetFaces();
    auto edges = mesh->GetEdges();

    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    //计算面积和
    double area_sum = 0;
    for (int i = 0; i < faceNum; i++) {

        auto f = mesh->GetFace(i);

        vcnt = faces->GetCellIds(i, vhs); //面的三个点id
        assert(vcnt == 3);
        Point p1 = mesh->GetPoint(vhs[0]);
        Point p2 = mesh->GetPoint(vhs[1]);
        Point p3 = mesh->GetPoint(vhs[2]);

        auto p12 = p2 - p1;
        auto p13 = p3 - p1;
        double area = 0.5 * p12.cross(p13).norm();
        area_sum += area;
    }

    //把边界点映射到circle上
    
    igIndex bound[1024];

    int boundaryNum = 0;
    int k = 0;
    for (int i = 0; i < vNum; i++) {
        if (mesh->IsBoundaryPoint(i)) { 
            boundaryNum++;
            k = i;
        }
    }
    std::cout << "boundaryNum:" << boundaryNum << std::endl;
 
    int startid = k;
    Point startp = mesh->GetPoint(k);
    Point endp = {0, 0, 0};
    int lastid = startid;
    int t = 0;
    while (!compareVectors(startp,endp)) {
        bound[t++] = startid;
        
        vcnt = mesh->GetPointToOneRingPoints(startid, vhs);
        for (int i = 0; i < vcnt; i++) {
            if (mesh->IsBoundaryPoint(vhs[i]) && vhs[i]!=lastid) { 
                
                lastid = startid;
                startid=vhs[i];
                endp = mesh->GetPoint(vhs[i]);
                break; 
            }
        }
        
    }
    std::cout << "t:" << t << std::endl;

    double delta_angle = 2 * PI / boundaryNum;
    double area_factor = sqrt(area_sum / PI);
    Eigen::VectorXd position_mesh;
    position_mesh.resize(2 * vNum);
    for (int i = 0; i < boundaryNum; i++) { 
        
        position_mesh[bound[i]] = area_factor * cos(i * delta_angle);
        position_mesh[bound[i] + vNum] = area_factor * sin(-i * delta_angle);
        
    }

    //计算矩阵
    typedef Eigen::Triplet<double> T;
    typedef Eigen::SparseMatrix<double> SMatrix;
    std::vector<T> tripletlist;
    Eigen::VectorXd bu = Eigen::VectorXd::Zero(vNum);
    Eigen::VectorXd bv = Eigen::VectorXd::Zero(vNum);

    for (int i = 0; i < vNum; i++) {

        if (mesh->IsBoundaryPoint(i)) { 
            
            tripletlist.push_back(T(i, i, 1));
            auto p = mesh->GetPoint(i);

            bu(i) = position_mesh[i];
            bv(i) = position_mesh[i + vNum];

        } else {
            vcnt = mesh->GetPointToOneRingPoints(i, vhs);
            for (int j = 0; j < vcnt; j++) {
                tripletlist.push_back(T(i, vhs[j], -1));
            }
            tripletlist.push_back(T(i, i, vcnt));
        }
    }
    SMatrix coff(vNum, vNum);
    coff.setFromTriplets(tripletlist.begin(), tripletlist.end());
    Eigen::SparseLU<SMatrix> solver;
    solver.compute(coff);
    Eigen::VectorXd xu = solver.solve(bu);
    Eigen::VectorXd xv = solver.solve(bv);

    Eigen::MatrixX2d uv;
    uv.resize(vNum, 2);
    uv.col(0) = xu;
    uv.col(1) = xv;
    for (int i = 0; i < vNum; i++) {
        Vector3f tmp = {static_cast<float>(uv(i, 0)),
                         static_cast<float>(uv(i, 1)), 0};
        mesh->SetPoint(i, tmp);
    }

    std::cout << "Tutte Parameterization end" << std::endl;
    mesh->GarbageCollection();

    SetOutput(mesh); // 设置输出数据
    return true; 

}

IGAME_NAMESPACE_END