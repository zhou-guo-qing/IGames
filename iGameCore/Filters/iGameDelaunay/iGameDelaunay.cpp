#include "iGameDelaunay.h"
#include "iGameVector.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#define PI 3.1415926
IGAME_NAMESPACE_BEGIN

#define M_PI 3.14159265358979323846


double Delaunay::get_TriFace_Area(int i) {

    auto faces = mesh->GetFaces();
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    vcnt = faces->GetCellIds(i, vhs);
    assert(vcnt == 3);
    Point v0 = mesh->GetPoint(vhs[0]);
    Point v1 = mesh->GetPoint(vhs[1]);
    Point v2 = mesh->GetPoint(vhs[2]);


    return 0.5*((v1 - v0).cross(v2 - v0)).norm();
}
//得到外心
Point Delaunay::get_TriFace_Circumcenter(int i) {
    auto faces = mesh->GetFaces();
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    vcnt = faces->GetCellIds(i, vhs);
    assert(vcnt == 3);
    Point v0 = mesh->GetPoint(vhs[0]);
    Point v1 = mesh->GetPoint(vhs[1]);
    Point v2 = mesh->GetPoint(vhs[2]);

    double x1, y1,z1, x2, y2,z2, x3, y3,z3;
    x1 = v0[0];
    y1 = v0[1];
    z1 = v0[2];
    x2 = v1[0];
    y2 = v1[1];
    z2 = v1[2];
    x3 = v2[0];
    y3 = v2[1];
    z3 = v2[2];

    double a1 = (y1 * z2 - y2 * z1 - y1 * z3 + y3 * z1 + y2 * z3 - y3 * z2);
    double b1 = -(x1 * z2 - x2 * z1 - x1 * z3 + x3 * z1 + x2 * z3 - x3 * z2);
    double c1 = (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
    double d1 = -(x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 +
                  x3 * y1 * z2 - x3 * y2 * z1);

    double a2 = 2 * (x2 - x1);
    double b2 = 2 * (y2 - y1);
    double c2 = 2 * (z2 - z1);
    double d2 = x1 * x1 + y1 * y1 + z1 * z1 - x2 * x2 - y2 * y2 - z2 * z2;

    double a3 = 2 * (x3 - x1);
    double b3 = 2 * (y3 - y1);
    double c3 = 2 * (z3 - z1);
    double d3 = x1 * x1 + y1 * y1 + z1 * z1 - x3 * x3 - y3 * y3 - z3 * z3;

    Point circumcenter(0.0, 0.0, 0.0);
    circumcenter[0] = -(b1 * c2 * d3 - b1 * c3 * d2 - b2 * c1 * d3 +
                       b2 * c3 * d1 + b3 * c1 * d2 - b3 * c2 * d1) /
                     (a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 +
                      a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1);
    circumcenter[1] = (a1 * c2 * d3 - a1 * c3 * d2 - a2 * c1 * d3 +
                      a2 * c3 * d1 + a3 * c1 * d2 - a3 * c2 * d1) /
                     (a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 +
                      a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1);
    circumcenter[2] = -(a1 * b2 * d3 - a1 * b3 * d2 - a2 * b1 * d3 +
                       a2 * b3 * d1 + a3 * b1 * d2 - a3 * b2 * d1) /
                     (a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 +
                      a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1);

    return circumcenter;
}

bool Delaunay::is_filp_ok_triangle(IGsize e)
{ 
    
    if (mesh->IsBoundaryEdge(e)) return false;

    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];                                                                                               
    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];
    igIndex vvcnt = 0;
    igIndex vvhs[IGAME_CELL_MAX_SIZE];   

    vcnt = mesh->GetEdgePointIds(e, vhs);
    assert(vcnt == 2);

    fcnt = mesh->GetEdgeToNeighborFaces(e, fhs);
    assert(fcnt == 2);

    int va[2];
    for (int i = 0; i < fcnt; i++) {

        igIndex v[3];
        mesh->GetFacePointIds(fhs[i], v);
        va[i] = v[0] + v[1] + v[2] - vhs[0] - vhs[1];

    }
    if (va[0] == va[1]) return false;

    vvcnt = mesh->GetPointToOneRingPoints(va[0], vvhs);
    for (int i = 0; i < vvcnt; i++) {

        if (vvhs[i] == va[1]) return false;
    }
    return true;

}

bool Delaunay::Execute()
{
    std::cout << "Delaunay start \n";
    UpdateProgress(0);
    mesh = DynamicCast<SurfaceMesh>(GetInput(0));

    if (mesh == nullptr) {
        std::cout << " mesh is null  \n";
        return false;
    }

    std::cout << "mesh gets \n";

    int iter_num = 100 ;
    int eNum = mesh->GetNumberOfEdges();

    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];

    igIndex vvcnt = 0;
    igIndex vvhs[IGAME_CELL_MAX_SIZE];

    igIndex vfcnt = 0;
    igIndex vfhs[IGAME_CELL_MAX_SIZE];
    for (int i = 0; i < iter_num; i++) {
        std::cout << "iter_num:" << i << std::endl;
        mesh->RequestEditStatus();
        for (int j = 0; j < eNum; j++) {

            auto edges = mesh->GetEdges();
            auto faces = mesh->GetFaces();
            if (mesh->IsBoundaryEdge(j) || !is_filp_ok_triangle(j)) continue;

            vcnt = edges->GetCellIds(j, vhs);

            Point v1 = mesh->GetPoint(vhs[0]);
            Point v2 = mesh->GetPoint(vhs[1]);

            
            int v[2];
            fcnt = mesh->GetEdgeToNeighborFaces(j, fhs);
            assert(fcnt == 2);
            for (int t = 0; t < fcnt; t++) {

                vvcnt = mesh->GetFacePointIds(fhs[t], vvhs);
                v[t] = vvhs[0] + vvhs[1] + vvhs[2] - vhs[0] - vhs[1];
                
            }
            Point v3 = mesh->GetPoint(v[0]);
            Point v4 = mesh->GetPoint(v[1]);

            double alpha(0.0), alpha1(0.0), alpha2(0.0);
            alpha1 = acos((pow((v1 - v3).norm(), 2) + pow((v2 - v3).norm(), 2) -
                           pow((v1 - v2).norm(), 2)) /
                          (2 * (v1 - v3).norm() * (v2 - v3).norm()));
            alpha2 = acos((pow((v1 - v4).norm(), 2) + pow((v2 - v4).norm(), 2) -
                           pow((v1 - v2).norm(), 2)) /
                          (2 * (v1 - v4).norm() * (v2 - v4).norm()));
            alpha = alpha1 + alpha2;
            if (alpha > M_PI) {

                igIndex face[3]{vhs[0], v[0], v[1]};
                igIndex face2[3]{vhs[1], v[0], v[1]};

                mesh->DeleteEdge(j);
                mesh->AddFace(face, 3);
                mesh->AddFace(face2, 3);
            }
        }
        int vNum = mesh->GetNumberOfPoints();
        for (int j = 0; j < vNum; j++) {

            if (mesh->IsBoundaryPoint(j)) continue;

            Point tmp(0.0, 0.0, 0.0);
            double area(0.0), sum_area(0.0);
            vfcnt = mesh->GetPointToNeighborFaces(j, vfhs);

            for (int k = 0; k < vfcnt; k++) {

                area = get_TriFace_Area(vfhs[k]);
                sum_area += area;

                Point center = get_TriFace_Circumcenter(vfhs[k]);
                tmp = tmp +  center*area;
            }
            mesh->SetPoint(j, tmp / sum_area);
        }
        mesh->GarbageCollection();
    }
    std::cout << "Delaunay end" << std::endl;

    SetOutput(mesh); // 设置输出数据
    return true; 

}

IGAME_NAMESPACE_END