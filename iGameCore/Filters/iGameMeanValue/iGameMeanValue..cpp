#include"iGameMeanValue.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "iGameVector.h"
IGAME_NAMESPACE_BEGIN

#define pi 3.1415926

bool compareVectorsW(const Vector3f& v1, const Vector3f& v2) {
    for (int i = 0; i < 3; ++i) {
        if (v1[i] != v2[i]) { return false; }
    }
    return true;
}
bool MeanValue::Execute() 
{
    std::cout << "MeanValue Parameterization start \n";
    UpdateProgress(0);
    mesh = DynamicCast<SurfaceMesh>(GetInput(0));

    if (mesh == nullptr) {
        std::cout << " mesh is null  \n";
        return false;
    }

    std::cout << "mesh gets \n";
    mesh->RequestEditStatus();

    igIndex vNum = mesh->GetNumberOfPoints();
    igIndex facesNum = mesh->GetNumberOfFaces();
    auto faces = mesh->GetFaces();

    std::vector<int> inner_id(vNum);
    std::vector<double> u(vNum);
    std::vector<double> v(vNum);
    
    
    //映射到圆的边上
    int iter = 0;
    int boundary_num, inner_num=0;
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    int k = 0;
    for (int i = 0; i < vNum; i++) {
        if (mesh->IsBoundaryPoint(i)) { 
            iter++;
            k = i;
        } else {
            inner_id[i] = inner_num;
            inner_num++;
        }
    }
    //inner_num = vNum - iter;
    boundary_num = iter;

    std::cout << "boundaryNum:" << boundary_num << std::endl;
    std::cout << "innerNum:" << inner_num << std::endl;

    igIndex bound[1024];
    int startid = k;
    Point startp = mesh->GetPoint(k);
    Point endp = {0, 0, 0};
    int lastid = startid;
    int t = 0;
    while (!compareVectorsW(startp, endp)) {
        bound[t++] = startid;

        vcnt = mesh->GetPointToOneRingPoints(startid, vhs);
        for (int i = 0; i < vcnt; i++) {
            if (mesh->IsBoundaryPoint(vhs[i]) && vhs[i] != lastid) {

                lastid = startid;
                startid = vhs[i];
                endp = mesh->GetPoint(vhs[i]);
                break;
            }
        }
    }
    for (int i = 0; i < boundary_num; i++) {

        u[bound[i]] = cos(double(2 * pi * i / boundary_num));
        v[bound[i]] = sin(double(2 * pi * i / boundary_num));
    }
    Eigen::SparseMatrix<double> weight(vNum, vNum);  //稀疏矩阵
    std::vector<Eigen::Triplet<double>> triplet; //结构体三元组 行，列，值

    for (int i = 0; i < facesNum; i++) {
        //面的三个点id
        vcnt = faces->GetCellIds(i, vhs);
        assert(vcnt == 3);

        Point p0 = mesh->GetPoint(vhs[0]);
        Point p1 = mesh->GetPoint(vhs[1]);
        Point p2 = mesh->GetPoint(vhs[2]);

        double l2 = (p1 - p0).norm();
        double l1 = (p2 - p0).norm();
        double l0 = (p2 - p1).norm();

        double angle = acos((p2-p0).dot(p1-p0) /(l1 * l2)); //l1和l2的夹角值
                            
        triplet.push_back(Eigen::Triplet<double>(vhs[0], vhs[1],tan(angle * 0.5) / l2));
        triplet.push_back(Eigen::Triplet<double>(vhs[0], vhs[2],tan(angle * 0.5) / l1));

        angle = acos((p0-p1).dot(p2-p1) /(l2 * l0));
                     
        triplet.push_back(Eigen::Triplet<double>(vhs[1], vhs[0],tan(angle * 0.5) / l2));
        triplet.push_back(Eigen::Triplet<double>(vhs[1], vhs[2],tan(angle * 0.5) / l0));

        angle = acos((p0-p2).dot(p1-p2) / (l0 * l1));
                    
        triplet.push_back(Eigen::Triplet<double>(vhs[2], vhs[0],tan(angle * 0.5) / l1));
        triplet.push_back(Eigen::Triplet<double>(vhs[2], vhs[1],tan(angle * 0.5) / l0));

    }
    weight.setFromTriplets(triplet.begin(), triplet.end());

    Eigen::SparseMatrix<double> inner(inner_num, inner_num);
    Eigen::VectorXd ub(inner_num), vb(inner_num);

    triplet.clear();

    for (int i = 0; i < vNum; i++) {

        if (!mesh->IsBoundaryPoint(i)) {

            double tempu = 0;
            double tempv = 0;

            vcnt = mesh->GetPointToOneRingPoints(i, vhs);
            for (int j = 0; j < vcnt; j++) {

                if (mesh->IsBoundaryPoint(vhs[j])) {

                    tempu += u[vhs[j]] * weight.coeff(i, vhs[j]);
                    tempv += v[vhs[j]] * weight.coeff(i, vhs[j]);

                }
                else
                {
                    triplet.push_back(Eigen::Triplet<double>(
                            inner_id[i],inner_id[vhs[j]],-weight.coeff(i, vhs[j])));

                }
                triplet.push_back(Eigen::Triplet<double>(inner_id[i], inner_id[i],
                                                  weight.coeff(i, vhs[j])));
            }
            ub[inner_id[i]] = tempu;
            vb[inner_id[i]] = tempv;
        }
    }
    inner.setFromTriplets(triplet.begin(), triplet.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(inner);
    solver.factorize(inner);
    Eigen::VectorXd result_u = solver.solve(ub);
    Eigen::VectorXd result_v = solver.solve(vb);

    for (int i = 0; i < vNum; i++) {

        if (!mesh->IsBoundaryPoint(i)) {
            
            Vector3f temp = {static_cast<float>(result_u[inner_id[i]]),
                             static_cast<float>(result_v[inner_id[i]]), 0};
            mesh->SetPoint(i, temp);
        } else {

            Vector3f temp = {static_cast<float>(u[i]), static_cast<float>(v[i]),
                             0};
            mesh->SetPoint(i, temp);
        }
    }

    std::cout << "MeanValue Parameterization end" << std::endl;
    mesh->GarbageCollection();

    SetOutput(mesh); // 设置输出数据
    return true; 
}
IGAME_NAMESPACE_END