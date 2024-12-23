#include"iGameARAPInterpolation.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "iGameVector.h"
#include<cmath>
#include <Eigen/SVD>

IGAME_NAMESPACE_BEGIN

#define precision 1e-9
bool ARAPInterpolation::Execute()
{

	std::cout << "ARAPInterpolation start \n";
    UpdateProgress(0);
    source = DynamicCast<SurfaceMesh>(GetInput(0));
    target = DynamicCast<SurfaceMesh>(GetInput(1));


    if (source == nullptr||target==nullptr) {
        std::cout << " mesh is null  \n";
        return false;
    }

    std::cout << "mesh gets \n";
    source->RequestEditStatus();
    target->RequestEditStatus();

    int vNum = source->GetNumberOfPoints();
    int facesNum = source->GetNumberOfFaces();
    auto faces = source->GetFaces();
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    std::vector<Eigen::Matrix2d> S(facesNum);
    std::vector<Eigen::Vector3d> xx(facesNum);
    std::vector<Eigen::Vector3d> yy(facesNum);
    std::vector<double> area(facesNum);
    std::vector<double> angle(facesNum);

    Eigen::Matrix2d I = Eigen::MatrixXd::Identity(2, 2);

    std::vector<std::vector<int>> v_id(facesNum);

    std::cout << "target" << target->GetPoint(0) << std::endl;
    std::cout << "source" << source->GetPoint(0) << std::endl;
    for (int i = 0; i < facesNum; i++) {

        vcnt = faces->GetCellIds(i, vhs);
        assert(vcnt == 3);
        Point p0 = source->GetPoint(vhs[0]);
        Point p1 = source->GetPoint(vhs[1]);
        Point p2 = source->GetPoint(vhs[2]);

        v_id[i].push_back(vhs[0]);
        v_id[i].push_back(vhs[1]);
        v_id[i].push_back(vhs[2]);

        

        area[i] = (p1 - p0).cross(p2 - p0).norm() / 2.0;

        xx[i] = Eigen::Vector3d(p2[0] - p1[0],p0[0] - p2[0],p1[0] - p0[0]);

        yy[i] = Eigen::Vector3d(p1[1] - p2[1],p2[1] - p0[1],p0[1] - p1[1]);
        
        Eigen::Vector3d u(target->GetPoint(vhs[0])[0], target->GetPoint(vhs[1])[0],
                   target->GetPoint(vhs[2])[0]);
        Eigen::Vector3d v(target->GetPoint(vhs[0])[1], target->GetPoint(vhs[1])[1],
                   target->GetPoint(vhs[2])[1]);

        Eigen::Matrix2d J;
        J << yy[i].dot(u) / (2 * area[i]), xx[i].dot(u) / (2 * area[i]),
                yy[i].dot(v) / (2 * area[i]), xx[i].dot(v) / (2 * area[i]);
        Eigen::JacobiSVD<Eigen::Matrix2d> svd(J, Eigen::ComputeFullV | Eigen::ComputeFullU);

        Eigen::Matrix2d V = svd.matrixV();
        Eigen::Matrix2d U = svd.matrixU();
        Eigen::Matrix2d R = U * V.transpose();
        Eigen::Matrix2d sigma;
        sigma(0, 1) = sigma(1, 0) = 0;
        sigma(0, 0) = svd.singularValues()[0];
        sigma(1, 1) = svd.singularValues()[1];
        S[i] = V * sigma * V.transpose();
        angle[i] = atan2(R(1, 0), R(1, 1));
    }
    // fix nv-1
    double u0 = source->GetPoint(vNum - 1)[0];
    double v0 = source->GetPoint(vNum - 1)[1];

    std::cout << "u0" << u0 << std::endl;
    std::cout << "v0" << v0 << std::endl;
    Eigen::SparseMatrix<double> A(2 * vNum - 2, 2 * vNum - 2);
    std::vector<Eigen::Triplet<double>> triplet;

    for (int d= 0; d < facesNum; d++) {


        for (int i = 0; i < 3; i++) {

            for (int j = 0; j < 3; j++) {

                if ((v_id[d][i] == vNum - 1) || (v_id[d][j] == vNum - 1)) {
                    continue;
                } else {
                    triplet.push_back(Eigen::Triplet<double>(
                            2 * v_id[d][i], 2 * v_id[d][j],
                            yy[d][i] * yy[d][j] / (area[d] * area[d] * 2)));
                    triplet.push_back(Eigen::Triplet<double>(
                            2 * v_id[d][i], 2 * v_id[d][j],
                            xx[d][i] * xx[d][j] / (area[d] * area[d] * 2)));
                    triplet.push_back(Eigen::Triplet<double>(
                            2 * v_id[d][i] + 1, 2 * v_id[d][j] + 1,
                            yy[d][i] * yy[d][j] / (area[d] * area[d] * 2)));
                    triplet.push_back(Eigen::Triplet<double>(
                            2 * v_id[d][i] + 1, 2 * v_id[d][j] + 1,
                            xx[d][i] * xx[d][j] / (area[d] * area[d] * 2)));
                }
            }
        }
    }
    A.setFromTriplets(triplet.begin(), triplet.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    double t=0;
    while (t<=1) 
    {
    
        std::cout << "iter:" << t << std::endl;

        Eigen::VectorXd b(2 * vNum - 2);
        b.setZero();

        for (int d = 0; d < facesNum; d++) {

            Eigen::Matrix2d A, R;
            double theta_t = angle[d] * t;

            R << cos(theta_t), -sin(theta_t), sin(theta_t), cos(theta_t);
            A = R * ((1.0 - t) * I + t * S[d]);

            for (int i = 0; i < 3; i++) {

                for (int j = 0; j < 3; j++) {

                    if ((v_id[d][i] == vNum - 1) || (v_id[d][j] == vNum - 1)) {

                        if ((v_id[d][i] != vNum - 1) &&
                            (v_id[d][j] == vNum - 1)) {
                                b[2 * v_id[d][i]] += -yy[d][i] * yy[d][j] * u0 /
                                                      (area[d] * area[d] * 2);
                                b[2 * v_id[d][i]] += -xx[d][i] * xx[d][j] * u0 /
                                                      (area[d] * area[d] * 2);
                                b[2 * v_id[d][i] + 1] += -yy[d][i] * yy[d][j] * v0 /
                                                          (area[d] * area[d] * 2);
                                b[2 * v_id[d][i] + 1] += -xx[d][i] * xx[d][j] * v0 /
                                                          (area[d] * area[d] * 2);
                        }
                     
                    }
                }
                if (v_id[d][i] != vNum - 1) {
                    b[2 * v_id[d][i]] += yy[d][i] * A(0, 0) / area[d];
                    b[2 * v_id[d][i]] += xx[d][i] * A(0, 1) / area[d];
                    b[2 * v_id[d][i] + 1] += yy[d][i] * A(1, 0) / area[d];
                    b[2 * v_id[d][i] + 1] += xx[d][i] * A(1, 1) / area[d];
                }
            }
            
        }
        Eigen::VectorXd result = solver.solve(b);

        for (int i = 0; i < vNum-1; i++) {
            
            if (fabs(t - m_T)< precision) {
                    Vector3f tmp = {static_cast<float>(result[2 * i]),
                                        static_cast<float> (result[2 * i + 1]),
                                        0};
                   source->SetPoint(i, tmp);
            }
            
        }
        t += 0.1;
    }

  
    std::cout << "ARAPInterpolation end" << std::endl;
    source->GarbageCollection();

    SetOutput(source); // 设置输出数据
    return true;
}


IGAME_NAMESPACE_END