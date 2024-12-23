#include"iGameCrossField.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "iGameVector.h"
#include <complex>
#include <cmath>
#include <cstring>
#define PI 3.1415926
using namespace std;
IGAME_NAMESPACE_BEGIN

void CrossField::calculateMeshFaceBase(SurfaceMesh* mesh, std::vector<Point>& f_base)
{
    
    f_base.resize(mesh->GetNumberOfFaces() * 2);

    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];
    auto faces = mesh->GetFaces();
    for (int i = 0; i < mesh->GetNumberOfFaces(); i++) {

        vcnt = faces->GetCellIds(i, vhs); //面的三个点id
        assert(vcnt == 3);
        Point p1 = mesh->GetPoint(vhs[0]);
        Point p2 = mesh->GetPoint(vhs[1]);
        Point p3 = mesh->GetPoint(vhs[2]);

        Point no = (p2 - p1).cross(p3 - p1);

        no.normalize();
        f_base[i * 2] = (p2 - p1).normalized();
        f_base[i * 2 + 1] = no.cross(f_base[i * 2]).normalized();

    }

}
//#include "IOManager.h"
//void writeCrossField(const std::string& filename,
//                     std::vector<Point>& crossfield) {
//    std::fstream ofile(filename.c_str(), std::ios_base::out);
//    int num = crossfield.size();
//    ofile << num << std::endl;
//    for (int i = 0; i < num; i++) {
//        ofile << crossfield[i][0] << " " << crossfield[i][1] << " "
//              << crossfield[i][2] << std::endl;
//    }
//    ofile.close();
//}
// 计算点到平面的投影点，点和平面上的点都通过索引表示
const int X_INDEX = 0;
const int Y_INDEX = 1;
const int Z_INDEX = 2;
// 计算两个三维向量的叉积，向量通过索引表示
std::vector<double> crossProduct(const std::vector<double>& v1,
                                 const std::vector<double>& v2) {
    std::vector<double> result(3, 0);
    result[X_INDEX] = v1[Y_INDEX] * v2[Z_INDEX] - v1[Z_INDEX] * v2[Y_INDEX];
    result[Y_INDEX] = v1[Z_INDEX] * v2[X_INDEX] - v1[X_INDEX] * v2[Z_INDEX];
    result[Z_INDEX] = v1[X_INDEX] * v2[Y_INDEX] - v1[Y_INDEX] * v2[X_INDEX];
    return result;
}
//Point projectPointToPlane(Point point, Point p1, Point p2, Point p3) {
//    // 计算平面的两个向量
//
//    Point AB = p2 - p1;
//    Point AC = p3 - p1;
//
//    // 计算平面的法向量
//    Point normal = AB.cross(AC);
//
//    // 计算平面方程的系数
//    double a = normal[X_INDEX];
//    double b = normal[Y_INDEX];
//    double c = normal[Z_INDEX];
//    double d = -(a * p1[X_INDEX] + b * p1[Y_INDEX] + c * p1[Z_INDEX]);
//
//    // 计算投影系数k
//    double k = -((a * point[X_INDEX] + b * point[Y_INDEX] + c * point[Z_INDEX] +
//                  d)) /
//               (a * a + b * b + c * c);
//
//    // 计算投影点
//    Point projection;
//    projection = point - normal * k;
//
//    return projection;
//}

// 判断点P是否在点P1和P2的连线上
bool isPointOnLine(const Point& P, const Point& P1, const Point& P2) {
    // 计算向量P1P和P1P2
    Point P1P = P - P1;
    Point P1P2 = P2 - P1;

    // 计算向量叉积的模
    double crossProductX = (P1P[1] * P1P2[2]) - (P1P[2] * P1P2[1]);
    double crossProductY = (P1P[2] * P1P2[0]) - (P1P[0] * P1P2[2]);
    double crossProductZ = (P1P[0] * P1P2[1]) - (P1P[1] * P1P2[0]);
    double crossProductMagnitude = std::sqrt(crossProductX * crossProductX +
                                             crossProductY * crossProductY +
                                             crossProductZ * crossProductZ);

    // 如果叉积的模为0，则三点共线
    if (crossProductMagnitude < 1e-9) {
        // 计算向量点积
        double dotProduct =
                (P1P[0] * P1P2[0]) + (P1P[1] * P1P2[1]) + (P1P[2] * P1P2[2]);
        // 计算P1P2的模的平方
        double P1P2MagnitudeSquared =
                (P1P2[0] * P1P2[0]) + (P1P2[1] * P1P2[1]) + (P1P2[2] * P1P2[2]);

        // 如果点积大于0且小于P1P2的模的平方，则点P在P1和P2之间
        return dotProduct >= 0 && dotProduct <= P1P2MagnitudeSquared;
    }

    return false;
}
bool CrossField::Execute()
{

	std::cout << "CrossField start \n";
    UpdateProgress(0);
    mesh = DynamicCast<SurfaceMesh>(GetInput(0));

    if (mesh == nullptr) {
        std::cout << " mesh is null  \n";
        return false;
    }

    std::cout << "mesh gets \n";
    mesh->RequestEditStatus();

    std::vector<int> cons_id(0);
    std::vector<Point> cons_vec(0), crossfield(0);

    cons_id.push_back(2746);
    cons_vec.push_back({0.004945, 0.0141348, - 0.004092});

    cons_id.push_back(8692);
    cons_vec.push_back({0.011338, - 0.0192, - 0.00949099999999999});

    typedef complex<double> COMPLEX;

    int fnum = mesh->GetNumberOfFaces();
    std::vector<int> status(fnum, 0);
    std::vector<COMPLEX> f_dir(fnum);

    crossfield.clear();
    crossfield.resize(fnum);
    std::vector<Point> f_base(fnum * 2);
    calculateMeshFaceBase(mesh, f_base); //在每个面上算一组基

    for (int i = 0; i < cons_id.size(); i++) {

        int fid = cons_id[i];
        status[fid] = 1;
        Point cf = cons_vec[i].normalized();
        f_dir[fid] = std::pow(COMPLEX(cf.dot(f_base[fid * 2]), cf.dot(f_base[fid * 2 + 1])),4);
    }

    std::vector<int> id2sln(fnum, -1);
    std::vector<int> sln2id(0);
    int count = 0;
    for (int i = 0; i < fnum; i++) {
        if (status[i] == 0) {
            sln2id.push_back(i);
            id2sln[i] = count;
            count++;
        }
    }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<COMPLEX>> slu;
    Eigen::SparseMatrix<COMPLEX> A;
    Eigen::VectorXcd b_pre(mesh->GetNumberOfEdges());
    Eigen::VectorXcd b;
    b_pre.setZero();
    vector<Eigen::Triplet<COMPLEX>> tris;

    igIndex vcnt=0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    igIndex ecnt=0;
    igIndex ehs[IGAME_CELL_MAX_SIZE];//面的边id

    igIndex fcnt=0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];

    auto faces=mesh->GetFaces();
    auto edges=mesh->GetEdges();

    count=0;
    for (int i = 0; i < fnum; i++) {

        int id_f = i;
        COMPLEX sum = 0;
        ecnt= mesh->GetFaceEdgeIds(i, ehs);
        for (int j = 0; j < ecnt; j++) {

            if (!mesh->IsBoundaryEdge(ehs[j])) {

                vcnt = edges->GetCellIds(ehs[j], vhs);
                assert(vcnt == 2);
                Point p1 = mesh->GetPoint(vhs[0]);
                Point p2 = mesh->GetPoint(vhs[1]);

                fcnt = mesh->GetEdgeToNeighborFaces(ehs[j], fhs);
                assert(fcnt == 2);
                int id_g = (fhs[0] == i)? fhs[1] : fhs[0];

                if (id_f < id_g) { 
                    
                    Point e = (p2 - p1).normalized();

                    COMPLEX e_f = COMPLEX(e.dot(f_base[id_f * 2]),e.dot(f_base[id_f * 2 + 1]));
                    COMPLEX e_g = COMPLEX(e.dot(f_base[id_g * 2]),e.dot(f_base[id_g * 2 + 1]));

                    COMPLEX e_f_c_4 = pow(conj(e_f), 4);
                    COMPLEX e_g_c_4 = pow(conj(e_g), 4);

                    if (status[id_f] == 1 && status[id_g] == 1) continue;
                    if (status[id_f] == 0) {
                        tris.push_back(Eigen::Triplet<COMPLEX>(
                                count, id2sln[id_f], e_f_c_4));
                    } else {
                        b_pre[count] += -e_f_c_4 * f_dir[id_f];
                    }
                    if (status[id_g] == 0) {
                        tris.push_back(Eigen::Triplet<COMPLEX>(
                                count, id2sln[id_g], -e_g_c_4));
                    } else {
                        b_pre[count] += e_g_c_4 * f_dir[id_g];
                    }

                    count++;
                }


            }
        }

    }
    A.resize(count, sln2id.size());
	b.resize(count);
	b = b_pre.head(count);
	A.setFromTriplets(tris.begin(), tris.end());
	Eigen::SparseMatrix<COMPLEX> AT = A.adjoint();
	slu.compute(AT * A);
	Eigen::VectorXcd x = slu.solve(AT * b);

    crossfield.resize(4 * fnum);
	for (int i = 0; i < fnum; i++)
	{
        if (status[i] == 0) { f_dir[i] = x(id2sln[i]); }
        double length = 1;
        double arg = std::arg(f_dir[i]) / 4;
        for (int j = 0; j < 4; j++) {
            crossfield[i * 4 + j] =f_base[i * 2] * length * cos(arg + j * PI / 2) +f_base[i * 2 + 1] * length * sin(arg + j * PI / 2);

        }
	}
    //std::cout << "crossfield.size():"<< crossfield.size()<<std::endl;

    //writeCrossField("D:/test.txt", crossfield);
    for (int i = 0; i < fnum; i++)
    {
        IGsize id[4] = {0};

        auto f = mesh->GetFace(i);
        vcnt = faces->GetCellIds(i, vhs);
        assert(vcnt == 3);
        Point p1 = mesh->GetPoint(vhs[0]);
        Point p2 = mesh->GetPoint(vhs[1]);
        Point p3 = mesh->GetPoint(vhs[2]);
        Point p4 = (p1 + p2 + p3) / 3;
        float Lmean =
                ((p4 - p1).norm() + (p4 - p2).norm() + (p4 - p3).norm()) / 6;
        //std::cout << p4 << ":" << crossfield[4 * i] << ":"
        //          << crossfield[4 * i + 2]<<endl;
        //std::cout << "alean"<< Lmean << std::endl;
        
        for (int j = 0; j < 4; j++) {

            float x = (crossfield[(4 * i + j) % 4] -
                       crossfield[(4 * i + j + 2) % 4])
                              .norm();
            /*std::cout << "x:"<< x << std::endl;*/
            Point pt = (crossfield[4 * i + j])*Lmean/x+p4;
            id[j] = mesh->AddPoint(pt);
            
        }
        
        IGsize a = mesh->AddEdge(id[0], id[2]);
        IGsize b = mesh->AddEdge(id[1], id[3]);

    }
    //mesh->SetPoints(points);
    

    std::cout << "CrossField end" << std::endl;
    mesh->GarbageCollection();

    SetOutput(mesh); // 设置输出数据
    return true; 


}



IGAME_NAMESPACE_END