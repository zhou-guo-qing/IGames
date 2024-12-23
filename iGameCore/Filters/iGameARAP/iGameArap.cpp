#include"iGameArap.h"
#include"iGameVector.h"
#include<Eigen/Dense>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<vector>
#include<unordered_map>
#include<Eigen/SVD>
#include<math.h>
#include<algorithm>
#include<omp.h>
IGAME_NAMESPACE_BEGIN
#define PI 3.1415926

bool compareVectorsT(const Vector3f& v1, const Vector3f& v2) {
     for (int i = 0; i < 3; ++i) {
         if (v1[i] != v2[i]) { return false; }
     }
     return true;
 }

struct Key {
     int id1;
     int id2;
 };

 // 重载哈希函数
 struct KeyHasher {
     std::size_t operator()(const Key& k) const {
         return std::hash<int>()(k.id1) ^ std::hash<int>()(k.id2);
     }
 };

 // 重载相等判断函数
 struct KeyEqual {
     bool operator()(const Key& lhs, const Key& rhs) const {
         return lhs.id1 == rhs.id1 && lhs.id2 == rhs.id2;
     }
 };


 void Tutte(SurfaceMesh* mesh,Eigen::MatrixX2d& uv) 
{

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
    while (!compareVectorsT(startp, endp)) {
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

    //inner
    uv.col(0) = xu;
    uv.col(1) = xv;

}
 int num = 0;
 bool Arap::Execute()
{

	std::cout << "ARAP Parameterization start \n";
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
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    //计算局部坐标
    Eigen::MatrixXd localcoord;
    localcoord.resize(faceNum, 6);

    for (int i = 0; i < faceNum; i++) {

        auto f = mesh->GetFace(i);
        auto n = f->GetNormal();
        //std::cout << "mainfaxiang"<< n << std::endl;
        vcnt = faces->GetCellIds(i, vhs); //面的三个点id
        assert(vcnt == 3);
        Point p1 = mesh->GetPoint(vhs[0]);
        Point p2 = mesh->GetPoint(vhs[1]);
        Point p3 = mesh->GetPoint(vhs[2]);

        auto p12 = p2 - p1;
        auto p13 = p3 - p1;

        auto x_ = p12.normalized();
        auto y_ = n.cross(x_);

        localcoord.row(i) << 0, 0, p12.norm(), 0, p13.dot(x_), p13.dot(y_);

    }

    //计算cot_weigth的拉普拉斯矩阵
    //std::vector<double> cots(mesh->GetNumberOfEdges()*3);
    std::unordered_map<Key, double, KeyHasher, KeyEqual> cottan;
    std::vector<Eigen::Triplet<double>> trivec;

    for (int i = 0; i < mesh->GetNumberOfFaces(); i++) {

        auto f = mesh->GetFace(i);

        vcnt = faces->GetCellIds(i, vhs); //面的三个点id
        assert(vcnt == 3);
        Point p0 = mesh->GetPoint(vhs[0]);
        Point p1 = mesh->GetPoint(vhs[1]);
        Point p2 = mesh->GetPoint(vhs[2]);


        //1,2边
        auto e0 = p1 - p0;
        auto e1 = p2 - p0;
        double cotangle = e0.dot(e1) / e1.cross(e0).norm();

        
		//cots[mesh->GetEdgeIdFormPointIds(vhs[1], vhs[2])] = cotangle;
        if (cottan.find({vhs[1], vhs[2]}) == cottan.end()) {
            cottan[{vhs[1], vhs[2]}] = cotangle;
            int vid0 = vhs[1];
            int vid1 = vhs[2];
            trivec.emplace_back(vid0, vid0, cotangle);
            trivec.emplace_back(vid1, vid1, cotangle);
            trivec.emplace_back(vid0, vid1, -cotangle);
            trivec.emplace_back(vid1, vid0, -cotangle);
        }
        else
        {
            std::cout << "键已存在，值为：" <<  num++ << "\n";
        }


        //0，1边
        e0 = p1 - p2; 
        e1 = p0 - p2;
        cotangle = e0.dot(e1) / e1.cross(e0).norm();


        //cots[mesh->GetEdgeIdFormPointIds(vhs[1], vhs[2])] = cotangle;
        if (cottan.find({ vhs[0], vhs[1] }) == cottan.end()) {
            cottan[{vhs[0], vhs[1]}] = cotangle;
            int vid0 = vhs[0];
            int vid1 = vhs[1];
            trivec.emplace_back(vid0, vid0, cotangle);
            trivec.emplace_back(vid1, vid1, cotangle);
            trivec.emplace_back(vid0, vid1, -cotangle);
            trivec.emplace_back(vid1, vid0, -cotangle);
        } else {
            std::cout << "键已存在，值为：" << num++ << "\n";
        }


        //2,0边
        e0 = p0 - p1;
        e1 = p2 - p1;
        cotangle = e0.dot(e1) / e1.cross(e0).norm();


        //cots[mesh->GetEdgeIdFormPointIds(vhs[1], vhs[2])] = cotangle;
        if (cottan.find({ vhs[2], vhs[0] }) == cottan.end()) {
            cottan[{vhs[2], vhs[0]}] = cotangle;
            int vid0 = vhs[2];
            int vid1 = vhs[0];
            trivec.emplace_back(vid0, vid0, cotangle);
            trivec.emplace_back(vid1, vid1, cotangle);
            trivec.emplace_back(vid0, vid1, -cotangle);
            trivec.emplace_back(vid1, vid0, -cotangle);
        } else {
            std::cout << "键已存在，值为：" << num++ << "\n";
        }

    }
    std::cout << trivec.size() << std::endl;
    std::cout << "cotan:" << cottan.size() << std::endl;
    Eigen::SparseMatrix<double> smat;
    smat.resize(vNum, vNum);
    smat.setFromTriplets(trivec.begin(), trivec.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(smat);


    Eigen::MatrixX2d uv;
    uv.resize(vNum, 2);

    //Tutte 
    Tutte(mesh,uv);

    std::vector<Eigen::Matrix2d> lts;
    lts.resize(faceNum);

    //local-global iteration
    for (int iter = 0; iter < 100; iter++) {

        //local calc Lt
#pragma omp parallel for
        for (int i = 0; i < faceNum; i++) {

            auto f = mesh->GetFace(i);

            vcnt = faces->GetCellIds(i, vhs); //面的三个点id
            assert(vcnt == 3);
            Point p1 = mesh->GetPoint(vhs[0]);
            Point p2 = mesh->GetPoint(vhs[1]);
            Point p3 = mesh->GetPoint(vhs[2]);

            Eigen::Matrix2d P, S, J;

            P << uv(vhs[1], 0) - uv(vhs[0], 0), uv(vhs[2], 0) - uv(vhs[0], 0),
                    uv(vhs[1], 1) - uv(vhs[0], 1), uv(vhs[2], 1) - uv(vhs[0], 1);
            S << localcoord(i, 2) - localcoord(i, 0),
                    localcoord(i, 4) - localcoord(i, 0),
                    localcoord(i, 3) - localcoord(i, 1),
                    localcoord(i, 5) - localcoord(i, 1);
            J = P * S.inverse();

            Eigen::JacobiSVD<Eigen::Matrix2d> svd(
                    J, Eigen::ComputeFullU | Eigen::ComputeFullV);

            Eigen::Matrix2d U = svd.matrixU();
            Eigen::Matrix2d V = svd.matrixV();

            Eigen::Matrix2d R = U * V.transpose();

            if (R.determinant() < 0) {
                U(0, 1) = -U(0, 1);
                U(1, 1) = -U(1, 1);
                R = U * V.transpose();
            }
            lts[i] = R;
        }

        //global calc b

		Eigen::VectorXd bu, bv;
        bu.setZero(vNum);
        bv.setZero(vNum);

        for (int i = 0; i < faceNum; i++) {

            auto f = mesh->GetFace(i);
            auto n = f->GetNormal();

            vcnt = faces->GetCellIds(i, vhs); //面的三个点id
            assert(vcnt == 3);
            Point p0 = mesh->GetPoint(vhs[0]);
            Point p1 = mesh->GetPoint(vhs[1]);
            Point p2 = mesh->GetPoint(vhs[2]);


            Eigen::Vector2d e0, e1, e2;
            e0 << localcoord(i, 2), localcoord(i, 3);
            e1 << localcoord(i, 4) - localcoord(i, 2),
                    localcoord(i, 5) - localcoord(i, 3);
            e2 << -localcoord(i, 4), -localcoord(i, 5);

            Eigen::Vector2d b0 = cottan[{vhs[0],vhs[1]}] * lts[i] * e0;
            bu[vhs[0]] -= b0[0];
            bv[vhs[0]] -= b0[1];

            bu[vhs[1]] += b0[0];
            bv[vhs[1]] += b0[1];
            Eigen::Vector2d b1 = cottan[{vhs[1], vhs[2]}] * lts[i] * e1;
            bu[vhs[1]] -= b1[0];
            bv[vhs[1]] -= b1[1];

            bu[vhs[2]] += b1[0];
            bv[vhs[2]] += b1[1];

            Eigen::Vector2d b2 = cottan[{vhs[2], vhs[0]}] * lts[i] * e2;
            bu[vhs[2]] -= b2[0];
            bv[vhs[2]] -= b2[1];

            bu[vhs[0]] += b2[0];
            bv[vhs[0]] += b2[1];
        }
        //global solve
        uv.col(0) = solver.solve(bu);
        uv.col(1) = solver.solve(bv);

    }
    for (int i = 0; i < vNum; i++) {
        Vector3f tmp = {static_cast<float>(uv(i, 0)),
                        static_cast<float>(uv(i, 1)), 0};
        mesh->SetPoint(i, tmp);
    }
    std::cout << "ARAP Parameterization end" << std::endl;
    mesh->GarbageCollection();

    SetOutput(mesh); // 设置输出数据
    return true; 

}

IGAME_NAMESPACE_END