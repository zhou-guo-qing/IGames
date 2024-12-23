#ifndef iGameARAPTest_h
#define iGameARAPTest_h

#include "iGameFilter.h"
#include "iGameSurfaceMesh.h"
#include "iGameModel.h"
#include <random>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iGameVolume.h>
#include <iGameVolume.h>
#include <Quadratic/Base/iGameQuadraticVolume.h>

IGAME_NAMESPACE_BEGIN
bool compareVectorsQ(const Vector3f& v1, const Vector3f& v2) {
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
class ARAPTest : public Filter {
public:
	I_OBJECT(ARAPTest);
	static Pointer New() { return new ARAPTest; }

	bool Initialize() {
		mesh = DynamicCast<SurfaceMesh>(GetInput(0));
		model = m_Model;
		if (mesh == nullptr) { return false; }

        std::cout << "mesh gets \n";
        mesh->RequestEditStatus();
        

		// 这里请求进行选点
		Points* ps = mesh->GetPoints();
		fixed = Selection::New();
		fixed->SetFilterEvent(&ARAPTest::CallbackEvent, this, std::placeholders::_1);
		model->RequestPointSelection(ps, fixed);

		// 执行算法初始化
		//auto painter = model->GetPainter();
		//painter->SetPen(Color::Red);
		//painter->SetPen(3);
		//painter->DrawLine(mesh->GetPoint(0), mesh->GetPoint(1));

		//FloatArray::Pointer color = FloatArray::New();
		//color->SetName("test");
		//for (int i = 0; i < mesh->GetNumberOfPoints(); i++) {
		//	color->AddValue(mesh->GetPoint(i)[0]);
		//}
		//mesh->GetAttributeSet()->AddAttribute(IG_SCALAR, IG_POINT, color);

		return true;
	}

	bool Begin() {
		// 这里请求拖动点
		Points* ps = mesh->GetPoints();
		moved = Selection::New();
        std::cout << "1" << std::endl;
		moved->SetFilterEvent(&ARAPTest::CallbackEvent, this, std::placeholders::_1);
        std::cout << "2" << std::endl;
        model->RequestDragPoint(ps, moved);
        std::cout << "3" << std::endl;

		return true;
	}

	void calc_cot_weight(
            std::unordered_map<Key, double, KeyHasher, KeyEqual> &cottan) {
        for (int i = 0; i < mesh->GetNumberOfFaces(); i++) {
            igIndex vcnt = 0;
            igIndex vhs[IGAME_CELL_MAX_SIZE];
            auto f = mesh->GetFace(i);
            auto faces = mesh->GetFaces();

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
                //cottan[{vhs[1], vhs[2]}] = cotangle;
                cottan[{vhs[1], vhs[2]}] = 1;
                //int vid0 = vhs[1];
                //int vid1 = vhs[2];
                //trivec.emplace_back(vid0, vid0, cotangle);
                //trivec.emplace_back(vid1, vid1, cotangle);
                //trivec.emplace_back(vid0, vid1, -cotangle);
                //trivec.emplace_back(vid1, vid0, -cotangle);
            } else {
                std::cout << "键已存在，值为：" << num++ << "\n";
            }


            //0，1边
            e0 = p1 - p2;
            e1 = p0 - p2;
            cotangle = e0.dot(e1) / e1.cross(e0).norm();


            //cots[mesh->GetEdgeIdFormPointIds(vhs[1], vhs[2])] = cotangle;
            if (cottan.find({vhs[0], vhs[1]}) == cottan.end()) {
                //cottan[{vhs[0], vhs[1]}] = cotangle;
                cottan[{vhs[0], vhs[1]}] = 1;
                //int vid0 = vhs[0];
                //int vid1 = vhs[1];
                //trivec.emplace_back(vid0, vid0, cotangle);
                //trivec.emplace_back(vid1, vid1, cotangle);
                //trivec.emplace_back(vid0, vid1, -cotangle);
                //trivec.emplace_back(vid1, vid0, -cotangle);
            } else {
                std::cout << "键已存在，值为：" << num++ << "\n";
            }


            //2,0边
            e0 = p0 - p1;
            e1 = p2 - p1;
            cotangle = e0.dot(e1) / e1.cross(e0).norm();


            //cots[mesh->GetEdgeIdFormPointIds(vhs[1], vhs[2])] = cotangle;
            if (cottan.find({vhs[2], vhs[0]}) == cottan.end()) {
                //cottan[{vhs[2], vhs[0]}] = cotangle;
                cottan[{vhs[2], vhs[0]}] = 1;
                //int vid0 = vhs[2];
                //int vid1 = vhs[0];
                //trivec.emplace_back(vid0, vid0, cotangle);
                //trivec.emplace_back(vid1, vid1, cotangle);
                //trivec.emplace_back(vid0, vid1, -cotangle);
                //trivec.emplace_back(vid1, vid0, -cotangle);
            } else {
                std::cout << "键已存在，值为：" << num++ << "\n";
            }
        }
    }
	bool Execute() override {
		// TODO: 执行算法
        std::cout << "start" << std::endl;
		//这是移动点的id和坐标
		//mesh->SetPoint(dragId, dragNew);
        //mesh->SetPoint(1, dragNew);

		std::vector<int> handles_m;
		std::vector<Vector3f> handles_m_pos;
        handles_m.push_back(dragId);
        handles_m_pos.push_back(dragNew);

		std::set<int> handles = handles_f;
        handles.insert(handles_m.begin(), handles_m.end());

		igIndex vNum = mesh->GetNumberOfPoints();
        auto edges = mesh->GetEdges();

		igIndex vcnt = 0;
        igIndex vhs[IGAME_CELL_MAX_SIZE];
        igIndex ecnt = 0;
        igIndex ehs[IGAME_CELL_MAX_SIZE];
        //position backup
        //Drag point id : 173 - 0.0625835 0.137891 0.036679
		//calc cot-weight laplacian matrix
        std::unordered_map<Key, double, KeyHasher, KeyEqual> cots;
        calc_cot_weight(cots);
        pos_mesh_ref.resize(mesh->GetNumberOfPoints());
        for (int i = 0; i < mesh->GetNumberOfPoints(); i++) {
            pos_mesh_ref[i] = mesh->GetPoint(i);
        }
        std::cout << "nei" << mesh->GetPoint(173) << std::endl;
        std::cout << "nei2" << pos_mesh_ref[173] << std::endl;
        std::vector<Eigen::Triplet<double>> trivec;
        for (int i = 0; i < vNum; i++) {
            if (handles.count(i) > 0) {
                trivec.emplace_back(i, i, 1.);
                continue;
            }
            auto v_h = mesh->GetPoint(i);
            double weight_sum = 0.;
            vcnt = mesh->GetPointToNeighborEdges(i, vhs);
            Point p = mesh->GetPoint(i);
            //std::cout << "p" << p<<std::endl;
            for (int j = 0; j < vcnt;j++) { 
                ecnt = edges->GetCellIds(vhs[j], ehs); //边的二个端点id
                Point p1 = mesh->GetPoint(ehs[0]);
                Point p2 = mesh->GetPoint(ehs[1]);
                assert(ecnt == 2);

                int a;
                int b;
                if (compareVectorsQ(p, p1)) { 
                    
                    a = i;
                    b = ehs[1];
                } else {
                
                    a = i;
                    b = ehs[0];
                }
                //auto v_to_h = ehs[1];
                double weight_ = cots[{a, b}] + cots[{b, a}];

                weight_sum += weight_;
                trivec.emplace_back(i, b, -weight_);
            }

            trivec.emplace_back(i, i, weight_sum);
        }
        Eigen::SparseMatrix<double> smat;
        smat.resize(vNum, vNum);
        smat.setFromTriplets(trivec.begin(), trivec.end());

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(smat);
		
        Eigen::MatrixX3d uv;
        uv.resize(vNum, 3);

        std::vector<Eigen::Matrix3d> Lts;
        Lts.resize(vNum);
        Eigen::MatrixX3d b;
        b.resize(vNum, 3);

       //local-global iteration
        for (int iter = 0; iter < 10; iter++) {

            //local calc Lt
            for (int i = 0; i < vNum; i++) {

                auto p = mesh->GetPoint(i);
                Eigen::Matrix3d J = Eigen::Matrix3d::Zero();

                vcnt = mesh->GetPointToNeighborEdges(i, vhs);
                for (int k = 0; k < vcnt; k++) {

                    ecnt = edges->GetCellIds(vhs[k], ehs); //边的二个端点id
                    assert(ecnt == 2);
                    Point p1 = mesh->GetPoint(ehs[0]);
                    Point p2 = mesh->GetPoint(ehs[1]);
                    int a;
                    int b;
                    if (compareVectorsQ(p, p1)) {

                        a = i;
                        b = ehs[1];
                    } else {

                        a = i;
                        b = ehs[0];
                    }
                    //auto v_to_h = ehs[1];
                    auto pp = mesh->GetPoint(b);
                    auto e_ = pos_mesh_ref[i] - pos_mesh_ref[b];
                    auto ep_ = p - pp;
                    double weight_ =
                            cots[{ehs[0], ehs[1]}] + cots[{ehs[1], ehs[0]}];

                    Eigen::Vector3d ep(ep_[0], ep_[1], ep_[2]);
                    Eigen::Vector3d e(e_[0], e_[1], e_[2]);
                    J += weight_ * (e * ep.transpose());

                }

                Eigen::JacobiSVD<Eigen::Matrix3d> svd(J, Eigen::ComputeFullU |
                                                        Eigen::ComputeFullV);
                Eigen::Matrix3d U = svd.matrixU();
                Eigen::Matrix3d V = svd.matrixV();

                Eigen::Matrix3d R = V * U.transpose();

                if (R.determinant() < 0) {
                    U(0, 2) *= -1;
                    U(1, 2) *= -1;
                    U(2, 2) *= -1;
                    R = V * U.transpose();
                }
                Lts[i] = R;
            }

            //global calc b
#pragma omp parallel for
            for (int i = 0; i < vNum; i++) {
                auto v_h = mesh->GetPoint(i);
                Eigen::Vector3d b_tmp(0., 0., 0.);
                auto p = mesh->GetPoint(i);

                vcnt = mesh->GetPointToNeighborEdges(i, vhs);
                for (int j = 0; j < vcnt;j++) {
                    ecnt = edges->GetCellIds(vhs[j], ehs); //边的二个端点id
                    assert(ecnt == 2);
                    Point p1 = mesh->GetPoint(ehs[0]);
                    Point p2 = mesh->GetPoint(ehs[1]);
                    int a;
                    int b;
                    if (compareVectorsQ(p, p1)) {

                        a = i;
                        b = ehs[1];
                    } else {

                        a = i;
                        b = ehs[0];
                    }
                    //auto v_to_h = ehs[1];
                    auto ep_ = pos_mesh_ref[i] - pos_mesh_ref[b];
                    Eigen::Vector3d ep(ep_[0], ep_[1], ep_[2]);

                    Eigen::Matrix3d JR = Lts[i] + Lts[b];
                    double weight_ =
                            (cots[{ehs[0], ehs[1]}] + cots[{ehs[1], ehs[0]}])/2;
                    b_tmp += weight_ * (JR * ep);
                }
                b(i, 0) = b_tmp[0];
                b(i, 1) = b_tmp[1];
                b(i, 2) = b_tmp[2];
            }

            //set handles
            for (int i: handles_f) {
                auto b_tmp = pos_mesh_ref[i];
                b(i, 0) = b_tmp[0];
                b(i, 1) = b_tmp[1];
                b(i, 2) = b_tmp[2];
            }

            for (int i = 0; i < handles_m.size(); i++) {
                auto b_tmp = handles_m_pos[i];
                b(handles_m[i], 0) = b_tmp[0];
                b(handles_m[i], 1) = b_tmp[1];
                b(handles_m[i], 2) = b_tmp[2];
            }

            //global solve
            uv.col(0) = solver.solve(b.col(0));
            uv.col(1) = solver.solve(b.col(1));
            uv.col(2) = solver.solve(b.col(2));
            
            for (int i = 0; i < vNum; i++) {
                    Vector3f tmp = {static_cast<float>(uv(i, 0)),
                                    static_cast<float>(uv(i, 1)),
                                    static_cast<float>(uv(i, 2))};
                    mesh->SetPoint(i, tmp);		
}                
        } 
        

        mesh->Modified();
		model->Update(); // 更新模型
        std::cout << "end" << std::endl;
		return true;
    }

	void CallbackEvent(Selection::Event _event) {
		switch (_event.type)
		{
		case Selection::Event::PickPoint:
			// 选几个固定点, 并保存下来
			std::cout << "Pick point id: " << _event.pickId << std::endl;
            handles_f.emplace(_event.pickId);
			break;
		case Selection::Event::PickFace:
			std::cout << "Pick face id: " << _event.pickId << std::endl;
			break;
		case Selection::Event::DragPoint:
			std::cout << "Drag point id: " << _event.pickId << " " << _event.pos << std::endl;
			dragId = _event.pickId;
			dragNew = _event.pos;
			this->Execute();
			break;
		default:
			break;
		}
	}

protected:
	ARAPTest()
	{
		SetNumberOfInputs(1);
	}
	~ARAPTest() override = default;

	igIndex dragId{ -1 };
	Vector3f dragNew{};
	
	int num = 0;
    std::set<int> handles_f;
    std::vector<Vector3f> pos_mesh_ref;

	Selection::Pointer fixed{};
	Selection::Pointer moved{};
	SurfaceMesh::Pointer mesh{};
	Model::Pointer model{};
};
IGAME_NAMESPACE_END
#endif