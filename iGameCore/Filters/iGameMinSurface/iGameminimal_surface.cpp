#include "iGameminimal_surface.h"
#include <numeric>
#include <Eigen/Dense>
#include "iGameFilter.h"
#include "iGameVolumeMesh.h"
#include "iGamePriorityQueue.h"
#include <Eigen/Sparse>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/src/SparseCore/SparseVector.h>

IGAME_NAMESPACE_BEGIN
#include <vector>


struct MyStruct {
	Vector3f position{ 0.f ,0.f,0.f};
	Vector3f normal{ 0.f,0.f ,0.f };
	int index{ -1 };//��������
	int Bound_index{ -1 };//�߽������
	Vector2f UV{ 0.f,0.f };//������UV����
};
std::vector<MyStruct> v(1000);

// �Զ���ȽϺ���
bool areVectorsApprox(const iGame::Vector<float, 3>& v1, const iGame::Vector<float, 3>& v2, float epsilon = 1e-9) {
	return (std::abs(v1[0] - v2[0]) < epsilon) &&
		(std::abs(v1[1] - v2[1]) < epsilon) &&
		(std::abs(v1[2] - v2[2]) < epsilon);
}

// ���Һ���
int findIndexByPosition(const iGame::Vector<float, 3>& pos, const std::vector<MyStruct>& vertices) {
	for (size_t i = 0; i < vertices.size(); ++i) {
		if (areVectorsApprox(vertices[i].position, pos)) {
            return static_cast<int>(vertices[i].index);
		}
	}
	return -1; // ���û���ҵ������� -1
}
int findBoundIndexByPosition(const iGame::Vector<float, 3>& pos, const std::vector<MyStruct>& vertices) {
	for (size_t i = 0; i < vertices.size(); ++i) {
		if (areVectorsApprox(vertices[i].position, pos)) {
            return static_cast<int>(vertices[i].Bound_index);
		}
	}
	return -1; // ���û���ҵ������� -1
}


bool minmalSurface::Execute()
{
	std::cout << "  start \n";
	UpdateProgress(0);
	mesh = DynamicCast<SurfaceMesh>(GetInput(0));

	if (mesh == nullptr) {
		std::cout << " mesh is null  \n";
		return false;
	}

	std::cout << "mesh gets \n";
	mesh->RequestEditStatus();


	//GlobalLaplacianSmoothing(mesh);


	Parameterization(mesh);


	return true;
}
void minmalSurface::GlobalLaplacianSmoothing(SurfaceMesh::Pointer mesh)
{
	// std::vector<MyStruct> v(1000);
	igIndex vcnt = 0;
	igIndex vhs[IGAME_CELL_MAX_SIZE];
	igIndex fcnt = 0;
	igIndex fhs[IGAME_CELL_MAX_SIZE];

	auto vec = mesh->GetPoints(); //��ȡ���ж���
	int vNum = mesh->GetNumberOfPoints(); //��������
	int vBoundary = 0; //�߽����
	int vInner = 0;  //�Ǳ߽����

	std::cout << vNum;

	for (int i = 0; i < vNum; i++) {
		//������������
		v[i].index = i;
		v[i].position = vec->GetPoint(i);
		if (mesh->IsBoundaryPoint(i)) {
			//�������λ�ڱ߽�㣬���ӱ߽����
			vBoundary++;
		}
	}
	vInner = vNum - vBoundary;


	std::cout << vInner << "�ڲ�";

	//���ȫ�ּ�С�������
	Eigen::SparseMatrix <float> A(3 * (vNum + vBoundary), 3 * vNum);
	Eigen::SparseVector<float> b(3 * (vNum + vBoundary));



	//��������
	auto Idx = [&](int idx, int k) {
		return k * vNum + idx;
		};
	auto B_Idx = [&](int idx, int k) {
		return 3 * vNum + idx * 3 + k;
		};

	int Boundary_index = 0;
	for (int i = 0; i < vNum; i++) {
		//�������ж��㣬����ϡ�����
		const auto P = mesh->GetPoint(i);//��ǰ��������

		int ind = findIndexByPosition(P, v);
		//A_ij & x_i
		if (mesh->IsBoundaryPoint(i)) {
			//����Ǳ߽綥��
			//�洢1 or 0��ʾ�Ƿ�Ϊ�߽綥��
			for (auto k : { 0,1,2 }) {
				A.coeffRef(B_Idx(Boundary_index, k), Idx(ind, k)) = 1;
				b.coeffRef(B_Idx(Boundary_index, k)) = P[k];
			}
			Boundary_index++;
		}
		//�Ե�ǰ����1-���򶥵�
		vcnt = mesh->GetPointToOneRingPoints(i, vhs);
		// const auto& adj_v = v->AdjVertices();
		for (auto k : { 0, 1, 2 }) {
			//�Զ�������
			//�洢1-���򶥵������
			A.coeffRef(Idx(ind, k), Idx(ind, k)) = vcnt;
		}
		for (int t = 0; t < vcnt; t++) {
			//�Ե�ǰ����� 1-���򶥵�����
			//�ǵ�ǰ����� 1-���򶥵���Ϊ-1
			Vector3f x = mesh->GetPoint(vhs[t]);
			int in = findIndexByPosition(x, v);
			for (auto k : { 0, 1, 2 }) {
				A.coeffRef(Idx(ind, k), Idx(in, k)) = -1;
			}
		}
	}
	//���󹹽��ɹ�
	// spdlog::info("Build Matrix successful");
	std::cout << "Build Matrix successful";
	//������
	Eigen::ConjugateGradient<Eigen::SparseMatrix<float> > solver;
	solver.setTolerance(1e-4);
	solver.compute(A.transpose() * A);
	Eigen::SparseVector<float> x(3 * vNum);
	x = solver.solve(A.transpose() * b);

	//������
	//spdlog::info("Solve Matrix successful");
	std::cout << "Solve Matrix successful";
	std::cout << x << std::endl;

	//�������ж�������
	for (int ii = 0; ii < vNum; ii++) {
		Vector3f tmp;
		for (auto k : { 0,1,2 }) {
			tmp[k] = x.coeffRef(Idx(v[ii].index, k));

		}
		mesh->SetPoint(ii, tmp);//���¶������꣬
	}
	std::cout << "end" << std::endl;
	mesh->GarbageCollection();

	SetOutput(mesh); // �����������
}
void minmalSurface::Parameterization(SurfaceMesh::Pointer mesh)
{
	
	
	igIndex vcnt = 0;
	igIndex vhs[IGAME_CELL_MAX_SIZE];
	igIndex fcnt = 0;
	igIndex fhs[IGAME_CELL_MAX_SIZE];

	auto vec = mesh->GetPoints(); //��ȡ���ж���
    // std::cout << vec << std::endl;
	int vNum = mesh->GetNumberOfPoints(); //��������
	int vBoundary = 0; //�߽����
	int vInner = 0;  //�Ǳ߽����

	std::cout << vNum;

	for (int i = 0; i < vNum; i++) {
		//������������
		v[i].index = i;   
		v[i].position = vec->GetPoint(i);
		if (mesh->IsBoundaryPoint(i)) {
			if (v[i].Bound_index < 0) {
				//�����߽������
				v[i].Bound_index = vBoundary;
			}
			//�������λ�ڱ߽�㣬���ӱ߽����
			vBoundary++;
		}
	}

	vInner = vNum - vBoundary;


	std::cout << vInner << "�ڲ�";

	//������������
	Eigen::SparseMatrix <float> A(2 * (vNum + vBoundary), 2 * vNum);
	Eigen::SparseVector<float> b(2 * (vNum + vBoundary));


	//��������
	auto Idx = [&](int idx, int k) {
		return k * vNum + idx;
		};
	auto B_Idx = [&](int idx, int k) {
		return 2 * vNum + idx * 2 + k;
		};

	int Boundary_index = 0;int bnum = 0;
	for (int i = 0; i < vNum; i++) {
		//�������ж��㣬����ϡ�����
		const auto P = mesh->GetPoints()->GetPoint(i);//��ǰ��������

		int ind = findIndexByPosition(P, v);

		int Bound_index = findBoundIndexByPosition(P,v);
        
		if (mesh->IsBoundaryPoint(i) ){
			//����Ǳ߽綥��
			//���߽綥��ӳ�䵽uv����߽���[0,1]^2
            bnum++;
				float t = 4 * float(Bound_index) / float(vBoundary); //����t��Ϊ�ı�
				float U = { 0.0f }, V = { 0.0f };//��ʼ��uv����
				if (t <= 1.0f) V = t;
				else if (t <= 2.0f) U = t - 1, V = 1.0f;
				else if (t <= 3.0f) U = 1.0f, V = 3.0f - t;
				else U = 4.0f - t;

				//���UV����
				// spdlog::info(std::to_string(U) + " " + std::to_string(V));
				std::cout<< std::to_string(U)<<" "<< std::to_string(V)<<std::endl;
				//uv��������ֵ
				b.coeffRef(B_Idx(Bound_index, 0)) = U;
				b.coeffRef(B_Idx(Bound_index, 1)) = V;
				v[i].UV = { U,V };
			
			for (auto k : { 0, 1 })
				//����ֵ
				A.coeffRef(B_Idx(Bound_index, k), Idx(ind, k)) = 1;
		}
        
		//�Ե�ǰ����1-���򶥵�
		vcnt = mesh->GetPointToOneRingPoints(i, vhs);
        std::cout << "vcnt:" << vcnt << std::endl;
		// const auto& adj_v = v->AdjVertices();
		for (auto k : { 0, 1}) {
			//�Զ�������
			//�洢1-���򶥵������
			A.coeffRef(Idx(ind, k), Idx(ind, k)) = vcnt;
		}

		for (int t = 0; t < vcnt; t++) {
			//�Ե�ǰ����� 1-���򶥵�����
			//�ǵ�ǰ����� 1-���򶥵���Ϊ-1
			Vector3f x = vec->GetPoint(vhs[t]);
			int in = findIndexByPosition(x, v);
			for (auto k : { 0, 1}) {
				A.coeffRef(Idx(ind, k), Idx(in, k)) = -1;
			}
		}
	}
    std::cout << "bnum:" << bnum << std::endl;
	//���󹹽��ɹ�
	// spdlog::info("Build Matrix successful");
	std::cout << "Build Matrix successful";
	//std::cout << A << std::endl;

	//������
	Eigen::ConjugateGradient<Eigen::SparseMatrix<float> > solver;
	solver.setTolerance(1e-2);
	solver.compute(A.transpose() * A);
	Eigen::SparseVector<float> x(2 * vNum);
	x = solver.solve(A.transpose() * b);

	//������
	//spdlog::info("Solve Matrix successful");
	std::cout << "Solve Matrix successful";
	//std::cout << x << std::endl;
	//�������ж�������
	for (int ii = 0; ii < vNum; ii++) {
		Vector3f tmp;
        if (mesh->IsBoundaryPoint(ii)) {
            tmp = {v[ii].UV[0], v[ii].UV[1], 0.f}; // ʹ���Ѽ���� UV
        } else {
            tmp = {x.coeffRef(Idx(v[ii].index, 0)),
                   x.coeffRef(Idx(v[ii].index, 1)), 0.f};
        }
        mesh->SetPoint(ii, tmp); // ���¶�������
        std::cout << "Vertex " << ii << ": (" << tmp[0] << ", " << tmp[1]
                  << ", " << tmp[2] << ")" << std::endl;

	}
	
	// spdlog::info("Data transform successful");

	std::cout << "end" << std::endl;
	mesh->GarbageCollection();

	SetOutput(mesh); // �����������
}

IGAME_NAMESPACE_END