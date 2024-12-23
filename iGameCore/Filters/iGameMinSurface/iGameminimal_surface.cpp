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
	int index{ -1 };//顶点索引
	int Bound_index{ -1 };//边界点索引
	Vector2f UV{ 0.f,0.f };//参数化UV坐标
};
std::vector<MyStruct> v(1000);

// 自定义比较函数
bool areVectorsApprox(const iGame::Vector<float, 3>& v1, const iGame::Vector<float, 3>& v2, float epsilon = 1e-9) {
	return (std::abs(v1[0] - v2[0]) < epsilon) &&
		(std::abs(v1[1] - v2[1]) < epsilon) &&
		(std::abs(v1[2] - v2[2]) < epsilon);
}

// 查找函数
int findIndexByPosition(const iGame::Vector<float, 3>& pos, const std::vector<MyStruct>& vertices) {
	for (size_t i = 0; i < vertices.size(); ++i) {
		if (areVectorsApprox(vertices[i].position, pos)) {
            return static_cast<int>(vertices[i].index);
		}
	}
	return -1; // 如果没有找到，返回 -1
}
int findBoundIndexByPosition(const iGame::Vector<float, 3>& pos, const std::vector<MyStruct>& vertices) {
	for (size_t i = 0; i < vertices.size(); ++i) {
		if (areVectorsApprox(vertices[i].position, pos)) {
            return static_cast<int>(vertices[i].Bound_index);
		}
	}
	return -1; // 如果没有找到，返回 -1
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

	auto vec = mesh->GetPoints(); //获取所有顶点
	int vNum = mesh->GetNumberOfPoints(); //顶点数量
	int vBoundary = 0; //边界点数
	int vInner = 0;  //非边界点数

	std::cout << vNum;

	for (int i = 0; i < vNum; i++) {
		//建立顶点索引
		v[i].index = i;
		v[i].position = vec->GetPoint(i);
		if (mesh->IsBoundaryPoint(i)) {
			//如果顶点位于边界点，增加边界点数
			vBoundary++;
		}
	}
	vInner = vNum - vBoundary;


	std::cout << vInner << "内部";

	//求解全局极小曲面矩阵
	Eigen::SparseMatrix <float> A(3 * (vNum + vBoundary), 3 * vNum);
	Eigen::SparseVector<float> b(3 * (vNum + vBoundary));



	//索引函数
	auto Idx = [&](int idx, int k) {
		return k * vNum + idx;
		};
	auto B_Idx = [&](int idx, int k) {
		return 3 * vNum + idx * 3 + k;
		};

	int Boundary_index = 0;
	for (int i = 0; i < vNum; i++) {
		//遍历所有顶点，构建稀疏矩阵
		const auto P = mesh->GetPoint(i);//当前顶点坐标

		int ind = findIndexByPosition(P, v);
		//A_ij & x_i
		if (mesh->IsBoundaryPoint(i)) {
			//如果是边界顶点
			//存储1 or 0表示是否为边界顶点
			for (auto k : { 0,1,2 }) {
				A.coeffRef(B_Idx(Boundary_index, k), Idx(ind, k)) = 1;
				b.coeffRef(B_Idx(Boundary_index, k)) = P[k];
			}
			Boundary_index++;
		}
		//对当前结点的1-邻域顶点
		vcnt = mesh->GetPointToOneRingPoints(i, vhs);
		// const auto& adj_v = v->AdjVertices();
		for (auto k : { 0, 1, 2 }) {
			//对顶点数据
			//存储1-邻域顶点的数量
			A.coeffRef(Idx(ind, k), Idx(ind, k)) = vcnt;
		}
		for (int t = 0; t < vcnt; t++) {
			//对当前顶点的 1-邻域顶点数据
			//是当前顶点的 1-邻域顶点置为-1
			Vector3f x = mesh->GetPoint(vhs[t]);
			int in = findIndexByPosition(x, v);
			for (auto k : { 0, 1, 2 }) {
				A.coeffRef(Idx(ind, k), Idx(in, k)) = -1;
			}
		}
	}
	//矩阵构建成功
	// spdlog::info("Build Matrix successful");
	std::cout << "Build Matrix successful";
	//求解矩阵
	Eigen::ConjugateGradient<Eigen::SparseMatrix<float> > solver;
	solver.setTolerance(1e-4);
	solver.compute(A.transpose() * A);
	Eigen::SparseVector<float> x(3 * vNum);
	x = solver.solve(A.transpose() * b);

	//求解完成
	//spdlog::info("Solve Matrix successful");
	std::cout << "Solve Matrix successful";
	std::cout << x << std::endl;

	//更新所有顶点坐标
	for (int ii = 0; ii < vNum; ii++) {
		Vector3f tmp;
		for (auto k : { 0,1,2 }) {
			tmp[k] = x.coeffRef(Idx(v[ii].index, k));

		}
		mesh->SetPoint(ii, tmp);//更新顶点坐标，
	}
	std::cout << "end" << std::endl;
	mesh->GarbageCollection();

	SetOutput(mesh); // 设置输出数据
}
void minmalSurface::Parameterization(SurfaceMesh::Pointer mesh)
{
	
	
	igIndex vcnt = 0;
	igIndex vhs[IGAME_CELL_MAX_SIZE];
	igIndex fcnt = 0;
	igIndex fhs[IGAME_CELL_MAX_SIZE];

	auto vec = mesh->GetPoints(); //获取所有顶点
    // std::cout << vec << std::endl;
	int vNum = mesh->GetNumberOfPoints(); //顶点数量
	int vBoundary = 0; //边界点数
	int vInner = 0;  //非边界点数

	std::cout << vNum;

	for (int i = 0; i < vNum; i++) {
		//建立顶点索引
		v[i].index = i;   
		v[i].position = vec->GetPoint(i);
		if (mesh->IsBoundaryPoint(i)) {
			if (v[i].Bound_index < 0) {
				//建立边界点索引
				v[i].Bound_index = vBoundary;
			}
			//如果顶点位于边界点，增加边界点数
			vBoundary++;
		}
	}

	vInner = vNum - vBoundary;


	std::cout << vInner << "内部";

	//求解参数化矩阵
	Eigen::SparseMatrix <float> A(2 * (vNum + vBoundary), 2 * vNum);
	Eigen::SparseVector<float> b(2 * (vNum + vBoundary));


	//索引函数
	auto Idx = [&](int idx, int k) {
		return k * vNum + idx;
		};
	auto B_Idx = [&](int idx, int k) {
		return 2 * vNum + idx * 2 + k;
		};

	int Boundary_index = 0;int bnum = 0;
	for (int i = 0; i < vNum; i++) {
		//遍历所有顶点，构建稀疏矩阵
		const auto P = mesh->GetPoints()->GetPoint(i);//当前顶点坐标

		int ind = findIndexByPosition(P, v);

		int Bound_index = findBoundIndexByPosition(P,v);
        
		if (mesh->IsBoundaryPoint(i) ){
			//如果是边界顶点
			//将边界顶点映射到uv坐标边界上[0,1]^2
            bnum++;
				float t = 4 * float(Bound_index) / float(vBoundary); //参数t分为四边
				float U = { 0.0f }, V = { 0.0f };//初始化uv坐标
				if (t <= 1.0f) V = t;
				else if (t <= 2.0f) U = t - 1, V = 1.0f;
				else if (t <= 3.0f) U = 1.0f, V = 3.0f - t;
				else U = 4.0f - t;

				//输出UV坐标
				// spdlog::info(std::to_string(U) + " " + std::to_string(V));
				std::cout<< std::to_string(U)<<" "<< std::to_string(V)<<std::endl;
				//uv坐标结果赋值
				b.coeffRef(B_Idx(Bound_index, 0)) = U;
				b.coeffRef(B_Idx(Bound_index, 1)) = V;
				v[i].UV = { U,V };
			
			for (auto k : { 0, 1 })
				//矩阵赋值
				A.coeffRef(B_Idx(Bound_index, k), Idx(ind, k)) = 1;
		}
        
		//对当前结点的1-邻域顶点
		vcnt = mesh->GetPointToOneRingPoints(i, vhs);
        std::cout << "vcnt:" << vcnt << std::endl;
		// const auto& adj_v = v->AdjVertices();
		for (auto k : { 0, 1}) {
			//对顶点数据
			//存储1-邻域顶点的数量
			A.coeffRef(Idx(ind, k), Idx(ind, k)) = vcnt;
		}

		for (int t = 0; t < vcnt; t++) {
			//对当前顶点的 1-邻域顶点数据
			//是当前顶点的 1-邻域顶点置为-1
			Vector3f x = vec->GetPoint(vhs[t]);
			int in = findIndexByPosition(x, v);
			for (auto k : { 0, 1}) {
				A.coeffRef(Idx(ind, k), Idx(in, k)) = -1;
			}
		}
	}
    std::cout << "bnum:" << bnum << std::endl;
	//矩阵构建成功
	// spdlog::info("Build Matrix successful");
	std::cout << "Build Matrix successful";
	//std::cout << A << std::endl;

	//求解矩阵
	Eigen::ConjugateGradient<Eigen::SparseMatrix<float> > solver;
	solver.setTolerance(1e-2);
	solver.compute(A.transpose() * A);
	Eigen::SparseVector<float> x(2 * vNum);
	x = solver.solve(A.transpose() * b);

	//求解完成
	//spdlog::info("Solve Matrix successful");
	std::cout << "Solve Matrix successful";
	//std::cout << x << std::endl;
	//更新所有顶点坐标
	for (int ii = 0; ii < vNum; ii++) {
		Vector3f tmp;
        if (mesh->IsBoundaryPoint(ii)) {
            tmp = {v[ii].UV[0], v[ii].UV[1], 0.f}; // 使用已计算的 UV
        } else {
            tmp = {x.coeffRef(Idx(v[ii].index, 0)),
                   x.coeffRef(Idx(v[ii].index, 1)), 0.f};
        }
        mesh->SetPoint(ii, tmp); // 更新顶点坐标
        std::cout << "Vertex " << ii << ": (" << tmp[0] << ", " << tmp[1]
                  << ", " << tmp[2] << ")" << std::endl;

	}
	
	// spdlog::info("Data transform successful");

	std::cout << "end" << std::endl;
	mesh->GarbageCollection();

	SetOutput(mesh); // 设置输出数据
}

IGAME_NAMESPACE_END