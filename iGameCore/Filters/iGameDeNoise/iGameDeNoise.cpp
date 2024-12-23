#include "iGameDeNoise.h"
#include <numeric>
#include <Eigen/Dense>
#include "iGameFilter.h"
#include "iGameVolumeMesh.h"
#include "iGamePriorityQueue.h"

IGAME_NAMESPACE_BEGIN


bool DeNoise::Execute()
{
	std::cout << " DeNoise start \n" ;
	UpdateProgress(0);
	mesh = DynamicCast<SurfaceMesh>(GetInput(0));

	if (mesh == nullptr) { 
		std::cout << " mesh is null  \n";
		return false; 
	}

	std::cout << "mesh gets \n";
	mesh->RequestEditStatus();
	//进行极小曲面法迭代

	igIndex vcnt = 0;
	igIndex vhs[IGAME_CELL_MAX_SIZE];
	igIndex fcnt = 0;
	igIndex fhs[IGAME_CELL_MAX_SIZE];

	for (int t = 0; t < m_Cnt; ++t) {
		//迭代iter_time次

		igIndex PointNum = mesh->GetNumberOfPoints();
		for (int i = 0; i < PointNum; i++) {
			if (mesh->IsBoundaryPoint(i))
				continue;

			const auto P = mesh->GetPoint(i);//当前顶点坐标
			Vector3f  Hn = { 0.f,0.f,0.f };//当前顶点的平均曲率

			//计算每个更新后每个三角形面积的和（曲面总面积）
			float A = 0.f;

			//fcnt = mesh->GetPointToNeighborFaces(i, fhs);

			vcnt = mesh->GetPointToOneRingPoints(i, vhs);

			for (int k = 0; k < vcnt; k++)
			{
				Vector3f p0 = mesh->GetPoint(vhs[k]);
				Vector3f p1 = ((k+1) != vcnt ? mesh->GetPoint(vhs[k + 1]): mesh->GetPoint(vhs[0]));

				auto d1 = p1 - P;
				auto d2 = p0 - P;

				A += 0.5f * d2.cross(d1).norm();
			}
			std::cout << A;

			for (int j = 0; j < vcnt; j++)
			{
				auto  p0 = (j!=0 ? mesh->GetPoint(vhs[j-1]):mesh->GetPoint(vhs[vcnt-1]));
				auto  p1 = mesh->GetPoint(vhs[j]);
				auto  p2 = (j+1<vcnt ? mesh->GetPoint(vhs[j + 1]):mesh->GetPoint(vhs[0]));

				float cos_a = (p1 - p0).dot(P - p0);
				float cos_b = (p1 - p2).dot(P - p2);

				cos_a = std::min(cos_a, 0.9f);
				cos_b = std::min(cos_b, 0.9f);

				float cot_a = cos_a / sqrt(1 - cos_a * cos_a);
				float cot_b = cos_b / sqrt(1 - cos_b * cos_b);

				//当前顶点所在面的曲率K(xi)（求和）
				Hn -= (cot_a + cot_b) * (P - p1);

			}
			Hn /= (4.0f * A);
			auto newP = P + (float)m_Lamdba * Hn;
			//Hn.normalize();//法向量归一化
			mesh->SetPoint(i, newP);//更新顶点坐标，继续迭代


		}

	}
	std::cout << "end" << std::endl;
	mesh->GarbageCollection();

	SetOutput(mesh); // 设置输出数据
	return true;
}
IGAME_NAMESPACE_END