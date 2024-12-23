#ifndef iGameDeNoise_h
#define iGameDeNoise_h

#include "iGameFilter.h"
#include "iGamePointSet.h"
#include "iGameSurfaceMesh.h"
#include "iGamePoints.h"
#include "iGameObject.h"
#include "iGameQuad.h"
#include "../CellModel/Quadratic/iGameQuadraticTriangle.h"
#include <random>

IGAME_NAMESPACE_BEGIN
class DeNoise : public Filter {
public:
	I_OBJECT(DeNoise);
	static Pointer New() { return new DeNoise; }

	//void SetFilterRate(double value) {
	//	m_FilterRate = value;
	//}
	bool Execute()override;

	void SetMesh(SurfaceMesh::Pointer mesh) {
		this->mesh = mesh;
	}
	void SetLamdaAndCnt(double l, double c)
	{
		this->m_Lamdba = l;
		this->m_Cnt = c;
	}

protected:
	SurfaceMesh::Pointer mesh{};
	DeNoise() {
		mesh = nullptr;
		this->SetNumberOfOutputs(1);
		this->SetNumberOfInputs(1);
		this->m_Lamdba= 0.1;
		this->m_Cnt = 10;
	};
	DeNoise(SurfaceMesh::Pointer _mesh) {
		this->mesh = _mesh;
	}
	// double m_FilterRate{ 0.5 };
	double m_Lamdba;
	double m_Cnt;
};
IGAME_NAMESPACE_END
#endif