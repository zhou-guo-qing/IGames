#ifndef MINMALSURFACE_h
#define MINMALSURFACE_h

#include "iGameFilter.h"
#include "iGamePointSet.h"
#include "iGameSurfaceMesh.h"
#include "iGamePoints.h"
#include "iGameObject.h"
#include "iGameQuad.h"
#include "../CellModel/Quadratic/iGameQuadraticTriangle.h"
#include <random>

IGAME_NAMESPACE_BEGIN
class minmalSurface : public Filter {
public:
	I_OBJECT(minmalSurface);
	static Pointer New() { return new minmalSurface; }

	//void SetFilterRate(double value) {
	//	m_FilterRate = value;
	//}
	bool Execute()override;

	void GlobalLaplacianSmoothing(SurfaceMesh::Pointer mesh);
	void Parameterization(SurfaceMesh::Pointer mesh);


	void SetMesh(SurfaceMesh::Pointer mesh) {
		this->mesh = mesh;
	}
	void SetGlobal(bool b)
	{
		m_Global = b;
	}
	void SetParzmeterization(bool b)
	{
		m_Parameterization = b;
	}
protected:
	SurfaceMesh::Pointer mesh{};
	minmalSurface() {
		mesh = nullptr;
		this->SetNumberOfOutputs(1);
		this->SetNumberOfInputs(1);
		this->m_Global = false;
		this->m_Parameterization = false; 
	};
	minmalSurface(SurfaceMesh::Pointer _mesh) {
		this->mesh = _mesh;
	}
	// double m_FilterRate{ 0.5 };
	bool m_Global;
	bool m_Parameterization;
};
IGAME_NAMESPACE_END
#endif