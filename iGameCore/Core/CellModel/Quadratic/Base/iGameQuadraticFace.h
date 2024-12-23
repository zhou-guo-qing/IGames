#ifndef iGameQuadraticFace_h
#define iGameQuadraticFace_h

#include "iGameQuadraticLine.h"

IGAME_NAMESPACE_BEGIN
class QuadraticFace : public Cell {
public:
	I_OBJECT(QuadraticFace);

	int GetNumberOfFaces() override { return 0; }
	Cell* GetEdge(const int edgeId) override
	{
		int edgeIdPlus = (edgeId + 1) % this->GetNumberOfEdges();
		int edgeIdCenter = (edgeId + this->GetNumberOfEdges()) % this->GetNumberOfEdges();
		m_Line->PointIds->SetId(0, this->PointIds->GetId(edgeId));
		m_Line->PointIds->SetId(1, this->PointIds->GetId(edgeIdPlus));
		m_Line->PointIds->SetId(2, this->PointIds->GetId(edgeIdCenter));
		m_Line->Points->SetPoint(0, this->Points->GetPoint(edgeId));
		m_Line->Points->SetPoint(1, this->Points->GetPoint(edgeIdPlus));
		m_Line->Points->SetPoint(2, this->Points->GetPoint(edgeIdCenter));
		return m_Line;
	}
	Cell* GetFace(const int) override { return nullptr; }

	virtual Vector3f GetNormal() = 0;

	igIndex GetEdgeId(const int id) { return this->EdgeIds->GetId(id); }

	IdArray::Pointer EdgeIds{};

protected:
	QuadraticFace()
	{
		m_Line = QuadraticLine::New();
		EdgeIds = IdArray::New();
	}
	~QuadraticFace() override = default;

	QuadraticLine::Pointer m_Line;
};
IGAME_NAMESPACE_END
#endif
