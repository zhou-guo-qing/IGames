#ifndef iGameTriangle_h
#define iGameTriangle_h

#include "iGameFace.h"

IGAME_NAMESPACE_BEGIN
class Triangle : public Face {
public:
	I_OBJECT(Triangle);
	static Pointer New() { return new Triangle; }

	IGenum GetCellType() const noexcept override { return IG_TRIANGLE; }
    int GetCellSize() const noexcept override { return 3; }
    int GetNumberOfEdges() override { return 3; }

	Vector3f GetNormal() override
	{
		const Vector3f& v0 = this->Points->GetPoint(0);
		const Vector3f& v1 = this->Points->GetPoint(1);
		const Vector3f& v2 = this->Points->GetPoint(2);
		return this->Self::ComputrNormal(v0, v1, v2);
	}

	static Vector3f ComputrNormal(const Vector3f& v0, const Vector3f& v1, const Vector3f& v2)
	{
		Vector3f d10 = v1 - v0;
		Vector3f d20 = v2 - v0;

		Vector3f normal = CrossProduct(d10, d20);
		normal.normalize();
		return normal;
	}

protected:
	Triangle()
	{
		this->Points->SetNumberOfPoints(3);
		this->PointIds->SetNumberOfIds(3);
		for (int i = 0; i < 3; i++)
		{
			this->Points->SetPoint(i, 0.0, 0.0, 0.0);
			this->PointIds->SetId(i, 0);
		}
	}
	~Triangle() override = default;
};
IGAME_NAMESPACE_END
#endif