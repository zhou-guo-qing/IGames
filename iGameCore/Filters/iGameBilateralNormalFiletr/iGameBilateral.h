#ifndef iGame_Bilateral_H
#define iGame_Bilateral_H

#include "iGameObject.h"
#include "iGameFilter.h"
#include"iGameSurfaceMesh.h"
#include"iGamePoints.h"

IGAME_NAMESPACE_BEGIN
class Bilateral :public Filter
{
public:
    I_OBJECT(Bilateral);
    static Pointer New() { return new Bilateral; }
    bool Execute() override;
    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }

    void SetNoise(bool b){this->noise = b;}
    void SetSigmaNormal(double s) { this->m_SigmaNormal = s;}
    void SetNormalIterNum(double n) { this->m_NormalIterNum = n; }
    void SetVertexIterNum(double v) { this->m_VertexIterNum = v; }

protected:
    SurfaceMesh::Pointer mesh{};
    Bilateral() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
        noise = false;
        m_SigmaNormal = 1;
        m_VertexIterNum = 1;
        m_NormalIterNum = 1;
    };
    Bilateral(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }

    bool noise;
    double m_NormalIterNum , m_VertexIterNum;
    double m_SigmaNormal;
};
IGAME_NAMESPACE_END
#endif