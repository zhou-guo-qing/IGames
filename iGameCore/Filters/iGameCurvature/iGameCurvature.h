#ifndef iGameCurvature_h
#define iGameCurvature_h

#include "iGameObject.h"
#include "iGameFilter.h"
#include"iGameSurfaceMesh.h"
#include"iGamePoints.h"
IGAME_NAMESPACE_BEGIN

class Curvature :public Filter{
public:
    I_OBJECT(Curvature);
    static Pointer New() { return new Curvature; }

    bool Execute() override;

    void cal_local_ave_region(SurfaceMesh* const mesh,
                              std::vector<double>& vertexLAR);
    void cal_mean_curvature(SurfaceMesh* const mesh,
                             std::vector<double>& vertexLAR);
    void cal_gaussian_curvature(SurfaceMesh* const mesh,
                                const std::vector<double>& vertexLAR);
    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }

    void setMean(const std::vector<double>& mean) { vertexmean = mean; }

    std::vector<double> getMean()  { return vertexmean; }
    void setGauss(const std::vector<double>& gauss) { vertexGauss = gauss; }

    std::vector<double> getGauss() { return vertexGauss; }

protected: 
    SurfaceMesh::Pointer mesh{};
    Curvature() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
        vertexmean.resize(50000); // 将向量大小设置为 500
        vertexGauss.resize(50000);
        for (int i = 0; i < 50000; ++i) {
            vertexmean[i] = 0.0; // 将所有元素初始化为 0.0
            vertexGauss[i] = 0.0;
        }
    };
    Curvature(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }

    std::vector<double> vertexmean; //平均曲率
    std::vector<double> vertexGauss; //高斯曲率

};

IGAME_NAMESPACE_END
#endif