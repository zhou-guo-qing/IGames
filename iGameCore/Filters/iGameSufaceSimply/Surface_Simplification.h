//#ifndef IGAMESUFACESIMPLY_h
//#define IGAMESUFACESIMPLY_h
//
//#include "iGameFilter.h"
//#include "iGameObject.h"
//#include "iGamePointSet.h"
//#include "iGamePoints.h"
//#include "iGameQuad.h"
//#include "iGameSurfaceMesh.h"
//#include "PolyMesh.h"
//#include "MVector3.h"
//#include "MPoint3.h"
//#include "PolyMeshIterators.h"
//#include "PolyMesh_Base.h"
//#include "Mrgb.h"
//#include "Massage.h"
//#include "MemoryPool.h"
//#include "IOManager.h"
//#include "unistd.h"
//#include <random>
//
//IGAME_NAMESPACE_BEGIN
//class SurFaceSimply : public Filter {
//public:
//    I_OBJECT(SurFaceSimply);
//    static Pointer New() { return new SurFaceSimply; }
//
//    //void SetFilterRate(double value) {
//    //	m_FilterRate = value;
//    //}
//    bool Execute() override;
//
//    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }
//
//protected:
//    SurfaceMesh::Pointer mesh{};
//    SurFaceSimply() {
//        mesh = nullptr;
//        this->SetNumberOfOutputs(1);
//        this->SetNumberOfInputs(1);
//
//    };
//    SurFaceSimply(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }
//    // double m_FilterRate{ 0.5 };
//};
//
//
//// t is a threshold parameter
//// void Surface_Simplification(MyMesh& mesh, float ratio);
//IGAME_NAMESPACE_END
//#endif