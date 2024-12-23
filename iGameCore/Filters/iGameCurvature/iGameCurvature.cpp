#include"iGameCurvature.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include"iGameVector.h"

#define PI 3.14159265359

IGAME_NAMESPACE_BEGIN



//计算a,b,c三角形的外接圆的圆心
Point cal_circum_enter(const Point& a, const Point& b, const Point& c) {
    Point ac = c - a, ab = b - a;
    Point abXac = ab.cross(ac), abXacXab = abXac.cross(ab),
             acXabXac = ac.cross(abXac);
    return a + (abXacXab * ac.dot(ac) + acXabXac * ab.dot(ab)) / (2.0 * abXac.dot(abXac));
}

 double vectorAngle(const Vector3f& a, const Vector3f& b) {
    double cosTheta = a.dot(b);
    double sinTheta =(a.cross(b)).norm();
    return atan2(sinTheta, cosTheta);
}
void Curvature::cal_local_ave_region(SurfaceMesh* const mesh,
                          std::vector<double>& vertexLAR) {

    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    auto Faces = mesh->GetFaces();
    igIndex FaceNum = mesh->GetNumberOfFaces();

    std::vector<std::vector<Point>> FaceP(FaceNum);
    for (int i = 0; i < FaceNum; i++) { 
        bool isObtuseAngle = false;
        Point obtuseVertexHandle;
        
        vcnt = Faces->GetCellIds(i, vhs);
        Point p[3];
        for (int i = 0; i < vcnt; i++) { 
            p[i] = mesh->GetPoint(vhs[i]);
            //std::cout << "vcnt=" << vcnt << "vhs=" << vhs[i] << std::endl;
            //std::cout << "p=" << p[i] << std::endl;
        }
        Point to01 = p[1] - p[0];
        Point to12 = p[2] - p[1];
        Point to20 = p[0] - p[2];
        to01.normalize(), to12.normalize(), to20.normalize();
        
        if (vectorAngle(to01, -to20) > PI / 2.0)
        { 
            isObtuseAngle = true;
            obtuseVertexHandle = p[0];
        }
        else if (vectorAngle(to12, -to01) > PI / 2.0)
        {
            isObtuseAngle = true;
            obtuseVertexHandle = p[1];
        }
        else if (vectorAngle(to20, -to12) > PI / 2.0)
        {
            isObtuseAngle = true;
            obtuseVertexHandle = p[2];
        }

        if (isObtuseAngle) {
            double faceArea = 0.5 * to01.cross(-to20).norm();
            
            for (int k = 0; k < vcnt; k++)
            { 
                if (mesh->GetPoint(vhs[k]) == obtuseVertexHandle)
                    vertexLAR[vhs[k]] += faceArea / 2.0;
                else
                    vertexLAR[vhs[k]] += faceArea / 4.0;
            }
        } else {
            Vector3f cc = cal_circum_enter(p[0], p[1], p[2]);

            for (int j = 0; j < 3; j++) {
                Vector3f edgeMidpoint = (p[j % 3] + p[(j + 1) % 3])/2;

                double edgeLength = (p[j % 3] - p[(j + 1) % 3]).norm();
                double partArea = 0.5 * edgeLength * (edgeMidpoint - cc).norm();
                vertexLAR[vhs[j]] += 0.5 * partArea;
                vertexLAR[vhs[(j + 1)%3]] += 0.5 * partArea;
            }

        }

    }
}
double cotangent(const Vector3f& u, const Vector3f& v) {
    double dotProduct = u.dot(v);
    double uNorm = u.norm();
    double vNorm = v.norm();
    double cosine = dotProduct / (uNorm * vNorm);
    double sine = std::sqrt(1 - cosine * cosine);
    return cosine / sine;
}
void Curvature::cal_mean_curvature(SurfaceMesh* const mesh,
                         std::vector<double>& vertexLAR) {

    int vNum = mesh->GetNumberOfPoints();

    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];

    for (int i = 0; i < vNum; i++)
    { 
        Point tmp = {0, 0, 0}; //拉普拉斯算子
        Point vert = mesh->GetPoint(i);
        Point t_normal = {0, 0, 0}; //点的法向量
        fcnt = mesh->GetPointToNeighborFaces(i,fhs);


        for (int j = 0; j < fcnt; j++) { 
            
            auto f = mesh->GetFace(fhs[j]); //面
            //std::cout << "fhs" << j << "fhs" << fhs[j] << std::endl;
            Point p[2];
            int t = 0;
            //std::cout << "face points"<< f->GetNumberOfPoints() << std::endl;
            for (int k = 0; k < f->GetNumberOfPoints(); k++)
            { 
                if (f->GetPointId(k) != i)
                { 
                    p[t++] = f->GetPoint(k);
                }
            }
            double cotAlpha = cotangent(vert - p[0], p[1] - p[0]);
            double cotBeta = cotangent(vert - p[1], p[0] - p[1]);

            Vector3f e_0 = p[0] - vert;
            Vector3f e_1 = p[1] - vert;

            tmp += e_1 * cotAlpha + e_0 * cotBeta;

            t_normal += e_0.cross(e_1);
        }
        tmp = tmp /(vertexLAR[i]*4); //拉普拉斯
        t_normal.normalize(); //点的法向
        if (tmp.dot(t_normal) > 0) { 
            vertexmean[i] = -tmp.norm();
        }else
            vertexmean[i] = tmp.norm();
    }
}

void Curvature::cal_gaussian_curvature(SurfaceMesh* const mesh,
                            const std::vector<double>& vertexLAR) {
    int vNum = mesh->GetNumberOfPoints();

    for (int i = 0; i < vNum; i++)
    { 
        double angle_temp = 2 * PI;
        Point vert = mesh->GetPoint(i);
        igIndex fcnt = 0;
        igIndex fhs[IGAME_CELL_MAX_SIZE];
        
        fcnt = mesh->GetPointToNeighborFaces(i, fhs);
        for (int j = 0; j < fcnt; j++)
        {
            auto f = mesh->GetFace(fhs[j]); //面
            //std::cout << "fhs" << j << "fhs" << fhs[j] << std::endl;
            Point p[2];
            int t = 0;
            //std::cout << "face points" << f->GetNumberOfPoints() << std::endl;
            for (int k = 0; k < f->GetNumberOfPoints(); k++) {
                if (f->GetPointId(k) != i) { 
                    p[t++] = f->GetPoint(k); }
            }
            Vector3f e_0 = p[0] - vert;
            Vector3f e_1 = p[1] - vert;
            
            double angle = vectorAngle(e_0, e_1);
            angle_temp -= angle;
        }
        angle_temp /= vertexLAR[i];
        vertexGauss[i] = angle_temp;
        //std::cout << "gauss:" << i << ":" << vertexGauss[i] << std::endl;
 
    }
}
bool Curvature::Execute()
{
    std::cout << " Curvature start \n";
    UpdateProgress(0);
    mesh = DynamicCast<SurfaceMesh>(GetInput(0));

    if (mesh == nullptr) {
        std::cout << " mesh is null  \n";
        return false;
    }

    std::cout << "mesh gets \n";
    mesh->RequestEditStatus();

    std::vector<double> vertexLAR(mesh->GetNumberOfPoints(), 0);
    cal_local_ave_region(mesh, vertexLAR);

    //for (int i = 0; i < mesh->GetNumberOfPoints(); i++)
    //{
    //    std::cout << "LAR"<< vertexLAR[i] << std::endl;
    //}
    //cal_mean_curvature(mesh, vertexLAR);
    cal_gaussian_curvature(mesh, vertexLAR);

    std::cout << "end" << std::endl;
    mesh->GarbageCollection();

    SetOutput(mesh); // 设置输出数据
    return true;
}




IGAME_NAMESPACE_END