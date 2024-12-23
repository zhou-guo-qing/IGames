#include"iGameBilateral.h"
#include"iGameVector.h"
#include<random>
IGAME_NAMESPACE_BEGIN


bool Bilateral::Execute()
{ 
	std::cout << " Bilateral start \n";
    UpdateProgress(0);
    mesh = DynamicCast<SurfaceMesh>(GetInput(0));

    if (mesh == nullptr) {
        std::cout << " mesh is null  \n";
        return false;
    }

    std::cout << "mesh gets \n";
    mesh->RequestEditStatus();
	


    int faceNum = mesh->GetNumberOfFaces();
    std::vector<Vector3f> faceNormal(faceNum); //面法向
    std::vector<double> faceArea(faceNum); //面的面积
    std::vector<Vector3f> faceCenter(faceNum); //面

    igIndex vcnt = 0 ;
    igIndex vhs[IGAME_CELL_MAX_SIZE]; //一个面的点id

    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE]; //面的一邻面的id

    igIndex vfcnt = 0;
    igIndex vfhs[IGAME_CELL_MAX_SIZE]; //点的一邻面id

    double SigmaCenter;

    int vNum = mesh->GetNumberOfPoints();

    auto faces = mesh->GetFaces();

    if (this->noise == true) {
        for (int i = 0; i < vNum; i++) { 
            double noiseStrength = 0.01;
            Point v = mesh->GetPoint(i);
            double noiseX =
                    ((double) rand() / RAND_MAX * 2 - 1) * noiseStrength;
            double noiseY =
                    ((double) rand() / RAND_MAX * 2 - 1) * noiseStrength;
            double noiseZ =
                    ((double) rand() / RAND_MAX * 2 - 1) * noiseStrength;
            v[0] += noiseX;
            v[1] += noiseY;
            v[2] += noiseZ;
            mesh->SetPoint(i, v);
        }
        this->noise = false;
        std::cout << "noise end" << std::endl;
        mesh->GarbageCollection();

        SetOutput(mesh); // 设置输出数据
        return true; 
    }

    for (int i = 0; i < faceNum; i++) {

        auto f = mesh->GetFace(i);
        faceNormal[i] = f->GetNormal();  //面法向

        vcnt = faces->GetCellIds(i, vhs);  //面的三个点id
        assert(vcnt == 3);
        Point p1 = mesh->GetPoint(vhs[0]);
        Point p2 = mesh->GetPoint(vhs[1]);
        Point p3 = mesh->GetPoint(vhs[2]);

        auto p12 = p2 - p1;
        auto p13 = p3 - p1;
        double area = 0.5 * p12.cross(p13).norm();
        faceArea[i] = area;

        faceCenter[i] = (p1 + p2 + p3) / 3;
    }

    SigmaCenter = 0;

    for (int i = 0; i < faceNum; i++) { 
        fcnt = mesh->GetFaceToOneRingFaces(i, fhs);
        for (int j = 0; j < fcnt; j++) {
            SigmaCenter += (faceCenter[i] - faceCenter[fhs[j]]).norm();
        }
    }
    SigmaCenter /= faceNum * 3;

    double NormalIterNum=this->m_NormalIterNum, VertexIterNum=this->m_VertexIterNum;
    double SigmaNormal =this->m_SigmaNormal;

    for (int i = 0; i < NormalIterNum; i++) {
        for (int j = 0; j < faceNum; j++) { 
            
            double kp = 0;
            Vector3f NewN = {0, 0, 0};
            fcnt = mesh->GetFaceToOneRingFaces(i, fhs);
            for (int k = 0; k < fcnt; k++)
            {
                double delta_center = (faceCenter[j] - faceCenter[fhs[k]]).norm();
                double delta_normal = (faceNormal[j] - faceNormal[fhs[k]]).norm();
                double Aj = faceArea[fhs[k]];
                double Ws = exp(-delta_center * delta_center /
                                (2 * SigmaCenter * SigmaCenter));
                double Wr = exp(-delta_normal * delta_normal /
                                (2 * SigmaNormal * SigmaNormal));
                NewN += faceNormal[fhs[k]] * Aj * Ws * Wr;
                kp += Aj * Ws * Wr;
            }
            faceNormal[j] = NewN / kp;
            faceNormal[j] /= faceNormal[j].norm();
        }
    }
    for (int i = 0; i < VertexIterNum; i++) {
        for (int j = 0; j < vNum; j++) { 
            
            Point xi = mesh->GetPoint(j);
            Point delta_xi = {0, 0, 0};
            int Nei_count = 0;
            vfcnt = mesh->GetPointToNeighborFaces(j, vfhs);
            for (int k = 0; k < vfcnt; k++) { 
                
                Nei_count++;
                Point cj = faceCenter[vfhs[k]];
                Point nj = faceNormal[vfhs[k]];
                delta_xi = delta_xi + nj * (nj[0]*(cj-xi)[0]+nj[1]*(cj-xi)[1]+nj[2]*(cj-xi)[2]);

            }
            xi = xi + delta_xi / Nei_count;
            mesh->SetPoint(j, xi);
        }
    }
    std::cout << "Bilateral end" << std::endl;
    mesh->GarbageCollection();

    SetOutput(mesh); // 设置输出数据
	return true; 
}

IGAME_NAMESPACE_END