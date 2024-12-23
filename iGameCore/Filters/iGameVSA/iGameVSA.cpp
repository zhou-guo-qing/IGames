#include "iGameVSA.h"
#include "iGameVector.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#define PI 3.1415926
IGAME_NAMESPACE_BEGIN
double VSA::cal_area(int i) {

    auto faces = mesh->GetFaces();

    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    vcnt = faces->GetCellIds(i, vhs); //面的三个点id
    assert(vcnt == 3);
    Point p1 = mesh->GetPoint(vhs[0]);
    Point p2 = mesh->GetPoint(vhs[1]);
    Point p3 = mesh->GetPoint(vhs[2]);

    auto p12 = p2 - p1;
    auto p13 = p3 - p1;
    double area = 0.5 * p12.cross(p13).norm();

    return area;
}

bool VSA::Execute()
{
	std::cout << "VSA start \n";
    UpdateProgress(0);
    mesh = DynamicCast<SurfaceMesh>(GetInput(0));

    if (mesh == nullptr) {
        std::cout << " mesh is null  \n";
        return false;
    }
    std::cout << "mesh gets \n";
    mesh->RequestEditStatus();

    int k_num = 10;
    std::vector<Point> plane(k_num);

    int f_num = mesh->GetNumberOfFaces();
    plane.resize(k_num);
    partition.clear();
    partition.resize(f_num, -1);
    srand(time(0));
    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];

    std::vector<int> seed_set(k_num);
    std::vector<bool> is_conq(f_num, false);

    for (int i = 0; i < k_num;) {
        int ran = rand() % f_num;
        if (partition[ran] == -1) {
            partition[ran] = i;
            plane[i] = mesh->GetFace(i)->GetNormal(); 
            seed_set[i] = ran;  //种子面，K个
            i++;
        }
    }
    double last_energy = 0;

    for (int itnum = 0;; itnum++) {
    
        std::vector<double> small_energy(k_num, 1e10);
        double total_energy = 0;
    
        for (int i = 0; i < f_num; i++) {
            if (partition[i] != -1) {
                double energy = (plane[partition[i]] - mesh->GetFace(i)->GetNormal()).norm();
                total_energy += energy * cal_area(i);
                if (energy < small_energy[partition[i]]) {
                    small_energy[partition[i]] = energy;
                    seed_set[partition[i]] = i;
                }
            }
        }
        if (itnum > 0)
            std::cout << "Iteration " << itnum
                      << " Total energy : " << total_energy << "\n";

        if (itnum > 1 &&
                    abs(total_energy - last_energy) / total_energy < 1e-4 ||
            itnum > 20)
            break;

        last_energy = total_energy;

        for (int i = 0; i < f_num; i++) is_conq[i] = false;
        std::priority_queue<Tri> que;

        for (int i = 0; i < k_num; i++) {
            partition[seed_set[i]] = i;
            is_conq[seed_set[i]] = true;

            fcnt = mesh->GetFaceToNeighborFaces(seed_set[i],fhs);

            for (int k = 0; k < fcnt; k++) {

                double energy =
                        (mesh->GetFace(fhs[k])->GetNormal() - plane[i]).norm();
                que.push(Tri(fhs[k], i, energy));

            }
        }
        while (que.size() != 0) {
            int fid = que.top().id;
            int pid = que.top().tag;
            que.pop();
            if (!is_conq[fid]) {
                partition[fid] = pid;
                is_conq[fid] = true;

                fcnt = mesh->GetFaceToNeighborFaces(fid,fhs);

                for (int k = 0; k < fcnt; k++) {

                    if (!is_conq[fhs[k]]) {
                        double energy =
                                (mesh->GetFace(fhs[k])->GetNormal() -
                                     plane[pid]).norm();
                        que.push(Tri(fhs[k], pid, energy));
                    }

                }
            }
        }
        std::vector<Point> new_normal(k_num, Point(0, 0, 0));
        for (int j = 0; j < f_num; j++) {
            double area = cal_area(j);
            new_normal[partition[j]] += mesh->GetFace(j)->GetNormal() * area;
        }
        for (int i = 0; i < k_num; i++) {
            new_normal[i].normalize();
            plane[i] = new_normal[i];
        }

    }
    std::cout << "VSA end" << std::endl;
    mesh->GarbageCollection();

    SetOutput(mesh); // 设置输出数据
    return true;
}

IGAME_NAMESPACE_END