#include"iGameRemesh.h"
//#include"iGameVector.h"
#include "AABB_Tree.h"
#include <vector>
#include"TinyVector.h"
AABB_Tree* abtree;
IGAME_NAMESPACE_BEGIN

bool Remesh::is_collapse_ok(int eid, double low, double high)
{

    if (mesh->IsEdgeDeleted(eid) || mesh->IsBoundaryEdge(eid)) return false;
    igIndex fids[2];
    igIndex vids[2];
    
    mesh->GetEdgePointIds(eid, vids);

    igIndex survivalVid = vids[1];
    igIndex extinctVid = vids[0];

    Point survivalPoint = mesh->GetPoint(vids[1]);
    Point extinctPoint = mesh->GetPoint(vids[0]);

    if (mesh->IsBoundaryPoint(survivalVid) || mesh->IsBoundaryPoint(extinctVid))
        return false;
    if ((extinctPoint - survivalPoint).length() > low) return false;

    igIndex OneRingVids[64];
    int vvnum = mesh->GetPointToOneRingPoints(extinctVid, OneRingVids);
    for (int i = 0; i < vvnum; i++)
        if ((mesh->GetPoint(OneRingVids[i]) - survivalPoint).length() > high)
            return false;

    int facenum = mesh->GetEdgeToNeighborFaces(eid, fids);
    for (int i = 0; i < facenum; i++) {
        igIndex feids[64];
        int feNum = mesh->GetFaceEdgeIds(fids[i], feids);
        assert(feNum == 3);
        int loc = 0;
        while (feids[loc] != eid) loc++;
        if (mesh->IsBoundaryEdge(feids[(loc + 1) % 3]) &&
            mesh->IsBoundaryEdge(feids[(loc + 2) % 3]))
            return false;
    }

    return true;

}


double Remesh::calculateTargetEdgeLength() {

    double target_edge_length = 0.0;
    int eNum = mesh->GetNumberOfEdges();
    auto edges = mesh->GetEdges();

    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    for (int i = 0; i < mesh->GetNumberOfEdges();i++){

        vcnt = edges->GetCellIds(i, vhs);
        assert(vcnt == 2);
        Point p1 = mesh->GetPoint(vhs[0]);
        Point p2 = mesh->GetPoint(vhs[1]);
        Point pp = p2 - p1;
        target_edge_length += pp.length();
    }
    target_edge_length /= eNum;
    return target_edge_length;
}
void Remesh::split_long_edges(double high) {
    std::cout << "split_long_edges start" << std::endl;

    int NE = mesh->GetNumberOfEdges();
    auto edges = mesh->GetEdges();
    auto faces = mesh->GetFaces();
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];
    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];
    igIndex ecnt = 0;
    igIndex ehs[IGAME_CELL_MAX_SIZE];
    for (int i = 0; i < NE; i++) {
        vcnt = edges->GetCellIds(i, vhs);
        assert(vcnt == 2);
        Point p1 = mesh->GetPoint(vhs[0]);
        Point p2 = mesh->GetPoint(vhs[1]);
        Point pp = p2 - p1;
        double len = pp.length();
        if (len > high) {
            //分裂长边，增加点和二条边
            Point P12center = (p1 + p2) / 2.0;
            
            fcnt = mesh->GetEdgeToNeighborFaces(i, fhs);
            
            assert(fcnt == 2 || fcnt == 1);

            IGsize newp = mesh->AddPoint(P12center);
            for (int j = 0; j < fcnt; j++) {

                ecnt=mesh->GetFacePointIds(fhs[j], ehs);
                assert(ecnt == 3);
                int fvid = ehs[0] + ehs[1] + ehs[2] - vhs[0] - vhs[1];

                igIndex face[3]{vhs[0], newp, fvid};
                igIndex face2[3]{vhs[1], newp, fvid};
                mesh->AddFace(face, 3);
                mesh->AddFace(face2, 3);
                
            }
            mesh->DeleteEdge(i);
        }
    }
    std::cout << "split_long_edges end" << std::endl;
}


void Remesh::collapse_short_edges(double high, double low) {
    std::cout << "collapse_short_edges start" << std::endl;
    int ppv[2];
    int NE = mesh->GetNumberOfEdges();
    auto edges = mesh->GetEdges();
    auto faces = mesh->GetFaces();
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    igIndex ecnt = 0;
    igIndex ehs[IGAME_CELL_MAX_SIZE]; 

    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];

    igIndex ffcnt = 0;
    igIndex ffhs[IGAME_CELL_MAX_SIZE];

    igIndex vvcnt = 0;
    igIndex vvhs[IGAME_CELL_MAX_SIZE];
    for (int i = NE - 1; i >= 0; i--) {
        if (!is_collapse_ok(i, low, high)) continue;
        if (i > mesh->GetNumberOfEdges() - 1) continue;

        auto e = mesh->GetEdge(i);
        vcnt = edges->GetCellIds(i, vhs);
        assert(vcnt == 2);
        Point p1 = mesh->GetPoint(vhs[0]);
        Point p2 = mesh->GetPoint(vhs[1]);

        
        if (mesh->IsBoundaryPoint(vhs[0]) || mesh->IsBoundaryPoint(vhs[1])) continue;  
        Point pp = p2 - p1;
        double len = pp.length();

        if (len < low) {
            
            bool is_collapse = true;
            ecnt = mesh->GetPointToOneRingPoints(vhs[0], ehs);
            for (int j = 0; j < ecnt; j++) {

                Point vv = mesh->GetPoint(ehs[j]);
                double len = (p1 - vv).norm();
                if (len > high)
                { 
                    is_collapse = false;
                    break;
                }
            }
            int ve[2];
            int num = 0;
            if (is_collapse)
            {
                //坍塌这条边
                Point P12center = (p1 + p2) / 2;

                igIndex fcnt = 0;
                igIndex fhs[IGAME_CELL_MAX_SIZE];

                igIndex fvcnt = 0;
                igIndex fvhs[IGAME_CELL_MAX_SIZE];

                fcnt = mesh->GetPointToNeighborFaces(vhs[0], fhs); //0的1领域面

                for (int i = 0; i < fcnt; i++) {

                    fvcnt = faces->GetCellIds(fhs[i], fvhs);

                    bool isflag = false;
                    assert(fvcnt == 3);
                    int num = 0;
                    for (int j = 0; j < fvcnt; j++) {

                        if (fvhs[j] == vhs[1]) {
                            isflag = true;
                            break;
                        }
                        if (fvhs[j] != vhs[0]) {
                            ppv[num++] = fvhs[j];
                            isflag = false;
                        }
                    }
                    
                    if (!isflag) {
                        assert(num == 2);
                        igIndex face[3]{vhs[1], ppv[0], ppv[1]};
                        mesh->AddFace(face, 3);
                    }
                }
                mesh->SetPoint(vhs[1], P12center);
                mesh->DeletePoint(vhs[0]);
            }
            
        }
    }
    std::cout << "collapse_short_edges end" << std::endl;
}
IGsize vv[2] = {0, 0};
void Remesh::equalize_valences() {
    std::cout << "equalize_valences start" << std::endl;
   
    std::vector<int> target_valence;
    int deviation_pre, deviation_post;
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    igIndex ecnt = 0;
    igIndex ehs[IGAME_CELL_MAX_SIZE];

    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];

    igIndex cnt = 0;
    igIndex chs[IGAME_CELL_MAX_SIZE];
    int vNum = mesh->GetNumberOfPoints();
    for (int i = 0; i < vNum; i++) { 

        if (mesh->IsBoundaryPoint(i)) target_valence.push_back(4);
        else
            target_valence.push_back(6);
    }
    int eNum = mesh->GetNumberOfEdges();
    auto edges = mesh->GetEdges();
    auto faces = mesh->GetFaces();

    for (int i = 0; i < eNum; i++) {
        if (mesh->IsBoundaryEdge(i)) continue;
        auto e = mesh->GetEdge(i);
        vcnt = edges->GetCellIds(i, vhs);
        assert(vcnt == 2);
        Point p1 = mesh->GetPoint(vhs[0]);
        Point p2 = mesh->GetPoint(vhs[1]);

        fcnt = mesh->GetEdgeToNeighborFaces(i, fhs);
        
        if (fcnt > 2) { std::cout << "error" << std::endl;
        }
        int num = 0;
        for (int j = 0; j < fcnt; j++) {

            ecnt=faces->GetCellIds(fhs[j], ehs);
            assert(ecnt == 3);
            
            for (int k = 0; k < ecnt; k++) { 

                if (ehs[k] != vhs[0] && ehs[k] != vhs[1]) { 
                    
                    vv[num++] = ehs[k];
                }
            }
        }
   
        int v1 = mesh->GetPointToOneRingPoints(vhs[0], chs);
        int v2 = mesh->GetPointToOneRingPoints(vhs[1], chs);
        int v3 = mesh->GetPointToOneRingPoints(vv[0], chs);
        int v4 = mesh->GetPointToOneRingPoints(vv[1], chs);

        deviation_pre = abs(int(v1 - target_valence[vhs[0]])) +
                        abs(int(v2 - target_valence[vhs[1]])) +
                        abs(int(v3 - target_valence[vv[0]])) +
                        abs(int(v4 - target_valence[vv[1]]));
        deviation_post = abs(int(v1 - 1 - target_valence[vhs[0]])) +
                         abs(int(v2 - 1 -target_valence[vhs[1]])) +
                         abs(int(v3 + 1 -target_valence[vv[0]])) +
                         abs(int(v4 + 1 -target_valence[vv[1]]));
        if (deviation_pre > deviation_post) {

            //翻转边
            igIndex face[3]{vhs[0], vv[0], vv[1]};
            igIndex face2[3]{vhs[1], vv[0], vv[1]};
            mesh->AddFace(face, 3);
            mesh->AddFace(face2, 3);
            mesh->DeleteEdge(i); 
        }
    }
    std::cout << "equalize_valences end" << std::endl;
}
void Remesh::tangential_relaxation() {
    std::cout << "tangential_relaxation satrt" << std::endl;
    
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];
    int vNum = mesh->GetNumberOfPoints();
    for (int i = 0; i < vNum; i++)
    {
        if (mesh->IsBoundaryPoint(i)) continue;
        double count = 0;
        Point p = mesh->GetPoint(i);
        Point q(0, 0, 0);
        vcnt = mesh->GetPointToOneRingPoints(i, vhs);
        for (int j = 0; j < vcnt; j++) {

            q += mesh->GetPoint(vhs[j]);
            count++;
        }
        q /= count;

        fcnt = mesh->GetPointToNeighborFaces(i, fhs);
        Point n(0, 0, 0);
        for (int k = 0; k < fcnt; k++) {
            
            n += mesh->GetFace(fhs[k])->GetNormal();
        }
        n.normalize();
        mesh->SetPoint(i, q + (n.dot(p-q)) * n);
    }
    std::cout << "tangential_relaxation end" << std::endl;
}
void Remesh::project_to_surface() {
    std::cout << "project_to_surface start" << std::endl;
    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];
    
    int vNum = mesh->GetNumberOfPoints();
    for (int i = 0; i < vNum; i++) {

        Point p = mesh->GetPoint(i);
        Vertex pp(p[0], p[1], p[2]);
        Vertex ab_nearst_point;
        abtree->findNearstPoint(pp, ab_nearst_point);
         
        Point new_point(ab_nearst_point[0], ab_nearst_point[1],
                        ab_nearst_point[2]);
        mesh->SetPoint(i, new_point);
    }
    std::cout << "project_to_surface end" << std::endl;
}
void Remesh::get_AABB_tree() {
    std::vector<Vertex> point_set;
    point_set.clear();

    igIndex fcnt = 0;
    igIndex fhs[IGAME_CELL_MAX_SIZE];
    auto faces = mesh->GetFaces();

    int fNum = mesh->GetNumberOfFaces();
    for (int i = 0; i < fNum; i++) {

        fcnt = faces->GetCellIds(i, fhs);
        for (int j = 0; j < fcnt; j++) {
            Vertex a;
            Point pp = mesh->GetPoint(fhs[j]);
            a[0] = pp[0];
            a[1] = pp[1];
            a[2] = pp[2];
            point_set.push_back(a);
        }
        
    }
    abtree = new AABB_Tree(point_set);
}
bool Remesh::Execute()
{
    std::cout << "Remesh start \n";
    UpdateProgress(0);
    mesh = DynamicCast<SurfaceMesh>(GetInput(0));
    if (mesh == nullptr) {
        std::cout << " mesh is null  \n";
        return false;
    }
    std::cout << "mesh gets \n";
    

    double target_edge_length, high, low;
    target_edge_length = calculateTargetEdgeLength() / 2.0;
    high = 4.0 / 3.0 * target_edge_length;
    low = 4.0 / 5.0 * target_edge_length;

    get_AABB_tree();
      

    for (int i = 0; i < 10; i++) {
        mesh->RequestEditStatus();
        split_long_edges(high);
        collapse_short_edges(high, low);
        equalize_valences();
        tangential_relaxation();
        project_to_surface();
        mesh->GarbageCollection();
    }


    std::cout << "Remesh end" << std::endl;
    SetOutput(mesh); // 设置输出数据
    return true; 

}

IGAME_NAMESPACE_END