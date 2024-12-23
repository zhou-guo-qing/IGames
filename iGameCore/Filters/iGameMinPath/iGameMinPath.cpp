#include "iGameMinPath.h"
#include "iGameFilter.h"
#include <queue>
#include <algorithm>
#include <unordered_map>
#include <numeric>
#define DBL_MAX 1.7976931348623158e+308 // double最大值
IGAME_NAMESPACE_BEGIN

double Length(Vector3f a,Vector3f b) {
    Vector3f t = a - b;
    return t.norm();
}
//计算数组lmk之间的点的路径
void MinPath::Dijkstra(SurfaceMesh& Mesh, int s_p, int e_p, std::vector<int>& path) {
    std::cout << "Dijkstra strat" << std::endl;
    std::priority_queue<node> que;

    std::vector<double> distance(Mesh.GetNumberOfPoints(), DBL_MAX);
    std::vector<int> v_p(Mesh.GetNumberOfPoints(), -1);
    std::vector<int> is_visited(Mesh.GetNumberOfPoints(), 0);

    igIndex vcnt = 0;
    igIndex vhs[IGAME_CELL_MAX_SIZE];

    v_p[s_p] = s_p;
    distance[s_p] = 0;
    que.push(node(s_p,0));

    while (!que.empty()) {

        node tmp = que.top();
        que.pop();
        if (is_visited[tmp.id]) continue;
        if (tmp.id == e_p) break;
        is_visited[tmp.id] = 1;

        auto v1 = Mesh.GetPoint(tmp.id);
        vcnt = Mesh.GetPointToOneRingPoints(tmp.id, vhs);

        for (int i = 0; i < vcnt; i++) {
            auto v2 = Mesh.GetPoint(vhs[i]);
            auto e12 = Mesh.GetEdgeIdFormPointIds(tmp.id, vhs[i]);

            if (distance[tmp.id] + Length(v1, v2) < distance[vhs[i]]) {
                distance[vhs[i]] = distance[tmp.id] + Length(v1, v2);
                que.push(node(vhs[i], distance[vhs[i]]));
                v_p[vhs[i]] = tmp.id;
            }
        }
    }
    path.clear();
    path.push_back(e_p);
    do { 
        if (e_p < Mesh.GetNumberOfPoints())
        {
            e_p = v_p[e_p];
            path.push_back(e_p);
        }
    } while (e_p < Mesh.GetNumberOfPoints() && e_p != v_p[e_p]);
}


void MinPath::Dijkstra_group(SurfaceMesh& Mesh, std::vector<int>& lmk,
                    std::vector<std::vector<int>>& path) {
    std::sort(lmk.begin(), lmk.end());
    std::priority_queue<PathInfo> complete_graph;
    for (int i = 0; i < lmk.size()-1; i++) {
        std::vector<int> is_lmk(Mesh.GetNumberOfPoints(), 0);
        for (int j = i; j < lmk.size(); j++) { is_lmk[lmk[j]] = 1; }
        int count = lmk.size() - i;
        int s_p = lmk[i];
        std::vector<int> is_visited(Mesh.GetNumberOfPoints(), 0);
        std::vector<double> distance(Mesh.GetNumberOfPoints(), DBL_MAX);
        std::priority_queue<node> que;
        std::vector<int> v_p(Mesh.GetNumberOfPoints(), -1);

        igIndex vcnt = 0;
        igIndex vhs[IGAME_CELL_MAX_SIZE];

        v_p[s_p] = s_p;
        distance[s_p] = 0;
        que.push(node(s_p,0));
        while (count != 0 && !que.empty()) {
            node tmp = que.top();
            que.pop();
            if (is_visited[tmp.id] == 1) continue;

            if (is_lmk[tmp.id] == 1) {
                count--;
                is_lmk[tmp.id] = 0;
            }

            is_visited[tmp.id] = 1;

            auto v1 = Mesh.GetPoint(tmp.id);
            vcnt = Mesh.GetPointToOneRingPoints(tmp.id, vhs);

            for (int i = 0; i < vcnt; i++) {
                auto v2 = Mesh.GetPoint(vhs[i]);
                auto e12 = Mesh.GetEdgeIdFormPointIds(tmp.id, vhs[i]);

                if (distance[tmp.id] + Length(v1, v2) < distance[vhs[i]]) {
                    distance[vhs[i]] = distance[tmp.id] + Length(v1, v2);
                    que.push(node(vhs[i], distance[vhs[i]]));
                    v_p[vhs[i]] = tmp.id;
                }
            }
        }
            for (int j = i + 1; j < lmk.size(); j++) {
                std::vector<int> p;
                int e_p = lmk[j];
                double l = distance[e_p];
                p.push_back(e_p);
                do {
                    e_p = v_p[e_p];
                    p.push_back(e_p);

                } while (e_p != v_p[e_p]);
                PathInfo p_info(s_p, lmk[j], l);
                p_info.path = p;
                complete_graph.push(p_info);
            }
        }
        std::vector<PathInfo> Mst;
        int nv = lmk.size();
        std::vector<int> spanning_tree_current(nv);
        for (int i = 0; i < nv; i++)
        {
            spanning_tree_current[i] = i;
        }
        int index = 0;
        int j = 0;
        while (index < nv - 1)
        {   
            PathInfo tmp = complete_graph.top();
            complete_graph.pop();
            int m_id = -1, n_id = -1;
            for (int u = 0; u < nv; u++) {
                if (lmk[u] == tmp.s_p) m_id = u;
                if (lmk[u] == tmp.e_p) n_id = u;
            }
            m_id = spanning_tree_current[m_id];
            n_id = spanning_tree_current[n_id];
            if (m_id < n_id) std::swap(m_id, n_id);
            if (m_id != n_id) {
                path.push_back(tmp.path);
                index++;
                for (int i = 0; i < nv; ++i) {
                    if (spanning_tree_current[i] == n_id)
                        spanning_tree_current[i] = m_id;
                }
            }
        }
}
bool MinPath::Execute() 
{ 
    std::cout << "MinPath start" << std::endl; 
    UpdateProgress(0);
    mesh = DynamicCast<SurfaceMesh>(GetInput(0));
    if (mesh == nullptr) {
        std::cout << " mesh is null  \n";
        return false;
    }

    std::cout << "mesh gets \n";
    mesh->RequestEditStatus();
    if (m_path) {
        std::vector<int> path(mesh->GetNumberOfPoints(),-1); 
        Dijkstra(*mesh, 1, 200, path);
        setPath(path);
    } else {
        std::vector<std::vector<int>> path; 
        std::vector<int> lmk = {1, 20,50,100, 200};
        Dijkstra_group(*mesh, lmk, path);
        setPath2(path);
    }
    std::cout << "end" << std::endl;
    mesh->GarbageCollection();

    SetOutput(mesh); // 设置输出数据
    return true;

}

IGAME_NAMESPACE_END