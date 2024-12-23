#ifndef iGameMinPath_h
#define iGameMinPath_h

#include"iGameObject.h"
#include"iGameSurfaceMesh.h"
#include"iGameFilter.h"
#include"iGamePoints.h"
#include"iGamePointSet.h"
#include<iostream>

IGAME_NAMESPACE_BEGIN
struct node {
    int id;
    double dis;
    node(int id, double d) { 
        this->id = id;
        this->dis = d;
    }
    bool operator<(const node& rhs) const { return dis > rhs.dis;
    }
};
struct PathInfo {
    int s_p, e_p;
    double length;
    std::vector<int> path;
    PathInfo(int a, int b, double c) {
        s_p = a;
        e_p = b;
        length = c;
    }
    bool operator<(const PathInfo& rhs) const { return length > rhs.length; }
};
class MinPath : public  Filter 
{
public: 
    I_OBJECT(MinPath);
    static Pointer New() { return new MinPath; }

    bool Execute() override;
    //求二点之间的最短距离
    void Dijkstra(SurfaceMesh& Mesh, int s_p, int e_p, std::vector<int>& path);

    void Dijkstra_group(SurfaceMesh& Mesh, std::vector<int>& lmk,
                        std::vector<std::vector<int>>& path);

    void SetMesh(SurfaceMesh::Pointer mesh) { this->mesh = mesh; }

    void setPath(const std::vector<int>& newPath) { path = newPath; }

    const std::vector<int>& getPath() const { return path; }

    void setPath2(const std::vector<std::vector<int>>& newPath) {
        path2 = newPath;
    }

    const std::vector<std::vector<int>>& getPath2() const { return path2; }

    void setMPath(bool a) { m_path = a; }

protected:
    SurfaceMesh::Pointer mesh{};
    MinPath() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
        m_path = false;
    };
    MinPath(SurfaceMesh::Pointer _mesh) { this->mesh = _mesh; }
    std::vector<int> path;
    std::vector<std::vector<int>>  path2;
    bool m_path;
};

IGAME_NAMESPACE_END
#endif // !1

