#pragma once
#ifndef iGameOctree_h
#define iGameOctree_h
#include<iGamePoints.h>
#include <assert.h>
#include <vector>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include<iGameVolumeMesh.h>
using namespace iGame;
#define COMPUTE_SIDE(i, bit, p, mid, newMin, newMax) \
if (p >= mid)         \
{                     \
    i |= bit;         \
    newMin = mid;     \
}                     \
else                  \
{                     \
    newMax = mid;     \
}
// (x,y,z)   index(i)  < =  p<mid
// < < <        0
// >= < <       1
// < >= <       2
// >= >= <      3
// < < >=       4
// >= < >=      5
// < >= >=      6
// >= >= >=     7
struct Particle
{
    int id;
    Particle(int _id) : id(_id) {}
};
struct OctreePoint
{
    float x;
    float y;
    float z;
    OctreePoint(const Point& p2) : x(p2[0]), y(p2[1]), z(p2[2]) {}
    OctreePoint(const OctreePoint& p2) : x(p2.x), y(p2.y), z(p2.z) {}
    OctreePoint& operator=(const OctreePoint& p2) { x = p2.x; y = p2.y; z = p2.z; return *this; }
    OctreePoint(float in_x, float in_y, float in_z) : x(in_x), y(in_y), z(in_z) {}
    OctreePoint(const float p2[3]) : x(p2[0]), y(p2[1]), z(p2[2]) {}
    operator float* () { return &x; }
    operator const float* () const { return &x; }
    OctreePoint operator+(const OctreePoint& p2) const { return OctreePoint(x + p2.x, y + p2.y, z + p2.z); }
    OctreePoint operator-(const OctreePoint& p2) const { return OctreePoint(x - p2.x, y - p2.y, z - p2.z); }
    OctreePoint operator*(float f) const { return OctreePoint(x * f, y * f, z * f); }
    bool operator< (const OctreePoint& p2) const { return x < p2.x && y < p2.y && z < p2.z; }
    bool operator>=(const OctreePoint& p2) const { return x >= p2.x && y >= p2.y && z >= p2.z; }
    float length() { return sqrt(x * x + y * y + z * z); }
};
// A node in the quadtree.
struct Node
{
    std::vector<Particle*> particles;
};

#define NUM_PARTICLES 10000000
#define R 0.01f
#define EPSILON 0.0001f
template <class N>
class Octree
{
protected:

    struct OctreeNode
    {
        N _nodeData;
        OctreeNode* _children[8];
        OctreeNode()
        {
            for (int i = 0; i < 8; i++)
                _children[i] = 0;
        }
        virtual ~OctreeNode()
        {
            for (int i = 0; i < 8; i++)
                if (_children[i])
                    delete _children[i];
        }
    };

    iGame::Point _min;
    iGame::Point _max;
    iGame::Point _cellSize;
    OctreeNode* _root;
    VolumeMesh::Pointer mesh{};
public:
    Octree(float min[3], float max[3], float cellSize[3]) : _min(min), _max(max), _cellSize(cellSize), _root(0) {}
    virtual ~Octree() { delete _root; }
    void SetMesh(VolumeMesh::Pointer mesh) {
        this->mesh = mesh;
    };
    void buildOctree() {
        if (mesh == nullptr) { 
            std::cout << "no mesh!" << std::endl;
            return; }
        int numOfPoints = mesh->GetNumberOfPoints();
    }
    class Callback
    {
    public:
        // Return value: true = continue; false = abort.
        virtual bool operator()(const Point min, const Point max, N& nodeData) = 0;
    };

    N& getCell(const int id, Callback* callback = NULL)
    {
        Point pos = mesh->GetPoint(id);
        iGame::Point ppos(pos);
        assert(ppos >= _min && ppos < _max);
        iGame::Point currMin(_min);
        iGame::Point currMax(_max);
        iGame::Point delta = _max - _min;
        if (!_root)
            _root = new OctreeNode();
        OctreeNode* currNode = _root;
        while (delta >= _cellSize)
        {
            bool shouldContinue = true;
            if (callback)
                shouldContinue = callback->operator()(currMin, currMax, currNode->_nodeData);
            if (!shouldContinue)
                break;
            Point mid = (delta * 0.5f) + currMin;
            Point newMin(currMin);
            Point newMax(currMax);
            int index = 0;
            COMPUTE_SIDE(index, 1, ppos.x, mid.x, newMin.x, newMax.x)
            COMPUTE_SIDE(index, 2, ppos.y, mid.y, newMin.y, newMax.y)
            COMPUTE_SIDE(index, 4, ppos.z, mid.z, newMin.z, newMax.z)
                if (!(currNode->_children[index]))
                    currNode->_children[index] = new OctreeNode();
            currNode = currNode->_children[index];
            currMin = newMin;
            currMax = newMax;
            delta = currMax - currMin;
        }
        return currNode->_nodeData;
    }

    void traverse(Callback* callback)
    {
        assert(callback);
        traverseRecursive(callback, _min, _max, _root);
    }

    void clear()
    {
        delete _root;
        _root = NULL;
    }

    class Iterator
    {
    public:
        Iterator getChild(int i)
        {
            return Iterator(_currNode->_children[i]);
        }
        N* getData()
        {
            if (_currNode)
                return &_currNode->_nodeData;
            else return NULL;
        }
    protected:
        OctreeNode* _currNode;
        Iterator(OctreeNode* node) : _currNode(node) {}
        friend class Octree;
    };

    Iterator getIterator()
    {
        return Iterator(_root);
    }

protected:
    void traverseRecursive(Callback* callback, const Point& currMin, const Point& currMax, OctreeNode* currNode)
    {
        if (!currNode)
            return;
        bool shouldContinue = callback->operator()(currMin, currMax, currNode->_nodeData);
        if (!shouldContinue)
            return;
        Point delta = currMax - currMin;
        Point mid = (delta * 0.5f) + currMin;
        traverseRecursive(callback, currMin, mid, currNode->_children[0]);
        traverseRecursive(callback, Point(mid.x, currMin.y, currMin.z),
            Point(currMax.x, mid.y, mid.z), currNode->_children[1]);
        traverseRecursive(callback, Point(currMin.x, mid.y, currMin.z),
            Point(mid.x, currMax.y, mid.z), currNode->_children[2]);
        traverseRecursive(callback, Point(mid.x, mid.y, currMin.z),
            Point(currMax.x, currMax.y, mid.z), currNode->_children[3]);
        traverseRecursive(callback, Point(currMin.x, currMin.y, mid.z),
            Point(mid.x, mid.y, currMax.z), currNode->_children[4]);
        traverseRecursive(callback, Point(mid.x, currMin.y, mid.z),
            Point(currMax.x, mid.y, currMax.z), currNode->_children[5]);
        traverseRecursive(callback, Point(currMin.x, mid.y, mid.z),
            Point(mid.x, currMax.y, currMax.z), currNode->_children[6]);
        traverseRecursive(callback, mid, currMax, currNode->_children[7]);
    }
};
class CallbackTraverse : public Octree<Node>::Callback
{
private:
    VolumeMesh::Pointer mesh{};
public:
    int count;
    Particle* p0;
    void SetMesh(VolumeMesh::Pointer mesh) {
        this->mesh = mesh;
    };
    virtual bool operator()(const Point min, const Point max, Node& n)
    {
        Point pmin(min), pmax(max);
        float cellSize = (pmax - pmin).length();
        float maxDist = ((cellSize) * 0.5f) + R + EPSILON;
        Point temPos0 = mesh->GetPoint(p0->id);
        Point center = (pmin + pmax) * 0.5f;
        Point vectCenter = center - temPos0;
        float distCenter = vectCenter.length();
        if (distCenter > maxDist * maxDist)
            return false; // Too far; don't subdivide cell.

        // Iterate through particles in this cell.
        std::vector<Particle*>::const_iterator it = n.particles.begin();
        for (; it != n.particles.end(); it++)
        {
            Particle* p = *it;
            if (p == p0)
                continue;
            Point temPos = mesh->GetPoint(p0->id);
            float dsq = (temPos - temPos0).length();
            // If particle is within the radius, increment counter.
            if (dsq <= R * R)
                count++;
        }
        // Subdivide cell if needed.
        return true;
    }
};
#endif