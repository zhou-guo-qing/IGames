#include "ImageSampler.h"
#include "iGameFilter.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

#include "iGameFilter.h"
#include "iGamePointSet.h"
#include "iGameSurfaceMesh.h"
#include "iGamePoints.h"
#include "iGameObject.h"
#include "iGameQuad.h"
#include "../CellModel/Quadratic/iGameQuadraticTriangle.h"
#include <random>
#define JCV_REAL_TYPE double
#define JC_VORONOI_IMPLEMENTATION


IGAME_NAMESPACE_BEGIN

#include "jc_voronoi.h"
class CentroidVoronoiDiagram : public Filter
{
public:

    I_OBJECT(CentroidVoronoiDiagram);

    static Pointer New() { return new CentroidVoronoiDiagram ; }

    explicit CentroidVoronoiDiagram(jcv_rect bbox, size_t n_sample, std::string density_image = "")
        : density_(density_image), bbox_(bbox)
    {
        if (density_.is_constant()) points_.resize(n_sample);
        else
            points_.resize(128 * 128 * density_.height() / density_.width());
        srand(0);
        for (size_t i = 0; i < points_.size(); i++)
        {
            points_[i].x = (float)(rand() / (1.0f + RAND_MAX));
            points_[i].y = (float)(rand() / (1.0f + RAND_MAX));
        }
    }
    bool Execute()override;

    ~CentroidVoronoiDiagram() { if (diagram_.internal) jcv_diagram_free(&diagram_); }

    jcv_point polygon_centroid(const jcv_graphedge* graph_edge);

    jcv_point varying_polygon_centroid(const jcv_graphedge* graph_edge, size_t n_sample = 30);

    double relax_points()
    {
        const jcv_site* sites = jcv_diagram_get_sites(&diagram_);
        double max_diff = 1e-9;
        for (int i = 0; i < diagram_.numsites; ++i)
        {
            const jcv_site* site = &sites[i];
            const jcv_point pre_p = points_[site->index];
            if (!density_.is_constant())
                points_[site->index] = varying_polygon_centroid(site->edges);
            else
                points_[site->index] = polygon_centroid(site->edges);

            // max_diff = std::max(max_diff, jcv_point_dist_sq(&pre_p, &points_[site->index]));
        }
        return max_diff;
    }

    double lloyd()
    {
        if (diagram_.internal)
            jcv_diagram_free(&diagram_);
        memset(&diagram_, 0, sizeof(jcv_diagram));
        jcv_diagram_generate(points_.size(), points_.data(), &bbox_, 0, &diagram_);
        return relax_points();
    }

    std::string export_svg(std::string file);
protected:
    SurfaceMesh::Pointer mesh{};
    CentroidVoronoiDiagram() {
        mesh = nullptr;
        this->SetNumberOfOutputs(1);
        this->SetNumberOfInputs(1);
    };
    CentroidVoronoiDiagram(SurfaceMesh::Pointer _mesh) {
        this->mesh = _mesh;
    }
private:
    jcv_diagram            diagram_{ 0 };
    std::vector<jcv_point> points_;
    ImageSampler           density_;
    jcv_rect               bbox_;
};

jcv_point CentroidVoronoiDiagram::polygon_centroid(const jcv_graphedge* graph_edge) {
    double total_det = 0;
    jcv_point center{ 0, 0 };
    for (; graph_edge; graph_edge = graph_edge->next)
    {
        jcv_point p1 = graph_edge->pos[0], p2 = graph_edge->pos[1];
        double det = p1.x * p2.y - p2.x * p1.y;
        total_det += det;
        center.x += (p1.x + p2.x) * det;
        center.y += (p1.y + p2.y) * det;
    }
    center.x /= 3 * total_det;
    center.y /= 3 * total_det;
    return center;
}

jcv_point CentroidVoronoiDiagram::varying_polygon_centroid(const jcv_graphedge* graph_edge,
    size_t n_sample)
{
    jcv_point pc = polygon_centroid(graph_edge);
    jcv_point center{ 0, 0 };
    double W = 0;
    for (; graph_edge; graph_edge = graph_edge->next)
    {
        jcv_point p1 = graph_edge->pos[0], p2 = graph_edge->pos[1];
        Eigen::RowVector3d centroid = density_.centric(p1.x, p1.y, p2.x, p2.y, pc.x, pc.y);
        center.x += centroid.x();
        center.y += centroid.y();
        W += centroid.z();
    }
    center.x /= W;
    center.y /= W;
    return center;
}

std::string  CentroidVoronoiDiagram::export_svg(std::string output_file) {
    std::stringstream svg_str;
    double w = 512;
    double h = density_.is_constant() ? 512 : density_.height() / double(density_.width()) * 512;
    svg_str << "<svg width=\"" << w << "\" height=\"" << h << "\" viewBox=\"0 0 " << w << ' ' << h
        << R"(" xmlns="http://www.w3.org/2000/svg" >)"
        << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    const jcv_site* sites = jcv_diagram_get_sites(&diagram_);
    for (size_t i = 0; i < diagram_.numsites && density_.is_constant(); i++)
    {
        svg_str << "<polygon points=\"";
        jcv_graphedge* graph_edge = sites[i].edges;
        while (graph_edge)
        {
            jcv_point p1 = graph_edge->pos[0], p2 = graph_edge->pos[1];
            svg_str << p1.x * w << ',' << p1.y * h << ' ' << p2.x * w << ',' << p2.y * h << ' ';
            graph_edge = graph_edge->next;
        }
        svg_str << R"(" fill="#1d1d9b" stroke="white" stroke-width="2" />)" << "\n";
    }

    for (int i = 0; i < diagram_.numsites; ++i)
    {
        const jcv_site* site = &sites[i];
        jcv_point p = site->p, pc = polygon_centroid(site->edges);
        if (density_.is_constant()) {
            svg_str << "<circle cx=\"" << pc.x * w << "\" cy=\"" << pc.y * h << R"(" r="1" fill="red"/>)" << '\n';
            svg_str << "<circle cx=\"" << p.x * w << "\" cy=\"" << p.y * h << R"(" r="1" fill="yellow"/>)" << '\n';
        }
        else
            svg_str << "<circle cx=\"" << p.x * w << "\" cy=\"" << p.y * h << R"(" r="1" fill="black"/>)" << '\n';
    }

    svg_str << "</svg>";

    if (output_file.size() != 0)
    {
        std::ofstream svg_file(output_file);
        svg_file << svg_str.str();
        svg_file.close();
    }


    return svg_str.str();
}
IGAME_NAMESPACE_END