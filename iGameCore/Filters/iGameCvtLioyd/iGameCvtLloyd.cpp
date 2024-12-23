#include <iGameCvtLioyd/iGameCvtLloyd.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include "iGameFilter.h"
#include "iGameVolumeMesh.h"
#include "iGamePriorityQueue.h"


IGAME_NAMESPACE_BEGIN

#include "jc_voronoi.h"
bool CentroidVoronoiDiagram::Execute()
{
    jcv_rect bounding_box = { {0.0f, 0.0f}, {1.0f, 1.0f} };
    std::string density_image = "empty ";
    CentroidVoronoiDiagram cvd(bounding_box, 100, density_image);
    for (int i = 0; i < 1000; i++)
    {
        double max_diff = cvd.lloyd();
        if ((i + 1) % 4 == 0)
            cvd.export_svg("D:\\igames\\pmp_pupa-games102\\data\\bmp\\" + std::to_string(int(i / 4)) + ".svg");
        if (max_diff < 1e-7) break;
        std::cout << "\r" << "lloyd: " << i << ' ' << max_diff << std::endl;
    }
    return true;
}
IGAME_NAMESPACE_END