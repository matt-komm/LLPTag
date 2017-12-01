#include "XTag/DataFormats/interface/DisplacedGenVertex.h"

namespace xtag
{

double DisplacedGenVertex::d3d() const
{
    if (motherVertex.isNull())
    {
        return -1;
    }
    return std::sqrt((motherVertex->vertex-vertex).mag2());
}

double DisplacedGenVertex::dx() const
{
    if (motherVertex.isNull())
    {
        return -1;
    }
    return std::fabs(motherVertex->vertex.x()-vertex.x());
}

double DisplacedGenVertex::dy() const
{
    if (motherVertex.isNull())
    {
        return -1;
    }
    return std::fabs(motherVertex->vertex.y()-vertex.y());
}

double DisplacedGenVertex::dz() const
{
    if (motherVertex.isNull())
    {
        return -1;
    }
    return std::fabs(motherVertex->vertex.z()-vertex.z());
}

double DisplacedGenVertex::dxy() const
{
    if (motherVertex.isNull())
    {
        return -1;
    }
    const double x = dx();
    const double y = dy();
    return std::sqrt(x*x+y*y);
}

}
