#ifndef SUPERVOXELS_HXX
#define SUPERVOXELS_HXX

#include <string>
#include <vigra/multi_array.hxx>

void gStatistics(const std::string& geomFile);
void supervoxelStatistics(const vigra::MultiArrayView<3, uint32_t>& seg);

#endif /* SUPERVOXELS_HXX */