#include "supervoxels.hxx"

#include <vigra/accumulator.hxx>

#include <cgp/GeometryReader.hxx>

typedef unsigned int label_type;
typedef short coordinate_type;
typedef cgp::hdf5::GeometryReader<label_type, coordinate_type> GeometryReader;

void statistic(const GeometryReader& g, int dimension) {
    std::vector<size_t> sizes(g.maxLabel(dimension));
    for(size_t i=1; i<g.maxLabel(dimension); ++i) {
        sizes[i-1] = g.size(dimension, i);
    }
    std::cout << "average size for dim=" << dimension << " : "
              << std::accumulate(sizes.begin(), sizes.end(), 0.0) / ((double)sizes.size())
              << std::endl;
}

void supervoxelStatistics(const vigra::MultiArrayView<3, uint32_t>& seg) {
    using namespace vigra;
    using namespace vigra::acc;
    AccumulatorChainArray<CoupledArrays<3, uint32_t, uint32_t>, 
                                        Select<DataArg<1>, LabelArg<2>,
                                        RegionCenter> > a;
    auto start = createCoupledIterator(seg, seg);
    auto end = start.getEndIterator();
    extractFeatures(start, end, a);
    double avg = 0.0;
    for(size_t i=1; i<=a.maxRegionLabel(); ++i) {
        avg += get<Count>(a,i);
    }
    avg /= ((double)a.maxRegionLabel());
    std::cout << "avg supervoxel size: " << avg << std::endl;
}

void gStatistics(const std::string& geomFile) {
    GeometryReader g(geomFile, GeometryReader::EnableZeroSetBoundsQuery |
        GeometryReader::EnableOneSetBoundsQuery |
        GeometryReader::EnableTwoSetBoundsQuery |
        GeometryReader::EnableOneSetBoundedByQuery |
        GeometryReader::EnableTwoSetBoundedByQuery |
        GeometryReader::EnableThreeSetBoundedByQuery);
    for(int d=1; d<3; ++d) {
        statistic(g, d);
    }
}