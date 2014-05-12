#include <boost/program_options.hpp>

#include <vigra/accumulator.hxx>
#include <vigra/hdf5impex.hxx>

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

namespace po = boost::program_options;

int main(int argc, char** argv) {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("geom", po::value<std::string>(),
         "geometry file")
        ("seg", po::value<std::string>(),
         "seg file")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);   
   
    std::string geomFile;
    std::string segFile;
    if (vm.count("help")) {
        cout << desc << endl;
        return 1;
    }
    if (vm.count("geom")) {
        geomFile = vm["geom"].as<std::string>();
    } 
    if (vm.count("seg")) {
        segFile = vm["seg"].as<std::string>();
    } 
    
    using namespace vigra;
    using namespace vigra::acc;
    
    HDF5File f(segFile, HDF5File::OpenReadOnly);
    MultiArray<3, uint32_t> seg;
    f.readAndResize("seg", seg);
    f.close();
    
    AccumulatorChainArray<CoupledArrays<3, uint32_t, uint32_t>, 
                                        Select<DataArg<1>, LabelArg<2>,
                                        RegionCenter> > a;
    auto start = createCoupledIterator(seg, seg);
    auto end = start.getEndIterator();
    extractFeatures(start, end, a);
    std::cout << get<Count>(a,1);
   
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
