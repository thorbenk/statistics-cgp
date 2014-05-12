#include <map>

#include <boost/program_options.hpp>

#include <vigra/accumulator.hxx>
#include <vigra/hdf5impex.hxx>
#include <vigra/compression.hxx>
#include <vigra/timing.hxx>

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


int main(int argc, char** argv) {
    namespace po = boost::program_options;
    using namespace vigra;
    using std::cout; using std::endl; using std::flush;
    
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("geom", po::value<std::string>(),
         "geometry file")
        ("seg", po::value<std::string>(),
         "seg file")
        ("tg", po::value<std::string>(),
         "tg file")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);   
   
    std::string geomFile;
    std::string segFile;
    std::string segGroup;
    std::string tgFile;
    if (vm.count("help")) {
        cout << desc << endl;
        return 1;
    }
    if (vm.count("geom")) {
        geomFile = vm["geom"].as<std::string>();
    } 
    if (vm.count("tg")) {
        tgFile = vm["tg"].as<std::string>();
    } 
    if (vm.count("seg")) {
        std::string s = vm["seg"].as<std::string>();
        auto pos = s.find_last_of("/");
        segFile = s.substr(0,pos);
        segGroup = s.substr(pos+1,s.size());
    } 
    if (geomFile.empty() && segFile.empty() && tgFile.empty()) {
        cout << "Error: Need at least one of --geom and --seg options!" << endl << endl;
        cout << desc << endl;
        return 1;
    }
    
    
    if(!segFile.empty()) {
        std::cout << "* Loading segmentation from " << segFile << "/" << segGroup << std::endl;
        HDF5File f(segFile, HDF5File::OpenReadOnly);
        MultiArray<3, uint32_t> seg;
        f.readAndResize(segGroup, seg);
        f.close();
        supervoxelStatistics(seg);
    }
    
    if(!geomFile.empty()) {
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
    
    if(!tgFile.empty()) {
        std::map<std::string, CompressionMethod> cm;
        cm["DEFAULT_COMPRESSION"] = DEFAULT_COMPRESSION;
        cm["NO_COMPRESSION     "] = NO_COMPRESSION;
        cm["ZLIB_NONE          "] = ZLIB_NONE;
        cm["ZLIB_FAST          "] = ZLIB_FAST;
        cm["ZLIB               "] = ZLIB;
        cm["ZLIB_BEST          "] = ZLIB_BEST;
        cm["LZ4                "] = LZ4;
        
        USETICTOC;
        
        MultiArray<3, uint32_t> tg;

        for(const auto& kv : cm) {
            const auto cname = kv.first;
            const auto cflag = kv.second;
            cout << "compressing with " << cname << " ... " << flush;
            
            ArrayVector<char> dest;
            double avg = 0.0;
        
            //void compress(char const * source, std::size_t size, ArrayVector<char> & dest, CompressionMethod method);
            HDF5File f(tgFile, HDF5File::OpenReadOnly);
            f.cd("blocks");
            auto ls = f.ls();
            for(const auto& x : ls) {
                f.cd(x);
                f.readAndResize("topological-grid", tg);
                TIC;
                compress(reinterpret_cast<const char*>(tg.data()), tg.size()*sizeof(uint32_t), dest, cflag);
                avg += TOCN;
                f.cd_up();
            }
            cout << " (" << ls.size() << " blocks)" << " took on average " << avg << endl;
        }
    }
}
