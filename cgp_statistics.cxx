#include <iomanip>
#include <map>
#include <thread>

#include <boost/program_options.hpp>

#include <vigra/hdf5impex.hxx>
#include <vigra/compression.hxx>
#include <vigra/timing.hxx>

#include "supervoxels.hxx"
#include "compressors.hxx"

int main(int argc, char** argv) {
    namespace po = boost::program_options;
    using namespace vigra;
    using std::cout; using std::endl; using std::flush; using std::setw;
    
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("geom", po::value<std::string>(),
         "geometry file")
        ("seg", po::value<std::string>(),
         "seg file")
        ("tg", po::value<std::string>(),
         "tg file")
        ("cwx", po::value<std::string>(),
         "cwx file")
        ("maxTgBlocks", po::value<int>(),
         "maximum number of tg blocks considered")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);   
   
    std::string geomFile;
    std::string segFile;
    std::string segGroup;
    std::string tgFile;
    std::string cwxFile;
    int maxTgBlocks = 10;
    
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
    if (vm.count("cwx")) {
        cwxFile = vm["cwx"].as<std::string>();
    } 
    if (vm.count("seg")) {
        std::string s = vm["seg"].as<std::string>();
        auto pos = s.find_last_of("/");
        segFile = s.substr(0,pos);
        segGroup = s.substr(pos+1,s.size());
    } 
    if (vm.count("maxTgBlocks")) {
        maxTgBlocks = vm["maxTgBlocks"].as<int>();
    }
    if (geomFile.empty() && segFile.empty() && tgFile.empty() && cwxFile.empty()) {
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
        gStatistics(geomFile);
    }
    
    if(!tgFile.empty()) {
        
        std::map<vigra::CompressionMethod, std::ofstream> files;
        
        std::map<std::string, vigra::CompressionMethod> cl = compressorList();
        for(const auto& kv : cl) {
            files[kv.second].open("stat_"+toString(kv.second)+".txt", std::ios::trunc);
        }
        
        HDF5File f(tgFile, HDF5File::OpenReadOnly);
        f.cd("blocks");
        auto ls = f.ls();
        std::sort(ls.begin(), ls.end());
        cout << " (" << ls.size() << " blocks) " << endl;
        int n = 0;
        for(const auto& x : ls) {
            if(n >= maxTgBlocks) { break; }
            
            MultiArray<3, uint32_t> tg;
            f.cd(x);
            f.readAndResize("topological-grid", tg);
            f.cd_up();
          
            std::cout << "compressing shape=" << tg.shape() << std::endl;
            auto stats = statCompressors(tg); 
            
            for(const auto& kv : stats) {
                vigra::CompressionMethod cflag = kv.first;
                const CompressionStatistics& stat = kv.second;
                
                files[cflag] /* 0 */ << stat.sizeBytesUncompressed << " "
                             /* 1 */ << stat.sizeBytesUncompressed << " "
                             /* 2 */ << stat.timeCompress << " "
                             /* 3 */ << stat.timeUncompress << " "
                             /* 4 */ << stat.msPerMB_compress() << " "
                             /* 5 */ << stat.msPerMB_uncompress() << " "
                             /* 6 */ << stat.compessionRatio()
                                     << std::endl;
            }
        }
        
#if 0
        std::map<std::string, CompressionStatistics> stats;
        
        std::map<std::string, CompressionMethod> cm = compressorList();
        
        int nthreads = std::thread::hardware_concurrency();
        cout << "nthreads = " << nthreads << endl;
        
        USETICTOC;
        
        MultiArray<3, uint32_t> tg;

        //go over all compression methods available
        for(const auto& kv : cm) {
            const auto cname = kv.first;
            const auto cflag = kv.second;
            cout << "compressing with " << cname << flush;
            
            ArrayVector<char> dest;
            CompressionStatistics stat;
        
            HDF5File f(tgFile, HDF5File::OpenReadOnly);
            f.cd("blocks");
            auto ls = f.ls();
            std::sort(ls.begin(), ls.end());
            cout << " (" << ls.size() << " blocks) " << flush;
            int n = 0;
            for(const auto& x : ls) {
                if(n >= maxTgBlocks) { break; }
                
                cout << "c" << flush;
                f.cd(x);
                f.readAndResize("topological-grid", tg);
                stat.totSize += tg.size()*sizeof(uint32_t);
                
                TIC;
                compress(reinterpret_cast<const char*>(tg.data()), tg.size()*sizeof(uint32_t),
                         dest, cflag, sizeof(uint32_t), nthreads);
                stat.avgTimeCompress += TOCN;
                stat.avgCompression += dest.size();
                
                cout << "u" << flush;
                
                TIC; 
                uncompress(dest.data(), dest.size(),
                           reinterpret_cast<char *>(tg.data()), tg.size()*sizeof(uint32_t),
                           cflag, nthreads);
                stat.avgTimeUncompress += TOCN;
                
                f.cd_up();
                ++n;
            }
            stat.avgTimeCompress   /= ((double)n);
            stat.avgTimeUncompress /= ((double)n);
            stat.avgCompression /= stat.totSize;
            cout << endl;
            cout << "  compress   " << stat.avgTimeCompress << " ms/block" << endl;
            cout << "  uncompress " << stat.avgTimeUncompress << " ms/block" << endl;
            cout << "  ratio      " << stat.avgCompression << endl;
            
            
            stats[cname] = stat;
        }
       
        std::cout << "# method | compress | uncompress | ratio" << std::endl;
        for(const auto& kv : stats) {
            const CompressionStatistics& stat = kv.second;
            std::cout << kv.first
                      << " | " << setw(10) << stat.avgTimeCompress
                      << " | " << setw(10) << stat.avgTimeUncompress
                      << " | " << setw(10) << stat.avgCompression
                      << std::endl;
        }
        
    }
    
    if(!cwxFile.empty()) {
        std::map<std::string, CompressionMethod> cm = compressorList();
        
        cout << "nthreads = " << nthreads << endl;
        
        USETICTOC;
        
        MultiArray<3, uint32_t> tg;

        for(const auto& kv : cm) {
            const auto cname = kv.first;
            const auto cflag = kv.second;
            cout << "compressing with " << cname << flush;
            
            ArrayVector<char> dest;
            double avgTimeCompress   = 0.0;
            double avgTimeUncompress = 0.0;
            double totSize = 0.0;
            double avgCompression = 0.0;
        
            HDF5File f(cwxFile, HDF5File::OpenReadOnly);
            cout << "c" << flush;
            f.readAndResize("cwx", tg);
            totSize += tg.size()*sizeof(uint32_t);
            
            TIC;
            compress(reinterpret_cast<const char*>(tg.data()), tg.size()*sizeof(uint32_t),
                        dest, cflag, sizeof(uint32_t), nthreads);
            avgTimeCompress += TOCN;
            avgCompression += dest.size();
            
            cout << "u" << flush;
            
            TIC; 
            uncompress(dest.data(), dest.size(),
                        reinterpret_cast<char *>(tg.data()), tg.size()*sizeof(uint32_t),
                        cflag, nthreads);
            avgTimeUncompress += TOCN;
            
            f.cd_up();
            
            avgCompression /= totSize;
            cout << endl;
            cout << "  compress   " << avgTimeCompress << " ms/block" << endl;
            cout << "  uncompress " << avgTimeUncompress << " ms/block" << endl;
            cout << "  ratio      " << avgCompression << endl;
        }
#endif
    }
}
