#include <map>
#include <thread>
#include <iomanip>

#include <vigra/compression.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/timing.hxx>

#include "compressors.hxx"

std::map<vigra::CompressionMethod, CompressionStatistics> stats;

std::string toString(const vigra::CompressionMethod m) {
    using namespace vigra;
    switch(m) {
        case NO_COMPRESSION:
        return "NO_COMPRESSION";
        case ZLIB_NONE:
        return "ZLIB_NONE";
        case ZLIB_FAST:
        return "ZLIB_FAST";
        case ZLIB:
        return "ZLIB";
        case ZLIB_BEST:
        return "ZLIB_BEST";
        case LZ4:
        return "LZ4";
        case BLOSC_BLOSCLZ_FAST:
        return "BLOSC_BLOSCLZ_FAST";
        case BLOSC_LZ4_FAST:
        return "BLOSC_LZ4_FAST";
        case BLOSC_LZ4HC_FAST:   
        return "BLOSC_LZ4HC_FAST"; 
        case BLOSC_SNAPPY_FAST:   
        return "BLOSC_SNAPPY_FAST";   
        case BLOSC_ZLIB_FAST:
        return "BLOSC_ZLIB_FAST";
        case BLOSC_BLOSCLZ_BEST:  
        return "BLOSC_BLOSCLZ_BEST";  
        case BLOSC_LZ4_BEST:
        return "BLOSC_LZ4_BEST";
        case BLOSC_LZ4HC_BEST:   
        return "BLOSC_LZ4HC_BEST";  
        case BLOSC_SNAPPY_BEST:  
        return "BLOSC_SNAPPY_BEST";  
        case BLOSC_ZLIB_BEST: 
        return "BLOSC_ZLIB_BEST";
    }
}

std::map<std::string, vigra::CompressionMethod> compressorList() {
    using namespace vigra;
    std::map<std::string, CompressionMethod> cm;
    cm["NO_COMPRESSION     "] = NO_COMPRESSION;
    cm["ZLIB_NONE          "] = ZLIB_NONE;
    cm["ZLIB_FAST          "] = ZLIB_FAST;
    cm["ZLIB               "] = ZLIB;
    cm["ZLIB_BEST          "] = ZLIB_BEST;
    cm["LZ4                "] = LZ4;
    cm["BLOSC_BLOSCLZ_FAST "] = BLOSC_BLOSCLZ_FAST;
    cm["BLOSC_LZ4_FAST     "] = BLOSC_LZ4_FAST; 
    cm["BLOSC_LZ4HC_FAST   "] = BLOSC_LZ4HC_FAST;   
    cm["BLOSC_SNAPPY_FAST  "] = BLOSC_SNAPPY_FAST;   
    cm["BLOSC_ZLIB_FAST    "] = BLOSC_ZLIB_FAST;  
    cm["BLOSC_BLOSCLZ_BEST "] = BLOSC_BLOSCLZ_BEST;  
    cm["BLOSC_LZ4_BEST     "] = BLOSC_LZ4_BEST; 
    cm["BLOSC_LZ4HC_BEST   "] = BLOSC_LZ4HC_BEST;   
    cm["BLOSC_SNAPPY_BEST  "] = BLOSC_SNAPPY_BEST;  
    cm["BLOSC_ZLIB_BEST    "] = BLOSC_ZLIB_BEST; 
    return cm;
}

Stats statCompressors(
    const vigra::MultiArrayView<3, uint32_t>& a
) {
    using std::cout; using std::endl; using std::flush; using std::setw;
    using namespace vigra;
    USETICTOC;
    
    std::map<std::string, CompressionMethod> cm = compressorList();
    int nthreads = std::thread::hardware_concurrency();
    
    Stats stats; 
    
    for(const auto& kv : cm) {
        const auto cname = kv.first;
        const auto cflag = kv.second;
        cout << "compressing with " << cname << flush;
        
        ArrayVector<char> dest;
    
        CompressionStatistics stat;
        
        cout << "c" << flush;
        stat.sizeBytesUncompressed = a.size()*sizeof(uint32_t);
        
        TIC;
        compress(reinterpret_cast<const char*>(a.data()), a.size()*sizeof(uint32_t),
                    dest, cflag, sizeof(uint32_t), nthreads);
        stat.timeCompress = TOCN;
        stat.sizeBytesCompressed = dest.size();
        
        cout << "u" << flush;
        
        TIC; 
        uncompress(dest.data(), dest.size(),
                    reinterpret_cast<char *>(a.data()), a.size()*sizeof(uint32_t),
                    cflag, nthreads);
        stat.timeUncompress += TOCN;
        cout << endl;
        cout << "  compress   " << stat.msPerMB_compress()   << " MB/ms" << endl;
        cout << "  uncompress " << stat.msPerMB_uncompress() << " MB/ms" << endl;
        cout << "  ratio      " << stat.compessionRatio()    << endl;
        
        stats[cflag] = stat;
    }
    
    /*
    for(const auto& kv : stats) {
        const CompressionStatistics& s = kv.second;
        cout << kv.first;
        cout << setw(10) << s.compessionRatio() << " | "
             << setw(10) << s.msPerMB_compress() << " | "
             << setw(10) << s.msPerMB_uncompress() << std::endl;
    }
    */
    
    return stats;
}