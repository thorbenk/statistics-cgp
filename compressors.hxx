#ifndef COMPRESSORS_HXX
#define COMPRESSORS_HXX

#include <vector>

#include <vigra/multi_array.hxx>

struct CompressionStatistics {
    CompressionStatistics(vigra::CompressionMethod m)
      : timeCompress(0)
      , timeUncompress(0)
      , sizeBytesUncompressed(0)
      , sizeBytesCompressed(0)
      , compressor(m) {}
      
    CompressionStatistics()
      : timeCompress(0)
      , timeUncompress(0)
      , sizeBytesUncompressed(0)
      , sizeBytesCompressed(0)
      {}
        
    double compessionRatio() const {
        return sizeBytesCompressed / sizeBytesUncompressed;
    }
    double msPerMB_compress() const {
        return timeCompress / (sizeBytesUncompressed/(1024*1024));
    }
    double msPerMB_uncompress() const {
        return timeUncompress / (sizeBytesUncompressed/(1024*1024));
    }
        
    double timeCompress;
    double timeUncompress;
    double sizeBytesUncompressed;
    double sizeBytesCompressed;
    vigra::CompressionMethod compressor;
};

std::map<std::string, vigra::CompressionMethod> compressorList();

std::string toString(const vigra::CompressionMethod m);

typedef std::map<vigra::CompressionMethod, CompressionStatistics> Stats;

Stats statCompressors(
    const vigra::MultiArrayView<3, uint32_t>& a,
    bool verbose = true
);

#endif /* COMPRESSORS_HXX */
