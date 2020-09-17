#ifndef UTIL_IO_PROTEIN_READER_H_
#define UTIL_IO_PROTEIN_READER_H_

#include <string>
#include <vector>
#include "../../model/protein/protein.h"

namespace util {
namespace io {

class ProteinReader
{
public:
    ProteinReader(std::string path): path_(path){}
    virtual ~ProteinReader(){}
    
    virtual std::vector<model::protein::Protein> Read()
    {
         std::vector<model::protein::Protein> proteins;
         return proteins;
    }

    std::string Path() const { return path_; }
    void set_path(std::string path) { path_ = path; }

protected:
    std::string trim(const std::string& str)
    {
        size_t first = str.find_first_not_of(" \t\n\r\f\v");
        if (std::string::npos == first)
        {
            return str;
        }
        size_t last = str.find_last_not_of(" \t\n\r\f\v");
        return str.substr(first, (last - first + 1));
    }

    std::string path_;

};

} // namespace io
} // namespace util


#endif