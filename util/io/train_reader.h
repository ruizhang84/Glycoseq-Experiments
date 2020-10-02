#ifndef UTIL_IO_TRAIN_READER_H_
#define UTIL_IO_TRAIN_READER_H_

#include <string>
#include <map> 
#include <vector>
#include <fstream>
#include <algorithm>

namespace util {
namespace io {

class TrainReader
{   
public:
    TrainReader(std::string path): path_(path) {}
    
    std::string Path(){ return path_; }
    void set_path(std::string path) { path_ = path; }

    void Init(bool header=true)
    {

        std::ifstream file(path_);
        std::string line;
        data_.clear();
        
        if (file.is_open()){
            while(std::getline(file, line)){
                if (header)
                {
                    header = false;
                    continue;
                }
                std::vector<int> scans = GetScans(line);
                std::sort(scans.begin(), scans.end());
                data_[GetFile(line)] = scans;
            }
        }
    }

    std::map<std::string, std::vector<int>> Dataset()
        { return data_;}

protected:
    std::string GetFile(const std::string& line)
    {
        std::string file_name;
        for(const char& c : line)
        {
            if (c == ',') break;
            file_name += c;
        }
        return file_name;
    }

    std::vector<int> GetScans(const std::string& line)
    {
        std::vector<int> scans;
        int idx = 0;
        for(; idx < (int) line.size(); idx++)
        {
            if (line[idx] == ',') break;
        }

        std::string temp = "";
        for(idx++; idx < (int) line.size(); idx++)
        {
            if (line[idx] == ' ')
            {
                if (temp.size() > 0)
                    scans.push_back(std::stoi(temp));
                temp.clear();
            }
            else
            {
                temp += line[idx];
            }
        }
        if (temp.size() > 0)
        {
            scans.push_back(std::stoi(temp));
        }

        return scans;
    }

    
private:
    std::string path_;
    std::map<std::string, std::vector<int>> data_;
};


} // namespace io
} // namespace util


#endif