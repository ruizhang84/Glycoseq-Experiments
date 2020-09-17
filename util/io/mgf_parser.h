#ifndef UTIL_IO_MGF_PARSER_H_
#define UTIL_IO_MGF_PARSER_H_

#include <string>
#include <map> 
#include <fstream>
#include <regex>
#include "spectrum_reader.h"

namespace util {
namespace io {

class MGFParser : public SpectrumParser
{   
public:
    MGFParser(std::string path, SpectrumType type): 
        type_(type){ path_ = path; }
        
    void Init() override
    {
        MGFData data;
        int scan_num = -1;

        std::ifstream file(path_);
        std::string line;

        std::smatch result;
        std::regex start("BEGIN\\s+IONS");
        std::regex end("END\\s+IONS");
        std::regex title("TITLE=(.*)");
        std::regex pepmass("PEPMASS=(\\d+\\.?\\d*)");
        std::regex charge("CHARGE=(\\d+)");
        std::regex rt_second("RTINSECONDS=(\\d+\\.?\\d*)");
        std::regex scan("SCANS=(\\d+)");
        std::regex mz_intensity("^(\\d+\\.?\\d*)\\s+(\\d+\\.?\\d*)");

        if (file.is_open()){
            while(std::getline(file, line)){
                if (std::regex_search(line, result, start))
                {
                    data = MGFData();
                    scan_num++;
                }else if (std::regex_search(line, result, mz_intensity))
                {
                    data.mz.push_back(std::stod(result[1]));
                    data.intensity.push_back(std::stod(result[2]));
                }
                else if (std::regex_search(line, result, pepmass))
                {
                    data.pep_mass = std::stod(result[1]);
                }
                else if (std::regex_search(line, result, charge)){
                    data.charge = std::stoi(result[1]);
                }
                else if (std::regex_search(line, result, scan))
                {
                    scan_num = std::stoi(result[1]);
                    data.scans = scan_num;
                }
                else if (std::regex_search(line, result, title))
                {
                    data.title = std::string(result[1]);
                }
                else if (std::regex_search(line, result, rt_second))
                {
                    data.rt_seconds = std::stod(result[1]);
                }
                else if (std::regex_search(line, result, end))
                {
                    data_set_.emplace(scan_num, data);
                } 
            }
        }
    }

    double ParentMZ(int scan_num) override 
    { 
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.pep_mass;
        }
        return 0;
    }
    int ParentCharge(int scan_num) override
    { 
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.charge;
        }
        return 0;
    }
    int GetFirstScan() override 
    {
        auto it = data_set_.begin();
        if (it != data_set_.end())
        {
            return it->first;
        }
        return -1;
    }
    int GetLastScan() override
    {
        auto it = data_set_.rbegin();
        if (it != data_set_.rend())
        {
            return it->first;
        }
        return -1;
    }
    std::vector<Peak> Peaks(int scan_num) override
    { 
        std::vector<Peak> peaks;
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            for(size_t i = 0; 
                i < it->second.intensity.size(); i++)
            {
                Peak pk(it->second.mz[i], 
                        it->second.intensity[i]);
                peaks.push_back(pk);
            }
        }
        return peaks;
    }
    SpectrumType GetSpectrumType(int scan_num) override 
        { return type_; };
    double RTFromScanNum(int scan_num) override
    {
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.rt_seconds;
        }
        return -1;
    }
    std::string GetScanInfo(int scan_num) override
    {
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.title;
        }
        return "";
    }
    bool Exist(int scan_num) override
    {
        return data_set_.find(scan_num) != data_set_.end();
    }
    
private:
    class MGFData
    {
    public:
        MGFData() = default;

        std::vector<double> mz;
        std::vector<double> intensity;
        double pep_mass;
        int charge;
        double rt_seconds;
        int scans;
        std::string title;
    };
    SpectrumType type_;
    std::map<int, MGFData> data_set_;
};


} // namespace io
} // namespace util


#endif