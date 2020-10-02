#define BOOST_TEST_MODULE MGFParserTest
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "train_reader.h"

namespace util {
namespace io {

BOOST_AUTO_TEST_CASE( read_training_data ) 
{
    TrainReader reader("/home/ruiz/Documents/Glycoseq-Experiments/data/identified_scan.csv");
    reader.Init();

    std::map<std::string, std::vector<int>> data = reader.Dataset();
    for(const auto& it : data)
    {
        std::cout << it.first << std::endl;
        for(const auto& scan: it.second)
        {
            std::cout << scan << "\t";
        }
        std::cout << std::endl;
    }

}


} // namespace io
} // namespace util
