#define BOOST_TEST_MODULE LSHTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <chrono> 

#include "spectrum_sim.h"
#include "calc.h"
#include "../../model/spectrum/spectrum.h"
#include "../../util/io/spectrum_reader.h"
#include "../../util/io/mgf_parser.h"

namespace util {
namespace calc {

using namespace model::spectrum;
using namespace std;
using namespace util::io;
using namespace std::chrono; 

BOOST_AUTO_TEST_CASE( calc_test ) 
{
    // Calc calculator;
    // std::vector<double> v1 {1, 2, 3, 4};
    // std::vector<double> v2 {2, 4, 6, 8};
   
    // cout << calculator.DotProduct(v1, v2) << endl; 
}

BOOST_AUTO_TEST_CASE( sim_test ) 
{
    string path = "/home/yu/Documents/MultiGlycan-Cpp/data/test_EThcD.mgf";
    std::unique_ptr<SpectrumParser> parser = 
        std::make_unique<MGFParser>(path, SpectrumType::EThcD);
    SpectrumReader spectrum_reader(path, std::move(parser));
    spectrum_reader.Init();

    Spectrum s1 = spectrum_reader.GetSpectrum(3);
    Spectrum s2 = spectrum_reader.GetSpectrum(64);

    auto start = high_resolution_clock::now(); 
    SpectrumSim sim;
    double cos = sim.ComputeCosine(s1, s2);
    std::cout << cos << std::endl;
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<microseconds>(stop - start); 
    std::cout << duration.count() << std::endl; 
}

} // namespace io
} // namespace calc


