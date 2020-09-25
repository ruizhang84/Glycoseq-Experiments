#define BOOST_TEST_MODULE SearchTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <random>
#include "neural_network.h"

namespace engine{
namespace learn{

BOOST_AUTO_TEST_CASE( neural_network_test ) 
{
    // data 3x1+2x2+1;
    const int nrolls=100;  // number of experiments

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0,0.1);

    std::vector<std::vector<double>> X;
    std::vector<double> Y;

    for (int i=0; i<nrolls; ++i) {
        for (int j=0; j < nrolls; ++j) 
        {
            double number = distribution(generator);
            std::vector<double> x;
            x.push_back(i);
            x.push_back(j);
            X.push_back(x);
            double y = 3 * i + 2 * j + number;
            Y.push_back(y);
        }
    }

    // neural network setup
    std::vector<double> weight; 
    double bias = 1.0;
    for (int i = 0; i < 2; i++)
    {
        weight.push_back(0.0);
    }

    Classifier classifer(weight, bias);
    classifer.set_debug(true);
    classifer.Train(X, Y, 100, 0.1);

    std::cout<<"Weight: ";
    for(const auto& it : classifer.Weight())
    {
        std::cout << " " << it << " ";
    }
    std::cout << std::endl;
    std::cout<<"Bias: " << classifer.Bias() << std::endl;
}

} // namespace engine
} // namespace learn