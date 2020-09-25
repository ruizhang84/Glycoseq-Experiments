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
    const int nrolls=10000;  // number of experiments

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(1.0,0.1);

    std::vector<std::vector<double>> X;
    std::vector<int> Y;

    for (int i=0; i<nrolls; ++i) {
        double a = distribution(generator);
        double b = distribution(generator);
        double c = distribution(generator);
        int y = 6 * c > (3 * a + 2 * b );

        std::vector<double> x;
        x.push_back(a);
        x.push_back(b);
        x.push_back(c);
        X.push_back(x);
        Y.push_back(y);
    }
    // training and testing 
    std::vector<int> Y_train = std::vector<int>(Y.begin(), Y.begin() + 7000);
    std::vector<int> Y_test = std::vector<int>(Y.begin()+7000, Y.end());

    // neural network setup
    Classifier classifer;
    classifer.set_debug(true);
    classifer.Train(X, Y_train, 10000, 0.1);

    // accuracy
    double accu = classifer.Test(X, Y_test);
    std::cout <<"Accuracy: " << accu << std::endl;

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