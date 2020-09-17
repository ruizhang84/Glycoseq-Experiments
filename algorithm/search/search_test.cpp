#define BOOST_TEST_MODULE ProteinTest
#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <iostream>
#include "search.h"
#include "bucket_search.h"
#include <unordered_map>


using namespace std;
namespace algorithm {
namespace search {

std::shared_ptr<Point<double>> CreatePoint(double value)
{
    return make_shared<Point<double>>(value, value);
}

BOOST_AUTO_TEST_CASE( Algorithm_test ) 
{
    BasicSearch<double> searcher(20.0, ToleranceBy::Dalton);
    std::vector<std::shared_ptr<Point<double>>> box; 

    for(int i=1; i<100; i++)
    {
        box.push_back(CreatePoint(i));
    }

    searcher.set_data(box);
    std::vector<double> res = searcher.Query(40);
    BOOST_CHECK(res.size() == 39);

    BucketSearch<double> bucket_searcher(20.0, ToleranceBy::Dalton);
    bucket_searcher.set_data(box);
    bucket_searcher.Init();
    res = bucket_searcher.Query(40);
    BOOST_CHECK(bucket_searcher.Search(50));
    BOOST_CHECK(res.size() == 39);
}


} // namespace algorithm
} // namespace search 