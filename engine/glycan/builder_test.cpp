#define BOOST_TEST_MODULE BuilderTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include "glycan_builder.h"


namespace engine{
namespace glycan {

BOOST_AUTO_TEST_CASE( glycan_builder_test ) 
{
    GlycanBuilder builder(2, 3, 1, 0, 0);
    builder.Build();
    // std::cout <<  builder.Isomer().Map().size() << std::endl;
    // for(auto& it :  builder.Isomer().Map()){
    //     std::cout << it.first << std::endl;
    //     for (auto& j: it.second)
    //     {
    //         std::cout << j << std::endl;
    //     }
    //     std::cout << std::endl;
    // }
    BOOST_CHECK(builder.Isomer().Map().size() > 10);
}


BOOST_AUTO_TEST_CASE( nglycan_builder_test ) 
{
    // GlycanBuilder builder2(4, 5, 1, 1, 0);
    // builder2.Build();


    NGlycanBuilder builder(4, 5, 0, 0, 0);
    builder.Build();

    for(auto& it :  builder.Core().Map()){
        std::cout << it.first << std::endl;
        for (auto& j: it.second)
        {
            std::cout << j << std::endl;
        }
        std::cout << std::endl;
    }
    // for(auto& it : builder2.Mass().Map())
    // {
    //     std::unordered_set<double> s = it.second;
    //     std::unordered_set<double> p;
    //     std::cout << it.first<< std::endl;
    //     for(auto j : builder.Core().Map()[it.first])
    //     {
    //         p.insert(j);
    //     }
    //     for(auto j : builder.Branch().Map()[it.first])
    //     {
    //         p.insert(j);
    //     }
    //     for(auto j : builder.Terminal().Map()[it.first])
    //     {
    //         p.insert(j);
    //     }        
        
    //     std::cout << s.size() << " " << p.size() << std::endl;
    //     BOOST_CHECK(s.size() == p.size());
    // }

    // BOOST_CHECK(builder.Isomer().Map().size() == builder2.Isomer().Map().size());
    // BOOST_CHECK(builder.Core().Map().size() == builder2.Mass().Map().size());



}

}
}