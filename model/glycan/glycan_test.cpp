#define BOOST_TEST_MODULE GlycanTest
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "moiety.h"
#include "nglycan_complex.h"

using namespace std;
using namespace model::glycan;
std::unique_ptr<Moiety> CreateRoot()
{
    std::unique_ptr<Moiety> a = 
        std::make_unique<Moiety>(Monosaccharide::GlcNAc);
    std::unique_ptr<Moiety> b = 
        std::make_unique<Moiety>(Monosaccharide::Gal);
    std::unique_ptr<Moiety> c = 
        std::make_unique<Moiety>(Monosaccharide::Man);
    std::unique_ptr<Moiety> d = 
        std::make_unique<Moiety>(Monosaccharide::Fuc);

    c->Children().push_back(std::move(d));
    b->Children().push_back(std::move(c));
    a->Children().push_back(std::move(b));
    
    std::unique_ptr<Moiety> e = a->Clone();
    return e;
}

BOOST_AUTO_TEST_CASE( Monosaccharide_test ) 
{

    std::unique_ptr<Moiety> e = CreateRoot();
    BOOST_CHECK( e->Name() ==  Monosaccharide::GlcNAc); 
    BOOST_CHECK( e->Children().front()->Name() ==  Monosaccharide::Gal); 
    BOOST_CHECK( e->Children().front()->Children().front()
                    ->Children().front()->Name() ==  Monosaccharide::Fuc); 
    BOOST_CHECK( e->Children().front()->Children().front()
                    ->Parent()->Name() ==  Monosaccharide::Gal); 
}

BOOST_AUTO_TEST_CASE( Glycan_test ) 
{
    NGlycanComplex nglycan;
    nglycan.set_table(0, 2); // core 
    nglycan.set_table(1, 3); // 
    nglycan.set_table(2, 1); // fuc
    nglycan.set_table(3, 1); // bisect
    nglycan.set_table(4, 1); // branch 
    nglycan.set_table(5, 1);
    nglycan.set_table(8, 1);
    nglycan.set_table(9, 1);  

    std::vector<std::unique_ptr<Glycan>> glycans = nglycan.Grow(Monosaccharide::GlcNAc);
    BOOST_CHECK(glycans.size() == 2);
    for (int i = 0; i < (int) glycans.size(); i++){
        std::cout << glycans[i]->Name() << std::endl;
        std::cout << glycans[i]->ID() <<std::endl;
    }

    std::string table_str = nglycan.Serialize();
    NGlycanComplex nglycan_dup;
    nglycan_dup.Deserialize(table_str);
    std::cout << table_str <<std::endl;
    BOOST_CHECK(table_str == nglycan_dup.Serialize());

    NGlycanComplex nglycan_1;
    std::map<Monosaccharide, int> composite;
    composite[Monosaccharide::GlcNAc] = 12;
    composite[Monosaccharide::Gal] = 12;
    composite[Monosaccharide::Fuc] = 1;
    composite[Monosaccharide::Man] = 3;
    composite[Monosaccharide::NeuAc] = 6;

    nglycan_1.set_composition(composite);
    std::string compos = nglycan_1.Name();
    NGlycanComplex nglycan_2;
    nglycan_2.set_composition(compos);

    std::cout << compos << std::endl;
    std::cout << nglycan_2.Name() << std::endl;
    BOOST_CHECK(compos == nglycan_2.Name());

}

BOOST_AUTO_TEST_CASE( glycan_row_test ) 
{
    NGlycanComplex nglycan;
    std::vector<std::unique_ptr<Glycan>> glycans = nglycan.Grow(Monosaccharide::GlcNAc);
    std::cout << glycans.front()->ID() << std::endl;
    std::vector<std::unique_ptr<Glycan>> glycans_2 = glycans.front()->Grow(Monosaccharide::GlcNAc);
    std::cout << glycans_2.front()->ID() << std::endl;
    std::vector<std::unique_ptr<Glycan>> glycans_3 = glycans_2.front()->Grow(Monosaccharide::Man);
    std::string name = glycans_3.front()->Name();
    std::cout << name << std::endl;
    BOOST_CHECK(glycans_3.front()->CompositionConst()[Monosaccharide::GlcNAc] == 2);
}

// int add( int i, int j ) { return i+j; }

// BOOST_AUTO_TEST_CASE( my_test )
// {
//     // seven ways to detect and report the same error:
//     BOOST_CHECK( add( 2,2 ) == 4 );        // #1 continues on error

//     BOOST_REQUIRE( add( 2,2 ) == 4 );      // #2 throws on error

//     if( add( 2,2 ) != 4 )
//       BOOST_ERROR( "Ouch..." );            // #3 continues on error

//     if( add( 2,2 ) != 4 )
//       BOOST_FAIL( "Ouch..." );             // #4 throws on error

//     if( add( 2,2 ) != 4 ) throw "Ouch..."; // #5 throws on error

//     BOOST_CHECK_MESSAGE( add( 2,2 ) == 4,  // #6 continues on error
//                          "add(..) result: " << add( 2,2 ) );

//     BOOST_CHECK_EQUAL( add( 2,2 ), 4 );	  // #7 continues on error
// }








