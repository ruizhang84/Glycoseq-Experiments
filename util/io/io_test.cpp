#define BOOST_TEST_MODULE MGFParserTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>

#include "mgf_parser.h"
#include "fasta_reader.h"

namespace util {
namespace io {

BOOST_AUTO_TEST_CASE( mgf_read_test ) 
{
    MGFParser parser("/home/yu/Documents/MultiGlycan-Cpp/data/test_EThcD_min.mgf", 
                    SpectrumType::EThcD);
    parser.Init();
    BOOST_CHECK( parser.GetFirstScan() == 3); 
    BOOST_CHECK( parser.ParentMZ(64) == 1289.262); 
    BOOST_CHECK( parser.ParentCharge(64) == 2); 
    BOOST_CHECK( parser.RTFromScanNum(64) == 24.422337857); 
    BOOST_CHECK( parser.GetScanInfo(64) == "C:\\Users\\iruiz\\Desktop\\app3\\ZC_20171218_H68_R1.raw"); 

    Peak pk = parser.Peaks(64).front();
    BOOST_CHECK( pk.MZ() == 113.3392); 
    BOOST_CHECK( pk.Intensity() == 238.3); 
}

BOOST_AUTO_TEST_CASE( fasta_read_test ) 
{
    FASTAReader fasta_reader("/home/yu/Documents/MultiGlycan-Cpp/data/test_fasta.fasta");
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
    BOOST_CHECK( proteins.front().Sequence() == 
        "MSALGAVIALLLWGQLFAVDSGNDVTDIADDGCPKPPEIAHGYVEHSVRYQCKNYYKLRTEGDGVYTLND"); 
}

} // namespace io
} // namespace util


