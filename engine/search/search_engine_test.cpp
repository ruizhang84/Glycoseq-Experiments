#define BOOST_TEST_MODULE SearchTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <iomanip>
#include "spectrum_search.h"
#include "../../util/io/mgf_parser.h"
#include "../../util/io/fasta_reader.h"
#include "../protein/protein_digest.h"
#include "../protein/protein_ptm.h"
#include "../glycan/glycan_builder.h"
#include "../spectrum/normalize.h"
#include <chrono> 

namespace engine{
namespace search {

BOOST_AUTO_TEST_CASE( search_engine_test ) 
{
    // read spectrum
    std::string path = "/home/ruiz/Documents/Glycoseq-Experiments/data/ZC_20171218_C22_R1.mgf";
    std::unique_ptr<util::io::SpectrumParser> parser = 
        std::make_unique<util::io::MGFParser>(path, util::io::SpectrumType::EThcD);
    util::io::SpectrumReader spectrum_reader(path, std::move(parser));
    spectrum_reader.Init();
    int start_scan = spectrum_reader.GetFirstScan();
    int last_scan = spectrum_reader.GetLastScan();
    BOOST_CHECK(start_scan < last_scan);

    // process spectrum by normalization
    model::spectrum::Spectrum spec = spectrum_reader.GetSpectrum(start_scan);
    engine::spectrum::Normalizer::Transform(spec);
    BOOST_CHECK(spec.Peaks().front().Intensity() < 1);

    // read fasta and build peptides
    util::io::FASTAReader fasta_reader("/home/ruiz/Documents/Glycoseq-Experiments/data/haptoglobin.fasta");
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
 
    engine::protein::Digestion digest;
    digest.SetProtease(engine::protein::Proteases::Trypsin);
    std::unordered_set<std::string> seqs = digest.Sequences(proteins.front().Sequence(),
         engine::protein::ProteinPTM::ContainsNGlycanSite);
    digest.SetProtease(engine::protein::Proteases::GluC);
    std::vector<std::string> peptides;
    for(auto& it : seqs)
    {
        std::unordered_set<std::string> seq = digest.Sequences(it,
         engine::protein::ProteinPTM::ContainsNGlycanSite);
        peptides.insert(peptides.end(), seq.begin(), seq.end());
    }
    BOOST_CHECK(std::find(peptides.begin(), peptides.end(), "MVSHHNLTTGATLINE") != peptides.end());



    // // build glycans
    int hexNAc = 12, hex = 12, Fuc = 5, NeuAc = 4, NeuGc = 0;
    std::unique_ptr<engine::glycan::NGlycanBuilder> builder =
        std::make_unique<engine::glycan::NGlycanBuilder>(hexNAc, hex, Fuc, NeuAc, NeuGc);
    builder->Build();

    model::glycan::NGlycanComplex glycan;
    std::map<model::glycan::Monosaccharide, int> composite;
    composite[model::glycan::Monosaccharide::GlcNAc] = 5;
    composite[model::glycan::Monosaccharide::Man] = 3;
    composite[model::glycan::Monosaccharide::Gal] = 3;
    composite[model::glycan::Monosaccharide::Fuc] = 1;
    composite[model::glycan::Monosaccharide::NeuAc] = 1;  
    glycan.set_composition(composite);
    std::string glycan_name = glycan.Name();
    BOOST_CHECK(builder->Isomer().QueryMass(glycan_name) == util::mass::GlycanMass::Compute(composite));
    BOOST_CHECK(util::mass::GlycanMass::Compute(composite) == 
        util::mass::GlycanMass::Compute(model::glycan::Glycan::Interpret(glycan_name)));
    std::vector<std::string> collection = builder->Isomer().Collection();
    BOOST_CHECK(std::find(collection.begin(), collection.end(), glycan.Name()) != collection.end());


    // spectrum matching
    int special_scan = 7678; // 2742; // 8729, 8778
    double ms1_tol = 10, ms2_tol = 0.01;
    int isotopic_count = 2;
    algorithm::search::ToleranceBy ms1_by = algorithm::search::ToleranceBy::PPM;
    algorithm::search::ToleranceBy ms2_by = algorithm::search::ToleranceBy::Dalton;

    PrecursorMatcher precursor_runner(ms1_tol, ms1_by, builder->Isomer());
    std::vector<std::string> glycans_str = builder->Isomer().Collection();
    precursor_runner.Init(peptides, glycans_str);

    SpectrumSearcher spectrum_runner(ms2_tol, ms2_by, 2, builder.get(), true);
    spectrum_runner.Init();

    auto special_spec = spectrum_reader.GetSpectrum(special_scan);
    double special_target = util::mass::SpectrumMass::Compute(special_spec.PrecursorMZ(), special_spec.PrecursorCharge());
    MatchResultStore special_r = precursor_runner.Match(special_target, special_spec.PrecursorCharge(), isotopic_count);    
    std::cout << special_spec.Scan() << " : " << std::endl;
    special_r.Add("NLFLNHSE", "GlcNAc-4-Man-3-Gal-2-NeuAc-2-");
    for(auto it : special_r.Map())
    {
        std::cout << it.first << std::endl;
        for(auto g: it.second)
        {
            std::cout << g << std::endl;
        }
    }
    BOOST_CHECK(!special_r.Empty());
    // GlcNAc-5-Man-3-Gal-3-Fuc-2-NeuAc-2-
    // MVSHHNLTTGATLINE
    // 202896

    

    std::cout << "scan start:\n"; 
    auto start = std::chrono::high_resolution_clock::now(); 
    special_spec = spectrum_reader.GetSpectrum(special_scan);
    engine::spectrum::Normalizer::Transform(special_spec);
    spectrum_runner.set_candidate(special_r);
    spectrum_runner.set_spectrum(special_spec);
    std::vector<SearchResult> special_res = spectrum_runner.Search();

    for (const auto& it : special_res)
    {
        std::cout << it.Sequence() << std::endl;
        std::cout << it.Glycan() << std::endl;
        std::cout << it.RawScore() << std::endl;
    }

    // compute score
    BOOST_CHECK(!special_res.empty());
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

}

} // namespace search
} // namespace engine