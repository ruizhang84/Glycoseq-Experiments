#include <map>
#include <unordered_map>
#include <fstream>

#include "../../util/io/fasta_reader.h"
#include "../../engine/protein/protein_digest.h"
#include "../../engine/protein/protein_ptm.h"
#include "../../engine/search/search_result.h"
#include "../../engine/score/extra_scorer.h"

// generate peptides by digestion
std::unordered_set<std::string> PeptidesDigestion
    (const std::string& fasta_path, SearchParameter parameter)
{
    util::io::FASTAReader fasta_reader(fasta_path);
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
   
    engine::protein::Digestion digest;
    digest.set_miss_cleavage(parameter.miss_cleavage);
    std::unordered_set<std::string> peptides;

    // digestion
    std::deque<engine::protein::Proteases> proteases(parameter.proteases);
    engine::protein::Proteases enzyme = proteases.front();
    digest.SetProtease(enzyme);
    proteases.pop_front();
    for(const auto& protein : proteins)
    {
        std::unordered_set<std::string> seqs = digest.Sequences(protein.Sequence(),
            engine::protein::ProteinPTM::ContainsNGlycanSite);
        peptides.insert(seqs.begin(), seqs.end());
    }
        
    // double digestion or more
    while (proteases.size() > 0)
    {
        std::unordered_set<std::string> double_seqs;
        enzyme = proteases.front();
        digest.SetProtease(enzyme);
        proteases.pop_front();
        for(const auto& seq : peptides)
        {
            std::unordered_set<std::string> seqs = digest.Sequences(seq,
                engine::protein::ProteinPTM::ContainsNGlycanSite);
            double_seqs.insert(seqs.begin(), seqs.end());
        }
        peptides.insert(double_seqs.begin(), double_seqs.end());
    }   

    return peptides;
}

// assign score to searching results
void ScoringWorker(std::vector<engine::search::SearchResult>& results)
{
    engine::score::ExtraScorer scorer;
    scorer.UpdateScore(results);
}

// report glycopeptide identification of spectrum
void ReportResults(const std::string& out_path,
    const std::vector<engine::search::SearchResult>&  results)
{
    std::ofstream outfile;
    outfile.open (out_path);
    outfile << "scan#,peptide,glycan,score\n";

    for(auto it : results)
    {
        outfile << it.Scan() << ",";
        outfile << it.Sequence() << ",";
        outfile << it.Glycan() << ",";
        outfile << it.RawScore() << "\n";
    }
    outfile.close();
}
