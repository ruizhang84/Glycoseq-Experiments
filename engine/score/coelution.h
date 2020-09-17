#ifndef ENGINE_SCORE_COELUTION_H
#define ENGINE_SCORE_COELUTION_H

#include <climits> 
#include "../search/search_result.h"
#include "../../util/io/spectrum_reader.h"

namespace engine{
namespace score {

class CoElution
{
public:
    CoElution() = default;

    void Update(std::vector<engine::search::SearchResult>& results)
    {
        // compute coelution of peptide sequence
        std::unordered_map<std::string, int> total;
        int start = INT_MAX,  end = 0;
        for(const auto& r : results)
        {
            int scan = r.Scan();
            start = std::min(start, scan);
            end = std::max(end, scan);
            if (total.find(r.Sequence()) == total.end())
            {
                total[r.Sequence()] = 0;
            }
            total[r.Sequence()] += 1;
        }

        // count peptide sequence in each range bucket
        std::vector<std::unordered_map<std::string, int>> range_count;
        range_count.assign(kRange + 1, std::unordered_map<std::string, int>());
        for(const auto& r : results)
        {
            int scan = r.Scan();
            int index = Index(scan, start, end);
            std::string peptide = r.Sequence();

            std::unordered_map<std::string, int>& count = range_count[index];
            
            if (count.find(peptide) == count.end())
            {
                count[peptide] = 0;
            }
            count[peptide] += 1;
        }

        // udpate score
        for(auto& it : results)
        {
            
            int scan = it.Scan();
            int index = Index(scan, start, end);
            std::string s = it.Sequence();
            int counts = range_count[index][s];
            // find neighbor counts
            if (index > 0)
            {
                std::unordered_map<std::string, int>& count = range_count[index-1];
                if (count.find(s) != count.end())
                    counts += count[s];
            }

            std::unordered_map<std::string, int>& count = range_count[index+1];
            if (count.find(s) != count.end())
                counts += count[s];
 
 
            double score = counts * 1.0 / total[s];
            score = score > kLimit ? 1.0 : score;
            it.set_extra(score, engine::search::ScoreType::Elution);
        }
    }


protected:
    const int kRange = 30;
    const int kLimit = 0.8;
    int Index(int scan, int start, int end)
    {
        return (scan - start) * kRange / (end + 1- start);
    }
};


}   // namespace score 
}   // namespace engine



#endif