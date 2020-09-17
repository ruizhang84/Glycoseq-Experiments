#ifndef ENGINE_SCORE_EXTRA_SCORER_H
#define ENGINE_SCORE_EXTRA_SCORER_H

#include "coelution.h"
#include "../search/search_result.h"
#include "../../util/io/spectrum_reader.h"

namespace engine{
namespace score {

class ExtraScorer
{
public:
    ExtraScorer() = default;
    void UpdateScore(std::vector<engine::search::SearchResult>& results)
    {
        UpdateElutionScore(results);
        for(auto& it : results)
        {
            double score = it.Score() *
                it.ExtraScore(engine::search::ScoreType::Elution);
            it.set_score(score);
        }
    }

protected:
    void UpdateElutionScore(std::vector<engine::search::SearchResult>& results)
    {
        elutor_.Update(results);
    }
    CoElution elutor_;
};


}   // namespace score 
}   // namespace engine



#endif