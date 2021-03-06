#ifndef ENGINE_ANALYSIS_FDR_PROB_FILTER_H
#define ENGINE_ANALYSIS_FDR_PROB_FILTER_H

#include <algorithm>
#include <numeric>
#include <boost/math/distributions/normal.hpp>
#include  "../search/search_result.h"

using boost::math::normal_distribution;

namespace engine{
namespace analysis{

class FDRProbFilter
{
public:
    FDRProbFilter(double fdr): fdr_(fdr){}

    std::vector<engine::search::SearchResult> Filter(
        const std::vector<engine::search::SearchResult>& targets,
        const std::vector<engine::search::SearchResult>& decoys)
    {
        std::vector<engine::search::SearchResult> results;

         // get score lists of each spectrum
        std::map<int, std::vector<double>> v;
        for(const auto& it : decoys)
        {
            int i = it.Scan();
            if (v.find(i) == v.end())
            {
                v[i] = std::vector<double>();
            }
            v[i].push_back(it.RawScore());
        }

        // compute p value
        std::map<int, engine::search::SearchResult> s_results;
        std::map<int, double> p_v;
        for(const auto& it : targets)
        {
            if (v.find(it.Scan()) == v.end()) // no decoys find
            {
                results.push_back(it);
            }
            else
            {
                std::vector<double>& score_list = v[it.Scan()];
                double avg = mean(score_list);
                if (avg > it.RawScore()) continue;

                double p = pValue(score_list, it.RawScore());

                p_v[it.Scan()] = p;
                s_results[it.Scan()] = it;
            }
        }

        //Bonferroni
        int size = (int) p_v.size();
        for(const auto& it : p_v)
        {
            double p = it.second;
            if (p * size < fdr_)
            {
                results.push_back(s_results[it.first]);
            }
        }
        if (!results.empty())
            std::sort(results.begin(), results.end(), OrderByScan);
        return results;
    }

protected:
    double fdr_;

    static double mean(std::vector<double> v)
    {
        int size = (int) v.size();
        double sum = 0;
        for(const auto& it : v)
        {
            sum += it;
        }
        return sum * 1.0 / size;
    }

    static double stdv(std::vector<double> v)
    {
        
        double m = 0;
        double avg = mean(v);
        int size = (int) v.size();
        for(const auto& it : v)
        {
            m += (it - avg) * (it - avg);
        }
        return std::sqrt( m * 1.0 / size);
    }

    static double pValue(std::vector<double> v, double q)
    {
        double m = mean(v);
        double s = stdv(v);
        return 0.5 * erfc( (q-m) * 1.0 / (s * std::sqrt(2)) );
    }
    
    static bool OrderByScan(const engine::search::SearchResult& r1, const engine::search::SearchResult& r2)
        { return r1.Scan() < r2.Scan(); }
};

}
}

#endif