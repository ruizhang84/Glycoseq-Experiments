#ifndef ALGORITHM_BASE_UNION_FIND_H
#define ALGORITHM_BASE_UNION_FIND_H

#include <vector>
#include <unordered_map> 

namespace algorithm {
namespace base {

class UnionFind
{
public: 
    UnionFind() = default;

    void Clear()
    {
        map_.clear();
        size_.clear();
        rank_.clear();
    }

    int Find(int i) 
    { 
        if (map_.find(i) == map_.end())
        {
            map_[i] = i;
            size_[i] = 1;
            rank_[i] = 0;
        }
        return (map_[i] == i) ? i : (map_[i] = Find(map_[i])); 
    }
    bool IsSameSet(int i, int j) { return Find(i) == Find(j); }
    virtual void Union(int i, int j) 
    { 
        if (!IsSameSet(i, j)) 
        {
            int x = Find(i), y = Find(j);
            // rank is used to keep the tree short
            if (rank_[x] > rank_[y]) 
            { 
                map_[y] = x; 
                size_[x] += size_[y]; 
            }
            else                   
            { 
                map_[x] = y; 
                size_[y] += size_[x];
                if (rank_[x] == rank_[y]) 
                rank_[y]++; 
            } 
        } 
    }

protected:
    std::unordered_map<int, int> map_; 
    std::unordered_map<int, int> size_; 
    std::unordered_map<int, int> rank_; 
};

} // namespace base
} // namespace algorithm


#endif