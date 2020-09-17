#ifndef ENGINE_GLYCAN_GLYCAN_STORE_H
#define ENGINE_GLYCAN_GLYCAN_STORE_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace engine {
namespace glycan {

typedef std::unordered_map<std::string, 
        std::unordered_set<std::string>> StringsMapping;

typedef std::unordered_map<std::string, 
        std::unordered_set<double>> DoublesMapping;

class GlycanStore
{
public:
    StringsMapping Map() const { return map_; }
    std::unordered_map<std::string, double> Mass() const { return mass_; }
    std::unordered_set<std::string> Query(const std::string item) const
    {
        if (map_.find(item) != map_.end())
        {
           return Map()[item];
        }
        std::unordered_set<std::string> result;
        return result;
    }
    double QueryMass(const std::string item) const
    {
        double mass = 0;
        if (map_.find(item) != map_.end())
        {
           return Mass()[item];
        }
        return mass;
    }
    std::vector<std::string> Collection() const
    {
        std::vector<std::string> collection;
        for(const auto& it : map_)
        {
            collection.push_back(it.first);
        }
        return collection;
    }
    bool Contains(const std::string item) const
    {
        return map_.find(item) != map_.end();
    }
    void Add(const std::string& name, const std::string& table_id)
    {
        if (map_.find(name) == map_.end())
        {
            map_[name] = std::unordered_set<std::string>();
        }
        map_[name].insert(table_id);
    }
    void Add(const std::string& name, const double mass)
        { mass_[name] = mass; }

    void AddSubset(const std::string& table_id, const std::string& subset_id)
    {
        Add(table_id, subset_id);
        if (map_.find(subset_id) != map_.end())
        {
            std::unordered_set<std::string> subset =  map_[subset_id];
            map_[table_id].insert(subset.begin(), subset.end());
        }
    }
    void Clear(){ map_.clear(); }

protected:
    // glycan composition_str(name) -> table_str(id) or mass, by isomer
    StringsMapping map_;
    std::unordered_map<std::string, double> mass_;
};

class GlycanMassStore
{
public:
    DoublesMapping Map() const
        { return map_; }

    std::unordered_set<double> Query(const std::string item) const
    {
        if (map_.find(item) != map_.end())
        {
           return Map()[item];
        }
        std::unordered_set<double> result;
        return result;
    }
    bool Contains(const std::string item) const
    {
        return map_.find(item) != map_.end();
    }
    void Add(const std::string& name, const double mass)
    {
        if (map_.find(name) == map_.end())
        {
            map_[name] = std::unordered_set<double>();
        }
        map_[name].insert(mass);
    }
    void AddSubset(const std::string& name, 
        const std::string& subset_id, const double mass)
    {
        if (mass > 0)
            Add(name, mass);
        if (map_.find(subset_id) != map_.end())
        {
            std::unordered_set<double> subset =  map_[subset_id];
            map_[name].insert(subset.begin(), subset.end());
        }
    }
    void Clear(){ map_.clear(); }

protected:
    // glycan_id -> mass of its subset, by biosynthesis
    DoublesMapping map_;
};


} // namespace glycan
} // namespace engine



#endif
