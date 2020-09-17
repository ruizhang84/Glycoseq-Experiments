#ifndef MODEL_GLYCAN_GLYCAN_H
#define MODEL_GLYCAN_GLYCAN_H

#include <string>
#include <vector>
#include <map> 
#include <memory>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <regex>
#include <iostream>

namespace model {
namespace glycan {

enum class Monosaccharide
{ GlcNAc, Man, Gal, Fuc, NeuAc, NeuGc};

class Glycan
{
public:
    Glycan() = default;
    virtual ~Glycan(){}
    
    std::string Name() const 
    { 
        std::string name = "";
        for (const auto& it : composite_)
        {
            switch (it.first)
            {
            case Monosaccharide::GlcNAc:
                name += "GlcNAc-" + std::to_string(it.second) + "-";
                break;
            case Monosaccharide::Man:
                name += "Man-" + std::to_string(it.second) + "-";
                break;
            case Monosaccharide::Gal:
                name += "Gal-" + std::to_string(it.second) + "-";
                break;
            case Monosaccharide::Fuc:
                name += "Fuc-" + std::to_string(it.second) + "-";
                break;    
            case Monosaccharide::NeuAc:
                name += "NeuAc-" + std::to_string(it.second) + "-";
                break;
            case Monosaccharide::NeuGc:
                name += "NeuGc-" + std::to_string(it.second) + "-";
                break;        
            default:
                break;
            }
        }
        return name;
    } // for print
    std::string ID() const { return Serialize(); }  // use as key

    void set_name(const std::string& name) 
        { name_ = name; }
    void set_id(const std::string& id) 
        { id_ = id; }

    std::vector<int>& Table() { return table_; }
    void set_table(const std::vector<int>& table) 
        { table_ = table; }
    void set_table(int index, int num)
    {
        if (index >= 0 && index < (int) table_.size())
            table_[index] = num;
    }

    std::string Serialize() const
    {
        std::stringstream result;
        std::copy(table_.begin(), table_.end(), 
            std::ostream_iterator<int>(result, " "));
        return result.str();
    }

    void Deserialize(std::string table_str)
    {
        std::istringstream iss(table_str);
        std::string item;
        std::vector<std::string> tokens 
        {
            std::istream_iterator<std::string>{iss}, 
            std::istream_iterator<std::string>{}
        };
        table_.clear();
        for (auto& s : tokens)
        {
            table_.push_back(std::stoi(s));
        }
    }

    std::map<Monosaccharide, int>&  Composition()
        { return composite_; }
    void set_composition(const std::map<Monosaccharide, int>& composite)
        { composite_ = composite; }
    
    static std::map<Monosaccharide, int> Interpret(const std::string& name)
    {
        std::map<Monosaccharide, int> composite;
        std::smatch result;
        std::regex rGlcNAc("GlcNAc-(\\d+)-");
        std::regex rGal("Gal-(\\d+)-");
        std::regex rMan("Man-(\\d+)-");
        std::regex rFuc("Fuc-(\\d+)-");
        std::regex rNeuAc("NeuAc-(\\d+)-");
        std::regex rNeuGc("NeuGc-(\\d+)-");

        if (std::regex_search(name.begin(), name.end(), result, rGlcNAc))
        {
            composite[Monosaccharide::GlcNAc] = std::stoi(result[1]);
        }

        if (std::regex_search(name.begin(), name.end(), result, rGal))
        {
            composite[Monosaccharide::Gal] = std::stoi(result[1]);
        }

        if (std::regex_search(name.begin(), name.end(), result, rMan))
        {
            composite[Monosaccharide::Man] = std::stoi(result[1]);
        }

        if (std::regex_search(name.begin(), name.end(), result, rFuc))
        {
            composite[Monosaccharide::Fuc] = std::stoi(result[1]);
        }

        if (std::regex_search(name.begin(), name.end(), result, rNeuAc))
        {
            composite[Monosaccharide::NeuAc] = std::stoi(result[1]);
        }

        if (std::regex_search(name.begin(), name.end(), result, rNeuGc))
        {
            composite[Monosaccharide::NeuGc] = std::stoi(result[1]);
        }
        return composite;
    }
    void set_composition(const std::string& name)
    {
       set_composition(Interpret(name));
    }

    const std::map<Monosaccharide, int>&  CompositionConst() const
        { return composite_; }

    virtual std::vector<std::unique_ptr<Glycan>> Grow(Monosaccharide suger)
    {
        std::vector<std::unique_ptr<Glycan>> result;
        return result;
    }

protected:
    std::string name_;
    std::string id_;
    std::vector<int> table_;
    std::map<Monosaccharide, int> composite_; 

};


}  //  namespace glycan
}  //  namespace model

#endif

