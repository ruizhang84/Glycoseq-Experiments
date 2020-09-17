#ifndef MODEL_GLYCAN_NGLYCAN_COMPLEX_H
#define MODEL_GLYCAN_NGLYCAN_COMPLEX_H
// The complex type of n-glycan

#include <sstream>
#include <iterator>
#include <typeinfo>
#include "glycan.h"

//GlcNAc(2) - Man(3) - Fuc(1) - GlcNAc(bisect,1) -0,1,2,3
//[GlcNAc(branch1) - GlcNAc(branch2) - GlcNAc(branch3) - GlcNAc(branch4)] -4,5,6,7
//[Gal(branch1) - Gal(branch2) - Gal(branch3) - Gal(branch4)] -8,9,10,11
//[Fuc(branch1) - Fuc(branch2) - Fuc(branch3) - Fuc(branch4)] -12,13,14,15
//[NeuAc(branch1) - NeuAc(branch2) - NeuAc(branch3) - NeuAc(branch4)] -16,17,18,19
//[NeuGc(branch1) - NeuGc(branch2) - NeuGc(branch3) - NeuGc(branch4)] -20,21,22,23

namespace model {
namespace glycan {
class NGlycanComplex : public Glycan 
{
public:
    NGlycanComplex()
    { 
        table_.assign(24, 0);
    }
    ~NGlycanComplex(){}
    
    std::vector<std::unique_ptr<Glycan>> Grow(Monosaccharide suger) override;

    static std::map<Monosaccharide, int> InterpretID(const std::string& table_str)
    {
        std::map<Monosaccharide, int> composite;
        std::vector<int> table;
        std::istringstream iss(table_str);
        std::string item;
        std::vector<std::string> tokens 
        {
            std::istream_iterator<std::string>{iss}, 
            std::istream_iterator<std::string>{}
        };
        for (auto& s : tokens)
        {
            table.push_back(std::stoi(s));
        }
        composite[Monosaccharide::GlcNAc] = 
            table[0] + table[3] + table[4] + table[5] + table[6] + table[7];;
        composite[Monosaccharide::Man] = table[1];
        composite[Monosaccharide::Gal] = table[8] + table[9] + table[10] + table[11];
        composite[Monosaccharide::Fuc] =  table[2] + table[12] + table[13] + table[14] + table[15];
        composite[Monosaccharide::NeuAc] = table[16] + table[17] + table[18] + table[19];
        composite[Monosaccharide::NeuGc] = table[20] + table[21] + table[22] + table[23];
        return composite;
    }


protected:
    void AddMonosaccharide(Monosaccharide suger)
    {
        auto it = composite_.find(suger);
        if (it != composite_.end())
        {
            composite_[suger] += 1;
        }
        else
        {
            composite_[suger] = 1;
        }
    }

    bool ValidAddGlcNAcCore();
    std::unique_ptr<NGlycanComplex> CreateByAddGlcNAcCore();
    bool ValidAddGlcNAc();
    bool ValidAddGlcNAcBisect();
    std::unique_ptr<NGlycanComplex> CreateByAddGlcNAcBisect();
    bool ValidAddGlcNAcBranch();
    std::vector<std::unique_ptr<NGlycanComplex>> CreateByAddGlcNAcBranch();

    bool ValidAddMan();
    std::unique_ptr<NGlycanComplex> CreateByAddMan();

    bool ValidAddGal();
    std::vector<std::unique_ptr<NGlycanComplex>> CreateByAddGal();

    bool ValidAddFucCore();
    std::unique_ptr<NGlycanComplex> CreateByAddFucCore();

    bool ValidAddFucTerminal();
    std::vector<std::unique_ptr<NGlycanComplex>> CreateByAddFucTerminal();

    bool ValidAddNeuAc();
    std::vector<std::unique_ptr<NGlycanComplex>> CreateByAddNeuAc();

    bool ValidAddNeuGc();
    std::vector<std::unique_ptr<NGlycanComplex>> CreateByAddNeuGc();

}; 

}  //  namespace glycan
}  //  namespace model


#endif