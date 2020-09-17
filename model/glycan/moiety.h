#ifndef MODEL_GLYCAN_MOIETY_H
#define MODEL_GLYCAN_MOIETY_H

#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <unordered_map> 
#include "glycan.h"

namespace model {
namespace glycan {

class Moiety
{
public:
    Moiety(Monosaccharide name): name_(name){}
    Moiety(const Moiety& other)
    {
        name_ = other.name_;
        parent_ = other.parent_;
        for(auto& child: other.children_)
        {
            std::unique_ptr<Moiety> new_child = child->Clone();
            new_child->set_parent(this);
            children_.push_back(std::move(new_child));
        }
    }

    Monosaccharide Name() const { return name_; }
    Moiety* Parent() { return parent_; }
    void set_parent(Moiety* parent) { parent_ = parent; }
    std::vector<std::unique_ptr<Moiety>>& Children() { return children_; } 
    std::vector<Monosaccharide> ChildrenName() const
    {
        std::vector<Monosaccharide> names;
        for(auto& child : children_)
        {
            names.push_back(child->Name());
        }
        return names;
    }
    std::vector<Monosaccharide> ParentChildrenName() const
    {
        return parent_->ChildrenName();
    }

    std::unique_ptr<Moiety> Clone()
    {
        std::unique_ptr<Moiety> root_ = 
            std::make_unique<Moiety>(name_);

        for(auto& child: children_)
        {
            std::unique_ptr<Moiety> new_child = child->Clone();
            new_child->set_parent(root_.get());
            root_->Children().push_back(std::move(new_child));
            
        }
        return root_;
    }

protected:
    Monosaccharide name_;
    Moiety* parent_;
    std::vector<std::unique_ptr<Moiety>> children_;
};


}  //  namespace glycan
}  //  namespace model

#endif

