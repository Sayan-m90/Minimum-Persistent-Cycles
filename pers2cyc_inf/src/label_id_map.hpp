// An efficient map which can return an increment id for an object of any type,
// as long as the type has implemented a hash functor class an equal functor class.
// The map has its own unique copy of each object.

#ifndef _LABEL_ID_MAP_HPP_
#define _LABEL_ID_MAP_HPP_

#include <unordered_set>

struct LabelIdMapElem {
    // the void pointer points to a 'vector<LabelType>' if id is non-negative;
    // it points to a 'LabelType' otherwise.
    void* p;
    int id;
};

template<class LabelType, class LabelHash>
struct LabelIdMapHash { 
    size_t operator()(const LabelIdMapElem& elem) const; 
};

template<class LabelType, class LabelHash>
size_t LabelIdMapHash<LabelType, LabelHash>
    ::operator()(const LabelIdMapElem& elem) const {

    LabelHash label_hash;

    if (elem.id < 0) {
        return label_hash( *((LabelType*)elem.p) );
    } else {
        auto p_label_pool = (vector<LabelType>*)elem.p;
        return label_hash( (*p_label_pool)[elem.id] );
    }
}

template<class LabelType, class LabelEqual>
struct LabelIdMapEqual { 
    bool operator()(const LabelIdMapElem& elem1, const LabelIdMapElem& elem2) const; 
};

template<class LabelType, class LabelEqual>
bool LabelIdMapEqual<LabelType, LabelEqual>
    ::operator()(const LabelIdMapElem& elem1, const LabelIdMapElem& elem2) const {

    LabelType *label1, *label2;

    if (elem1.id < 0) {
        label1 = (LabelType*)elem1.p;
    } else {
        auto p_label_pool = (vector<LabelType>*)elem1.p;
        label1 = &((*p_label_pool)[elem1.id]);
    }

    if (elem2.id < 0) {
        label2 = (LabelType*)elem2.p;
    } else {
        auto p_label_pool = (vector<LabelType>*)elem2.p;
        label2 = &((*p_label_pool)[elem2.id]);
    }

    LabelEqual label_equal;
    return label_equal(*label1, *label2);
}

template <
    class LabelType, class NodeIdType, 
    class LabelHash, class LabelEqual>
class LabelIdMap {

    typedef std::unordered_set<LabelIdMapElem, 
        LabelIdMapHash<LabelType, LabelHash>, 
        LabelIdMapEqual<LabelType, LabelEqual>> LabelMap;

public:
    LabelIdMap() { label_map_ = new LabelMap; }
    ~LabelIdMap() { 
        if (label_map_ != NULL) {
            delete label_map_;
        }
    }

    NodeIdType getId(const LabelType &label);
    void clearLabelMap() { delete label_map_; label_map_ = NULL; }

// private:
public:
    LabelMap* label_map_;
    vector<LabelType> label_pool_;
};

template <
    class LabelType, class NodeIdType, 
    class LabelHash, class LabelEqual>
NodeIdType LabelIdMap<LabelType, NodeIdType, LabelHash, LabelEqual>
    ::getId(const LabelType &label) {

    LabelIdMapElem elem;
    elem.p = (void*)&label;
    elem.id = -1;
    
    auto iter = label_map_->find(elem);
    if (iter == label_map_->end()) {
        NodeIdType new_id = label_pool_.size();
        label_pool_.push_back(label);
        elem.p = (void*)&(this->label_pool_);
        elem.id = new_id;
        label_map_->insert(elem);
        return new_id;
    } else {
        return iter->id;
    }
}

#endif
