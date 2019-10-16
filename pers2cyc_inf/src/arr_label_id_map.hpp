#ifndef _ARR_LABEL_ID_MAP_HPP_
#define _ARR_LABEL_ID_MAP_HPP_

#include <unordered_set>
#include <boost/functional/hash.hpp>

// IdType should support negative values, i.e., it should be 
// 'int' rather than 'unsigned int'.

template <class IdType>
struct ArrLabelIdMapElem {
    // the void pointer points to a 'vector<ElemType>' if id is non-negative;
    // it points to a 'ElemType' otherwise.
    void* p;
    IdType id;
};

template <class ElemType, class IdType>
class ArrLabelIdMapHash { 
public:
    ArrLabelIdMapHash(int length) : length_(length) {}

    size_t operator()(const ArrLabelIdMapElem<IdType>& elem) const; 

private:
    int length_;
};

template <class ElemType, class IdType>
size_t ArrLabelIdMapHash<ElemType, IdType>
    ::operator()(const ArrLabelIdMapElem<IdType>& elem) const {

    // cout << "ArrLabelIdMapHash: " << this->length_ << endl;

    std::size_t seed = 0;

    for (auto i = 0; i < this->length_; i ++) {
        if (elem.id < 0) {
            boost::hash_combine( seed, ((ElemType*)elem.p)[i] );
        } else {
            auto pool = (vector<ElemType>*)elem.p;
            boost::hash_combine( seed, (*pool)[elem.id * this->length_ + i] );
        }
    }

    return seed;
}

template <class ElemType, class IdType>
class ArrLabelIdMapEqual { 
public:
    ArrLabelIdMapEqual(int length) : length_(length) {}

    bool operator()(const ArrLabelIdMapElem<IdType>& elem1, 
        const ArrLabelIdMapElem<IdType>& elem2) const; 

private:
    int length_;
};

template <class ElemType, class IdType>
bool ArrLabelIdMapEqual<ElemType, IdType>
    ::operator()(const ArrLabelIdMapElem<IdType>& elem1, 
        const ArrLabelIdMapElem<IdType>& elem2) const {

    // cout << "ArrLabelIdMapEqual: " << this->length_ << endl;

    for (auto i = 0; i < this->length_; i ++) {
        ElemType *label1, *label2;

        if (elem1.id < 0) {
            label1 = (ElemType*)elem1.p + i;
        } else {
            auto pool = (vector<ElemType>*)elem1.p;
            label1 = &( (*pool)[elem1.id * this->length_ + i] );
        }

        if (elem2.id < 0) {
            label2 = (ElemType*)elem2.p + i;
        } else {
            auto pool = (vector<ElemType>*)elem2.p;
            label2 = &( (*pool)[elem2.id * this->length_ + i] );
        }

        if (*label1 != *label2) {
            return false;
        }
    }

    return true;
}

template <class ElemType, class IdType>
class ArrayLabelIdMap {

    typedef std::unordered_set< ArrLabelIdMapElem<IdType>, 
        ArrLabelIdMapHash<ElemType, IdType>, 
        ArrLabelIdMapEqual<ElemType, IdType> > ArrayLabelMap;

public:
    ArrayLabelIdMap() { 
        arr_label_map_ = NULL; 
        elem_cnt_ = 0;
    }

    ~ArrayLabelIdMap() { 
        clearLabelMap();
    }

    void init(const int arr_label_length, bool need_query_map = true) {
        arr_label_length_ = arr_label_length;

        if (need_query_map) {
            ArrLabelIdMapHash<ElemType, IdType> hash(arr_label_length);
            ArrLabelIdMapEqual<ElemType, IdType> equal(arr_label_length);
            arr_label_map_ = new ArrayLabelMap(10, hash, equal);
        }
    }

    void reserve(int size) {
        if (arr_label_map_ != NULL) {
            arr_label_map_->reserve(size);
        }
        arr_label_pool_.reserve(size * arr_label_length_);
    }

    // the size is the total number of array labels the pool used to hold,
    // though some may be deleted and marked null
    int size() { return arr_label_pool_.size() / arr_label_length_; }

    int elemCount() { return elem_cnt_; }

    IdType getId(const ElemType* arr);

    IdType getId(const vector<ElemType>& vec) {
        if (vec.size() != arr_label_length_) { return -1; }
        return this->getId(vec.data());
    }

    // no duplication check is done, invoker should self-check
    IdType addArrayLabel(const ElemType* arr);

    IdType addArrayLabel(const vector<ElemType>& vec) {
        if (vec.size() != arr_label_length_) { return -1; }
        return this->addArrayLabel(vec.data());
    }

    void getArrayLabel(const IdType id, vector<ElemType>* vec);

    void deleteArrayLabel(const IdType id, ElemType null_val);

    void clearLabelMap() { 
        if (arr_label_map_ != NULL) {
            delete arr_label_map_; 
            arr_label_map_ = NULL; 
        }
    }

    void printLoadStat();

public:

    ArrayLabelMap* arr_label_map_;

    vector<ElemType> arr_label_pool_;

    int arr_label_length_;

    int elem_cnt_;
};

template <class ElemType, class IdType>
IdType ArrayLabelIdMap<ElemType, IdType>
    ::getId(const ElemType* arr) {
    
    ArrLabelIdMapElem<IdType> elem;
    elem.p = (void*)arr;
    elem.id = -1;
    
    auto iter = arr_label_map_->find(elem);
    if (iter == arr_label_map_->end()) {
        return -1;
    } else {
        return iter->id;
    }
}

template <class ElemType, class IdType>
IdType ArrayLabelIdMap<ElemType, IdType>
    ::addArrayLabel(const ElemType* arr) {

    IdType new_id = arr_label_pool_.size() / arr_label_length_;
    for (auto i = 0; i < arr_label_length_; i ++) {
        arr_label_pool_.push_back(arr[i]);
    }
    
    if (arr_label_map_ != NULL) {
        ArrLabelIdMapElem<IdType> elem;
        elem.p = (void*)&arr_label_pool_;
        elem.id = new_id;
        arr_label_map_->insert(elem);
    }

    elem_cnt_ ++;

    return new_id;
}

template <class ElemType, class IdType>
void ArrayLabelIdMap<ElemType, IdType>
    ::getArrayLabel(const IdType id, vector<ElemType>* vec) {

    vec->clear();
    vec->reserve(arr_label_length_);

    for (auto i = 0; i < arr_label_length_; i ++) {
        vec->push_back(arr_label_pool_[id * arr_label_length_ + i]);
    }
}

template <class ElemType, class IdType>
void ArrayLabelIdMap<ElemType, IdType>
    ::deleteArrayLabel(const IdType id, ElemType null_val) {

    if (arr_label_map_ != NULL) {
        ArrLabelIdMapElem<IdType> elem;
        elem.p = (void*)&arr_label_pool_;
        elem.id = id;

        auto iter = arr_label_map_->find(elem);
        if (iter != arr_label_map_->end()) {
            arr_label_map_->erase(iter);

            arr_label_pool_[id * arr_label_length_] = null_val;
            elem_cnt_ --;
        }
    }
}

template <class ElemType, class IdType>
void ArrayLabelIdMap<ElemType, IdType>
    ::printLoadStat() {

    cout << "label_pool_: pool_capacity=" << arr_label_pool_.capacity() 
        << ",pool_size=" << arr_label_pool_.size() 
        << ",size=" << arr_label_pool_.size() / arr_label_length_ << endl;

    if (arr_label_map_ != NULL) {
        cout << "label_map_: ";
        printHashMapLoad(*arr_label_map_);
    } else {
        cout << "no arr_label_map_" << endl;
    }
    cout << "elem_cnt_=" << elem_cnt_ << endl;
}

// template <class ElemType, class LenType>
// class WrapperArray {
// public:
//     WrapperArray(const ElemType* data, const LenType length) :
//        data_(data), length_(length) {}

//     // WrapperArray(const WrapperArray<ElemType, LenType> &arr) :
//     //    data_(arr.data_), length_(arr.length_) {
//        // cout << "WrapperArray constructor: " << arr << ", " << *this << endl;
//     // }

//     const ElemType* p() const { return data_; }
//     void wrap(const ElemType* data, const LenType length);

//     const ElemType operator[](const LenType i) const { return data_[i]; }
//     LenType length() const { return length_; }

// private:
//     const ElemType* data_;
//     LenType length_;
// };

// template <class ElemType, class LenType>
// void WrapperArray<ElemType, LenType>::wrap(const ElemType* data, const LenType length) {
//     data_ = data;
//     length_ = length;
// }

// template <class ElemType, class LenType>
// ostream& operator<<(ostream& os, const WrapperArray<ElemType,LenType>& a) {
//     if (a.length() == 0) {
//         os << "()";
//         return os;
//     }

//     os << '(' << a[0];
//     for (auto i = 1; i < a.length(); i ++) {
//         os << ',' << a[i];
//     }
//     os << ')';

//     return os;
// }

// struct IntWrapArrayHash { 
//     size_t operator()(const WrapperArray<int,char>& arr) const; 
// };

// struct IntWrapArrayEqual { 
//     bool operator()(const WrapperArray<int,char>& arr1, 
//         const WrapperArray<int,char>& arr2) const; 
// };

// size_t IntWrapArrayHash::operator()(const WrapperArray<int,char>& arr) const {

//     std::size_t seed = 0;
//     for (auto i = 0; i < arr.length(); i ++) {
//         boost::hash_combine(seed, arr[i]);
//         // cout << "seed" << i << ": " << seed << endl;
//     }

//     cout << "IntWrapArrayHash: " << arr << " s: " << seed << endl;

//     return seed;
// }

// bool IntWrapArrayEqual::operator()(const WrapperArray<int,char>& arr1, 
//     const WrapperArray<int,char>& arr2) const {

//     cout << "IntWrapArrayEqual: " << arr1 << ", " << arr2 << ", ";

//     if (arr1.length() != arr2.length()) {
//         cout << "false 1" << endl;
//         return false;
//     }

//     for (auto i = 0; i < arr1.length(); i ++) {
//         if (arr1[i] != arr2[i]) {
//             cout << "false 2" << endl;
//             return false;
//         }
//     }

//     cout << "true" << endl;
//     return true;
// }

#endif
