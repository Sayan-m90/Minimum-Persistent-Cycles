#ifndef _MESH_WRITER_H_
#define _MESH_WRITER_H_

#include <unordered_map>
#include "cell_complex.h"

class MeshWriter {
public:
    MeshWriter(CellComplex* complex,  
        const int face_vert_cnt) : 
        complex_(complex),
        face_vert_cnt_(face_vert_cnt) {}

    void addFace(const vector<int>& face);

    void write(const std::string& filename);

private:
    CellComplex* complex_;
    std::unordered_map<int,int> vert_ind_map_;
    vector<int> face_pool_;
    vector<int> vert_pool_;
    const int face_vert_cnt_;
};

#endif
