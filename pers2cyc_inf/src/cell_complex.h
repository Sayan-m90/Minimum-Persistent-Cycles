// - cell adding needs not be a filtration, however, 
//   whenever a cell is added, none of its face can be added after this.
// - top dimension doesn't store cofaces and vertices query map by defaulst

#ifndef _CELL_COMPLEX_H_
#define _CELL_COMPLEX_H_

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "utils.h"
#include "arr_label_id_map.hpp"

#define _FUNC_DESC virtual
// #define _FUNC_DESC

enum ComplexType {
    SIMP_CMPLX, CUBE_CMPLX
};

typedef ArrayLabelIdMap<int,int> VertArrIdMap;

typedef vector<vector<int>> CofaceContainer;

class CellComplex {

public:

    virtual ~CellComplex();

    _FUNC_DESC
    void init(const int min_dim, const int max_dim, const bool top_dim_query = false);

    _FUNC_DESC
    void setCofaceCountHint(const int d, const int cnt);

    _FUNC_DESC
    void reserveCellSize(const int d, const int cnt);

    // the cell size may be greater than the cell count if some cells are deleted.
    // cell size is used to traverse cells. validity check of cell id is needed 
    // if ever some cells are deleted.
    _FUNC_DESC
    // !! for bm cube complex this should alwasy return all possible cells !!
    int getCellSize(const int d);

    _FUNC_DESC
    int getCellCount(const int d);

    virtual void toCanonVerts(const int d, vector<int> *p_verts) = 0;


    // !! getId should return -1 if the cell does not reside in the complex !!
    _FUNC_DESC
    int getCellIdNoConvCanon(const int d, const vector<int> &verts);

    int getCellId(const int d, const vector<int>& verts) {
        vector<int> verts_c(verts);
        this->toCanonVerts(d, &verts_c);
        return this->getCellIdNoConvCanon(d, verts_c);
    }

    // !! the returned 'verts' should be in canon order !!
    _FUNC_DESC
    void getCellVerts(const int d, const int cell_id, vector<int>* verts);

    _FUNC_DESC
    const vector<int>* getCellCofaces(const int d, const int cell_id);

    _FUNC_DESC
    int addCellNoConvCanon(const int d, const vector<int> &verts);

    int addCell(const int d, const vector<int> &verts) {
        vector<int> verts_c(verts);
        this->toCanonVerts(d, &verts_c);
        return this->addCellNoConvCanon(d, verts_c);
    }

    // bool getCofaces(const vector<Cell*> &cells, vector<Cell*>* cofaces);

    // return the number of faces of a d-cell
    int getFaceCnt(const int d) { return face_cnt_[d]; }

    // return the number of vertices of a d-cell
    int getVertCnt(const int d) { return vert_cnt_[d]; }

    // get the spannin vertices of the i-th face of a d-cell,
    // the returning vertices is not necessarily in canonical order
    virtual void getFaceVerts(const int d, const vector<int>& verts, 
        const int i, vector<int>* face_verts) = 0;

    _FUNC_DESC
    // CAUTION!!: the deletion does not guarantee the validity of the complex
    void deleteCell(const int d, const int cell_id);

    _FUNC_DESC
    // CAUTION!!: the deletion does not guarantee the validity of the complex
    void deleteCoface(const int d, const int cell_id, const int coface_id);

    _FUNC_DESC
    // check if the cell valid (not deleted)
    bool isCellValid(const int d, const int cell_id);


    /* CAUTION!!: the following are just for memory efficiency, 
       the adjacency of the complex will be ruined */

    _FUNC_DESC
    void deleteDimCellData(const int d);

    _FUNC_DESC
    void deleteDimCofaces(const int d);
    
    _FUNC_DESC
    void deleteDimVertQueryMap(const int d);


    virtual void getVertPos(const int v, Vector3d* pos) = 0;

    virtual double getWeight(const int d, const int cell_id) = 0;

    void writeMesh(const std::string& filename);

    _FUNC_DESC
    void print();

    _FUNC_DESC
    void printLoadStat();

    // check whether every face of each cell reside in the complex
    // and every coface stored for each cell reside in the complex
    // bool checkComplexValid();

protected:

    // spanning vertices maps for cells of each dimension
    vector<VertArrIdMap*> vert_arr_maps_;

    vector<CofaceContainer*> cofaces_pool_;

    // vector<bool> has_vert_query_;

    int min_dim_, max_dim_;

    // number of faces for cells of each dimension
    vector<int> face_cnt_; 

    // number of vertices for cells of each dimension
    vector<int> vert_cnt_; 

    vector<int> coface_cnt_hint_;
};

#endif














