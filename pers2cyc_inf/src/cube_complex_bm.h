#ifndef _CUBE_COMPLEX_BM_H_
#define _CUBE_COMPLEX_BM_H_

#include "cell_complex.h"

/*

 The vertices for the 3d-cell (cube) should be specified in the following order:

          7----------6
         /|         /|
        / |        / |
       /  |       /  |
      /   3------/---2
     4----------5   /
     |  /       |  /
     | /        | /
     |/         |/  
     0----------1

 */

class CubeComplex : public CellComplex {
public:
    CubeComplex();

    virtual ~CubeComplex();


    /* implemented virtual methods */

    virtual void toCanonVerts(const int d, vector<int> *p_verts);

    virtual void getFaceVerts(
        const int d, const vector<int>& verts, 
        const int i, vector<int>* face_verts);

    virtual void getVertPos(const int v, Vector3d* pos);

    virtual double getWeight(const int d, const int cell_id);



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

    // !! getId should return -1 if the cell does not reside in the complex !!
    _FUNC_DESC
    int getCellIdNoConvCanon(const int d, const vector<int> &verts);

    // !! the returned 'verts' should be in canon order !!
    _FUNC_DESC
    void getCellVerts(const int d, const int cell_id, vector<int>* verts);

    _FUNC_DESC
    const vector<int>* getCellCofaces(const int d, const int cell_id);

    _FUNC_DESC
    int addCellNoConvCanon(const int d, const vector<int> &verts);

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

    _FUNC_DESC
    void print();

    _FUNC_DESC
    void printLoadStat();



    // resolution of the original volume, 
    // number of slices on each dimension (vertex resolution)
    void setResolution(int x_res_vert, int y_res_vert, int z_res_ver);

    int getXRes() { return x_res_vert_; }

    int getYRes() { return y_res_vert_; }
    
    int getZRes() { return z_res_vert_; }

    void getCellInternalVerts(const int d, const int cell_id, vector<int>* verts);

    int getInternalIdNoConvCanon(const int d, const vector<int> &verts);

    int getCell1dId_X(int x, int y, int z);

    int getCell1dId_Y(int x, int y, int z);

    int getCell1dId_Z(int x, int y, int z);

    int getCell2dId_XY(int x, int y, int z);

    int getCell2dId_XZ(int x, int y, int z);

    int getCell2dId_YZ(int x, int y, int z);

    int getCell3dId(int x, int y, int z);

    void getCell1dCoord_X(int cell_id, int* x, int* y, int* z);

    void getCell1dCoord_Y(int cell_id, int* x, int* y, int* z);

    void getCell1dCoord_Z(int cell_id, int* x, int* y, int* z);

    void getCell2dCoord_XY(int cell_id, int* x, int* y, int* z);

    void getCell2dCoord_XZ(int cell_id, int* x, int* y, int* z);

    void getCell2dCoord_YZ(int cell_id, int* x, int* y, int* z);

    void getCell3dCoord(int cell_id, int* x, int* y, int* z);

    vector<int>* getIntVec() {
        auto vec = &(int_vec_pool_.at(int_vec_pool_index_));
        vec->clear();
        int_vec_pool_index_ ++;
        int_vec_pool_index_ %= int_vec_pool_.size();

        return vec;
    }

    void writeEdges(const std::string& filename);
    

    static void toCanonVertsStatic(const int d, vector<int> *p_verts);

    static const int cube_faces_lookup_[6][4];

public:
    int x_res_vert_, y_res_vert_, z_res_vert_;
    int x_res_cube_, y_res_cube_, z_res_cube_;
    // int x_cnt_, y_cnt_, z_cnt_;
    Decimal x_width_, y_width_, z_width_;

    int cell_1d_x_size_, cell_1d_y_size_, cell_1d_z_size_;
    int cell_2d_xy_size_, cell_2d_xz_size_, cell_2d_yz_size_;

    vector<vector<bool>> cell_in_complex_;
    vector<int> cell_cnt_;
    vector<int> cell_size_;

    vector<vector<int>> int_vec_pool_;
    int int_vec_pool_index_;
};

// bool loadCubeComplex(const std::string& filename, 
//     CubeComplex* complex, vector<int>* pos_cell_2d_verts);

#endif
