#include "cube_complex_bm.h"
#include <iostream>
#include <algorithm>
#include "utils.h"

/*
 id numbering:

 1-cells: 
    cells parallel to x axis first,
    cells parallel to y axis second,
    cells parallel to z axis third.

 2-cells:
    cells parallel to x-y plane first,
    cells parallel to x-z plane second,
    cells parallel to y-z plane third.
 */

const int CubeComplex::cube_faces_lookup_[6][4] = {
    {0,1,2,3},
    {4,5,6,7},
    {0,1,5,4},
    {3,2,6,7},
    {0,3,7,4},
    {1,2,6,5}
};

CubeComplex::CubeComplex() {
    vert_cnt_ = {1,2,4,8};
    face_cnt_ = {0,2,4,6};
}

CubeComplex::~CubeComplex() {}

void CubeComplex::toCanonVertsStatic(const int d, vector<int> *p_verts) {
    if (d == 1) {
        std::sort(p_verts->begin(), p_verts->end());
    } else if (d == 2) {
        vector<int> verts_c(*p_verts);
        auto min_ind = 0;
        for (auto i = 1; i < verts_c.size(); i ++) {
            if (verts_c[i] < verts_c[min_ind]) {
                min_ind = i;
            }
        }

        p_verts->clear();
        if (verts_c[mod(min_ind+1, verts_c.size())] < verts_c[mod(min_ind-1, verts_c.size())]) {
            for (auto i = 0; i < verts_c.size(); i ++) {
                p_verts->push_back(verts_c[mod(min_ind+i, verts_c.size())]);
            }
        } else {
            for (auto i = 0; i < verts_c.size(); i ++) {
                p_verts->push_back(verts_c[mod(min_ind-i, verts_c.size())]);
            }
        }
    }
#ifdef DEBUG_OPT
    else if (d == 3) {
        auto min_iter = std::min_element(std::begin(*p_verts), std::end(*p_verts));
        if (min_iter != std::begin(*p_verts)) {
            cout << "FATAL: 3-cell verts not start with min vert in "
                "CubeComplex::toCanonVertsStatic" << endl;
            exit(-1);
        }
    }
#endif
}

void CubeComplex::init(const int min_dim, const int max_dim, const bool top_dim_query) {
    min_dim_ = min_dim;
    max_dim_ = max_dim;
}

void CubeComplex::toCanonVerts(const int d, vector<int> *p_verts) {
    CubeComplex::toCanonVertsStatic(d, p_verts);
}

void CubeComplex::getFaceVerts(
    const int d, const vector<int>& verts, 
    const int i, vector<int>* face_verts) {

    if (d == 2) {
        *face_verts = { verts[i], verts[mod(i+1, verts.size())] };
    } else if (d == 3) {
        face_verts->clear();

        for (int j = 0; j < 4; j ++) {
            auto ind = CubeComplex::cube_faces_lookup_[i][j];
            face_verts->push_back(verts[ind]);
        }
    }
}

void CubeComplex::getVertPos(const int v, Vector3d* pos) {
    (*pos)[2] = (v / (y_res_vert_*x_res_vert_));
    auto left = v % (y_res_vert_*x_res_vert_);
    (*pos)[1] = (left / x_res_vert_);
    (*pos)[0] = (left % x_res_vert_);
}

double CubeComplex::getWeight(const int d, const int cell_id) {
    return 1.0;
}

void CubeComplex::setResolution(int x_res_vert, int y_res_vert, int z_res_vert) {
    x_res_vert_ = x_res_vert;
    y_res_vert_ = y_res_vert;
    z_res_vert_ = z_res_vert;

    x_res_cube_ = x_res_vert_ - 1;
    y_res_cube_ = y_res_vert_ - 1;
    z_res_cube_ = z_res_vert_ - 1;

    cell_cnt_.resize(4, 0);

    // number of 1-cells parallel to x axis
    cell_1d_x_size_ = x_res_cube_ * (y_res_cube_ + 1) * (z_res_cube_ + 1);
    // number of 1-cells parallel to y axis
    cell_1d_y_size_ = (x_res_cube_ + 1) * y_res_cube_ * (z_res_cube_ + 1);
    // number of 1-cells parallel to z axis
    cell_1d_z_size_ = (x_res_cube_ + 1) * (y_res_cube_ + 1) * z_res_cube_;

    // number of 2-cells parallel to x-y plane
    cell_2d_xy_size_ = x_res_cube_ * y_res_cube_ * (z_res_cube_ + 1);
    // number of 2-cells parallel to x-z plane
    cell_2d_xz_size_ = x_res_cube_ * (y_res_cube_ + 1) * z_res_cube_;
    // number of 2-cells parallel to y-z plane
    cell_2d_yz_size_ = (x_res_cube_ + 1) * y_res_cube_ * z_res_cube_;

    cell_size_.resize(4);
    cell_size_[1] = cell_1d_x_size_ + cell_1d_y_size_ + cell_1d_z_size_;
    cell_size_[2] = cell_2d_xy_size_ + cell_2d_xz_size_ + cell_2d_yz_size_;
    cell_size_[3] = x_res_cube_ * y_res_cube_ * z_res_cube_;

    cell_in_complex_.resize(4);
    cell_in_complex_[1].resize(cell_size_[1], false);
    cell_in_complex_[2].resize(cell_size_[2], false);
    cell_in_complex_[3].resize(cell_size_[3], false);

    int_vec_pool_.resize(20);
    int_vec_pool_index_ = 0;
}

int CubeComplex::getCell1dId_X(int x, int y, int z) {
    int cell_local_id;
    packDimId(x, y, z, 
        x_res_cube_, y_res_cube_ + 1, z_res_cube_ + 1, 
        &cell_local_id);

    return cell_local_id;
}

int CubeComplex::getCell1dId_Y(int x, int y, int z) {
    int cell_local_id;
    packDimId(x, y, z, 
        x_res_cube_ + 1, y_res_cube_, z_res_cube_ + 1, 
        &cell_local_id);

    return cell_local_id + cell_1d_x_size_;
}

int CubeComplex::getCell1dId_Z(int x, int y, int z) {
    int cell_local_id;
    packDimId(x, y, z, 
        x_res_cube_ + 1, y_res_cube_ + 1, z_res_cube_, 
        &cell_local_id);

    return cell_local_id + cell_1d_x_size_ + cell_1d_y_size_;
}

int CubeComplex::getCell2dId_XY(int x, int y, int z) {
    if (x >= 0 && x < x_res_cube_ &&
        y >= 0 && y < y_res_cube_ &&
        z >= 0 && z < z_res_cube_ + 1) {

        int cell_local_id;
        packDimId(x, y, z, 
            x_res_cube_, y_res_cube_, z_res_cube_ + 1, 
            &cell_local_id);

        return cell_local_id;
    } else {
        return -1;
    }
}

int CubeComplex::getCell2dId_XZ(int x, int y, int z) {
    if (x >= 0 && x < x_res_cube_ &&
        y >= 0 && y < y_res_cube_ + 1 &&
        z >= 0 && z < z_res_cube_) {

        int cell_local_id;
        packDimId(x, y, z, 
            x_res_cube_, y_res_cube_ + 1, z_res_cube_, 
            &cell_local_id);

        return cell_local_id + cell_2d_xy_size_;
    } else {
        return -1;
    }
}

int CubeComplex::getCell2dId_YZ(int x, int y, int z) {
    if (x >= 0 && x < x_res_cube_ + 1 &&
        y >= 0 && y < y_res_cube_ &&
        z >= 0 && z < z_res_cube_) {

        int cell_local_id;
        packDimId(x, y, z, 
            x_res_cube_ + 1, y_res_cube_, z_res_cube_, 
            &cell_local_id);

        return cell_local_id + cell_2d_xy_size_ + cell_2d_xz_size_;
    } else {
        return -1;
    }
}

int CubeComplex::getCell3dId(int x, int y, int z) {
    if (x >= 0 && x < x_res_cube_ &&
        y >= 0 && y < y_res_cube_ &&
        z >= 0 && z < z_res_cube_) {

        int cell_local_id;
        packDimId(x, y, z, 
            x_res_cube_, y_res_cube_, z_res_cube_, 
            &cell_local_id);

        return cell_local_id;
    } else {
        return -1;
    }
}

void CubeComplex::getCell1dCoord_X(int cell_id, int* x, int* y, int* z) {
    unpackDimId(cell_id, 
        x_res_cube_, y_res_cube_ + 1, z_res_cube_ + 1, 
        x, y, z);
}

void CubeComplex::getCell1dCoord_Y(int cell_id, int* x, int* y, int* z) {
    unpackDimId(cell_id - cell_1d_x_size_, 
        x_res_cube_ + 1, y_res_cube_, z_res_cube_ + 1, 
        x, y, z);
}

void CubeComplex::getCell1dCoord_Z(int cell_id, int* x, int* y, int* z) {
    unpackDimId(cell_id - cell_1d_x_size_ - cell_1d_y_size_, 
        x_res_cube_ + 1, y_res_cube_ + 1, z_res_cube_, 
        x, y, z);
}

void CubeComplex::getCell2dCoord_XY(int cell_id, int* x, int* y, int* z) {
    unpackDimId(cell_id, 
        x_res_cube_, y_res_cube_, z_res_cube_ + 1, 
        x, y, z);
}

void CubeComplex::getCell2dCoord_XZ(int cell_id, int* x, int* y, int* z) {
    unpackDimId(cell_id - cell_2d_xy_size_, 
        x_res_cube_, y_res_cube_ + 1, z_res_cube_, 
        x, y, z);
}

void CubeComplex::getCell2dCoord_YZ(int cell_id, int* x, int* y, int* z) {
    unpackDimId(cell_id - cell_2d_xy_size_ - cell_2d_xz_size_, 
        x_res_cube_ + 1, y_res_cube_, z_res_cube_, 
        x, y, z);
}

void CubeComplex::getCell3dCoord(int cell_id, int* x, int* y, int* z) {
    unpackDimId(cell_id, 
        x_res_cube_, y_res_cube_, z_res_cube_, 
        x, y, z);
}

int CubeComplex::getInternalIdNoConvCanon(const int d, const vector<int> &verts) {
    int x_start, y_start, z_start;
    unpackDimId(verts.at(0), 
        x_res_vert_, y_res_vert_, z_res_vert_, 
        &x_start, &y_start, &z_start);

    if (d == 1 || d == 2) {
        int op_v;
        if (d == 1) { op_v = 1; } 
        else        { op_v = 2; }

        int x, y, z;
        unpackDimId(verts.at(op_v), 
            x_res_vert_, y_res_vert_, z_res_vert_, 
            &x, &y, &z);

        int x_vari = 0, y_vari = 0, z_vari = 0;
        if (x > x_start) { x_vari = 1; }
        if (y > y_start) { y_vari = 1; }
        if (z > z_start) { z_vari = 1; }

        if (x_vari + y_vari + z_vari != d) {
            cout << "FATAL: x_vari + y_vari + z_vari != d in "
                "CubeComplex::getInternalIdNoConvCanon" << endl;
            exit(-1);
        }

        if (d == 1) {
            if (x_vari == 1) {
                return getCell1dId_X(x_start, y_start, z_start);
            } else if (y_vari == 1) {
                return getCell1dId_Y(x_start, y_start, z_start);
            } else {
                return getCell1dId_Z(x_start, y_start, z_start);
            }
        } else {
            if (x_vari == 1 && y_vari == 1) {
                // parallel x-y plane
                return getCell2dId_XY(x_start, y_start, z_start);
            } else if (x_vari == 1 && z_vari == 1) {
                // parallel x-z plane
                return getCell2dId_XZ(x_start, y_start, z_start);
            } else {
                // parallel y-z plane
                return getCell2dId_YZ(x_start, y_start, z_start);
            }
        }
    } else {  // d == 3
        return getCell3dId(x_start, y_start, z_start);
    }
}

// the cell size may be greater than the cell count if some cells are deleted.
// cell size is used to traverse cells. validity check of cell id is needed 
// if ever some cells are deleted.
// !! for bm cube complex this should alwasy return all possible cells !!
int CubeComplex::getCellSize(const int d) {
    return cell_size_[d];
}

int CubeComplex::getCellCount(const int d) {
    return cell_cnt_[d];
}

// !! getId should return -1 if the cell does not reside in the complex !!
int CubeComplex::getCellIdNoConvCanon(const int d, const vector<int> &verts) {
    int cell_id = getInternalIdNoConvCanon(d, verts);

    if (cell_in_complex_[d][cell_id]) {
        return cell_id;
    } else {
        return -1;
    }
}

void CubeComplex::getCellVerts(const int d, const int cell_id, vector<int>* verts) {
    verts->clear();
    if (cell_in_complex_[d][cell_id]) {
        getCellInternalVerts(d, cell_id, verts);
    }
}

// !! the returned 'verts' should be in canon order !!
void CubeComplex::getCellInternalVerts(const int d, const int cell_id, vector<int>* verts) {
    int x_start, y_start, z_start;
    int x_local_cnt, y_local_cnt, z_local_cnt;

    if (d == 1) {
        if (cell_id < cell_1d_x_size_) {
            // parallel to x axis

            getCell1dCoord_X(cell_id, &x_start, &y_start, &z_start);
            x_local_cnt = 2;
            y_local_cnt = 1;
            z_local_cnt = 1;
        } else if (cell_id < cell_1d_x_size_ + cell_1d_y_size_) {
            // parallel to y axis

            getCell1dCoord_Y(cell_id, &x_start, &y_start, &z_start);
            x_local_cnt = 1;
            y_local_cnt = 2;
            z_local_cnt = 1;
        } else if (cell_id < cell_size_[1]) {
            // parallel to z axis

            getCell1dCoord_Z(cell_id, &x_start, &y_start, &z_start);
            x_local_cnt = 1;
            y_local_cnt = 1;
            z_local_cnt = 2;
        } else {
            cout << "FATAL: invalid 1-cell id: " << cell_id << endl;
            exit(-1);
        }
    } else if (d == 2) {
        if (cell_id < cell_2d_xy_size_) {
            // parallel x-y plane

            getCell2dCoord_XY(cell_id, &x_start, &y_start, &z_start);
            x_local_cnt = 2;
            y_local_cnt = 2;
            z_local_cnt = 1;
        } else if (cell_id < cell_2d_xy_size_ + cell_2d_xz_size_) {
            // parallel x-z plane

            getCell2dCoord_XZ(cell_id, &x_start, &y_start, &z_start);
            x_local_cnt = 2;
            y_local_cnt = 1;
            z_local_cnt = 2;
        } else if (cell_id < cell_size_[2]) {
            // parallel y-z plane

            getCell2dCoord_YZ(cell_id, &x_start, &y_start, &z_start);
            x_local_cnt = 1;
            y_local_cnt = 2;
            z_local_cnt = 2;
        } else {
            cout << "FATAL: invalid 2-cell id: " << cell_id << endl;
            exit(-1);
        }
    } else if (d == 3) {
        getCell3dCoord(cell_id, &x_start, &y_start, &z_start);

        x_local_cnt = 2;
        y_local_cnt = 2;
        z_local_cnt = 2;
    }

    vector<int> temp_verts;
    temp_verts.reserve(8);

    for (int x = x_start; x < x_start + x_local_cnt; x ++) {
        for (int y = y_start; y < y_start + y_local_cnt; y ++) {
            for (int z = z_start; z < z_start + z_local_cnt; z ++) {
                int v_id;
                packDimId(x, y, z, x_res_vert_, y_res_vert_, z_res_vert_, &v_id);
                temp_verts.push_back(v_id);
            }
        }
    }

    if (d == 1) {
        *verts = temp_verts;
    } else if (d == 2) {
        *verts = {temp_verts[0], temp_verts[1], temp_verts[3], temp_verts[2]};
        this->toCanonVerts(d, verts);
    } else if (d == 3) {
        *verts = {temp_verts[0], temp_verts[1], temp_verts[3], temp_verts[2],
            temp_verts[4], temp_verts[5], temp_verts[7], temp_verts[6]};
    }
}

const vector<int>* CubeComplex::getCellCofaces(const int d, const int cell_id) {
    vector<int> cofaces;
    cofaces.reserve(4);
    int x_start, y_start, z_start;
    int x, y, z;

    if (d == 1) {
        if (cell_id < cell_1d_x_size_) {
            // parallel to x axis

            getCell1dCoord_X(cell_id, &x_start, &y_start, &z_start);

            // cofaces parallel to x-y plane
            x = x_start;
            z = z_start;
            for (y = y_start - 1; y <= y_start; y ++) {
                auto cof_id = getCell2dId_XY(x, y, z);
                if (cof_id >= 0) {
                    cofaces.push_back(cof_id);
                }
            }

            // cofaces parallel to x-z plane
            x = x_start;
            y = y_start;
            for (z = z_start - 1; z <= z_start; z ++) {
                auto cof_id = getCell2dId_XZ(x, y, z);
                if (cof_id >= 0) {
                    cofaces.push_back(cof_id);
                }
            }
        } else if (cell_id < cell_1d_x_size_ + cell_1d_y_size_) {
            // parallel to y axis

            getCell1dCoord_Y(cell_id, &x_start, &y_start, &z_start);

            // cofaces parallel to x-y plane
            y = y_start;
            z = z_start;
            for (x = x_start - 1; x <= x_start; x ++) {
                auto cof_id = getCell2dId_XY(x, y, z);
                if (cof_id >= 0) {
                    cofaces.push_back(cof_id);
                }
            }

            // cofaces parallel to y-z plane
            x = x_start;
            y = y_start;
            for (z = z_start - 1; z <= z_start; z ++) {
                auto cof_id = getCell2dId_YZ(x, y, z);
                if (cof_id >= 0) {
                    cofaces.push_back(cof_id);
                }
            }
        } else if (cell_id < cell_size_[1]) {
            // parallel to z axis

            getCell1dCoord_Z(cell_id, &x_start, &y_start, &z_start);

            // cofaces parallel to x-z plane
            y = y_start;
            z = z_start;
            for (x = x_start - 1; x <= x_start; x ++) {
                auto cof_id = getCell2dId_XZ(x, y, z);
                if (cof_id >= 0) {
                    cofaces.push_back(cof_id);
                }
            }

            // cofaces parallel to y-z plane
            x = x_start;
            z = z_start;
            for (y = y_start - 1; y <= y_start; y ++) {
                auto cof_id = getCell2dId_YZ(x, y, z);
                if (cof_id >= 0) {
                    cofaces.push_back(cof_id);
                }
            }        
        } else {
            cout << "FATAL: invalid 1-cell id: " << cell_id 
                << " in CubeComplex::getCellCofaces" << endl;
            exit(-1);
        }
    } else if (d == 2) {
        if (cell_id < cell_2d_xy_size_) {
            // parallel x-y plane

            getCell2dCoord_XY(cell_id, &x_start, &y_start, &z_start);

            x = x_start;
            y = y_start;
            for (z = z_start - 1; z <= z_start; z ++) {
                auto cof_id = getCell3dId(x, y, z);
                if (cof_id >= 0) {
                    cofaces.push_back(cof_id);
                }
            }
        } else if (cell_id < cell_2d_xy_size_ + cell_2d_xz_size_) {
            // parallel x-z plane

            getCell2dCoord_XZ(cell_id, &x_start, &y_start, &z_start);

            x = x_start;
            z = z_start;
            for (y = y_start - 1; y <= y_start; y ++) {
                auto cof_id = getCell3dId(x, y, z);
                if (cof_id >= 0) {
                    cofaces.push_back(cof_id);
                }
            }
        } else if (cell_id < cell_size_[2]) {
            // parallel y-z plane

            getCell2dCoord_YZ(cell_id, &x_start, &y_start, &z_start);

            y = y_start;
            z = z_start;
            for (x = x_start - 1; x <= x_start; x ++) {
                auto cof_id = getCell3dId(x, y, z);
                if (cof_id >= 0) {
                    cofaces.push_back(cof_id);
                }
            }
        } else {
            cout << "FATAL: invalid 2-cell id: " << cell_id 
                << " in CubeComplex::getCellCofaces" << endl;
            exit(-1);
        }
    } else {
        cout << "FATAL: d not right in CubeComplex::getCellCofaces" << endl;
        exit(-1);
    }


    // check if in complex
    auto true_cofaces = getIntVec();
    for (auto cof_id : cofaces) {
        if (cell_in_complex_[d + 1][cof_id]) {
            true_cofaces->push_back(cof_id);
        }
    }

    // cout << "true_cofaces addr: " << true_cofaces << endl;

    return true_cofaces;
}

int CubeComplex::addCellNoConvCanon(const int d, const vector<int> &verts) {
    int cell_id = this->getInternalIdNoConvCanon(d, verts);
    if (cell_in_complex_[d][cell_id]) {
        return -1;
    }

    cell_in_complex_[d][cell_id] = true;
    cell_cnt_[d] ++;

    if (d == min_dim_) {
        return cell_id;
    }

    vector<int> f_verts;
    f_verts.reserve(this->getVertCnt(d - 1));

    for (auto i = 0; i < this->getFaceCnt(d); i ++) {
        this->getFaceVerts(d, verts, i, &f_verts);
        this->toCanonVerts(d - 1, &f_verts);
        this->addCellNoConvCanon(d - 1, f_verts);
    }

    return cell_id;
}

// CAUTION!!: the deletion does not guarantee the validity of the complex
void CubeComplex::deleteCell(const int d, const int cell_id) {
    if (cell_in_complex_[d][cell_id]) {
        cell_in_complex_[d][cell_id] = false;
        cell_cnt_[d] --;
    }
}

// check if the cell valid (not deleted)
bool CubeComplex::isCellValid(const int d, const int cell_id) {
    return cell_in_complex_[d][cell_id];
}

void CubeComplex::writeEdges(const std::string& filename) {
    std::ofstream fout(filename);

    std::unordered_map<int,int> vert_ind_map;
    vector<int> vert_pool;
    vector<int> verts;

    for (auto cell_1d_id = 0; cell_1d_id < cell_size_[1]; cell_1d_id ++) {
        if (cell_in_complex_[1][cell_1d_id]) {
            this->getCellVerts(1, cell_1d_id, &verts);

            for (auto v : verts) {
                auto iter = vert_ind_map.find(v);
                if (iter == vert_ind_map.end()) {
                    auto new_v = vert_ind_map.size();
                    vert_ind_map.insert({v, new_v});
                    vert_pool.push_back(v);
                }
            }
        }
    }

    fout << "ply" << endl
        << "format ascii 1.0" << endl
        << "element vertex " << vert_ind_map.size() << endl
        << "property float x" << endl
        << "property float y" << endl
        << "property float z" << endl
        << "element edge " << cell_cnt_[1] << endl
        << "property int vertex1" << endl
        << "property int vertex2" << endl
        << "end_header" << endl;

    Vector3d pos;
    for (auto v : vert_pool) {
        getVertPos(v, &pos);
        fout << pos[0] << ' ' << pos[1] << ' ' << pos[2] << endl;
    }

    for (auto cell_1d_id = 0; cell_1d_id < cell_size_[1]; cell_1d_id ++) {
        if (cell_in_complex_[1][cell_1d_id]) {
            this->getCellVerts(1, cell_1d_id, &verts);
            fout << vert_ind_map[verts[0]] << ' ' << vert_ind_map[verts[1]] << endl;
        }
    }
}

// bool loadCubeComplex(const std::string& filename, 
//     CubeComplex* complex, vector<int>* pos_cell_2d_verts) {

//     std::ifstream fin(filename);

//     int x_res_vert, y_res_vert, z_res_vert;
//     fin >> x_res_vert >> y_res_vert >> z_res_vert;
//     // cout << x_res_vert << " " << y_res_vert << " " << z_res_vert << endl;
//     complex->setResolution(x_res_vert, y_res_vert, z_res_vert);

//     int cell_1d_cnt, cell_2d_cnt, cell_3d_cnt;
//     fin >> cell_1d_cnt >> cell_2d_cnt >> cell_3d_cnt;
//     // cout << cell_1d_cnt << " " << cell_2d_cnt << " " << cell_3d_cnt << endl;
//     complex->reserveCellSize(1, cell_1d_cnt);
//     complex->reserveCellSize(2, cell_2d_cnt);
//     complex->reserveCellSize(3, cell_3d_cnt);

//     loadComplex(fin, complex, pos_cell_2d_verts);

//     return true;
// }

void CubeComplex::setCofaceCountHint(const int d, const int cnt) {}

void CubeComplex::reserveCellSize(const int d, const int cnt) {}

void CubeComplex::deleteCoface(const int d, const int cell_id, const int coface_id) {}

void CubeComplex::deleteDimCellData(const int d) {}

void CubeComplex::deleteDimCofaces(const int d) {}

void CubeComplex::deleteDimVertQueryMap(const int d) {}

void CubeComplex::print() {}

void CubeComplex::printLoadStat() {}
