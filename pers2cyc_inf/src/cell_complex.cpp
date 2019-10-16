#include "cell_complex.h"
#include <iostream>
#include <string>
#include "mesh_writer.h"

CellComplex::~CellComplex() {
    for (auto i = 0; i < vert_arr_maps_.size(); i ++) {
        if (vert_arr_maps_[i] != NULL) {
            delete vert_arr_maps_[i];
        }
        if (cofaces_pool_[i] != NULL) {
            delete cofaces_pool_[i];
        }
    }
}

void CellComplex::init(
    const int min_dim, 
    const int max_dim, 
    const bool top_dim_query) {

    min_dim_ = min_dim;
    max_dim_ = max_dim;

    vert_arr_maps_.resize(max_dim_ - min_dim_ + 1);
    cofaces_pool_.resize(max_dim_ - min_dim_ + 1);
    coface_cnt_hint_.resize(max_dim_ - min_dim_ + 1, 0);

    for (auto i = 0; i < max_dim_ - min_dim_;  i ++) {
        vert_arr_maps_[i] = new VertArrIdMap;
        vert_arr_maps_[i]->init(vert_cnt_[min_dim_+i], true);
        cofaces_pool_[i] = new CofaceContainer;
    }

    vert_arr_maps_[max_dim_-min_dim_] = new VertArrIdMap;
    vert_arr_maps_[max_dim_-min_dim_]->init(vert_cnt_[max_dim_], top_dim_query);
    cofaces_pool_[max_dim_-min_dim_] = NULL;
}

void CellComplex::setCofaceCountHint(const int d, const int cnt) {
    coface_cnt_hint_[d - min_dim_] = cnt; 
}

void CellComplex::reserveCellSize(const int d, const int cnt) {
    if (vert_arr_maps_.at(d - min_dim_) != NULL) {
        vert_arr_maps_[d - min_dim_]->reserve(cnt);
        // cout << "vert_arr_maps_ reserve " << d << " " << cnt << endl;
    }
    if (cofaces_pool_.at(d - min_dim_) != NULL) {
        cofaces_pool_[d - min_dim_]->reserve(cnt);
        // cout << "cofaces_pool_ reserve " << d << " " << cnt << endl;
    }
}

int CellComplex::getCellSize(const int d) {
    return vert_arr_maps_[d - min_dim_]->size();
}

int CellComplex::getCellCount(const int d) {
    return vert_arr_maps_[d - min_dim_]->elemCount();
}

int CellComplex::getCellIdNoConvCanon(const int d, const vector<int> &verts) {
    return vert_arr_maps_[d - min_dim_]->getId(verts);
}

void CellComplex::getCellVerts(const int d, const int cell_id, vector<int>* verts) {
    vert_arr_maps_[d - min_dim_]->getArrayLabel(cell_id, verts);
}

const vector<int>* CellComplex::getCellCofaces(const int d, const int cell_id) {
    return &(cofaces_pool_[d - min_dim_]->at(cell_id));
}

int CellComplex::addCellNoConvCanon(const int d, const vector<int>& verts) {

    int cell_id = vert_arr_maps_[d - min_dim_]->addArrayLabel(verts);
    if (d != max_dim_) {
        cofaces_pool_[d - min_dim_]->emplace_back();

        if (coface_cnt_hint_[d - min_dim_] != 0) {
            cofaces_pool_[d - min_dim_]->at(cell_id)
                .reserve(coface_cnt_hint_[d - min_dim_]);
        }
    }

    if (d == min_dim_) {
        return cell_id;
    }

    vector<int> f_verts;
    f_verts.reserve(this->getVertCnt(d-1));

    for (auto i = 0; i < this->getFaceCnt(d); i ++) {
        this->getFaceVerts(d, verts, i, &f_verts);
        this->toCanonVerts(d-1, &f_verts);
        // this->addFaceNoConvCanon(cell_id, d-1, f_verts);

        auto face_id = vert_arr_maps_[d-1-min_dim_]->getId(f_verts);
        if (face_id < 0) {
            face_id = this->addCellNoConvCanon(d-1, f_verts);
        }
        cofaces_pool_[d-1-min_dim_]->at(face_id).push_back(cell_id);
    }

    return cell_id;
}

void CellComplex::deleteCell(const int d, const int cell_id) {
    vert_arr_maps_[d-min_dim_]->deleteArrayLabel(cell_id, -1);
}

void CellComplex::deleteCoface(
    const int d, const int cell_id, const int coface_id) {

    vector<int>* cofaces = &(cofaces_pool_[d-min_dim_]->at(cell_id));

    vector<int>::iterator iter;
    for (iter = cofaces->begin(); iter != cofaces->end(); iter ++) {
        if (*iter == coface_id) {
            break;
        }
    }

    if (iter != cofaces->end()) {
        cofaces->erase(iter);
    }
}

bool CellComplex::isCellValid(const int d, const int cell_id) {
    int ind = vert_arr_maps_[d-min_dim_]->
        arr_label_pool_.at(this->getVertCnt(d) * cell_id);
    return ind >= 0;
}

void CellComplex::deleteDimCellData(const int d) {
    if (vert_arr_maps_[d - min_dim_] != NULL) {
        delete vert_arr_maps_[d - min_dim_];
        vert_arr_maps_[d - min_dim_] = NULL;
    }
    if (cofaces_pool_[d - min_dim_] != NULL) {
        delete cofaces_pool_[d - min_dim_];
        cofaces_pool_[d - min_dim_] = NULL;
    }
}

void CellComplex::deleteDimCofaces(const int d) {
    if (cofaces_pool_[d - min_dim_] != NULL) {
        delete cofaces_pool_[d - min_dim_];
        cofaces_pool_[d - min_dim_] = NULL;
    }
}

void CellComplex::deleteDimVertQueryMap(const int d) {
    if (vert_arr_maps_[d - min_dim_] != NULL) {
        vert_arr_maps_[d - min_dim_]->clearLabelMap();
    }
}

void CellComplex::print() {
    vector<int> verts;

    for (auto i = 0; i <= max_dim_ - min_dim_; i ++) {
        cout << "[" << i+min_dim_ << "-cells:" 
            << vert_arr_maps_[i]->elemCount() << "]" << endl;

        vector<int> cell_order;
        for (auto j = 0; j < vert_arr_maps_[i]->size(); j ++) {
            if (this->isCellValid(i+min_dim_, j)) {
                cell_order.push_back(j);
            }
        }

        struct Comp {
            Comp(CellComplex* _complex, int _d) : complex(_complex), d(_d) {}

            bool operator()(const int cell_id1, const int cell_id2) {
                vector<int> verts1, verts2;
                complex->getCellVerts(d, cell_id1, &verts1);
                complex->getCellVerts(d, cell_id2, &verts2);

                for (auto i = 0; i < verts1.size(); i ++) {
                    if (verts1[i] < 0 || verts2[i] < 0) {
                        cout << "FATAL: vert ind < 0 in CellComplex::print Comp" << endl;
                        exit(-1);
                    }

                    if (verts1[i] < verts2[i]) {
                        return true;
                    } else if (verts1[i] > verts2[i]) {
                        return false;
                    }
                }

                return false;
            }

            CellComplex* complex;
            int d;
        };

        Comp comp(this, i+min_dim_);
        std::sort(cell_order.begin(), cell_order.end(), comp);

        for (auto cell_id : cell_order) {
            this->getCellVerts(i+min_dim_, cell_id, &verts);
            cout << verts << endl;

            if (i != max_dim_ - min_dim_) {
                cout << "  cofaces:" << endl;

                vector<int> coface_order( cofaces_pool_[i]->at(cell_id) );
                Comp comp2(this, i+1+min_dim_);
                std::sort(coface_order.begin(), coface_order.end(), comp2);

                for (auto cof_id : coface_order) {
                    this->getCellVerts(i+1+min_dim_, cof_id, &verts);
                    cout << "   " << verts << endl;
                }
            }
        }

        cout << endl;
    }
}

void CellComplex::printLoadStat() {
    for (auto i = 0; i <= max_dim_ - min_dim_; i ++) {
        cout << "[map" << min_dim_ + i << "]" << endl;
        vert_arr_maps_[i]->printLoadStat();

        if (i != max_dim_ - min_dim_) {
            cout << "[cofaces_pool" << min_dim_ + i << "]" << endl 
                << "capacity=" << cofaces_pool_[i]->capacity()
                << ",size=" << cofaces_pool_[i]->size() << endl;
        }
    }
}

void CellComplex::writeMesh(const std::string& filename) {
    MeshWriter mesh_writer(this, this->getVertCnt(2));
    vector<int> cell_2d_verts;

    for (auto cell_2d_id = 0; cell_2d_id < this->getCellSize(2); cell_2d_id ++) {
        if (this->isCellValid(2, cell_2d_id)) {
            this->getCellVerts(2, cell_2d_id, &cell_2d_verts);
            mesh_writer.addFace(cell_2d_verts);
        }
    }

    mesh_writer.write(filename);
}

// TODO: 
// bool CellComplex::checkComplexValid() {
//     return true;
// }

// bool CellComplex::getCofaces(const vector<Cell*> &cells, vector<Cell*>* cofaces) {
//     CellPtSet coface_set;

//     for (auto cell: cells) {
//         for (auto cof : cell->cofaces) {
//             coface_set.insert(cof);
//         }
//     }

//     for (auto cof : coface_set) {
//         cofaces->push_back(cof);
//     }

//     return true;
// }







