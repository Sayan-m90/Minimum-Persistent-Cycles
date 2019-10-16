#include "mesh_writer.h"
#include <fstream>

void MeshWriter::addFace(const vector<int>& face) {
    for (auto v : face) {
        int new_v;

        auto iter = vert_ind_map_.find(v);
        if (iter == vert_ind_map_.end()) {
            new_v = vert_ind_map_.size();
            vert_ind_map_.insert({v, new_v});
            vert_pool_.push_back(v);
        } else {
            new_v = iter->second;
        }

        face_pool_.push_back(new_v);

        // cout << "#" << v << "," << new_v << endl;
    }
}

void MeshWriter::write(const std::string& filename) {
    std::ofstream fout(filename);

    fout << "OFF" << endl;
    fout << vert_ind_map_.size() << ' ' << 
        face_pool_.size() / face_vert_cnt_ << " 0" << endl;

    Vector3d pos;
    for (auto v : vert_pool_) {
        complex_->getVertPos(v, &pos);
        fout << pos[0] << ' ' << pos[1] << ' ' << pos[2] << endl;
    }

    for (auto i = 0; i < face_pool_.size(); i += face_vert_cnt_) {
        fout << face_vert_cnt_;
        for (auto j = 0; j < face_vert_cnt_; j ++) {
            fout << ' ' << face_pool_[i + j];
        }
        fout << endl;
    }
}
