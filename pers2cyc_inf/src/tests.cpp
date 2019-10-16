#include "utils.h"
#include "label_id_map.hpp"
#include "arr_label_id_map.hpp"
#include "inf_pers2cyc.h"

// void testCubeComplex() {
//     CubeComplex* complex = new CubeComplex();
//     complex->init(1,3);
//     complex->setCofaceCountHint(2, 2);

//     complex->addCell(3, {0,1,2,3,4,5,6,7});
//     complex->addCell(1, {7,4});
//     complex->addCell(2, {4,7,17,16});
//     complex->addCell(1, {4,12});
//     complex->addCell(1, {5,10});
//     complex->addCell(2, {5,10,11,6});
//     complex->addCell(3, {1,8,9,2,5,10,11,6});
//     complex->addCell(2, {5,10,11,6});
//     complex->addCell(2, {4,5,13,12});
//     complex->addCell(3, {4,5,6,7,12,13,14,15});
//     complex->addCell(2, {1,2,6,5});
//     complex->addCell(2, {5,10,11,6});
//     complex->addCell(3, {1,8,9,2,5,10,11,6});
//     complex->addCell(2, {8,9,19,18});

//     complex->print();

//     delete complex;
// }

void printCells3d(CubeComplex* complex, const vector<int>& cell_ids) {
    int x = -1, y = -1, z = -1;

    for (auto cell_id : cell_ids) {
        complex->getCell3dCoord(cell_id, &x, &y, &z);
        cout << "    " << x << ',' << y << ',' << z << endl;
    }

    cout << endl;
}

void printCells2d(CubeComplex* complex, const vector<int>& cell_ids) {
    int x = -1, y = -1, z = -1;

    for (auto cell_id : cell_ids) {
        if (cell_id < complex->cell_2d_xy_size_) {
            complex->getCell2dCoord_XY(cell_id, &x, &y, &z);
            cout << "    xy " << x << ',' << y << ',' << z << endl;
        } else if (cell_id < complex->cell_2d_xy_size_ + complex->cell_2d_xz_size_) {
            complex->getCell2dCoord_XZ(cell_id, &x, &y, &z);
            cout << "    xz " << x << ',' << y << ',' << z << endl;
        } else if (cell_id < complex->cell_size_[2]) {
            complex->getCell2dCoord_YZ(cell_id, &x, &y, &z);
            cout << "    yz " << x << ',' << y << ',' << z << endl;
        }
    
    }

    cout << endl;
}

void testBMCubeComplex() {
    vector<int> res = {5, 5, 5};

    {
        CubeComplex *cube_complex = new CubeComplex;
        cube_complex->setResolution(res[0], res[1], res[2]);
    
        for (int id = 0; id < cube_complex->cell_1d_x_size_; id ++) {
            cube_complex->cell_in_complex_[1][id] = true;
            cube_complex->cell_cnt_[1] ++;
        }
    
        cube_complex->writeEdges(std::string("x_edges.ply"));
        delete cube_complex;
    }

    {
        CubeComplex *cube_complex = new CubeComplex;
        cube_complex->setResolution(res[0], res[1], res[2]);
    
        for (int id = cube_complex->cell_1d_x_size_; 
            id < cube_complex->cell_1d_x_size_ + cube_complex->cell_1d_y_size_; 
            id ++) {

            cube_complex->cell_in_complex_[1][id] = true;
            cube_complex->cell_cnt_[1] ++;
        }
    
        cube_complex->writeEdges(std::string("y_edges.ply"));
        delete cube_complex;
    }

    {
        CubeComplex *cube_complex = new CubeComplex;
        cube_complex->setResolution(res[0], res[1], res[2]);
    
        for (int id = cube_complex->cell_1d_x_size_ + cube_complex->cell_1d_y_size_; 
            id < cube_complex->cell_size_[1]; 
            id ++) {

            cube_complex->cell_in_complex_[1][id] = true;
            cube_complex->cell_cnt_[1] ++;
        }
    
        cube_complex->writeEdges(std::string("z_edges.ply"));
        delete cube_complex;
    }

    {
        CubeComplex *cube_complex = new CubeComplex;
        cube_complex->setResolution(res[0], res[1], res[2]);
    
        for (int id = 0; id < cube_complex->cell_2d_xy_size_; id ++) {
            cube_complex->cell_in_complex_[2][id] = true;
            cube_complex->cell_cnt_[2] ++;
        }
    
        cube_complex->writeMesh(std::string("xy_faces.off"));
        delete cube_complex;
    }

    {
        CubeComplex *cube_complex = new CubeComplex;
        cube_complex->setResolution(res[0], res[1], res[2]);
    
        for (int id = cube_complex->cell_2d_xy_size_; 
            id < cube_complex->cell_2d_xy_size_+cube_complex->cell_2d_xz_size_; 
            id ++) {
            cube_complex->cell_in_complex_[2][id] = true;
            cube_complex->cell_cnt_[2] ++;
        }
    
        cube_complex->writeMesh(std::string("xz_faces.off"));
        delete cube_complex;
    }

    {
        CubeComplex *cube_complex = new CubeComplex;
        cube_complex->setResolution(res[0], res[1], res[2]);
    
        for (int id = cube_complex->cell_2d_xy_size_+cube_complex->cell_2d_xz_size_; 
            id < cube_complex->cell_size_[2]; id ++) {
            cube_complex->cell_in_complex_[2][id] = true;
            cube_complex->cell_cnt_[2] ++;
        }
    
        cube_complex->writeMesh(std::string("yz_faces.off"));
        delete cube_complex;
    }

    {
        CubeComplex *cube_complex = new CubeComplex;
        cube_complex->init(1, 3);
        cube_complex->setResolution(res[0], res[1], res[2]);
    
        for (int id = 0; id < cube_complex->cell_size_[3]; id ++) {
            vector<int> verts;
            cube_complex->getCellInternalVerts(3, id, &verts);
            cube_complex->addCellNoConvCanon(3, verts);
        }
    
        cube_complex->writeMesh(std::string("test_faces.off"));

        if (false) {
            vector<int> coord = {2,2,2};
            cout << "  " << coord << endl;
            auto cell_id = cube_complex->getCell1dId_Z(coord[0], coord[1], coord[2]);
            // cout << cell_id << endl;
            printCells2d(cube_complex, *(cube_complex->getCellCofaces(1, cell_id)));
        }

        {
            vector<int> coord = {4,2,2};
            cout << "  " << coord << endl;
            auto cell_id = cube_complex->getCell2dId_YZ(coord[0], coord[1], coord[2]);
            // cout << cell_id << endl;
            printCells3d(cube_complex, *(cube_complex->getCellCofaces(2, cell_id)));
        }

        delete cube_complex;
    }
}

void testLabelIdMap() {
    int verts[][3] = {
        {1,2,3}, {4,5,6},
        {7,8,9}, {10,11,12},
        {4,5,6}, {1,2,3}, 
        {10,11,12}, {7,8,9}, 
        {1,2,3}, {4,5,6},
        {1,2,3}, {4,5,6},
        {7,8,9}, {4,5,6},
        {10,11,12}, {7,8,9},
        {1000,22221,7878}, {7,8,9},
    };

    typedef LabelIdMap<array<int,3>, int, IntArrayHash<3>, IntArrayEqual<3>> IntArrIdMap__;

    IntArrIdMap__ int_arr_id_map;

    for (int i = 0; i < sizeof(verts) / (3*sizeof(int)) ; i ++) {
        array<int,3> v1 = {verts[i][0], verts[i][1], verts[i][2]};
        int_arr_id_map.getId(v1);
    }

    // cout << orien_cell_2d_graph.label_map_.size() << endl;
    // cout << orien_cell_2d_graph.label_pool_.size() << endl;

    auto label_map = int_arr_id_map.label_map_;
    for (auto iter = label_map->begin(); iter != label_map->end(); iter ++) {
        auto p_label_pool = (vector<array<int,3>>*)iter->p;
        cout << iter->id << " " << (*p_label_pool)[iter->id] << endl;
    }

    cout << endl;

    for (auto i = 0; i < int_arr_id_map.label_pool_.size(); i ++) {
        cout << i << " " << int_arr_id_map.label_pool_[i] << endl;
    }

    int_arr_id_map.clearLabelMap();
}

void testGraph() {
    int g[][2] = {
        {1,2},
        {1,3},
        {1,4},
        {1,5},
        {5,4},
        {4,2},
    };

    Graph<int> graph;

    for (int i = 0; i < sizeof(g) / (2*sizeof(int)) ; i ++) {
        graph.addEdge(g[i][0], g[i][1]);
    }

    cout << endl;
    graph.print();
}

void testArrLabelIdMap() {
    vector<int> test_v;
    test_v.reserve(100);
    cout << test_v.capacity() << endl << endl;

    int verts[][3] = {
        {1,2,3}, {4,5,6},
        {7,8,9}, {10,11,12},
        {4,5,6}, {1,2,3}, 
        {10,11,12}, {7,8,9}, 
        {1,2,3}, {4,5,6},
        {1,2,3}, {4,5,6},
        {7,8,9}, {4,5,6},
        {10,11,12}, {7,8,9},
        {1000,22221,7878}, {7,8,9},
        {1001,211,781},
        {1002,211,781},
        {1003,211,781},
        {1004,211,781},
        {1005,211,781},
        {1006,211,781},
        {1007,211,781},
        {1008,211,781},
        {1009,211,781},
        {10010,211,781},
        {10011,211,781},
        {10012,211,781},
        {10013,211,781},
        {1,2,3}, {4,5,6},
        {7,8,9}, {10,11,12},
        {4,5,6}, {1,2,3}, 
        {10,11,12}, {7,8,9}, 
        {1,2,3}, {4,5,6},
        {1,2,3}, {4,5,6},
        {7,8,9}, {4,5,6},
        {10,11,12}, {7,8,9},
        {1000,22221,7878}, {7,8,9},
    };

    typedef ArrayLabelIdMap<int,int> IntArrIdMap;

    IntArrIdMap int_arr_id_map;
    int_arr_id_map.init(3);
    int_arr_id_map.reserve(20);

    for (int i = 0; i < sizeof(verts) / (3*sizeof(int)) ; i ++) {
        vector<int> v(verts[i], verts[i]+3);
        auto id = int_arr_id_map.getId(v);

        if (i == 0) {
            cout << "id size: " << sizeof(id) << endl << endl;
        }

        if (id < 0) {
            id = int_arr_id_map.addArrayLabel(v);
            cout << v << " added with id " << id << endl;
        } else {
            cout << v << " got id " << id << endl;
        }
        
        // int_arr_id_map.printLoadStat();
    }

    cout << endl;
    int_arr_id_map.printLoadStat();
    cout << endl;

    auto map = int_arr_id_map.arr_label_map_;
    for (auto iter = map->begin(); iter != map->end(); iter ++) {
        cout << iter->p << " " << iter->id << endl;
    }

    cout << endl;

    for (auto i = 0; i < int_arr_id_map.arr_label_pool_.size() / 3; i ++) {
        vector<int> v;
        int_arr_id_map.getArrayLabel(i, &v);

        cout << i << " " << v << endl;
    }

    cout << endl << int_arr_id_map.arr_label_pool_ << endl;

    // delete some elements

    int_arr_id_map.deleteArrayLabel(0, -1);
    int_arr_id_map.deleteArrayLabel(10, -1);
    int_arr_id_map.deleteArrayLabel(15, -1);  
    int_arr_id_map.deleteArrayLabel(0, -1);
    int_arr_id_map.deleteArrayLabel(10, -1);
    int_arr_id_map.deleteArrayLabel(15, -1);
    int_arr_id_map.deleteArrayLabel(0, -1);
    int_arr_id_map.deleteArrayLabel(10, -1);
    int_arr_id_map.deleteArrayLabel(15, -1);

    cout << endl << endl << "AFTER DELETION:" << endl << endl;

    int_arr_id_map.printLoadStat();
    cout << endl;

    for (auto iter = map->begin(); iter != map->end(); iter ++) {
        cout << iter->p << " " << iter->id << endl;
    }

    cout << endl;

    for (auto i = 0; i < int_arr_id_map.arr_label_pool_.size() / 3; i ++) {
        vector<int> v;
        int_arr_id_map.getArrayLabel(i, &v);

        cout << i << " " << v << endl;
    }

    cout << endl << int_arr_id_map.arr_label_pool_ << endl;
}

void loadToyCubeComplex(CubeComplex* complex, vector<int>* pos_cell_2d_verts) {
    
    complex->setResolution(3,2,3);
    complex->reserveCellSize(1, 60);
    complex->reserveCellSize(2, 30);
    complex->reserveCellSize(3, 5);

    // Vector3d p;
    // complex->getVertPos(17, &p);
    // cout << p << endl;
    // complex->getVertPos(7, &p);
    // std::cout << p << std::endl;
    // complex->getVertPos(10, &p);
    // std::cout << p << std::endl;
    // return 0;

    complex->addCell(1, {1,29});

    int cell_2d_verts[][4] {
        // faces of {0,6,9,3,1,7,10,4}
        {0,3,9,6}, {1,4,10,7}, {0,1,7,6}, {3,4,10,9}, {0,1,4,3}, {6,7,10,9},
        // faces of {1,7,10,4,2,8,11,5}
        {1,4,10,7}, {2,5,11,8}, {1,2,8,7}, {4,5,11,10}, {1,2,5,4}, {7,8,11,10},
        // faces of {6,12,15,9,7,13,16,10}
        {6,9,15,12}, {7,10,16,13}, {6,7,13,12}, {9,10,16,15}, {6,7,10,9}, {12,13,16,15},

        {1,4,18,19},
        {12,20,21,15},
        {19,18,30,31},
        {33,19,31,32},
        {20,21,25,34},
        {12,20,35,36},
        {15,21,24,37},
        {20,21,39,40},
    };

    for (int i = 0; i < sizeof(cell_2d_verts) / (4*sizeof(int)) ; i ++) {
        vector<int> v(cell_2d_verts[i], cell_2d_verts[i]+4);
        if (complex->getCellId(2, v) < 0) {
            complex->addCell(2, v);
        }
    }

    // complex->addCell(3, {0,6,9,3,1,7,10,4});
    // complex->addCell(3, {1,7,10,4,2,8,11,5});
    // complex->addCell(3, {6,12,15,9,7,13,16,10});
    complex->addCell(3, {21,22,23,24,25,26,27,28});

    *pos_cell_2d_verts = {0,3,9,6};
}










