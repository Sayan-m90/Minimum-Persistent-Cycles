#include "inf_pers2cyc.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <tuple>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/graph_utility.hpp>
#include "utils.h"
#include "mesh_writer.h"
#include "tests.h"
#include "flow_graph.hpp"

const size_t MAX_N_INTERVALS = 1000;

struct CompPersPairsByFiltIndex {
    CompPersPairsByFiltIndex(Bitmap_cubical_complex* _complex) : complex(_complex) {}

    bool operator()(const std::pair<int,int>& p1, const std::pair<int,int>& p2) {
        return complex->key(p1.second) - complex->key(p1.first) >
               complex->key(p2.second) - complex->key(p2.first);
    }

    Bitmap_cubical_complex* complex;
};

struct CompPersPairsByFuncVal {
    CompPersPairsByFuncVal(Bitmap_cubical_complex* _complex) : complex(_complex) {}

    bool operator()(const std::pair<int,int>& p1, const std::pair<int,int>& p2) {
        return complex->filtration(p1.second) - complex->filtration(p1.first) >
               complex->filtration(p2.second) - complex->filtration(p2.first);
    }

    Bitmap_cubical_complex* complex;
};

void computePersistenceCube(
    const std::string& perseus_fname, 
    const IntervSortType interv_sort_type, 
    const ExecOptions& ops,
    std::string* filt_fname) {

    Bitmap_cubical_complex bm_cube_cmplx(perseus_fname.c_str());
    cout << endl << "bitmap complex load done" << endl;

    Persistent_cohomology pcoh(bm_cube_cmplx);
    pcoh.init_coefficients(2);
    pcoh.compute_persistent_cohomology(0.0);
    auto pers_pairs = pcoh.get_persistent_pairs();
    
    vector<std::pair<int,int>> pers_pairs_dim2;
    for (auto p : pers_pairs) {
        if (bm_cube_cmplx.dimension(get<0>(p)) == 2) {
            pers_pairs_dim2.emplace_back(get<0>(p), get<1>(p));
        }
    }

    if (interv_sort_type == BY_FILT_INDEX) {
        CompPersPairsByFiltIndex comp(&bm_cube_cmplx);
        std::sort(pers_pairs_dim2.begin(), pers_pairs_dim2.end(), comp);
    } else {
        CompPersPairsByFuncVal comp(&bm_cube_cmplx);
        std::sort(pers_pairs_dim2.begin(), pers_pairs_dim2.end(), comp);
    }

    cout << endl << "computing persistence done, total intervals: " 
        << pers_pairs_dim2.size() << endl;

    int dim, x_res, y_res, z_res;

    std::ifstream fin(perseus_fname.c_str());
    fin >> dim;
    if (dim != 3) {
        cout << "FATAL: perseus file dim not 3" << endl;
        exit(-1);
    }

    fin >> x_res >> y_res >> z_res;

    std::string purename;
    getFilePurename(perseus_fname, &purename);

    *filt_fname = purename
        + (interv_sort_type == BY_FILT_INDEX ? "_IX" : "_FV")
        + ".filt";

    cout << endl << "writing filt file '" 
        << *filt_fname << "'" << endl;

    writeFiltCube(
        &bm_cube_cmplx, x_res, y_res, z_res, 
        pers_pairs_dim2, *filt_fname);
}

void writeFiltCube(Bitmap_cubical_complex* bm_cube_cmplx, 
    const int x_res, const int y_res, const int z_res, 
    vector<std::pair<int,int>> pers_pairs_dim2,
    const std::string& filename) {

    int num_intv = std::min(pers_pairs_dim2.size(), MAX_N_INTERVALS);

    // cell count of every dimension for each pair.
    // indexed the same as 'pers_pairs_dim2'.
    vector<array<int,4>> dim_cell_count_all(num_intv, {-1, -1, -1, -1});
    vector<array<int,2>> pairs_key_index(num_intv);

    for (auto i = 0; i < num_intv; i ++) {
        // i-th pair in 'pers_pairs_dim2'
        auto p = pers_pairs_dim2[i];
        pairs_key_index[i][0] = bm_cube_cmplx->key(p.first);
        pairs_key_index[i][1] = i;
    }

    struct  {
        bool operator()(const array<int,2>& ki1, const array<int,2>& ki2) {
            return ki1[0] < ki2[0];
        }
    } comp;

    std::sort(pairs_key_index.begin(), pairs_key_index.end(), comp);
    // for (auto ki : pairs_key_index) {
    //     cout << ki[0] << " " << ki[1] << endl;
    // }

    // calculate the cell count of every dimension for all pairs
    array<int,4> dim_cell_count {0, 0, 0, 0};
    int cur_ind = 0;

    for (auto key = 0; 
        key < bm_cube_cmplx->num_simplices() &&
        cur_ind < num_intv; 
        key ++) {

        auto cell_id = bm_cube_cmplx->simplex(key);
        dim_cell_count[bm_cube_cmplx->dimension(cell_id)] ++;

        if (key == pairs_key_index.at(cur_ind).at(0)) {
            auto intv_ind = pairs_key_index.at(cur_ind).at(1);
            dim_cell_count_all.at(intv_ind) = dim_cell_count;
            cur_ind ++;
        }
    }

    if (cur_ind != num_intv) {
        cout << "FATAL: cur_ind != num_intv in writeFiltCube" << endl;
        exit(-1);
    }
 
    std::ofstream fout(filename);

    fout << num_intv << endl;

    for (auto i = 0; i < num_intv; i ++) {
        auto p = pers_pairs_dim2[i];

        fout << i << " "                              // 0: seq of interval
            << dim_cell_count_all.at(i).at(2)         // 1: count of 2 and 3-cells
            + dim_cell_count_all.at(i).at(3) << " "   
            << bm_cube_cmplx->key(p.first) << " "     // 2: start key
            << bm_cube_cmplx->key(p.second) << " "    // 3: end key
            << p.first << " "                         // 4: start bm id 
            << bm_cube_cmplx->key(p.second)           // 5: index length
            - bm_cube_cmplx->key(p.first) << " "
            << bm_cube_cmplx->filtration(p.second)    // 6: function value length
            - bm_cube_cmplx->filtration(p.first) << " "
            << dim_cell_count_all.at(i).at(1) << " "  // 7: 1-cell tount
            << dim_cell_count_all.at(i).at(2) << " "  // 8: 2-cell tount
            << dim_cell_count_all.at(i).at(3) << endl;// 9: 3-cell tount
    }
    
    fout << endl;

    // the given resolution is in terms of 3-d cells (cubes), while
    // the resolution for the filt file as well as for CubeComplex 
    // is in terms of 0-d cells (vertices).
    fout << x_res+1 << " " << y_res+1 << " " << z_res+1 << endl << endl;

    vector<int> verts;
    for (auto key = 0; key <= bm_cube_cmplx->num_simplices(); key ++) {
        auto cell_id = bm_cube_cmplx->simplex(key);
        auto dim = bm_cube_cmplx->dimension(cell_id);

        if (dim == 2 || dim == 3) {
            getBitmapComplexCellVerts(cell_id, dim, 
                x_res, y_res, z_res, &verts);

            fout << verts[0];
            for (auto i = 1; i < verts.size(); i ++) {
                fout << ' ' << verts[i];
            }
            fout << endl;
        }
    }
}

void infPers2CycsFromFile(
    const std::string& filt_fname, 
    const ComplexType cmplx_type, 
    const int start_interval, 
    const int num_intervals, 
    const ExecOptions& ops) {

    int num_all_intervals;

    {
        std::ifstream fin(filt_fname);
        fin >> num_all_intervals;
    }

    cout << endl << "filt file num of intervals: " << num_all_intervals << endl;

    for (int i = start_interval; 
        i < num_all_intervals && i < start_interval + num_intervals; 
        i ++) {

        cout << endl << endl << "---- " << i << "-th interval ----" << endl;

        // this complex will be deleted in infPers2Cyc
        CellComplex* complex;
        if (cmplx_type == CUBE_CMPLX) {
            complex = new CubeComplex();
        }
        complex->init(1,3);

        vector<int> pos_cell_2d_verts;
        int start_key = -1, start_id = -1;

        // loadToyCubeComplex((CubeComplex*)complex, &pos_cell_2d_verts);
        // std::string purename = "hardcode_toy";

        loadComplex(filt_fname, cmplx_type, complex, 
            i, &start_key, &start_id, &pos_cell_2d_verts);

        char buf[50];
        snprintf(buf, 50, "%03d_", i);

        std::string purename;
        getFilePurename(filt_fname, &purename);
        purename = purename + buf + std::to_string(start_key) 
            + "_" + std::to_string(start_id);
        
        if (ops.verbose) {
            cout << endl << "positive cell: " << pos_cell_2d_verts << endl;
        }
        cout << endl << "load complex done" << endl;

        if (ops.write_intem) {
            complex->writeMesh(purename + "_orig.off");
            // ((CubeComplex*)complex)->writeEdges(purename + "_orig_edges.ply");
        }

        infPers2Cyc(cmplx_type, complex, pos_cell_2d_verts, ops, purename);
    }
}

void loadComplex(
    const std::string& filt_fname, 
    const ComplexType cmplx_type, 
    CellComplex* complex, 
    const int interval_id,
    int* interv_start_key,
    int* interv_start_id,
    vector<int>* pos_cell_2d_verts) {

    std::ifstream fin(filt_fname);

    int num_all_intervals, interv_cell_count;
    array<int,4> interv_dim_cell_count;

    fin >> num_all_intervals;

    for (int i = 0; i < num_all_intervals; i ++) {
        int interv_seq, cell_count, start_key, end_key, len_idx, start_id;
        double len_fv;
        array<int,4> dim_cell_count;

        fin >> interv_seq >> cell_count >> start_key >> end_key >> start_id 
            >> len_idx >> len_fv
            >> dim_cell_count[1] >> dim_cell_count[2] >> dim_cell_count[3];

        if (i == interval_id) {
            interv_cell_count = cell_count;
            interv_dim_cell_count = dim_cell_count;
            *interv_start_key = start_key;
            *interv_start_id = start_id;
        }
    }

    if (cmplx_type == CUBE_CMPLX) {
        int x_res, y_res, z_res;
        fin >> x_res >> y_res >> z_res;
        ((CubeComplex*)complex)->setResolution(x_res, y_res, z_res);
    }

    complex->reserveCellSize(1, interv_dim_cell_count[1]);
    complex->reserveCellSize(2, interv_dim_cell_count[2]);
    complex->reserveCellSize(3, interv_dim_cell_count[3]);

    std::string line, s;
    vector<int> verts;

    for (int k = 0; k < interv_cell_count;) {
        if (!std::getline(fin, line)) {
            cout << "FATAL: reach end of file but not finish reading in loadComplex" << endl;
            exit(-1);
        }

        if (line.compare("") == 0) {
            continue;
        }

        std::istringstream iss(line);
        verts.clear();
        while (getline(iss, s, ' ')) {
            verts.push_back(std::stoi(s));
        }

        // cout << verts << endl;

        if (verts.size() == complex->getVertCnt(2)) {
            complex->toCanonVerts(2, &verts);
            complex->addCellNoConvCanon(2, verts);
        } else if (verts.size() == complex->getVertCnt(3)) {
            complex->toCanonVerts(3, &verts);
            complex->addCellNoConvCanon(3, verts);
        } else {
            cout << "FATAL: filt file has cells other than 2d or 3d" << endl;
        }

        k ++;
    }

    if (verts.size() != complex->getVertCnt(2)) {
        cout << "FATAL: positive cell vert count not right in loadComplex" << endl;
        exit(-1);
    }

    *pos_cell_2d_verts = verts;
}

void infPers2Cyc(
    const ComplexType cmplx_type, CellComplex* orig_complex, 
    const vector<int>& pos_cell_2d_verts, const ExecOptions& ops, 
    const std::string& file_prefix) {

    /* pruning */

    pruneComplex(orig_complex);

    // orig_complex->print();
    cout << endl << "prune done" << endl;

    if (ops.verbose) {
        cout << endl << "orig_complex map load:" << endl;
        orig_complex->printLoadStat();
    }

    if (ops.write_intem) {
        orig_complex->writeMesh(file_prefix + "_pruned.off");
    }


    /* get 2-connected component */

    auto start_cell_2d_id = orig_complex->getCellIdNoConvCanon(2, pos_cell_2d_verts);
    if (start_cell_2d_id < 0) {
        cout << "FATAL: start_cell_2d_id < 0 in infPers2Cyc" << endl;
        exit(-1);
    }

    auto orig_complex_cell_1d_cnt = orig_complex->getCellCount(1);
    auto orig_complex_cell_2d_cnt = orig_complex->getCellCount(2);
    auto orig_complex_cell_3d_cnt = orig_complex->getCellCount(3);

    orig_complex->deleteDimCofaces(2);
    orig_complex->deleteDimVertQueryMap(2);

    CellComplex* conn_complex;

    if (cmplx_type == CUBE_CMPLX) {
        conn_complex = new CubeComplex();
        ((CubeComplex*)conn_complex)->setResolution(
            ((CubeComplex*)orig_complex)->getXRes(),
            ((CubeComplex*)orig_complex)->getYRes(),
            ((CubeComplex*)orig_complex)->getZRes());
    }

    conn_complex->init(1,3);
    conn_complex->setCofaceCountHint(2, 2);

    // just some hints
    conn_complex->reserveCellSize(1, orig_complex_cell_1d_cnt/4+1);
    conn_complex->reserveCellSize(2, orig_complex_cell_2d_cnt/4+1);
    conn_complex->reserveCellSize(3, orig_complex_cell_3d_cnt/4+1);
 
    cellConnectedComponent(orig_complex, start_cell_2d_id, conn_complex);

    int conn_compnt_cell_2d_cnt = conn_complex->getCellCount(2);
    
    if (ops.write_intem) {
        conn_complex->writeMesh(file_prefix + "_conn.off");
    }


    /* add necessary 3-cells */

    vector<int> cell_3d_verts;
    vector<int> f_2d_verts;

    for (auto cell_3d_id = 0; 
        cell_3d_id < orig_complex->getCellSize(3); 
        cell_3d_id ++) {

        if (orig_complex->isCellValid(3, cell_3d_id)) {
            orig_complex->getCellVerts(3, cell_3d_id, &cell_3d_verts);
            orig_complex->getFaceVerts(3, cell_3d_verts, 0, &f_2d_verts);

            if (conn_complex->getCellId(2, f_2d_verts) >= 0) {
                conn_complex->addCellNoConvCanon(3, cell_3d_verts);
            }
        }
    }

    if (conn_compnt_cell_2d_cnt != conn_complex->getCellCount(2)) {
        cout << "FATAL: 2-cell count not equal after 3-cells adding" << endl;
        exit(-1);
    }

    cout << endl << "2-conn component done" << endl;

    if (ops.verbose) {
        cout << endl <<  "conn_complex map load:" << endl;
        conn_complex->printLoadStat();
    }

    delete orig_complex;
    orig_complex = NULL;


    /* reconstruct void boundaries */

    vector<array<int,2>>* cell_2d_void_map = 
        new vector<array<int,2>>(conn_complex->getCellSize(2), {-1, -1});
    int void_cnt;

    reconVoidBound(cmplx_type, conn_complex, cell_2d_void_map, ops, file_prefix, &void_cnt);

    cout << endl << "void boundary reconstruction done" << endl;

    // free some memory space

    int src, sink;
    {
        auto pos_cell_2d_id = conn_complex->
            getCellIdNoConvCanon(2, pos_cell_2d_verts);

        vector<int> void_ids;
        getCell2dVoidIds(pos_cell_2d_id, conn_complex, cell_2d_void_map, &void_ids);
        src = void_ids.at(0);
        sink = void_ids.at(1);

        // auto cofaces = conn_complex->getCellCofaces(2, pos_cell_2d_id);
        // src = cofaces->at(0);
        // sink = cofaces->at(1);
    }

    if (src == sink) {
        cout << "FATAL: src sink the same in infPers2Cyc" << endl;
        exit(-1);
    }

    conn_complex->deleteDimCellData(1);
    conn_complex->deleteDimCellData(3);
    conn_complex->deleteDimVertQueryMap(2);


    /* record the corresponding set of 2-cells for each graph edge */

    GEdgeCellMap* gedge_cell_map = new GEdgeCellMap(conn_complex->getCellCount(2)/4+1);

    std::pair<int,int> gedge;
    const vector<int> empty_vector;
    for (auto cell_2d_id = 0;
        cell_2d_id < conn_complex->getCellSize(2);
        cell_2d_id ++) {

        if (conn_complex->isCellValid(2, cell_2d_id)) {
            vector<int> void_ids;
            getCell2dVoidIds(cell_2d_id, conn_complex, cell_2d_void_map, &void_ids);

            if (void_ids.at(0) != void_ids.at(1)) {
                gedge.first = void_ids.at(0);
                gedge.second = void_ids.at(1);
                toCanonEdge(&gedge);

                auto iter = gedge_cell_map->find(gedge);
                if (iter == gedge_cell_map->end()) {
                    std::tie(iter, std::ignore) = 
                        gedge_cell_map->insert({gedge, empty_vector});
                    iter->second.reserve(1);
                    iter->second.push_back(cell_2d_id);
                } else {
                    iter->second.push_back(cell_2d_id);
                }
            }
        }
    }

    if (ops.verbose) {
        cout << endl << "gedge_cell_map: ";
        printHashMapLoad(*gedge_cell_map);
    }
    // cout << endl;
    // printGEdgeCellMap(gedge_cell_map, conn_complex);


    /* compute the min cyc by flow network */

    delete cell_2d_void_map;
    conn_complex->deleteDimCofaces(2);

    if (cmplx_type == CUBE_CMPLX) {
        long max_w = 0;

        for (auto iter = gedge_cell_map->begin(); 
            iter != gedge_cell_map->end(); iter ++) {

            if (iter->second.size() > max_w) {
                max_w = iter->second.size();
            }
        }

        if (ops.verbose) {
            cout << endl << "max_w: " << max_w << endl;
        }

        // TODO: use cpp limits instead
        if (max_w <= 127/2) {
            if (ops.verbose) {
                cout << "use 'char' as weight type" << endl;
            }

            computeMinCyc<char>(
                conn_complex, gedge_cell_map, void_cnt,
                src, sink, ops, file_prefix);
        } else if (max_w <= 32767/2) {
            if (ops.verbose) {
                cout << "use 'short' as weight type" << endl;
            }

            computeMinCyc<short>(
                conn_complex, gedge_cell_map, void_cnt,
                src, sink, ops, file_prefix);
        } else {
            if (ops.verbose) {
                cout << "use 'int' as weight type" << endl;
            }
            
            computeMinCyc<int>(
                conn_complex, gedge_cell_map, void_cnt,
                src, sink, ops, file_prefix);
        }
    }

    cout << endl << "compute min cycle done" << endl;

    delete gedge_cell_map;
    delete conn_complex;
}

void cellConnectedComponent(
    CellComplex* in_complex, const int start_cell_2d_id, CellComplex* out_complex) {

    vector<bool> cell_2d_visited(in_complex->getCellSize(2), false);

    cell_2d_visited[start_cell_2d_id] = true;
    vector<int> cell_stack { start_cell_2d_id };

    // vertices container for the 2-cell being deleted
    vector<int> cell_2d_verts;
    // vertices container for 1-faces of the 2-cell being deleted
    vector<int> f_verts;

    while (!cell_stack.empty()) {
        auto cell_2d_id = cell_stack[cell_stack.size() - 1];
        cell_stack.pop_back();
        in_complex->getCellVerts(2, cell_2d_id, &cell_2d_verts);

#ifdef DEBUG_OPT
        if (out_complex->getCellIdNoConvCanon(2, cell_2d_verts) >= 0) {
            cout << "FATAL: cell in out_complex" << endl;
            exit(-1);
        }
#endif

        out_complex->addCellNoConvCanon(2, cell_2d_verts);

        for (auto i = 0; i < in_complex->getFaceCnt(2); i ++) {
            in_complex->getFaceVerts(2, cell_2d_verts, i, &f_verts);
            auto face_id = in_complex->getCellId(1, f_verts);
            const vector<int>* cofaces = in_complex->getCellCofaces(1, face_id);

            for (auto j = 0; j < cofaces->size(); j ++) {
#ifdef DEBUG_OPT
                if ( !in_complex->isCellValid( 2, cofaces->at(j) ) ) {
                    cout << "FATAL: coface invalid in cellConnectedComponent" << endl;
                    exit(-1);
                }
#endif

                if (!cell_2d_visited[ cofaces->at(j) ]) {
                    cell_stack.push_back( cofaces->at(j) );
                    cell_2d_visited[ cofaces->at(j) ] = true;
                }
            }
        }
    }
}

void pruneComplex(CellComplex* complex) {

    std::unordered_set<int> one_cells_w1cof;
    for (auto cell_1d_id = 0; cell_1d_id < complex->getCellSize(1); cell_1d_id ++) {
        if (complex->isCellValid(1, cell_1d_id)) {
            if (complex->getCellCofaces(1, cell_1d_id)->size() == 1) {
                one_cells_w1cof.insert(cell_1d_id);
            }
        }
    }

    // vertices container for the 2-cell being deleted
    vector<int> cell_2d_verts;
    // vertices container for 1-faces of the 2-cell being deleted
    vector<int> f_verts;

    while (!one_cells_w1cof.empty()) {
        auto cell_1d_id = *(one_cells_w1cof.begin());
        // the 2-cell being deleted
        auto cell_2d_id = complex->getCellCofaces(1, cell_1d_id)->at(0);

#ifdef DEBUG_OPT
        if (complex->getCellCofaces(1, cell_1d_id)->size() != 1) {
            cout << "FATAL: 1d cell has cofaces not 1" << endl;
            exit(-1);
        }
        if (complex->getCellCofaces(2, cell_2d_id)->size() != 0) {
            cout << "FATAL: 2d cell has coface in prune" << endl;
            exit(-1);
        }
#endif
 
        complex->getCellVerts(2, cell_2d_id, &cell_2d_verts);
        complex->deleteCell(2, cell_2d_id);

        // enumerate all faces of the 2-cell being deleted
        for (auto i = 0; i < complex->getFaceCnt(2); i ++) {
            complex->getFaceVerts(2, cell_2d_verts, i, &f_verts);
            auto face_id = complex->getCellId(1, f_verts);

            complex->deleteCoface(1, face_id, cell_2d_id);

            if (complex->getCellCofaces(1, face_id)->size() == 0) {
                one_cells_w1cof.erase(face_id);
            } else if (complex->getCellCofaces(1, face_id)->size() == 1) {                 
                one_cells_w1cof.insert(face_id);
            }
        }
    }
}

void get2DCellPtCubic(
    CellComplex* complex, const int cell_2d_id, 
    const vector<int>& cell_1d_verts, Vector3d* pt) {

    vector<int> cell_2d_verts;
    complex->getCellVerts(2, cell_2d_id, &cell_2d_verts);

    int verts[2];
    int j = 0;
    for (auto v : cell_2d_verts) {
        if (v != cell_1d_verts[0] && v != cell_1d_verts[1]) {
            verts[j] = v;
            j ++;
        }
    }

#ifdef DEBUG_OPT
    if (j != 2) {
        cout << "FATAL: j!=2 in get2DCellPtCubic" << endl;
        exit(-1);
    }
#endif

    Vector3d pt1, pt2;
    complex->getVertPos(verts[0], &pt1);
    complex->getVertPos(verts[1], &pt2);
    *pt = (pt1+pt2) / 2;
}

// get the canonical oriented 2d cell vertex sequence where the orientation is
// from the from_vert to the to_vert. the canonical vertex sequence then starts
// from the smallest vertex following the orientation.
// the input 'verts' is assumed to be in canon order.
void getCanonOrien2DCellVertsCubic(const vector<int>& verts, 
    const int& from_vert, const int& to_vert, vector<int>* canon_verts) {

    canon_verts->clear();

    int inc;
    for (auto i = 0; i < verts.size(); i ++) {
        if (verts[i] == from_vert) {
            if (verts[ mod(i+1, verts.size()) ] == to_vert) {
                inc = 1;
            } else {
                inc = -1;
            }
            
            break;
        }
    }

    // int min_ind = 0;
    // int min_vert = verts[0];

    // for (auto i = 1; i < verts.size(); i ++) {
    //     if (verts[i] < min_vert) {
    //         min_ind = i;
    //         min_vert = verts[i];
    //     }
    // }

    for (auto i = 0; i < verts.size(); i ++) {
        // canon_verts->push_back( verts[ mod(min_ind+i*inc, verts.size()) ] );
        canon_verts->push_back( verts[ mod(i*inc, verts.size()) ] );
    }
}

void addOrien2DCellPairs(
    const ComplexType cmplx_type, 
    CellComplex* complex, 
    const int cell_1d_id,
    const vector<int>& cell_1d_verts, 
    VertArrIdMap* orien_cell_2d_id_map, 
    Graph<int>* orien_cell_2d_graph) {

    // if (cell_1d->cofaces.size() == 2) {
    //     int dosomething = 0;
    // }

    const vector<int>* cell_1d_cofaces = complex->getCellCofaces(1, cell_1d_id);

    Vector3d cell_1d_pt0, cell_1d_pt1;
    complex->getVertPos(cell_1d_verts[0], &cell_1d_pt0);
    complex->getVertPos(cell_1d_verts[1], &cell_1d_pt1);

    // cout << endl << cell_1d_verts << endl;
    // cout << cell_1d_pt0 << endl;
    // cout << cell_1d_pt1 << endl;

    Vector3d cell_1d_mid_pt = (cell_1d_pt0+cell_1d_pt1) / 2;
    // cout << cell_1d_mid_pt << endl;


    /* get the plane equation of the first coface */

    auto first_cell_2d_id = cell_1d_cofaces->at(0);
    vector<int> first_cell_2d_verts;
    complex->getCellVerts(2, first_cell_2d_id, &first_cell_2d_verts);

    Vector3d first_cell_2d_pt0, first_cell_2d_pt1, first_cell_2d_pt2;
    complex->getVertPos(first_cell_2d_verts[0], &first_cell_2d_pt0);
    complex->getVertPos(first_cell_2d_verts[1], &first_cell_2d_pt1);
    complex->getVertPos(first_cell_2d_verts[2], &first_cell_2d_pt2);

    Vector4d first_cell_plane_eq;
    getPlaneEquation(first_cell_2d_pt0, first_cell_2d_pt1, first_cell_2d_pt2, &first_cell_plane_eq);
    // first_cell_plane_eq = -first_cell_plane_eq;
    // cout << first_cell_plane_eq << endl;
    // cout << evalPlaneEquation(first_cell_plane_eq, first_cell_2d_pt0) << endl;
    // cout << evalPlaneEquation(first_cell_plane_eq, first_cell_2d_pt1) << endl;
    // cout << evalPlaneEquation(first_cell_plane_eq, first_cell_2d_pt2) << endl;


    /* get the angles of other cofaces w.r.t the first coface */

    // TODO: this vector seems useless
    vector<Vector3d> cell_2d_pts(cell_1d_cofaces->size());
    vector<std::pair<Decimal,int>> cell_2d_order(cell_1d_cofaces->size());

    Vector3d first_cell_2d_vec;
    for (auto i = 0; i < cell_1d_cofaces->size(); i ++) {
        if (cmplx_type == CUBE_CMPLX) {
            get2DCellPtCubic(complex, cell_1d_cofaces->at(i), cell_1d_verts, &cell_2d_pts[i]);
        } else if (cmplx_type == SIMP_CMPLX) {
            int dosomething = 0;
        }

        if (i == 0) {
            cell_2d_order[i].first = 0;
            first_cell_2d_vec = cell_2d_pts[i] - cell_1d_mid_pt;
            first_cell_2d_vec.normalize();
        } else {
            Vector3d cell_2d_vec = cell_2d_pts[i] - cell_1d_mid_pt;
            cell_2d_vec.normalize();

            cell_2d_order[i].first = acos(first_cell_2d_vec.dot(cell_2d_vec));
            if (evalPlaneEquation(first_cell_plane_eq, cell_2d_pts[i]) < 0) {
                cell_2d_order[i].first = 2*M_PI - cell_2d_order[i].first;
            }
        }

        cell_2d_order[i].second = i;
    }

    struct {
        bool operator()(const std::pair<Decimal,int> &a, 
            const std::pair<Decimal,int> &b) const {   
            return a.first < b.first;
        }
    } comp;

    std::sort(cell_2d_order.begin(), cell_2d_order.end(), comp);

#ifdef DEBUG_OPT
    if (cell_2d_order[0].second != 0) {
        cout << "FATAL: cell_2d_order[0].second != 0 in addOrien2DCellPairs" << endl;
        exit(-1);
    }
#endif

    // for (auto pair : cell_2d_order) {
    //     cout << radian2Degree(pair.first) << ' ';
    // }
    // cout << endl;

    Vector3d cell_1d_vec = cell_1d_pt1 - cell_1d_pt0;
    cell_1d_vec.normalize();
    Vector3d cell_1d_ortho_pt = cell_1d_mid_pt + cell_1d_vec.cross(first_cell_2d_vec);

    bool cell_1d_0_to_1 = true;
    if (evalPlaneEquation(first_cell_plane_eq, cell_1d_ortho_pt) < 0) {
        cell_1d_0_to_1 = false;
    }

    // cout << "cell_1d_0_to_1: " << cell_1d_0_to_1 << endl;
    
    vector<int> cell_2d_verts1, cell_2d_verts2, canon_verts1, canon_verts2;
    canon_verts1.reserve(complex->getVertCnt(2));
    canon_verts2.reserve(complex->getVertCnt(2));

    for (auto i = 0; i < cell_2d_order.size(); i ++) {
        if (cmplx_type == CUBE_CMPLX) {
            auto cell_2d_order1 = cell_2d_order[i].second;
            auto cell_2d_order2 = cell_2d_order[ mod(i+1, cell_2d_order.size()) ].second;

            auto cell_2d_id1 = cell_1d_cofaces->at(cell_2d_order1);
            auto cell_2d_id2 = cell_1d_cofaces->at(cell_2d_order2);

            auto cell_2d_cofaces1 = complex->getCellCofaces(2, cell_2d_id1);
            auto cell_2d_cofaces2 = complex->getCellCofaces(2, cell_2d_id2);

            // check whether the two oriented 2-cells enclose a 3-cell
            if (vecIntersect(*cell_2d_cofaces1, *cell_2d_cofaces2)) {
                if (cell_1d_cofaces->size() == 2) {
                    auto degree1 = cell_2d_order[i].first;
                    auto degree2 = cell_2d_order[ mod(i+1, cell_2d_order.size()) ].first;

                    // cout << "deg:" << radian2Degree(degree1) << " " << radian2Degree(degree2) 
                        // << " " << radian2Degree(radianSub(degree2, degree1)) << endl;

                    if (radianSub(degree2, degree1) <= M_PI) {
                        continue;
                    }
                } else {
                    continue;
                }
            }

            complex->getCellVerts(2, cell_2d_id1, &cell_2d_verts1);
            complex->getCellVerts(2, cell_2d_id2, &cell_2d_verts2);

            if (cell_1d_0_to_1) {
                // CAUTION: should guarantee the input vertices in canon order
                getCanonOrien2DCellVertsCubic(
                    cell_2d_verts1, cell_1d_verts[0], cell_1d_verts[1], &canon_verts1);
                getCanonOrien2DCellVertsCubic(
                    cell_2d_verts2, cell_1d_verts[1], cell_1d_verts[0], &canon_verts2);
            } else {
                getCanonOrien2DCellVertsCubic(
                    cell_2d_verts1, cell_1d_verts[1], cell_1d_verts[0], &canon_verts1);
                getCanonOrien2DCellVertsCubic(
                    cell_2d_verts2, cell_1d_verts[0], cell_1d_verts[1], &canon_verts2);
            }

            // cout << cell_2d_verts1 << endl;
            // cout << "paired 2-cells: " << canon_verts1 << ' ' << canon_verts2 << endl;

            auto id1 = orien_cell_2d_id_map->getId(canon_verts1);
            if (id1 < 0) {
                id1 = orien_cell_2d_id_map->addArrayLabel(canon_verts1);
            }
            
            auto id2 = orien_cell_2d_id_map->getId(canon_verts2);
            if (id2 < 0) {
                id2 = orien_cell_2d_id_map->addArrayLabel(canon_verts2);
            }
            
            orien_cell_2d_graph->addEdge(id1, id2);
        }
    }
}

void reconVoidBound(
    const ComplexType cmplx_type, 
    CellComplex* complex, 
    vector<array<int,2>>* cell_2d_void_map,
    const ExecOptions& ops, 
    const std::string& file_prefix,
    int* void_cnt) {

    VertArrIdMap orien_cell_2d_id_map;
    orien_cell_2d_id_map.init(complex->getVertCnt(2));
    // cout << "complex->getVertCnt(2): " << complex->getVertCnt(2) << endl;

    Graph<int> orien_cell_2d_graph;


    /* get the number of boundary oriented 2-cells 
       and reserve the size for the graph and the id map */

    int orien_cell_2d_cnt = 0;
    for (auto cell_2d_id = 0; cell_2d_id < complex->getCellSize(2); cell_2d_id ++) {
        if (complex->isCellValid(2, cell_2d_id)) {
            orien_cell_2d_cnt += 2 - complex->getCellCofaces(2, cell_2d_id)->size();
        }
    }

    orien_cell_2d_id_map.reserve(orien_cell_2d_cnt);
    orien_cell_2d_graph.reserve(orien_cell_2d_cnt);


    /* traverse all boundary 2-cells and their 1-faces 
       to add all pairs of orented 2-cells */
    // TODO: traverse all 1-cells instead

    std::unordered_set<int>* visited_cells_1d = new std::unordered_set<int>;
    vector<int> cell_1d_verts, cell_2d_verts;

    for (auto cell_2d_id = 0; cell_2d_id < complex->getCellSize(2); cell_2d_id ++) {
        if (complex->isCellValid(2, cell_2d_id)) {
            if (complex->getCellCofaces(2, cell_2d_id)->size() < 2) {
                complex->getCellVerts(2, cell_2d_id, &cell_2d_verts);

                for (auto i = 0; i < complex->getFaceCnt(2); i ++) {
                    complex->getFaceVerts(2, cell_2d_verts, i, &cell_1d_verts);
                    complex->toCanonVerts(1, &cell_1d_verts);
                    // cout << "  " << cell_1d_verts << endl;
                    auto cell_1d_id = complex->getCellIdNoConvCanon(1, cell_1d_verts);

                    if (visited_cells_1d->find(cell_1d_id) == visited_cells_1d->end()) {
                        visited_cells_1d->insert(cell_1d_id);
                        addOrien2DCellPairs(cmplx_type, complex, cell_1d_id, cell_1d_verts,
                            &orien_cell_2d_id_map, &orien_cell_2d_graph);
                    }
                }
            }
        }
    }

    delete visited_cells_1d;
    visited_cells_1d = NULL;

    if (ops.verbose) {
        cout << endl << "orien_cell_2d_id_map size stat:" << endl;
        orien_cell_2d_id_map.printLoadStat();
        cout << endl << "orien_cell_2d_graph: capacity=" 
            << orien_cell_2d_graph.adj_nodes_.capacity() 
            << ",size=" << orien_cell_2d_graph.adj_nodes_.size() << endl;
    }

    // TODO: maybe can add a degree check for the graph: 4 for cubical complex?
    // orien_cell_2d_graph.print();
    // cout << endl;


    /* get connected components of the graph of oriented 2-cells by DFS */

    int void_id = 0;
    vector<bool> orien_cell_2d_visited;
    orien_cell_2d_visited.resize(orien_cell_2d_graph.size(), false);
    vector<int> node_stack;
    vector<int> orien_cell_2d;
    
    for (auto i = 0; i < orien_cell_2d_graph.size(); i ++) {
        if (orien_cell_2d_visited[i] == false) {
            node_stack = { i };
            orien_cell_2d_visited[i] = true;

            MeshWriter* mesh_writer;
            if (ops.write_intem) {
                mesh_writer = new MeshWriter(complex, complex->getVertCnt(2));
            }
            
            // cout << endl << "node_stack.capacity(): " << node_stack.capacity() << endl << endl;
            // cout << "void " << void_id << ": (";
            // cout << "void " << void_id << ":" << endl;

            while (!node_stack.empty()) {
                auto node = node_stack[node_stack.size() - 1];
                node_stack.pop_back();

                // cout << node << ",";
                orien_cell_2d_id_map.getArrayLabel(node, &orien_cell_2d);
                // cout << orien_cell_2d << endl;

                if (ops.write_intem) {
                    mesh_writer->addFace(orien_cell_2d);
                }

                complex->toCanonVerts(2, &orien_cell_2d);
                auto cell_2d_id = complex->getCellIdNoConvCanon(2, orien_cell_2d);

                if (cell_2d_void_map->at(cell_2d_id).at(0) < 0) {
                    cell_2d_void_map->at(cell_2d_id).at(0) = void_id + complex->getCellSize(3);
                } else if (cell_2d_void_map->at(cell_2d_id).at(1) < 0) {
                    cell_2d_void_map->at(cell_2d_id).at(1) = void_id + complex->getCellSize(3);
                } else {
                    cout << "FATAL cell_2d_void_map[cell_2d_id] both >= 0 in reconVoidBound" << endl;
                    exit(-1);
                }

                // auto cofaces = const_cast<vector<int>*>(complex->getCellCofaces(2, cell_2d_id));
                // cofaces->push_back(void_id + complex->getCellSize(3));

                for (auto j = 0; j < orien_cell_2d_graph.adj_nodes_[node].size(); j ++) {
                    auto adj_node = orien_cell_2d_graph.adj_nodes_[node][j];

                    if (orien_cell_2d_visited[adj_node] == false) {
                        node_stack.push_back(adj_node);
                        orien_cell_2d_visited[adj_node] = true;
                    }
                }
            }

            if (ops.write_intem) {
                mesh_writer->write(file_prefix + "_void" + std::to_string(void_id) + ".off");
                delete mesh_writer;
            }

            void_id ++;
            // cout << ")" << endl;
            // cout << endl;
        }
    }

    *void_cnt = complex->getCellSize(3) + void_id;

    if (ops.verbose) {
        cout << endl << "void count: " << void_id << endl;
    }
}

template<typename WeightType>
void computeMinCyc(
    CellComplex* complex, 
    const GEdgeCellMap* gedge_cell_map, 
    const int void_cnt,
    const int src,
    const int sink, 
    const ExecOptions& ops, 
    const std::string& file_prefix) {

    // if (ops.verbose) {
    //     cout << "vertex size:" << void_cnt + gedge_cell_map->size() << endl;
    // }

    // FlowGraph<int,WeightType> _graph(void_cnt + gedge_cell_map->size());

    // int edge_id = 0;
    // for (auto iter = gedge_cell_map->begin(); 
    //     iter != gedge_cell_map->end(); iter ++) {

    //     const std::pair<int,int>& edge = iter->first;

    //     WeightType w = 0;
    //     for (auto cell_2d_id : iter->second) {
    //         w += (WeightType)(complex->getWeight(2, cell_2d_id));
    //     }

    //     _graph.addEdge(edge.first, edge.second, w);
    //     _graph.addEdge(edge.second, void_cnt + edge_id, w);
    //     _graph.addEdge(void_cnt + edge_id, edge.first, w);

    //     edge_id ++;
    // }

    // return;


    /* compute the maximal flow */

    typedef boost::adjacency_list_traits < boost::vecS, boost::vecS, boost::directedS > Traits;

    typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS,
    // boost::property < boost::vertex_name_t, std::string,
    boost::property < boost::vertex_index_t, int, // originally it's 'long'
    boost::property < boost::vertex_color_t, boost::default_color_type,
    boost::property < boost::vertex_distance_t, int, // originally it's 'long'
    boost::property < boost::vertex_predecessor_t, Traits::edge_descriptor > > > >,

    boost::property < boost::edge_capacity_t, WeightType,
    boost::property < boost::edge_residual_capacity_t, WeightType,
    boost::property < boost::edge_reverse_t, Traits::edge_descriptor > > > > FlowGraphBoost;

    if (ops.verbose) {
        cout << endl << "graph elem sizes: " 
            << sizeof(Traits::vertex_descriptor) << " " 
            << sizeof(Traits::edge_descriptor) << endl;
    }

    FlowGraphBoost* graph = new FlowGraphBoost;
    typename boost::property_map<FlowGraphBoost, boost::edge_capacity_t>::type
        capacity = get(boost::edge_capacity, *graph);
    typename boost::property_map<FlowGraphBoost, boost::edge_residual_capacity_t>::type
        residual_capacity = get(boost::edge_residual_capacity, *graph);
    typename boost::property_map<FlowGraphBoost, boost::edge_reverse_t>::type 
        reverse_edge = get(boost::edge_reverse, *graph);

    for (auto iter = gedge_cell_map->begin(); 
        iter != gedge_cell_map->end(); iter ++) {

        const std::pair<int,int>& edge = iter->first;
        Traits::edge_descriptor e1, e2;

        boost::tie(e1, boost::tuples::ignore) 
            = add_edge(edge.first, edge.second, *graph);
        boost::tie(e2, boost::tuples::ignore) 
            = add_edge(edge.second, edge.first, *graph);

        WeightType w = 0;
        for (auto cell_2d_id : iter->second) {
            w += (WeightType)(complex->getWeight(2, cell_2d_id));
        }

        // cout << edge << " w=" << w << endl;

        capacity[e1] = w;
        capacity[e2] = w;
        reverse_edge[e1] = e2;
        reverse_edge[e2] = e1;
    }
    
    cout << endl << "build flow network done" << endl;
 
    auto max_flow = boykov_kolmogorov_max_flow(*graph ,src, sink);
    if (ops.verbose) {
        cout << endl << "max_flow: " << (double)max_flow << endl;
    }


    /* do a DFS to get the min-cut */

    if (ops.verbose) {
        cout << "graph num_vertices: " << num_vertices(*graph) << endl;
    }

    vector<bool> visited;

    {
        visited.resize(num_vertices(*graph), false);
        visited[src] = true;
        vector<int> vert_stack = { src };

        while (!vert_stack.empty()) {
            auto v = vert_stack.back();
            vert_stack.pop_back();

            typename boost::graph_traits<FlowGraphBoost>::out_edge_iterator e_iter, e_end;
            boost::tie(e_iter, e_end) = out_edges(v, *graph);
            for (; e_iter != e_end; e_iter ++) {
                if (residual_capacity[*e_iter] > 0) {
                    int adj_v = target(*e_iter, *graph);
                    if (!visited[adj_v]) {
                        vert_stack.push_back(adj_v);
                        visited[adj_v] = true;
                    }
                }
            }
        }
    }

    delete graph;

    
    /* collect edges across the min-cut */

    MeshWriter mesh_writer(complex, complex->getVertCnt(2));
    vector<int> cell_2d_verts;

    for (auto iter = gedge_cell_map->begin(); 
        iter != gedge_cell_map->end(); iter ++) {

        const std::pair<int,int>& edge = iter->first;
        if (visited[edge.first] != visited[edge.second]) {
            for (auto cell_2d_id : iter->second) {
                complex->getCellVerts(2, cell_2d_id, &cell_2d_verts);
                mesh_writer.addFace(cell_2d_verts);
            }
        }
    }

    mesh_writer.write(file_prefix + "_mincyc.off");
}

void getBitmapComplexCellVerts(const int cell_handle, const int cell_dim, 
    const int x_res, const int y_res, const int z_res, vector<int>* verts) {

    const int x_cnt = 2*x_res + 1;
    const int y_cnt = 2*y_res + 1;

    const int z_coord = cell_handle / (x_cnt*y_cnt);
    const int left = cell_handle % (x_cnt*y_cnt);
    const int y_coord = left / x_cnt;
    const int x_coord = left % x_cnt;

    int vari_cnt = 0;

    int x_vari = 1;
    if (x_coord % 2 == 1) {
        x_vari = 2;
        vari_cnt ++;
    }

    int y_vari = 1;
    if (y_coord % 2 == 1) {
        y_vari = 2;
        vari_cnt ++;
    }

    int z_vari = 1;
    if (z_coord % 2 == 1) {
        z_vari = 2;
        vari_cnt ++;
    }

#ifdef DEBUG_OPT
    if (vari_cnt != cell_dim) {
        cout << "FATAL: getBitmapCC2dCellVerts dim not equal" << endl;
        exit(-1);
    }
#endif

    vector<int> temp_verts;
    for (auto x = x_coord/2; x < x_coord/2 + x_vari; x ++) {
        for (auto y = y_coord/2; y < y_coord/2 + y_vari; y ++) {
            for (auto z = z_coord/2; z < z_coord/2 + z_vari; z ++) {
                temp_verts.push_back(x + y*(x_res+1) + z*(x_res+1)*(y_res+1));
            }
        }
    }

    if (cell_dim == 2) {
        *verts = {temp_verts[0], temp_verts[1], temp_verts[3], temp_verts[2]};
    } else if (cell_dim == 3) {
        *verts = {temp_verts[0], temp_verts[1], temp_verts[3], temp_verts[2],
            temp_verts[4], temp_verts[5], temp_verts[7], temp_verts[6]};
    }
}

void printGEdgeCellMap(const GEdgeCellMap* gedge_cell_map, CellComplex* complex) {
    cout << "3-cells [" << complex->getCellCount(3) << "]:" << endl;
    for (int cell_3d_id = 0; cell_3d_id < complex->getCellSize(3); cell_3d_id ++) {
        if (complex->isCellValid(3, cell_3d_id)) {
            vector<int> verts;
            complex->getCellVerts(3, cell_3d_id, &verts);
            cout << cell_3d_id << ": " << verts << endl;
        }
    }
    cout << endl;

    for (auto iter = gedge_cell_map->begin(); 
        iter != gedge_cell_map->end(); iter ++) {

        cout << iter->first << ":" << endl;
        for (auto cell_2d_id : iter->second) {
            vector<int> verts;
            complex->getCellVerts(2, cell_2d_id, &verts);
            cout << "  " << verts << endl;
        }
    }
}

void getCell2dVoidIds(
    const int cell_2d_id, CellComplex* complex, 
    vector<array<int,2>>* cell_2d_void_map, vector<int>* void_ids) {

    void_ids->clear();

    auto cofaces = complex->getCellCofaces(2, cell_2d_id);
    for (auto cell_3d_id : *cofaces) {
        void_ids->push_back(cell_3d_id);
    }

    for (auto void_id : cell_2d_void_map->at(cell_2d_id)) {
        if (void_id >= 0) {
            void_ids->push_back(void_id);
        }
    }

    if (void_ids->size() != 2) {
        cout << "FATAL: void_ids.size() != 2 in getCell2dVoidIds" << endl;
        exit(-1);
    }
}
  










