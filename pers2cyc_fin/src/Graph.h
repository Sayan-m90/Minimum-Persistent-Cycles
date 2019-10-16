#include <iostream>
#include <vector>
#include <queue>
#include <limits.h>
#include <string.h>
#include <algorithm>
#include "read_boost_nr.h"
#include <queue>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
using namespace std;
typedef Gudhi::cubical_complex::Bitmap_cubical_complex_base<double> Bitmap_cubical_complex_base;
typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base> Bitmap_cubical_complex;
typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;

using Simplex_tree = Gudhi::Simplex_tree<>;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Filtration_value = Simplex_tree::Filtration_value;
using typeVectorVertex = std::vector<Vertex_handle>;
using typePairSimplexBool = std::pair<Simplex_tree::Simplex_handle, bool>;

extern map <int,int> threesimplex_to_node; // maps index of tetrahedrons from simplex->first to dual-graph-nodes->second
extern map <int,int> twosimplex_to_tri; // maps index of traingles from simplex->first to dual-graph-edges->second
extern map <int,int> threesimplex_to_node_rev; // maps index of tetrahedrons from dual-graph-nodes->first to simplex->second
extern map <int,int> twosimplex_to_tri_rev; // maps index of traingles from dual-graph-edges->first to simplex->second
extern map<size_t, vector<int>>  sizescale;//0- edge, 1- vertex, for each deathindex, count of vertex and edges

//vector<int> graphNode;//negative weight: node is red
//vector<vector<int>> adjMatrix;  // negative weight: not present, inf= turned red
vector<double> x;
vector<double> y;
vector<double> z;
map<unsigned long,unsigned long> vertind; // vertex index in cubical complex(first) to actual index in x,y,z, filter
//Size at each scale


vector<vector<int>> boostCut(Bitmap_cubical_complex &cub,
                             size_t birthi, size_t deathi,bool verbose){
    
    using namespace boost;
    typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
    typedef adjacency_list < vecS, vecS, directedS,
    property < vertex_name_t, std::string,
    property < vertex_index_t, long,
    property < vertex_color_t, boost::default_color_type,
    property < vertex_distance_t, long,
    property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,
    
    property < edge_capacity_t, long,
    property < edge_residual_capacity_t, long,
    property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;
    
    Graph g;
    property_map < Graph, edge_capacity_t >::type
    capacity = get(edge_capacity, g);
    property_map < Graph, vertex_color_t >::type
    color = get(vertex_color, g);
    property_map < Graph, edge_residual_capacity_t >::type
    residual_capacity = get(edge_residual_capacity, g);
    property_map < Graph, edge_reverse_t >::type rev = get(edge_reverse, g);
    Traits::vertex_descriptor s, t;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
    std::vector<vertex_descriptor> verts;
    map<int,int> org2collapse;
    
    Graph gshort;
    property_map < Graph, edge_capacity_t >::type
    capacitysh = get(edge_capacity, gshort);
    property_map < Graph, vertex_color_t >::type
    colorsh = get(vertex_color, gshort);
    property_map < Graph, edge_residual_capacity_t >::type
    residual_capacitysh = get(edge_residual_capacity, gshort);
    property_map < Graph, edge_reverse_t >::type revsh = get(edge_reverse, gshort);
    
    
    int job, edge = 0, node = 0, buff = 0;
            
    vector<vector<int>> vectri;
    map<edge_descriptor,vector<int>> edging;
    //    std::filebuf fb;
    //    fb.open ("output",std::ios::in);
    //    std::istream is (&fb);
    //    return detail::read_dimacs_max_flow_internal(g, capacity, reverse_edge, src, sink, in1, true, "max");
    
    
    read_my_dimacs_max_flow(g, capacity, rev, s, t, cub, birthi, deathi);
    //    graph_traits < Graph >::out_edge_iterator ei, e_end;
//    cout<<"Read"; getchar();
    graph_traits < Graph >::out_edge_iterator eibfs, eiibfs, e_endbfs, ee_endbfs;
    graph_traits < Graph >::vertex_iterator u_iter, u_end;
    graph_traits < Graph >::out_edge_iterator ei, e_end;
    // Mark all the vertices as not visited
    bool *visited = new bool[ sizescale[deathi][1]];
    bool *selected = new bool[ sizescale[deathi][1]];
    for(int i = 0; i <  sizescale[deathi][1]; i++){
        visited[i] = false;
        selected[i] = false;
        verts.push_back(add_vertex(gshort));
    }
//    cout<<"INIT"; getchar();
    // Create a queue for BFS
    list<int> queue;
    // Mark the current node as visited and enqueue it
    visited[s] = true;
    queue.push_back(threesimplex_to_node[deathi]);
    if(verbose){
    cout<<"Getting in to copying graph for boost"; getchar();
    }
    
    int countcollapsebfs = 0, countcollapsebfs1 = 0, countcollapsebfs2 = 0, countcollapsebfs3 = 0;
    
    while(!queue.empty())
    {
        // Dequeue a vertex from queue and print it
        auto qf = queue.front();
//        cout << " qf: " << qf << " ";
        queue.pop_front();
        
        for (boost::tie(eibfs, e_endbfs) = out_edges(qf, g); eibfs != e_endbfs; ++eibfs){
            auto targ = target(*eibfs,g);
//            cout<<" targ: "<<targ<<" c: "<<capacity[*eibfs]<<" v: "<<visited[targ]<<" s: "<<selected[targ]<<" last: "<<sizescale[deathi][1]-1<<"\n";
            if (visited[targ] || capacity[*eibfs]!=INT_MAX || selected[targ] || targ == sizescale[deathi][1]-1)
                continue;
            bool cont = false;
            for (boost::tie(eiibfs, ee_endbfs) = out_edges(targ, g); eiibfs != ee_endbfs; ++eiibfs){
//                cout<<" inside "<<*eiibfs<<" c: "<<capacity[*eiibfs]<<" |";
                if(capacity[*eiibfs]!=INT_MAX && capacity[*eiibfs]!=0){
                    cont = true;
                    break;
                }
            }
            visited[targ] = true;
            if(cont==true)
                continue;
//            cout<<"came here ";
            selected[targ] = true;// will be collapsed
            countcollapsebfs++;
            queue.push_back(targ);
            
        }
//        getchar();
        
    }
    if(verbose){
        cout<<"Got out of copying graph. vertices collapsed from bfs:"<<countcollapsebfs; getchar();
    }
    
        verts.push_back(add_vertex(gshort));
    countcollapsebfs = 0;
    
    for (boost::tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter){
        for (boost::tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei){
            if(selected[*u_iter] && selected[target(*ei,g)]
               && *u_iter!=threesimplex_to_node[deathi] && target(*ei,g)!=threesimplex_to_node[deathi]
               && *u_iter!=sizescale[deathi][1]-1 && target(*ei,g)!=sizescale[deathi][1]-1)    // both true and neither are sink or source
            {
                countcollapsebfs++;
                continue;
                
            }   // dont need this edge
            if(selected[*u_iter]==true && *u_iter!=threesimplex_to_node[deathi] && selected[target(*ei,g)]==false && capacity[*ei]!=0){// transfer this node between green and initial
                edge_descriptor e1, e2;
                bool in1, in2;
                int tail, head;
                tail = threesimplex_to_node[deathi];
                head = target(*ei,g);
                boost::tie(e1, in1) = add_edge(verts[tail], verts[head], gshort);
                boost::tie(e2, in2) = add_edge(verts[head], verts[tail], gshort);
                capacitysh[e1] = capacity[*ei];
                capacitysh[e2] = 0;
                revsh[e1] = e2;
                revsh[e2] = e1;
                countcollapsebfs1++;
                edging[e1] = {tail,head};
            }
            else if(selected[*u_iter]==false && target(*ei,g)!=threesimplex_to_node[deathi] && selected[target(*ei,g)]==true && capacity[*ei]!=0){// transfer this node between green and initial
                edge_descriptor e1, e2;
                bool in1, in2;
                int tail, head;
                tail = *u_iter;
                head = threesimplex_to_node[deathi];
                boost::tie(e1, in1) = add_edge(verts[tail], verts[head], gshort);
                boost::tie(e2, in2) = add_edge(verts[head], verts[tail], gshort);
                capacitysh[e1] = capacity[*ei];
                capacitysh[e2] = 0;
                revsh[e1] = e2;
                revsh[e2] = e1;
                countcollapsebfs2++;
                edging[e1] = {tail,head};
            }
            else if(capacity[*ei]!=0){// keep original vertex
                edge_descriptor e1, e2;
                bool in1, in2;
                int tail, head;
                tail = *u_iter;
                head = target(*ei,g);
                boost::tie(e1, in1) = add_edge(verts[tail], verts[head], gshort);
                boost::tie(e2, in2) = add_edge(verts[head], verts[tail], gshort);
                capacitysh[e1] = capacity[*ei];
                capacitysh[e2] = 0;
                revsh[e1] = e2;
                revsh[e2] = e1;
                countcollapsebfs3++;
                edging[e1] = {tail,head};
            }
        }
    }
    cout<<"Got out of removing edges "<<countcollapsebfs<<" 1: "<<countcollapsebfs1<<" 2: "<<countcollapsebfs2<<" 3: "<<countcollapsebfs3 ;
    if(verbose)
        getchar();
//    for (boost::tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
//        if(selected[*u_iter])
//            remove_vertex(*u_iter,g);
//
//    cout<<"Got out of removing vertices "<<countcollapsebfs; getchar();
//    std::vector<default_color_type> color(num_vertices(g));
//    std::vector<long> distance(num_vertices(g));
    
    
    
    long flow = boykov_kolmogorov_max_flow(gshort ,s, t);
    
    cout<<"Got out of cut ";if(verbose) getchar();
    
    std::cout << "\nThe total flow:" << std::endl;
    std::cout << "s " << flow << std::endl << std::endl;
    
    for (boost::tie(u_iter, u_end) = vertices(gshort); u_iter != u_end; ++u_iter)
        for (boost::tie(ei, e_end) = out_edges(*u_iter, gshort); ei != e_end; ++ei)
            if (colorsh[*u_iter] != colorsh[target(*ei,gshort)] && (capacitysh[*ei] - residual_capacitysh[*ei])>0){
//                std::cout << "f " << (*u_iter) << " " << (target(*ei, g)) << " "
//                << (capacity[*ei] - residual_capacity[*ei]) << " COLOR: ("<<color[*u_iter]<<","<<color[target(*ei,g)]<<")"
//                x << "\t";//filtgrgrrev[birthi][ *u_iter]<<" "<< filtgrgrrev[birthi][target(*ei, g) ]<<" "<< std::endl;
                if(edging.find(*ei)!=edging.end()){
                    vectri.push_back({edging[*ei][0],edging[*ei][1]});
                        continue;
                }
                vectri.push_back({int(*u_iter),int(target(*ei, gshort)) });
            }
    cout<<"cut size: "<<vectri.size(); //getchar();
    verts.clear();
    gshort.clear();
    return vectri;
    
}


vector<unsigned long> returnSquare(Bitmap_cubical_complex &b, int k1 , int k2 ){
    vector<vector<vector<unsigned long>>> vert1, vert2;
    vector<vector<unsigned long>> vert;

    ////////////// DO IT FOR k1 /////////////////////////
    auto c3r = b.boundary_range(k1);
    if(c3r.size()!=6){cout<<"mistake 1 rs: "<<c3r.size(); exit(0);}
    
    for(size_t c2c=0; c2c<c3r.size(); c2c++){//going inside square
        auto c2r = b.boundary_range(c3r[c2c]);
        if(c2r.size()!=4){cout<<"mistake 2: "<<c2r.size(); exit(0);}
        vector<vector<unsigned long>> buff1;
        for(size_t c1c=0; c1c<c2r.size(); c1c++){// going inside edge
            auto c1r = b.boundary_range(c2r[c1c]);
            if(c1r.size()!=2){cout<<"mistake 3"; exit(0);}
            vector<unsigned long> buff2{vertind[c1r[0]],vertind[c1r[1]]};
            buff1.push_back(buff2);
//            vert1.push_back(c1r[1]);
//            cout<<"("<<vertind[(c1r[0])]<<","<<vertind[(c1r[1])]<<")";
        }
        vert1.push_back(buff1);
//        cout<<"\n";
    }
//    cout<<"v2:";
    ////////////// DO IT FOR k2 /////////////////////////
    c3r = (b.boundary_range(k2));
    if(c3r.size()!=6){cout<<"mistake 1 k2: "<<c3r.size(); exit(0);}
    
    for(size_t c2c=0; c2c<c3r.size(); c2c++){//going inside square
        auto c2r = b.boundary_range(c3r[c2c]);
        if(c2r.size()!=4){cout<<"mistake 2: "<<c2r.size(); exit(0);}
        
        vector<vector<unsigned long>> buff1;
        for(size_t c1c=0; c1c<c2r.size(); c1c++){// going inside edge
            auto c1r = b.boundary_range(c2r[c1c]);
            if(c1r.size()!=2){cout<<"mistake 3"; exit(0);}
            vector<unsigned long> buff2{vertind[c1r[0]],vertind[c1r[1]]};
            buff1.push_back(buff2);
//             cout<<"("<<vertind[(c1r[0])]<<","<<vertind[c1r[1]]<<")";
//           vert2.push_back(c1r[0]);
//            vert2.push_back(c1r[1]);
        }
        vert2.push_back(buff1);
//        cout<<"\n";
    }
//    cout<<"came till part 2"; getchar();
    bool found = false;
    int ob1,ob2;
    for(ob1=0;ob1<vert1.size();ob1++){
        for(ob2=0;ob2<vert2.size();ob2++){
            if(vert1[ob1]==vert2[ob2]){
                found = true;
                break;
            }

        }
        if(found==true)
            break;
    }
    
    if(found==false){
//        cout<<"wrong algo";
//         for(int ob1=0;ob1<vert1.size();ob1++)
//            cout<<vert1[ob1][0][0]<<" "<<vert1[ob1][1][0]<<" "<<vert1[ob1][2][0]<<" "<<vert1[ob1][3][0]<<"|";
//        cout<<"\n";
//        for(int ob2=0;ob2<vert2.size();ob2++)
//            cout<<vert2[ob2][0][0]<<" "<<vert2[ob2][1][0]<<" "<<vert2[ob2][2][0]<<" "<<vert2[ob2][3][0]<<"|";
//        getchar();
        return{};
        //exit(0);
    }
    

    vector<unsigned long> temp{0,0};

    temp[0] = vert1[ob1][1][0];
    temp[1] = vert1[ob1][1][1];
    vert1[ob1][1][0] = vert1[ob1][2][0];
    vert1[ob1][1][1] = vert1[ob1][2][1];
    vert1[ob1][2][0] = temp[0];
    vert1[ob1][2][1] = temp[1];

    if((vert1[ob1][0][0]==vert1[ob1][1][0])||(vert1[ob1][0][0]==vert1[ob1][1][1])){
        auto tempx=vert1[ob1][0][0];
        vert1[ob1][0][0]=vert1[ob1][0][1];
        vert1[ob1][0][1]=tempx;
    }

    for(int ii=1;ii<vert1[ob1].size();ii++){
        if(vert1[ob1][ii-1][1]==vert1[ob1][ii][0]){
            continue;
        }
        else if(vert1[ob1][ii-1][1]==vert1[ob1][ii][1]){
            auto tempx=vert1[ob1][ii][1];
        vert1[ob1][ii][1] = vert1[ob1][ii][0];
        vert1[ob1][ii][0] = tempx;

        }
        else{
            cout<<"wrong this part"; exit(0);
        }
    }

    return{vert1[ob1][0][0],vert1[ob1][1][0],vert1[ob1][2][0],vert1[ob1][3][0]};

}


vector<vector<unsigned long>> returnSquaresing(Bitmap_cubical_complex &b, int k , float death){
    auto c3r = (b.boundary_range(k));// triangles
    if(c3r.size()!=6){
        cout<<"mistake 0"; exit(0);
    }
//    cout<<"c3r: "<<c3r.size(); getchar();
    vector<vector<unsigned long>> vert1;
    
    //c2c: iterate over
    for(size_t c2c=0; c2c<c3r.size(); c2c++){//for each triangle
        auto c2r = b.coboundary_range(c3r[c2c]); // tetrahedrons: 2 or 1
        if(c2r.size()>2){
            cout<<"mistake 0.5"; exit(0);
        }
        
        vector<unsigned long> buff1;
        if(c2r.size()==1){// vertices of c3r has to be pushed in
            auto c2r = b.boundary_range(c3r[c2c]);
            if(c2r.size()!=4){
                cout<<"mistake 1"; exit(0);
            }
            for(size_t c1c=0; c1c<c2r.size(); c1c++){// going inside edge, push edges
                auto c1r = b.boundary_range(c2r[c1c]);
                if(c1r.size()!=2){
                    cout<<"mistake 2"; exit(0);
                }
                buff1.push_back(vertind[c1r[0]]);
                buff1.push_back(vertind[c1r[1]]);
            }
            if(buff1.size()!=8){
                cout<<"mistake 3"; exit(0);
            }
            
            sort(buff1.begin(),buff1.end());
            buff1.erase(std::unique(buff1.begin(), buff1.end()), buff1.end());
            auto temp = buff1[0]; buff1[0] = buff1[1]; buff1[1] = temp;
            if(buff1.size()!=4){
                cout<<"mistake 4"; exit(0);
            }

        }
        else if(b.filtration(c2r[0])>death && b.filtration(c2r[1])>death)
        {cout<<"\n\nExiting\n\n"; exit(0);}
        else if(b.filtration(c2r[0])>death || b.filtration(c2r[1])>death ){
            buff1 = returnSquare(b, c2r[0] , c2r[1] );
            
        }
        else
        continue;
        if(buff1.size()!=0)
        vert1.push_back(buff1);
    }
    return vert1;
    
}


void writeToFileOFF(Bitmap_cubical_complex &cub, vector<vector<int>> vectri,  size_t birthi, size_t deathi, string tr_file){
    ofstream ft(tr_file.c_str());
    
    int count = 0, ind = -1;
    for(int i=0;i<vectri.size();i++){
        
        auto findinf1 = threesimplex_to_node_rev.find(vectri[i][0]);
        auto findinf2 = threesimplex_to_node_rev.find(vectri[i][1]);
        if(findinf1 == threesimplex_to_node_rev.end() || findinf2 == threesimplex_to_node_rev.end()){
            if(findinf1 == threesimplex_to_node_rev.end()){
                int getindex1 = threesimplex_to_node_rev[vectri[i][1]];
                vector<vector<unsigned long>> tri =  returnSquaresing(cub, getindex1, cub.filtration(deathi));
                count+=tri.size();
            }
            if(findinf2 == threesimplex_to_node_rev.end()){
                int getindex2 = threesimplex_to_node_rev[vectri[i][0]];
                vector<vector<unsigned long>> tri = returnSquaresing(cub,getindex2, cub.filtration(deathi));
                count+=tri.size();
            }
        }
        
        else{
            int getindex1 = threesimplex_to_node_rev[vectri[i][0]];
            int getindex2 = threesimplex_to_node_rev[vectri[i][1]]; //indices in the filtration
            
            if(cub.dimension(getindex2)!=3 || cub.dimension(getindex1)!=3 ){
                cout<<"Should be 3-simplices.Error! "<<cub.dimension(getindex1)<<"->"<<vectri[i][0]<<" "<<cub.dimension(getindex2)<<"->"<<vectri[i][1]<<"\n"; exit(0);
            }
            vector<unsigned long> tri = returnSquare(cub, getindex1 , getindex2 );
            
            if(tri.size()==0){
//                 vector<vector<unsigned long>> tri2 = returnSquaresing(cub, getindex1, cub.filtration(deathi));
//                count+=tri2.size();
//                tri2 = returnSquaresing(cub, getindex2, cub.filtration(deathi));
//                count+=tri2.size();
            continue;
            }
            count++;
            
        }
    }
    
    
    ft<<"OFF\n"<<(count*4)<<" "<<count<<" 0\n";
    for(int i=0;i<vectri.size();i++){
        
        auto findinf1 = threesimplex_to_node_rev.find(vectri[i][0]);
        auto findinf2 = threesimplex_to_node_rev.find(vectri[i][1]);
        if(findinf1 == threesimplex_to_node_rev.end() && findinf2 == threesimplex_to_node_rev.end()){
            cout<<"both are zero\n"; exit(0);
        }
        if(findinf1 == threesimplex_to_node_rev.end() || findinf2 == threesimplex_to_node_rev.end()){
            //            cout<<"one zero\n"; //continue;
            if(findinf1 == threesimplex_to_node_rev.end()){
                int getindex1 = threesimplex_to_node_rev[vectri[i][1]];
                vector<vector<unsigned long>> tri =  returnSquaresing(cub, getindex1, cub.filtration(deathi));
                //                cout<<"tri size: "<<tri.size()<<" "<<tri[0].size(); getchar();
                for(int ix=0;ix<tri.size();ix++){
                    for(int ii=0;ii<4;ii++)
                        ft<<x[int(tri[ix][ii])]<<" "<<y[int(tri[ix][ii])]<<" "<<z[int(tri[ix][ii])]<<"\n";
//                    ft<<"\n";
                }
            }
            
            if(findinf2 == threesimplex_to_node_rev.end()){
                int getindex2 = threesimplex_to_node_rev[vectri[i][0]];
                vector<vector<unsigned long>> tri = returnSquaresing(cub,getindex2, cub.filtration(deathi));
                //                cout<<"tri size: "<<tri.size()<<" "<<tri[0].size(); getchar();
                for(int ix=0;ix<tri.size();ix++){
                    for(int ii=0;ii<4;ii++)
                        ft<<x[int(tri[ix][ii])]<<" "<<y[int(tri[ix][ii])]<<" "<<z[int(tri[ix][ii])]<<"\n";
//                    ft<<"\n";
                }
            }
        }
        
        
        else{
            int getindex1 = threesimplex_to_node_rev[vectri[i][0]];
            int getindex2 = threesimplex_to_node_rev[vectri[i][1]]; //indices in the filtration
            
            if(cub.dimension(getindex2)!=3 || cub.dimension(getindex1)!=3 ){
                cout<<"Should be 3-simplices.Error! "<<cub.dimension(getindex1)<<"->"<<vectri[i][0]<<" "<<cub.dimension(getindex2)<<"->"<<vectri[i][1]<<"\n"; exit(0);
            }
            vector<unsigned long> tri = returnSquare(cub, getindex1 , getindex2 );
            if(tri.size()==0){
                continue;
            }
            for(int ix=0;ix<4;ix++){
                ft<<x[int(tri[ix])]<<" "<<y[int(tri[ix])]<<" "<<z[int(tri[ix])]<<"\n";
            }
//            ft<<"\n";
            
        }
    }
    
    for(int i=0;i<count;i++){
        ft<<"4 "<<(ind+1)<<" "<<(ind+2)<<" "<<(ind+3)<<" "<<(ind+4)<<endl;
        ind+=4;
    }
    
    
}
void writeToFile(Bitmap_cubical_complex &cub, vector<vector<int>> vectri, size_t birthi, size_t deathi, string tr_file){//, ){
    
    
    ofstream ft(tr_file.c_str());
    
    ft<<vectri.size()<<"\n";
    
    for(int i=0;i<vectri.size();i++){

        auto findinf1 = threesimplex_to_node_rev.find(vectri[i][0]);
        auto findinf2 = threesimplex_to_node_rev.find(vectri[i][1]);
        if(findinf1 == threesimplex_to_node_rev.end() && findinf2 == threesimplex_to_node_rev.end()){
              cout<<"both are zero\n"; exit(0);
        }
        if(findinf1 == threesimplex_to_node_rev.end() || findinf2 == threesimplex_to_node_rev.end()){
//            cout<<"one zero\n"; //continue;
            if(findinf1 == threesimplex_to_node_rev.end()){
                int getindex1 = threesimplex_to_node_rev[vectri[i][1]];
                vector<vector<unsigned long>> tri =  returnSquaresing(cub, getindex1, cub.filtration(deathi));
//                cout<<"tri size: "<<tri.size()<<" "<<tri[0].size(); getchar();
                for(int ix=0;ix<tri.size();ix++){
                    for(int ii=0;ii<4;ii++)
                    ft<<x[int(tri[ix][ii])]<<" "<<y[int(tri[ix][ii])]<<" "<<z[int(tri[ix][ii])]<<" ";
                    ft<<"\n";
                }
            }
            
            if(findinf2 == threesimplex_to_node_rev.end()){
                int getindex2 = threesimplex_to_node_rev[vectri[i][0]];
                vector<vector<unsigned long>> tri = returnSquaresing(cub,getindex2, cub.filtration(deathi));
//                cout<<"tri size: "<<tri.size()<<" "<<tri[0].size(); getchar();
                for(int ix=0;ix<tri.size();ix++){
                    for(int ii=0;ii<4;ii++)
                        ft<<x[int(tri[ix][ii])]<<" "<<y[int(tri[ix][ii])]<<" "<<z[int(tri[ix][ii])]<<" ";
                    ft<<"\n";
                }
            }
        }
        
        
        else{
        int getindex1 = threesimplex_to_node_rev[vectri[i][0]];
        int getindex2 = threesimplex_to_node_rev[vectri[i][1]]; //indices in the filtration

        if(cub.dimension(getindex2)!=3 || cub.dimension(getindex1)!=3 ){
            cout<<"Should be 3-simplices.Error! "<<cub.dimension(getindex1)<<"->"<<vectri[i][0]<<" "<<cub.dimension(getindex2)<<"->"<<vectri[i][1]<<"\n"; exit(0);
        }
            vector<unsigned long> tri = returnSquare(cub, getindex1 , getindex2 );
            if(tri.size()==0)
                continue;
            for(int ix=0;ix<4;ix++){
                ft<<x[int(tri[ix])]<<" "<<y[int(tri[ix])]<<" "<<z[int(tri[ix])]<<" ";
            }
            ft<<"\n";
        
        }
    }
    
}

void buildFullCubicSetting(Bitmap_cubical_complex cub, string point_file, vector<size_t> deathiv){

    int threesimplex = 0, twosimplex = 0, countvert = 0;
    int xrange, yrange, zrange;
    char sLine[256] = "";
    char inf[] = "inf";
    ifstream openpt(point_file.c_str());
    if(openpt.good()==false)
    {
        cout<<"Point file "<<point_file<<" does not exist.";
        exit(0);
    }
    openpt.getline(sLine,256);
    if(sLine[0]!='3'){
        cout<<"Dimension must be 3.Error!"; exit(0);
    }
    
    openpt.getline(sLine,256);
    stringstream ss1(sLine);
    
    ss1 >> xrange;
    openpt.getline(sLine,256);
    stringstream ss2(sLine);
    ss2 >> yrange;
    openpt.getline(sLine,256);
    stringstream ss3(sLine);
    ss3 >> zrange;
    ++xrange;   ++yrange;   ++zrange;
    vector<int> simpcount = {xrange*yrange*zrange,0,0,0};
    cout<<"xrange: "<<xrange<<" yrange: "<<yrange<<" zrange: "<<zrange;
   // getchar();
    //Changed here !!! Rev it
    for(int k=0;k<zrange;k++)
        for(int j=0;j<yrange;j++)
            for(int i=0;i<xrange;i++){
                x.push_back(i);
                y.push_back(j);
                z.push_back(k);
                
            }
    
    for(auto iter: deathiv)
        sizescale[iter] = {0,0};
    
    for(size_t k=0;k<cub.num_simplices();k++){

        if(cub.dimension(k)==1){
            simpcount[1]++;

        }


        else if(cub.dimension(k)==2){
            simpcount[2]++;
            twosimplex_to_tri[k] = twosimplex;
            twosimplex_to_tri_rev[twosimplex] = k;
            twosimplex++;
        }

        else if(cub.dimension(k)==3){
            simpcount[3]++;
            threesimplex_to_node[k] = threesimplex;
            threesimplex_to_node_rev[threesimplex] = k;
            threesimplex++;
        }

        else if(cub.dimension(k)==0){
            simpcount[0]++;
            vertind[k]=countvert;
            countvert++;

        }
        else{
            cout<<"Unknown simplex error"; exit(0);
        }
        for(auto iter: deathiv){
            if(cub.dimension(k)==2 && cub.filtration(k)<=cub.filtration(iter))//Triangle:Push purple or green
                sizescale[iter][0]++;
            else if(cub.dimension(k)==3 && cub.filtration(k)<=cub.filtration(iter) && k>sizescale[iter][1])
                sizescale[iter][1]=k;
        }
    }
    if(countvert!=x.size()){
        cout<<"size mismatch "<<countvert<<" "<<simpcount.size()<<" "<<x.size()<<endl;
        exit(0);
    }

    
    

}

