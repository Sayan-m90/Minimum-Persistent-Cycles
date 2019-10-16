#ifndef _FLOW_GRAPH_HPP_
#define _FLOW_GRAPH_HPP_

#include <algorithm>
#include "utils.h"

template<class NodeIdType, class WeightType>
struct FlowGraphEdge {
    NodeIdType target_v;
    WeightType capacity;
    WeightType flow;
};

template <class NodeIdType, class WeightType>
class FlowGraph {
public:
    FlowGraph(long vertex_size) {
        adj_list_.resize(vertex_size);
    }

    void addEdge(const NodeIdType from_v, const NodeIdType to_v, const WeightType w) {
        this->adj_list_[from_v].push_back({to_v, w, 0});
    }

    // bool adjNodesHasDup();

    int size() { return adj_list_.size(); }

    void reserve(const int& size) { adj_list_.reserve(size); }

    // void print();

private:
    // void addAdjNode(const NodeIdType node_id1, const NodeIdType node_id2);

public:
    vector<vector<FlowGraphEdge<NodeIdType, WeightType>>> adj_list_;
};

// template <class NodeIdType, class WeightType> inline 
// void FlowGraph<NodeIdType, WeightType>::
//     addAdjNode(const NodeIdType node_id1, const NodeIdType node_id2) {

//     // if (std::find(adj_list_[node_id1].begin(), adj_list_[node_id1].end(), node_id2) 
//     //     == adj_list_[node_id1].end()) {
//     //     adj_list_[node_id1].push_back(node_id2);
//     // }

//     this->adj_list_[node_id1].push_back(node_id2);
// }

// template <class NodeIdType, class WeightType> 
// bool FlowGraph<NodeIdType, WeightType>::adjNodesHasDup() {
//     for (auto i = 0; i < this->adj_list_.size(); i ++) {
//         vector<NodeIdType> adj_nodes = this->adj_list_[i];
//         if (adj_nodes.size() < 2) {
//             continue;
//         }

//         std::sort(adj_nodes.begin(), adj_nodes.end());

//         for (auto j = 0; j < adj_nodes.size() - 1; j ++) {
//             if (adj_nodes[j] == adj_nodes[j+1]) {
//                 return true;
//             }
//         }
//     }

//     return false;
// }

// template <class NodeIdType, class WeightType> 
// void FlowGraph<NodeIdType>::print() {
//     for (auto i = 0; i < this->adj_list_.size(); i ++) {
//         std::sort(this->adj_list_[i].begin(), this->adj_list_[i].end());
//         cout << i << ": " << this->adj_list_[i] << endl;
//     }
// }

#endif
