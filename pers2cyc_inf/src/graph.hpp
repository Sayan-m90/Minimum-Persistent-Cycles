#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_

#include <algorithm>
#include "utils.h"

template <class NodeIdType>
class Graph {
public:
    void addEdge(const NodeIdType node_id1, const NodeIdType node_id2);
    bool adjNodesHasDup();
    int size() { return adj_nodes_.size(); }
    void reserve(const int& size) { adj_nodes_.reserve(size); }
    void print();

private:
    void addAdjNode(const NodeIdType node_id1, const NodeIdType node_id2);

public:
    vector<vector<NodeIdType>> adj_nodes_;
};

template <class NodeIdType> void Graph<NodeIdType>::
    addEdge(const NodeIdType node_id1, const NodeIdType node_id2) {

    auto max_node_id = std::max(node_id1, node_id2);
    if (max_node_id >= this->adj_nodes_.size()) {
        this->adj_nodes_.resize(max_node_id+1);
    }

    this->addAdjNode(node_id1, node_id2);
    this->addAdjNode(node_id2, node_id1);
}

template <class NodeIdType> inline void Graph<NodeIdType>::
    addAdjNode(const NodeIdType node_id1, const NodeIdType node_id2) {

    // if (std::find(adj_nodes_[node_id1].begin(), adj_nodes_[node_id1].end(), node_id2) 
    //     == adj_nodes_[node_id1].end()) {
    //     adj_nodes_[node_id1].push_back(node_id2);
    // }

    this->adj_nodes_[node_id1].push_back(node_id2);
}

template <class NodeIdType> bool Graph<NodeIdType>::adjNodesHasDup() {
    for (auto i = 0; i < this->adj_nodes_.size(); i ++) {
        vector<NodeIdType> adj_nodes = this->adj_nodes_[i];
        if (adj_nodes.size() < 2) {
            continue;
        }

        std::sort(adj_nodes.begin(), adj_nodes.end());

        for (auto j = 0; j < adj_nodes.size() - 1; j ++) {
            if (adj_nodes[j] == adj_nodes[j+1]) {
                return true;
            }
        }
    }

    return false;
}

template <class NodeIdType> void Graph<NodeIdType>::print() {
    for (auto i = 0; i < this->adj_nodes_.size(); i ++) {
        std::sort(this->adj_nodes_[i].begin(), this->adj_nodes_[i].end());
        cout << i << ": " << this->adj_nodes_[i] << endl;
    }
}

#endif
