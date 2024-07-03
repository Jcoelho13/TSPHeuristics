#ifndef DA_TP_CLASSES_GRAPH
#define DA_TP_CLASSES_GRAPH

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <list>
#include "vertexEdge.h"

using namespace std;

typedef vector<double> vd;
typedef vector<vd> vvd;

class Graph {
public:
    ~Graph();
    /*
    * Auxiliary function to find a vertex with a given ID.
    */
    Vertex *findVertex(const int &id) const;
    /*
     *  Adds a vertex with a given content or info (in) to a graph (this).
     *  Returns true if successful, and false if a vertex with that content already exists.
     */
    bool addVertex(const int &id);

    /*
     * Adds an edge to a graph (this), given the contents of the source and
     * destination vertices and the edge weight (w).
     * Returns true if successful, and false if the source or destination vertex does not exist.
     */
    bool addEdge(const int &sourc, const int &dest, double w);
    bool addBidirectionalEdge(const int &sourc, const int &dest, double w);

    int getNumVertex() const;
    bool removeEdge(const int &source, const int &dest);
    vector<Vertex *> getVertexSet() const;

    /**
    * @brief Implementation of the Prim algorithm
    * @return vector with the vertices of the MST
    * Complexity: O((V+E) log V), where E is the number of edges and V is the number of vertices
     */
    vector<Vertex *> prim();

    void dfsVisit(Vertex *v, std::vector<int> & res) const;
    vector<int> dfs(const int & source) const;

    /**
    * @brief Using backtracking, calculates the shortest path that visits all the nodes and returns to the starting node.
    * @return cost of the path
    * Complexity: O(2^n * n^2), where n is the number of vertices
    */
    double tsp_bt();

    /**
    * @brief Calculates the shortest path that visits all the nodes and returns to the starting node.
    * @param node current node
    * @param edge current edge
    * @param dists auxiliary matrix
    * @return cost of the path
    * Complexity: O(2^n * n^2), where n is the number of vertices
    */
    double calculate_tsp(int node, int edge, double** dists);

    vvd dist;

protected:
    vector<Vertex *> vertexSet;    // vertex set

    /*
     * Finds the index of the vertex with a given content.
     */
    int findVertexIdx(const int &id) const;
};

void deleteMatrix(int **m, int n);
void deleteMatrix(double **m, int n);

#endif /* DA_TP_CLASSES_GRAPH */