#include <map>
#include <cmath>
#include "graph.h"
#include "MutablePriorityQueue.h"

int Graph::getNumVertex() const {
    return vertexSet.size();
}

vector<Vertex *> Graph::getVertexSet() const {
    return vertexSet;
}

/*
 * Auxiliary function to find a vertex with a given content.
 */
Vertex * Graph::findVertex(const int &id) const {
    for (auto v : vertexSet)
        if (v->getId() == id)
            return v;
    return nullptr;
}

/*
 * Finds the index of the vertex with a given content.
 */
int Graph::findVertexIdx(const int &id) const {
    for (unsigned i = 0; i < vertexSet.size(); i++)
        if (vertexSet[i]->getId() == id)
            return i;
    return -1;
}

/*
 *  Adds a vertex with a given content or info (in) to a graph (this).
 *  Returns true if successful, and false if a vertex with that content already exists.
 */
bool Graph::addVertex(const int &id) {
    if (findVertex(id) != nullptr)
        return false;
    vertexSet.push_back(new Vertex(id));
    return true;
}

/*
 * Adds an edge to a graph (this), given the contents of the source and
 * destination vertices and the edge weight (w).
 * Returns true if successful, and false if the source or destination vertex does not exist.
 */
bool Graph::addEdge(const int &sourc, const int &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    v1->addEdge(v2, w);
    v2->addEdge(v1, w);
    dist[sourc][dest] = dist[dest][sourc] = w;
    return true;
}

bool Graph::addBidirectionalEdge(const int &sourc, const int &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    auto e1 = v1->addEdge(v2, w);
    auto e2 = v2->addEdge(v1, w);
    e1->setReverse(e2);
    e2->setReverse(e1);

    dist[sourc][dest] = dist[dest][sourc] = w;
    return true;
}



Graph::~Graph() {}



bool Graph::removeEdge(const int &source, const int &dest) {
    Vertex * srcVertex = findVertex(source);
    if (srcVertex == nullptr) {
        return false;
    }
    return srcVertex->removeEdge(dest);
}

void Graph::dfsVisit(Vertex *v, std::vector<int> & res) const {
    v->setVisited(true);
    res.push_back(v->getId());
    for (auto & e : v->getAdj()) {
        auto w = e->getDest();
        if (!w->isVisited()) {
            dfsVisit(w, res);
        }
    }
}

vector<int> Graph::dfs(const int & source) const {
    std::vector<int> res;
    // Get the source vertex
    auto s = findVertex(source);
    if (s == nullptr) {
        return res;
    }

    // Set that no vertex has been visited yet
    for (auto v : vertexSet) {
        v->setVisited(false);
    }

    // Perform the actual DFS using recursion
    dfsVisit(s, res);

    return res;
}

double Graph::tsp_bt() {
    auto** dists = new double*[this->getNumVertex()];

    for(int i = 0; i < this->getNumVertex(); i++) {
        dists[i] = new double[1 << this->getNumVertex()];

        for(int j = 0; j < (1 << this->getNumVertex()); j++) {
            dists[i][j] = -1;
        }
    }

    double res = calculate_tsp(0, 1, dists);

    return res;
}

double Graph::calculate_tsp(int node, int edge, double** dists) {
    if(dists[node][edge] != -1) return dists[node][edge];

    double ans = INF, tmp;

    if(edge == pow(2, this->getNumVertex()) - 1) {
        return dist[node][0] > 0 ? dist[node][0] : INF;
    }

    for(int i = 0; i < this->getNumVertex(); i++) {
        if(dist[node][i] == 0) continue;

        if(!(edge & (int) pow(2, i))) {
            tmp = dist[node][i] + calculate_tsp(i, edge | (int) pow(2, i), dists);
            ans = min(ans, tmp);
        }
    }

    return dists[node][edge] = ans;
}

vector<Vertex *> Graph::prim() {
    if (vertexSet.empty()) {return this->vertexSet;}

    // Reset auxiliary info
    for(auto v : vertexSet) {
        v->setDist(INF);
        v->setPath(nullptr);
        v->setVisited(false);
    }

    // start with an arbitrary vertex
    Vertex* s = vertexSet.front();
    s->setDist(0);

    // initialize priority queue
    MutablePriorityQueue<Vertex> q;
    q.insert(s);
    // process vertices in the priority queue
    while( ! q.empty() ) {
        auto v = q.extractMin();
        v->setVisited(true);
        for(auto &e : v->getAdj()) {
            Vertex* w = e->getDest();
            if (!w->isVisited()) {
                auto oldDist = w->getDist();
                if(e->getWeight() < oldDist) {
                    w->setDist(e->getWeight());
                    w->setPath(e);
                    //cout << v->getId() << " - " << e->getWeight() << " - " << w->getId() << endl;
                    if (oldDist == INF) {
                        q.insert(w);
                    }
                    else {
                        q.decreaseKey(w);
                    }
                }
            }
        }
    }

    return this->vertexSet;
}
