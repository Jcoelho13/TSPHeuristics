/*! \file */

/**
 * @brief Main file
 */

#include <iostream>
#include <limits>
#include <sstream>
#include <fstream>
#include "DataStructures/graph.h"
#include "DataStructures/MutablePriorityQueue.h"
#include "data/data.h"
#include "DataStructures/Node.h"
#include <cmath>
#include <unordered_map>
#include <stack>

using namespace std;

Graph tg1, tg2, tg3, rg1, rg2, rg3, fcg25, fcg50, fcg75, fcg100, fcg200, fcg300, fcg400, fcg500, fcg600, fcg700, fcg800, fcg900;

/**
 * @brief Reads a toy graph from a file, and saves it to a graph.
 * Complexity: O(n), where n is the number of lines in the file.
 */
void read_toy_graph(const string& file, Graph &g);

/**
 * @brief Reads the nodes of a real world graph from a file and saves it to a graph.
 * Complexity: O(n), where n is the number of lines in the file.
 */
void read_real_world_graph_nodes(const string& file, Graph &g);

/**
 * @brief Reads the edges of a real world graph from a file.
 * Complexity: O(n), where n is the number of lines in the file.
 */
void read_real_world_graph_edges(const string& file, Graph &g);

/**
 * @brief Reads a fully connected graph from a file, and saves it to a graph.
 * Complexity: O(n), where n is the number of lines in the file.
 */
void read_fully_connected_graph(const string& file, Graph &g);

/** @brief Calculates the Euler Path for a given graph
 * @param g Graph for which the Euler Path is calculated
 * @return The calculated Euler Path
 * Complexity: Hierholzer's algorithm is O(|E|), E being the number of edges
 */
vector<int> euler_path(Graph &g);
/**
 * @brief Calculates the minimum-weight perfect matching for the given vertices of the given graph
 * @param odd Vertices to match
 * @param g Graph to which the vertices belong
 * @return The minimum-weight perfect matches
 * Complexity: O(V^2|E|) where V i sth enumber of vertices and E is the number of edges
 */
vector<pair<int, int>> matching(vector<int> &odd, Graph &g);

vector<Vertex *> res;

/**
 * @brief Calculates the cost of a spanning tree.
 * @param ans Vector of vertices that represent the spanning tree.
 * @return The cost of the spanning tree.
 * Complexity: O(|V|*|E|), where V is the number of vertices in the Tree and E the number of edges.
 */
double spanningTreeCost(const std::vector<Vertex *> &ans);
vector<Node> nodes;


/**
 * @brief prints a menu that lets the user choose the toygraph to use, and after selecting it, calls the function that reads it and prints a menu with the options for that graph.
 * Complexity: O(1)
 */
void print_menu_1();

/**
 * @brief prints a menu that lets the user choose the real world graph to use, and after selecting it, calls the functions that read it and prints a menu with the options for that graph.
 * Complexity: O(1)
 */
void print_menu_2();

/**
 * @brief prints a menu that lets the user choose the fully connected graph to use, and after selecting it, calls the function that reads it and prints a menu with the options for that graph.
 * Complexity: O(1)
 */
void print_menu_3();

/**
 * @brief prints a menu with the options for the toy graph.
 * Complexity: O(1).
 * @param g graph to use.
 */
void print_menu_options(Graph &g);

/**
 * @brief prints a menu with the options for the real world graph / fully connected graph - No backtracking.
 * Complexity: O(1).
 * @param g graph to use.
 */
void print_menu_options_SBT(Graph &g);

/**
 * @brief clears the screen.
 * Complexity: O(1).
 */
void clear() {for (int i = 0; i < 50; i++) cout << endl;}

/**
 * @brief waits for the user to press enter.
 * Complexity: O(1).
 */
void wait() {cout << endl << "Press Enter to continue" << endl; cin.ignore(numeric_limits<streamsize>::max(), '\n'); cin.get();}

int main()  {
    int choice;
    do{
        clear();
        cout << "------------------------------------------" << endl;
        cout << "|            Choose a dataset            |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 1. Toy Graphs                          |" << endl;
        cout << "| 2. Real World Graphs                   |" << endl;
        cout << "| 3. Fully Connected Graphs              |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 0. Exit                                |" << endl;
        cout << "|                                        |" << endl;
        cout << "------------------------------------------" << endl;

        cout << "Choose an option:" << endl;
        cout << endl;
        cin >> choice;


        switch (choice) {
            case 0:
                choice = 0;
                break;
            case 1:
                print_menu_1();
                break;
            case 2:
                print_menu_2();
                break;
            case 3:
                print_menu_3();
                break;
            default:
                cout << "Invalid option! Try again" << endl;
                wait();
                choice = 1;
                break;
        }
    }while (choice != 0);

    return 0;
}

//  Print Menu Functions

void print_menu_1() {
    int choice;
    do{
        clear();
        cout << "------------------------------------------" << endl;
        cout << "|               Toy Graphs               |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 1. Shipping                            |" << endl;
        cout << "| 2. Stadiums                            |" << endl;
        cout << "| 3. Tourism                             |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 0. Go Back                             |" << endl;
        cout << "------------------------------------------" << endl;

        cout << "Choose an option:" << endl;
        cout << endl;
        cin >> choice;


        switch (choice) {
            case 1:
                read_toy_graph(shipping, tg1);
                print_menu_options(tg1);
                break;

            case 2:
                read_toy_graph(stadiums, tg2);
                print_menu_options(tg2);
                break;

            case 3:
                read_toy_graph(tourism, tg3);
                print_menu_options(tg3);
                break;

            case 0:
                choice = 0;
                break;

            default:
                cout << "Invalid option! Try again" << endl;
                wait();
                choice = 1;
                break;
        }
    }while (choice != 0);
}

void print_menu_2(){
    int choice;
    do{
        clear();
        cout << "------------------------------------------" << endl;
        cout << "|            Real World Graphs           |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 1. Graph 1                             |" << endl;
        cout << "| 2. Graph 2                             |" << endl;
        cout << "| 3. Graph 3                             |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 0. Go Back                             |" << endl;
        cout << "------------------------------------------" << endl;

        cout << "Choose an option:" << endl;
        cout << endl;
        cin >> choice;

        switch (choice) {
            case 1:
                read_real_world_graph_nodes(g1_nodes, rg1);
                read_real_world_graph_edges(g1_edges, rg1);
                print_menu_options_SBT(rg1);


                break;
            case 2:
                read_real_world_graph_nodes(g2_nodes, rg2);
                read_real_world_graph_edges(g2_edges, rg2);
                print_menu_options_SBT(rg2);

                break;
            case 3:
                read_real_world_graph_nodes(g3_nodes, rg3);
                read_real_world_graph_edges(g3_edges, rg3);
                print_menu_options_SBT(rg3);

                break;

            case 0:
                choice = 0;
                break;

            default:
                cout << "Invalid option! Try again" << endl;
                cout << "Press enter to continue..." << endl;
                wait();
                choice = 1;
                break;
        }

    }while (choice != 0);
}



void print_menu_3() {
    int choice;
    do{
        clear();
        cout << "------------------------------------------" << endl;
        cout << "|         Fully Connected Graphs         |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 1. Graph 1                             |" << endl;
        cout << "| 2. Graph 2                             |" << endl;
        cout << "| 3. Graph 3                             |" << endl;
        cout << "| 4. Graph 4                             |" << endl;
        cout << "| 5. Graph 5                             |" << endl;
        cout << "| 6. Graph 6                             |" << endl;
        cout << "| 7. Graph 7                             |" << endl;
        cout << "| 8. Graph 8                             |" << endl;
        cout << "| 9. Graph 9                             |" << endl;
        cout << "| 10. Graph 10                           |" << endl;
        cout << "| 11. Graph 11                           |" << endl;
        cout << "| 12. Graph 12                           |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 0. Go Back                             |" << endl;
        cout << "------------------------------------------" << endl;

        cout << "Choose an option:" << endl;
        cout << endl;
        cin >> choice;

        switch (choice) {
            case 1:
                read_fully_connected_graph(fcg25_data, fcg25);
                print_menu_options_SBT(fcg25);
                wait();
                break;
            case 2:
                read_fully_connected_graph(fcg50_data, fcg50);
                print_menu_options_SBT(fcg50);
                wait();
                break;
            case 3:
                read_fully_connected_graph(fcg75_data, fcg75);
                print_menu_options_SBT(fcg75);
                wait();
                break;
            case 4:
                read_fully_connected_graph(fcg100_data, fcg100);
                print_menu_options_SBT(fcg100);
                wait();
                break;

            case 5:
                read_fully_connected_graph(fcg200_data, fcg200);
                print_menu_options_SBT(fcg200);
                wait();
                break;

            case 6:
                read_fully_connected_graph(fcg300_data, fcg300);
                print_menu_options_SBT(fcg300);
                wait();
                break;

            case 7:
                read_fully_connected_graph(fcg400_data, fcg400);
                print_menu_options_SBT(fcg400);
                wait();
                break;

            case 8:
                read_fully_connected_graph(fcg500_data, fcg500);
                print_menu_options_SBT(fcg500);
                wait();
                break;

            case 9:
                read_fully_connected_graph(fcg600_data, fcg600);
                print_menu_options_SBT(fcg600);
                wait();
                break;

            case 10:
                read_fully_connected_graph(fcg700_data, fcg700);
                print_menu_options_SBT(fcg700);
                wait();
                break;

            case 11:
                read_fully_connected_graph(fcg800_data, fcg800);
                print_menu_options_SBT(fcg800);
                wait();
                break;

            case 12:
                read_fully_connected_graph(fcg900_data, fcg900);
                print_menu_options_SBT(fcg900);
                wait();
                break;

            case 0:
                choice = 0;
                break;

            default:
                cout << "Invalid option! Try again" << endl;
                cout << "Press enter to continue..." << endl;
                wait();
                choice = 1;
                break;
        }

    }while (choice != 0);
}

void print_menu_options(Graph &g) {
    int choice;
    do{
        clear();
        cout << "------------------------------------------" << endl;
        cout << "|            Choose an approach          |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 1. Backtracking Algorithm              |" << endl;
        cout << "| 2. Triangular Approximation Heuristic  |" << endl;
        cout << "| 3. Other Heuristics                    |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 0. Go Back                             |" << endl;
        cout << "|                                        |" << endl;
        cout << "------------------------------------------" << endl;

        cout << "Choose an option:" << endl;
        cout << endl;
        cin >> choice;


        switch (choice) {
            case 1:

                cout << g.tsp_bt() << endl;

                wait();

                break;
            case 2:

                res = g.prim();

                for(const auto v : res) {
                    cout << v->getId() << "<-";
                    if ( v->getPath() != nullptr ) {
                        cout << v->getPath()->getOrig()->getId();
                    }
                    cout << "|";
                }

                cout << endl;
                cout << spanningTreeCost(res) << endl;


                wait();


                break;
            case 3:{
                vector<Vertex *> mst = g.prim();

                for (Vertex *v : mst) v->setIndegree(0);
                for (Vertex *v : mst) {
                    if (v->getPath() == nullptr) continue;
                    v->setIndegree(v->getIndegree() + 1);
                }

                vector<int> odd;
                for (Vertex *v : mst) if (v->getIndegree() % 2 != 0) odd.push_back(v->getId());


                Graph final;
                final.dist = vvd(g.getNumVertex(), vd(g.getNumVertex(), 0));
                for (Vertex *v : mst) final.addVertex(v->getId());
                for (Vertex *v : mst) if (v->getPath() != nullptr) final.addEdge(v->getPath()->getOrig()->getId(), v->getId(), 1);

                auto matches = matching(odd, final);

                for (const auto &p : matches) {
                    final.addEdge(p.first, p.second, 1);
                }


                for (Vertex *v : final.getVertexSet()) {
                    for (Edge *e : v->getAdj()) {
                        e->setSelected(false);
                    }
                }

                vector<int> path = euler_path(final);


                unordered_map<int, bool> existing;
                vector<int> final_path;

                for (int id : path) {
                    if (existing.count(id) == 0) {
                        existing[id] = true;
                        final_path.push_back(id);
                    }
                }


                cout << "Path: ";
                for (int i : final_path) {
                    std::cout << i << (i == final_path.back() ? "\n" : " -> ");
                }

                double cost = 0;
                for (int i = 1; i < final_path.size(); i++) {
                    Vertex *v1 = g.findVertex(final_path[i - 1]);
                    for (Edge *e : v1->getAdj()) {
                        if (e->getDest()->getId() == final_path[i]) cost += e->getWeight();
                    }
                }
                cout << "Distance: " << cost << endl;

            }

                wait();

                break;

            case 0:
                choice = 0;
                break;

            default:
                cout << "Invalid option! Try again" << endl;
                wait();
                choice = 1;
                break;
        }
    }while (choice != 0);
}

void print_menu_options_SBT(Graph &g) {
    clock_t start, end;
    int choice;
    do{
        clear();
        cout << "------------------------------------------" << endl;
        cout << "|            Choose an approach          |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 1. Triangular Approximation Heuristic  |" << endl;
        cout << "| 2. Other Heuristics                    |" << endl;
        cout << "|                                        |" << endl;
        cout << "| 0. Go Back                             |" << endl;
        cout << "|                                        |" << endl;
        cout << "------------------------------------------" << endl;

        cout << "Choose an option:" << endl;
        cout << endl;
        cin >> choice;


        switch (choice) {
            case 1:
                start = clock();

                res = g.prim();
                cout << endl << spanningTreeCost(res) << endl;

                end = clock();
                cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;


                wait();

                break;

            case 0:
                choice = 0;
                break;

            case 2: {
                start = clock();
                vector<Vertex *> mst = g.prim();

                for (Vertex *v : mst) v->setIndegree(0);
                for (Vertex *v : mst) {
                    if (v->getPath() == nullptr) continue;
                    v->setIndegree(v->getIndegree() + 1);
                }

                vector<int> odd;
                for (Vertex *v : mst) if (v->getIndegree() % 2 != 0) odd.push_back(v->getId());


                Graph final;
                final.dist = vvd(g.getNumVertex(), vd(g.getNumVertex(), 0));
                for (Vertex *v : mst) final.addVertex(v->getId());
                for (Vertex *v : mst) if (v->getPath() != nullptr) final.addEdge(v->getPath()->getOrig()->getId(), v->getId(), 1);

                auto matches = matching(odd, final);

                for (const auto &p : matches) {
                    final.addEdge(p.first, p.second, 1);
                }


                for (Vertex *v : final.getVertexSet()) {
                    for (Edge *e : v->getAdj()) {
                        e->setSelected(false);
                    }
                }

                vector<int> path = euler_path(final);

                
                unordered_map<int, bool> existing;
                vector<int> final_path;

                for (int id : path) {
                    if (existing.count(id) == 0) {
                        existing[id] = true;
                        final_path.push_back(id);
                    }
                }
                

                cout << "Path: ";
                for (int i : final_path) {
                    std::cout << i << (i == final_path.back() ? "\n" : " -> ");
                }

                double cost = 0;
                for (int i = 1; i < final_path.size(); i++) {
                    Vertex *v1 = g.findVertex(final_path[i - 1]);
                    for (Edge *e : v1->getAdj()) {
                        if (e->getDest()->getId() == final_path[i]) cost += e->getWeight();
                    }
                }
                cout << "Distance: " << cost << endl;

                }
                end = clock();
                cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
                wait();
                break;
            default:
                cout << "Invalid option! Try again" << endl;
                wait();
                choice = 1;
                break;
        }
    }while (choice != 0);
}







void read_toy_graph(const string& file, Graph &g){
    ifstream toy_graph_file(file);

    g.dist = vvd(15, vd(15, 0));

    string line;
    getline(toy_graph_file, line);
    while(getline(toy_graph_file, line)){

        stringstream ss(line);
        string origem, destino, distancia;

        getline(ss, origem, ',');
        getline(ss, destino, ',');
        getline(ss, distancia, ',');

        if(!g.findVertex(stoi(origem))) {
            g.addVertex(stoi(origem));
        }
        if(!g.findVertex(stoi(destino))) {
            g.addVertex(stoi(destino));
        }

        g.addEdge(stoi(origem), stoi(destino), stod(distancia));
    }
}

void read_real_world_graph_nodes(const string& file, Graph &g){
    ifstream real_world_file(file);

    string line;

    getline(real_world_file, line);
    while(getline(real_world_file, line)){

        stringstream ss(line);
        string id, longitude, latitude;

        getline(ss, id, ',');
        getline(ss, longitude, ',');
        getline(ss, latitude, ',');

        g.addVertex(stoi(id));
        nodes.emplace_back(stoi(id), stod(latitude), stod(longitude));
    }
    g.dist = vvd(g.getNumVertex(), vd(g.getNumVertex(), 0));
}

void read_real_world_graph_edges(const string& file, Graph &g){
    ifstream real_world_file(file);

    string line;

    getline(real_world_file, line);
    while(getline(real_world_file, line)){

        stringstream ss(line);
        string origem, destino, distancia;

        getline(ss, origem, ',');
        getline(ss, destino, ',');
        getline(ss, distancia, ',');

        g.addEdge(stoi(origem), stoi(destino), stod(distancia));
    }
}

void read_fully_connected_graph(const string& file, Graph &g){
    ifstream fully_connected_graph_file(file);
    int n;

    if(file == "../data/Extra_Fully_Connected_Graphs/edges_25.csv" || file == "../data/Extra_Fully_Connected_Graphs/edges_50.csv" || file == "../data/Extra_Fully_Connected_Graphs/edges_75.csv") {
        n = stoi(file.substr(43, 2));
    }
    else {
        n = stoi(file.substr(43, 3));
    }

    g.dist = vvd(n, vd(n, 0));

    string line;
    while(getline(fully_connected_graph_file, line)){

        stringstream ss(line);
        string origem, destino, distancia;

        getline(ss, origem, ',');
        getline(ss, destino, ',');
        getline(ss, distancia, ',');

        if (!g.findVertex(stoi(origem))) {
            g.addVertex(stoi(origem));
        }
        if (!g.findVertex(stoi(destino))) {
            g.addVertex(stoi(destino));
        }

        g.addEdge(stoi(origem), stoi(destino), stod(distancia));
    }
}

double spanningTreeCost(const std::vector<Vertex *> &ans) {
    double ret = 0;
    for(const Vertex *v: ans){
        if(v->getPath() == nullptr) continue;
        const Vertex *u = v->getPath()->getOrig();
        for(const auto e: u->getAdj()){
            if(e->getDest()->getId() == v->getId()){
                ret += e->getWeight();
                break;
            }
        }
    }
    return ret;
}

double haversine(double lat1, double lon1, double lat2, double lon2) {
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;

    lat1 = lat1 * M_PI / 180.0;
    lat2 = lat2 * M_PI / 180.0;

    double a = pow(sin(dLat / 2), 2) +
               pow(sin(dLon / 2), 2) * cos(lat1) * cos(lat2);
    double rad = 6371000;
    double c = 2 * asin(sqrt(a));
    return rad * c;
}

double dist(int i, int j, Graph &g) {
    for(auto edge : g.findVertex(i)->getAdj()) {
        if(edge->getDest()->getId() == j) {
            return edge->getWeight();
        }
    }
    return haversine(nodes[i].lat, nodes[i].lon, nodes[j].lat, nodes[j].lon);
}

vector<int> euler_path(Graph &g) {
    vector<int> euler_path;

    queue<int> q;
    q.push(0);


    while (!q.empty()) {
        int v = q.front();
        Vertex *u = g.findVertex(v);
        Edge *e = nullptr;

        for (Edge *edge : u->getAdj()) {
            if (!edge->isSelected()) {
                e = edge;
                break;
            }
        }
        if (e != nullptr) {
            e->setSelected(true);
            q.push(e->getDest()->getId());
        } else {
            q.pop();
            euler_path.push_back(v);
        }
    }
    
    return euler_path;
}

vector<pair<int, int>> matching(vector<int> &odd, Graph &g) {
    std::vector<pair<int, int>> matches;

    for (int i : odd) g.findVertex(i)->setVisited(false);

    double dist = INF;
    int dest = -1;
    for (int i = 0; i < odd.size(); i++) {
        Vertex *v1 = g.findVertex(odd[i]);
        for (int j = 0; j < odd.size(); j++) {
            Vertex *v2 = g.findVertex(odd[j]);
            if (i == j || v2->isVisited()) continue;
            for (Edge *e: v1->getAdj()) {
                if (e->getDest()->getId() == v2->getId()) {
                    if (e->getWeight() < dist) {
                        dist = e->getWeight();
                        dest = v2->getId();
                        v2->setVisited(true);
                    }
                }
            }
        }
        matches.emplace_back(v1->getId(), dest);
    }

    return matches;
}