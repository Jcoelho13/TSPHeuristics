#ifndef DATP2_NODE_H
#define DATP2_NODE_H

/**
 * Class that represents a node in the graph.
 */

class Node {
public:
    /**
     * Node's ID.
     */
    int id;
    /**
     * Node's latitude and longitude.
     */
    double lat, lon;

    /**
     * Node's constructor.
     * @param id Node's ID.
     * @param lat Node's latitude.
     * @param lon Node's longitude.
     */
    Node(int id, double lat, double lon);
};

#endif //DATP2_NODE_H
