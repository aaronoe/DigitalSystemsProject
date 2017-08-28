
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>

typedef uint32_t Int;
typedef uint64_t Long;
typedef int64_t LongSigned;

// Matrikelnummer 582595 ; Aaron Oertel <oertelaa@informatik.hu-berlin.de>

/**
 * Used to temporarily store an input line of the kind:
 * 0 1 3924456639
 * All numbers are between 0 and 4,000,000,000 and hence we need at least 32 bit to store them.
 * The from and to values are stored as unsigned int and the weight as an unsigned long
 */
struct InputLine {
    Int from, to;
    Long weight;
    struct InputLine* next;
};


/**
 * Used to store the lines read from stdin
 */
struct InputLinkedList {
    struct InputLine *head;
    Int count;
};

/**
 * one element of the AdjList, is used for a node's adjacency list and each element
 * contains information on an edge from that node to a different one
 */
struct AdjListNode {
    Int destination;
    Long weight;
    struct AdjListNode* nextNode;
};

/**
 * Used by the dijkstra algorithm to compute the distances
 */
struct AdjList {
    struct AdjListNode* head;  // pointer to head node of list
};

/**
 * An Adjacency List to store the camps from the input
 */
struct CampList {
    struct CampElement* head;
    Int counter;
};

/**
 * One element of the CampList. It is used to store a single node number from the input
 */
struct CampElement {
    Int nodeNumber;
    struct CampElement* nextCamp;
};


/**
 * Data structure to represent a graph used by Dijkstra's algorithm.
 * For each node we store an AdjList of this node's outgoing edges
 */
struct Graph {
    Int numberOfNodes;
    struct AdjList* adjListArray;
};


/**
 * Creates a new AdjListNode representing an edge from a node
 * @param destination of the edge
 * @param weight of the edge
 * @return the AdjListNode containing the edge's information
 */
struct AdjListNode* createAdjListNode(Int destination, Long weight) {
    struct AdjListNode* newNode = (struct AdjListNode*) malloc(sizeof(struct AdjListNode));
    newNode->destination = destination;
    newNode->weight = weight;
    newNode->nextNode = NULL;
    return newNode;
}


/**
 * Creates a graph of nodeCount edges, this number is needed to create the array
 * of AdjLists used to store the edge information
 * @param nodeCount number of node's the graph contains
 * @return a pointer to a Graph structure
 */
struct Graph* createGraph(Int nodeCount) {

    struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));

    graph->numberOfNodes = nodeCount;

    // allocate memory for nodeCount AdjLists
    graph->adjListArray = (struct AdjList*) malloc(nodeCount * sizeof(struct AdjList));

    // Make the first element of each AdjList null
    for (Int i = 0; i < nodeCount; ++i) {
        graph->adjListArray[i].head = NULL;
    }

    return graph;
}


/**
 * Adds an edge to the graph at the given source's adjacency list's head
 * @param graph where the edge gets added to
 * @param source of the edge
 * @param destination of the edge
 * @param weight of the edge
 */
void addEdgeToGraph(struct Graph *graph, Int source, Int destination, Long weight) {
    // Create a new node and add it as the head of source's AdjList
    struct AdjListNode* newNode = createAdjListNode(destination, weight);
    newNode->nextNode = graph->adjListArray[source].head;
    graph->adjListArray[source].head = newNode;
}


/**
 * Represents a node in a MinHeap
 */
struct MinHeapNode {
    Int node;
    Long distance;
};


/**
 * Structure to represent a MinHeap
 */
struct MinHeap {
    Int size;
    Int* pos;
    struct MinHeapNode** array;
};

/**
 * Creates a node for the MinHeap
 * @param node number
 * @param distance to the node from source
 * @return the created node
 */
struct MinHeapNode* newMinHeapNode(Int node, Long distance) {
    struct MinHeapNode* minHeapNode = (struct MinHeapNode*) malloc(sizeof(struct MinHeapNode));
    minHeapNode->node = node;
    minHeapNode->distance = distance;
    return minHeapNode;
}

/**
 * Creates new MinHeap structure
 * @param capacity which is needed to allocate memory for storing the nodes
 * @return
 */
struct MinHeap* createMinHeap(Int capacity) {
    struct MinHeap* minHeap = (struct MinHeap*) malloc(sizeof(struct MinHeap));
    minHeap->pos = (Int*) malloc(capacity * sizeof(Int));
    minHeap->size = 0;
    minHeap->array = (struct MinHeapNode**) malloc(capacity * sizeof(struct MinHeapNode*));
    return minHeap;
}

/**
 * Swaps to Nodes in the MinHeap, used for Heapifying
 * @param one first node to be swapped
 * @param two second node to be swapped
 */
void swapMinHeapNode(struct MinHeapNode** one, struct MinHeapNode** two) {
    struct MinHeapNode* t = *one;
    *one = *two;
    *two = t;
}


/**
 * This utility method is responsible for recursively updating the heap to make sure
 * that both heap constraints are met
 * @param minHeap to be heapified
 * @param pos position to start
 */
void minHeapify(struct MinHeap* minHeap, Int pos) {

    Int smalledNode, leftChild, rightChild;
    smalledNode = pos;
    // Get the children of pos
    leftChild = 2 * pos + 1;
    rightChild = 2 * pos + 2;

    // Pick smalledNode node out of pos and it's children
    if (leftChild < minHeap->size && minHeap->array[leftChild]->distance < minHeap->array[smalledNode]->distance ) {
        smalledNode = leftChild;
    }

    if (rightChild < minHeap->size && minHeap->array[rightChild]->distance < minHeap->array[smalledNode]->distance ) {
        smalledNode = rightChild;
    }

    // Obviously we only need to swap nodes if pos is actually greater than it's children
    if (smalledNode != pos) {
        // The nodes to be swapped in min heap
        struct MinHeapNode *smallestNode = minHeap->array[smalledNode];
        struct MinHeapNode *idxNode = minHeap->array[pos];

        // swap positions of the smalledNode node and the input node
        minHeap->pos[smallestNode->node] = pos;
        minHeap->pos[idxNode->node] = smalledNode;

        // swap nodes
        swapMinHeapNode(&minHeap->array[smalledNode], &minHeap->array[pos]);

        minHeapify(minHeap, smalledNode);
    }
}


/**
 * Utility method which checks if heap is empty
 * @param minHeap
 * @return
 */
int isHeapEmpty(struct MinHeap *minHeap) {
    return (minHeap->size == 0);
}


/**
 * Method which extract minimum from heap and then updates the heap
 * using minHeapify() to make sure heap constraints are still met
 * @param minHeap
 * @return minimum element of the given MinHeap
 */
struct MinHeapNode* extractMin(struct MinHeap* minHeap) {

    if (isHeapEmpty(minHeap)) {
        return NULL;
    }

    struct MinHeapNode* rootNode = minHeap->array[0];

    // use last node of the heap-array to remove the first item
    struct MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    // we need to decrease the size of the heap and
    // update the new first item's position as 0
    minHeap->pos[rootNode->node] = minHeap->size-1;
    minHeap->pos[lastNode->node] = 0;
    --minHeap->size;
    // update the heap from the top to make sure that heap constraints are met
    minHeapify(minHeap, 0);

    return rootNode;
}


/**
 * Method to decrease distance to node from source
 * Updating the distance might require moving the node up in the heap
 * @param minHeap
 * @param nodeNumber for node which distance is updated
 * @param newDistance the new distance of the node
 */
void decreaseKey(struct MinHeap* minHeap, Int nodeNumber, Long newDistance) {

    // get the position of the Node in the heap-array
    Int index = minHeap->pos[nodeNumber];

    // update the node's distance in the heap-array
    minHeap->array[index]->distance = newDistance;

    // since we changed the node's distance from source we might have to move the node
    // up in the heap using swap operations if the parent of the element is grater
    while (index && minHeap->array[index]->distance < minHeap->array[(index - 1) / 2]->distance) {

        minHeap->pos[minHeap->array[index]->node] = (index-1)/2;
        minHeap->pos[minHeap->array[(index-1)/2]->node] = index;
        swapMinHeapNode(&minHeap->array[index],  &minHeap->array[(index - 1) / 2]);

        // this is the node's parent index
        index = (index - 1) / 2;
    }
}

/**
 * Checks if node is in MinHeap
 * @param minHeap
 * @param node
 * @return true or false
 */
bool isNodeInMinHeap(struct MinHeap *minHeap, Int node) {

    if (minHeap->pos[node] < minHeap->size) {
        return true;
    }

    return false;
}


/**
 * Main method to compute distances to all nodes from given source node
 * @param graph
 * @param sourceNode
 * @return
 */
Long *dijkstra(struct Graph *graph, Int sourceNode) {

    Int V = graph->numberOfNodes;

    // Array to store the distances to all V nodes
    Long* dist = malloc(V * sizeof(Int*));

    // create heap to be used as a priority queue to get the
    // closest node by distance
    struct MinHeap* minHeap = createMinHeap(V);

    // we initialize the distance for all nodes as ULONG_MAX
    for (Int v = 0; v < V; ++v) {
        dist[v] = ULONG_MAX;
        minHeap->array[v] = newMinHeapNode(v, dist[v]);
        minHeap->pos[v] = v;
    }

    // sourceNode's distance needs to be 0 because that's where we start
    dist[sourceNode] = 0;
    decreaseKey(minHeap, sourceNode, dist[sourceNode]);

    minHeap->size = V;

    // In this loop we extract the node with the smallest distance from the heap
    // It is removed from the heap and therefore finalized in dist
    while (!isHeapEmpty(minHeap)) {

        struct MinHeapNode* minHeapNode = extractMin(minHeap);
        Int nodeNumber = minHeapNode->node;

        // From the extracted minimum element, get all outgoing edges using
        // it's adjacency list and then update their distance values
        struct AdjListNode* outgoingVertex = graph->adjListArray[nodeNumber].head;
        while (outgoingVertex != NULL) {

            Int v = outgoingVertex->destination;

            // We only update the distance if the node has not been finalized and
            // if the new distance is actually smaller than the current one
            if (isNodeInMinHeap(minHeap, v) && dist[nodeNumber] != ULONG_MAX &&
                outgoingVertex->weight + dist[nodeNumber] < dist[v]) {

                dist[v] = dist[nodeNumber] + outgoingVertex->weight;

                // values in heap needs to be updated as well
                decreaseKey(minHeap, v, dist[v]);

            }

            // move to the next node of the AdjList for the next iteration
            outgoingVertex = outgoingVertex->nextNode;
        }
        // Since we lose the pointer of the MinHeapNode here (since we removed it from the heap)
        // we need to free the allocated memory here
        free(minHeapNode);
    }

    // Free all MinHeap related allocated memory since it's only used in this scope
    free(minHeap->array);
    free(minHeap->pos);
    free(minHeap);

    return dist;
}


/**
 * Utility method to save input line in linked list
 * Since the element is saved as the first element this operation is in O(1)
 * @param linkedList
 * @param from
 * @param to
 * @param weight
 */
void insertLineIntoLinkedList(struct InputLinkedList* linkedList, Int from, Int to, Long weight) {
    struct InputLine* line = (struct InputLine*) malloc(sizeof(struct InputLine));
    line->from = from;
    line->to = to;
    line->weight = weight;

    line->next = linkedList->head;
    linkedList->head = line;
    linkedList->count++;
}

/**
 * Inserts node into the linked list of camps
 * @param campList
 * @param nodeNumber
 */
void insertCampIntoList(struct CampList* campList, Int nodeNumber) {
    struct CampElement* sh = (struct CampElement*) malloc(sizeof(struct CampElement));
    sh->nextCamp = campList->head;
    sh->nodeNumber = nodeNumber;
    campList->head = sh;
    campList->counter++;
}


/**
 * Search if node is in list of camps
 * @param campList
 * @param nodeNumber
 * @return true or false
 */
bool searchCampInList(struct CampList* campList, Int nodeNumber) {

    struct CampElement* currentSh = campList->head;

    while (currentSh != NULL) {
        if (nodeNumber == currentSh->nodeNumber) {
            return true;
        }
        currentSh = currentSh->nextCamp;
    }

    return false;
}


/**
 * Frees each line of the input saved as an InputLinkedList
 * @param inputLinkedList
 * @return
 */
void freeInput(struct InputLinkedList* inputLinkedList) {

    struct InputLine* currentLine = inputLinkedList->head;
    struct InputLine* tempLine;

    while (currentLine != NULL) {
        tempLine = currentLine->next;
        free(currentLine);
        currentLine = tempLine;
    }

    free(inputLinkedList);
}


/**
 * Free the graphs including the array of linked lists and their items
 * @param graph - containing the graph info as an adjacency list (array of linked lists) -> hence we need to iterate through the
 * linked list for each of the V nodes and free each element of the adjacency list
 * @param V - this is the number of nodes in the graph -> used to iterate over the array of linked lists
 */
void freeGraph(struct Graph* graph, Int V) {

    for (Int j = 0; j < V; ++j) {

        struct AdjListNode* adjListNode = graph->adjListArray[j].head;
        struct AdjListNode* tempAdjListNode;

        while (adjListNode != NULL) {
            tempAdjListNode = adjListNode->nextNode;
            free(adjListNode);
            adjListNode = tempAdjListNode;
        }

    }

    free(graph->adjListArray);
    free(graph);
}


/**
 * Free the list of camps saved from the input
 * @param campList - read from the input
 */
void freeCampList(struct CampList *campList) {
    struct CampElement* elementPointer = campList->head;
    struct CampElement* tempElement;

    while (elementPointer != NULL) {
        tempElement = elementPointer->nextCamp;
        free(elementPointer);
        elementPointer = tempElement;
    }

    free(campList);
}


/**
 * This is used if we need to exit the program early if the input is invalid
 * @param inputLinkedList - which is the input to the program
 * @param campList - which is also input to the program
 */
void exitEarly(struct InputLinkedList* inputLinkedList, struct CampList* campList) {
    freeInput(inputLinkedList);
    freeCampList(campList);
    exit(-1);
}

/**
 *
 * Driver program which has the following functionality :
 * 1. Read the input from stdin and save it in linked lists
 * 2. Validate the input and exit if necessary
 * 3. Compute the distances from the start node to all other nodes using dijkstra's algorithm
 * 4. Compute the distances from the end node to all other nodes using dijkstra's algorithm and the transposed graph
 * 5. Outputs the nodes that fulfill the following criteria :
 *      a) Distance from start node is smaller than max weight
 *      b) Distance from end node is smaller than max weight
 *      c) Is in the linked list of camps stored from the input
 * 6. Frees the memory allocated by malloc
 *
 */
int main() {


    Int nodeCount = 0;
    Int startNode = 0;
    Int endNode = 0;
    Long maxWeight = 0;
    bool firstLine = true;
    Int maxNumber = 4000000000;

    // for reading from stdin
    LongSigned fromNode, toNode;
    LongSigned weight;
    // buffer
    char line[64];
    char eol;

    struct InputLinkedList* inputLinkedList = (struct InputLinkedList*) malloc(sizeof(struct InputLinkedList));
    inputLinkedList->count = 0;
    inputLinkedList->head = NULL;

    struct CampList* campList = (struct CampList*) malloc(sizeof(struct CampList));
    campList->counter = 0;
    campList->head = NULL;


    /*
     * read until linefeed is encountered or the first 63 characters
     * this is more than enough because numbers are at max 10 characters long
     */
    while ((fgets(line, sizeof line, stdin) != NULL) && (line[0] != '\n')) {

        /*
         * Parse the line and extract it's values
         * There are 4 legal cases for what lines can look like
         * 1. 3 numbers and a line feed
         * 2. 3 numbers without a line feed
         * 3. 1 number and a line feed
         * 4. 1 number without a line feed
         *
         * Therefore we need to catch the right case and act accordingly
         * to save the input data
         */
        if (sscanf(line, "%li%li%li", &fromNode, &toNode, &weight) == 3) {

            // Check if the input consists of negative or too large numbers
            if (fromNode >= maxNumber || toNode >= maxNumber || weight >= maxNumber
                || fromNode < 0 || toNode < 0 || weight < 0) {
                printf("Bad input\n");
                exitEarly(inputLinkedList, campList);
            }

            if (sscanf(line, "%li%li%li%c", &fromNode, &toNode, &weight, &eol) == 4) {

                if (eol != '\n') {
                    printf("Bad input\n");
                    exitEarly(inputLinkedList, campList);
                }

            }

            if (firstLine) {
                startNode = (Int) fromNode;
                endNode = (Int) toNode;
                maxWeight = (Long) weight;
                firstLine = false;

            } else {

                insertLineIntoLinkedList(inputLinkedList, (Int) fromNode, (Int) toNode, (Long) weight);

                if (fromNode > nodeCount) {
                    nodeCount = (Int) fromNode;
                }
                if (toNode > nodeCount) {
                    nodeCount = (Int) toNode;
                }
            }

        } else if (sscanf(line, "%li", &fromNode) == 1) {

            if (fromNode >= maxNumber || fromNode < 0) {
                printf("Bad input\n");
                exitEarly(inputLinkedList, campList);
            }

            if (sscanf(line, "%li%c", &fromNode, &eol) == 2) {

                // Case : One Int type number and a character which should be the linefeed

                if (eol != '\n') {
                    printf("Falscher input\n");
                    exitEarly(inputLinkedList, campList);
                }

            }


            if (firstLine) {
                printf("No start node, end node or max weight specified\n");
                exitEarly(inputLinkedList, campList);
            }

            insertCampIntoList(campList, fromNode);

        }
    }


    nodeCount++;

    // if no line has been read the input is faulty
    if (firstLine == true) {
        printf("No start node, end node or max weight specified\n");
        exitEarly(inputLinkedList, campList);
    }

    // TODO: look into this case
    if (startNode == endNode) {
        printf("%u\n", startNode);
        exitEarly(inputLinkedList, campList);
    }

    if (endNode > nodeCount) {
        exitEarly(inputLinkedList, campList);
    }

    struct Graph* graph = createGraph(nodeCount);
    struct Graph* reversedGraph = createGraph(nodeCount);

    struct InputLine* currentLine = inputLinkedList->head;
    struct InputLine* tempLine;

    // Construct graph from saved input and free input line
    while (currentLine != NULL) {
        addEdgeToGraph(graph, currentLine->from, currentLine->to, currentLine->weight);
        addEdgeToGraph(reversedGraph, currentLine->to, currentLine->from, currentLine->weight);
        tempLine = currentLine->next;
        free(currentLine);
        currentLine = tempLine;
    }

    free(inputLinkedList);

    Long *dist = dijkstra(graph, (Int) startNode);

    freeGraph(graph, nodeCount);

    Long *reverseDist = dijkstra(reversedGraph, (Int) endNode);

    freeGraph(reversedGraph, nodeCount);


    /*
     *  for every node if it's distance from source in the normal and transposed
     *  is smaller than maxWeight and it's a camp it needs to be printed
     */
    for (Int i = 0; i < nodeCount; ++i) {

        if (dist[i] <= maxWeight && reverseDist[i] <= maxWeight) {
            if (searchCampInList(campList, i)) {
                printf("%d\n", i);
            }
        }
    }

    free(dist);
    free(reverseDist);

    freeCampList(campList);

    return 0;
}
