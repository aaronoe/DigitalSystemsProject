
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>

u_int32_t INT;

// Matrikelnummer 582595 ; Aaron Oertel <oertelaa@informatik.hu-berlin.de>

struct InputLine {
    int from, to;
    long weight;
    struct InputLine* next;
};


struct InputLinkedList {
    struct InputLine *head;
    int count;
};

// A structure to represent a node in adjacency list
struct AdjListNode {
    int dest;
    long weight;
    struct AdjListNode* next;
};

// A structure to represent an adjacency liat
struct AdjList {
    struct AdjListNode *head;  // pointer to head node of list
};

struct CampList {
    struct CampElement* head;
    int counter;
};

struct CampElement {
    int nodeNumber;
    struct CampElement* nextCamp;
};

// A structure to represent a graph. A graph is an array of adjacency lists.
// Size of array will be V (number of vertices in graph)
struct Graph {
    int V;
    struct AdjList* array;
};

// A utility function to create a new adjacency list node
struct AdjListNode* newAdjListNode(int dest, long weight) {
    struct AdjListNode* newNode =
            (struct AdjListNode*) malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->weight = weight;
    newNode->next = NULL;
    return newNode;
}

// A utility function that creates a graph of V vertices
struct Graph* createGraph(int V) {

    struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
    graph->V = V;

    // Create an array of adjacency lists.  Size of array will be V
    graph->array = (struct AdjList*) malloc(V * sizeof(struct AdjList));

    // Initialize each adjacency list as empty by making head as NULL
    for (int i = 0; i < V; ++i)
        graph->array[i].head = NULL;

    return graph;
}

// Adds an edge to an undirected graph
void addEdge(struct Graph* graph, int src, int dest, long weight) {

    // Add an edge from src to dest.  A new node is added to the adjacency
    // list of src.  The node is added at the begining
    struct AdjListNode* newNode = newAdjListNode(dest, weight);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
}

// Structure to represent a min heap node
struct MinHeapNode {
    int  v;
    long dist;
};

// Structure to represent a min heap
struct MinHeap {
    int size;      // Number of heap nodes present currently
    int capacity;  // Capacity of min heap
    int *pos;     // This is needed for decreaseKey()
    struct MinHeapNode **array;
};

// A utility function to create a new Min Heap Node
struct MinHeapNode* newMinHeapNode(int v, long dist) {
    struct MinHeapNode* minHeapNode =
            (struct MinHeapNode*) malloc(sizeof(struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    return minHeapNode;
}

// A utility function to create a Min Heap
struct MinHeap* createMinHeap(int capacity) {
    struct MinHeap* minHeap =
            (struct MinHeap*) malloc(sizeof(struct MinHeap));
    minHeap->pos = (int *)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array =
            (struct MinHeapNode**) malloc(capacity * sizeof(struct MinHeapNode*));
    return minHeap;
}

// A utility function to swap two nodes of min heap. Needed for min heapify
void swapMinHeapNode(struct MinHeapNode** a, struct MinHeapNode** b) {
    struct MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

// A standard function to heapify at given idx
// This function also updates position of nodes when they are swapped.
// Position is needed for decreaseKey()
void minHeapify(struct MinHeap* minHeap, int idx) {
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < minHeap->size &&
        minHeap->array[left]->dist < minHeap->array[smallest]->dist )
        smallest = left;

    if (right < minHeap->size &&
        minHeap->array[right]->dist < minHeap->array[smallest]->dist )
        smallest = right;

    if (smallest != idx) {
        // The nodes to be swapped in min heap
        struct MinHeapNode *smallestNode = minHeap->array[smallest];
        struct MinHeapNode *idxNode = minHeap->array[idx];

        // Swap positions
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;

        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);

        minHeapify(minHeap, smallest);
    }
}

// A utility function to check if the given minHeap is ampty or not
int isEmpty(struct MinHeap* minHeap) {
    return minHeap->size == 0;
}

// Standard function to extract minimum node from heap
struct MinHeapNode* extractMin(struct MinHeap* minHeap) {

    if (isEmpty(minHeap))
        return NULL;

    // Store the root node
    struct MinHeapNode* root = minHeap->array[0];

    // Replace root node with last node
    struct MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    // Update position of last node
    minHeap->pos[root->v] = minHeap->size-1;
    minHeap->pos[lastNode->v] = 0;

    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}

// Function to decreasy dist value of a given vertex v. This function
// uses pos[] of min heap to get the current index of node in min heap
void decreaseKey(struct MinHeap* minHeap, int v, long dist) {

    // Get the index of v in  heap array
    int i = minHeap->pos[v];

    // Get the node and update its dist value
    minHeap->array[i]->dist = dist;

    // Travel up while the complete tree is not hepified.
    // This is a O(Logn) loop
    while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist) {

        // Swap this node with its parent
        minHeap->pos[minHeap->array[i]->v] = (i-1)/2;
        minHeap->pos[minHeap->array[(i-1)/2]->v] = i;
        swapMinHeapNode(&minHeap->array[i],  &minHeap->array[(i - 1) / 2]);

        // move to parent index
        i = (i - 1) / 2;
    }
}

// A utility function to check if a given vertex
// 'v' is in min heap or not
bool isInMinHeap(struct MinHeap *minHeap, int v) {

    if (minHeap->pos[v] < minHeap->size)
        return true;
    return false;
}


// The main function that calulates distances of shortest paths from src to all
// vertices. It is a O(ELogV) function
long *dijkstra(struct Graph *graph, int src, long maxWeight) {

    int V = graph->V;// Get the number of vertices in graph
    long *dist = malloc(V * sizeof(long*));      // dist values used to pick minimum weight edge in cut

    // minHeap represents set E
    struct MinHeap* minHeap = createMinHeap(V);

    // Initialize min heap with all vertices. dist value of all vertices
    for (int v = 0; v < V; ++v) {
        dist[v] = LONG_MAX;
        minHeap->array[v] = newMinHeapNode(v, dist[v]);
        minHeap->pos[v] = v;
    }

    // Make dist value of src vertex as 0 so that it is extracted first
    dist[src] = 0;
    decreaseKey(minHeap, src, dist[src]);

    // Initially size of min heap is equal to V
    minHeap->size = V;

    // In the followin loop, min heap contains all nodes
    // whose shortest distance is not yet finalized.
    while (!isEmpty(minHeap)) {
        // Extract the vertex with minimum distance value
        struct MinHeapNode* minHeapNode = extractMin(minHeap);
        int u = minHeapNode->v;


        // Store the extracted vertex number
        // Traverse through all adjacent vertices of u (the extracted
        // vertex) and update their distance values
        struct AdjListNode* pCrawl = graph->array[u].head;
        while (pCrawl != NULL) {

            int v = pCrawl->dest;

            // If shortest distance to v is not finalized yet, and distance to v
            // through u is less than its previously calculated distance
            if (isInMinHeap(minHeap, v) && dist[u] != LONG_MAX &&
                pCrawl->weight + dist[u] < dist[v]) {

                dist[v] = dist[u] + pCrawl->weight;

                // update distance value in min heap also
                decreaseKey(minHeap, v, dist[v]);

            }

            pCrawl = pCrawl->next;

        }

        free(minHeapNode);

    }


    free(minHeap->array);
    free(minHeap->pos);
    free(minHeap);
    // print the calculated shortest distances
    return dist;
}


void insertLineIntoLinkedList(struct InputLinkedList* linkedList, int from, int to, long weight) {
    struct InputLine* line = (struct InputLine*) malloc(sizeof(struct InputLine));
    line->from = from;
    line->to = to;
    line->weight = weight;

    line->next = linkedList->head;
    linkedList->head = line;
    linkedList->count++;
}

void insertShIntoShList(struct CampList *shList, int nodeNumber) {
    struct CampElement* sh = (struct CampElement*) malloc(sizeof(struct CampElement));
    sh->nextCamp = shList->head;
    sh->nodeNumber = nodeNumber;
    shList->head = sh;
    shList->counter++;
}


bool searchShInList(struct CampList* shList, int nodeNumber) {

    struct CampElement* currentSh = shList->head;

    // Construct graph from saved input
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
int freeInput(struct InputLinkedList* inputLinkedList) {

    struct InputLine* currentLine = inputLinkedList->head;
    struct InputLine* tempLine;

    // Construct graph from saved input and free input line
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
void freeGraph(struct Graph* graph, int V) {

    for (int j = 0; j < V; ++j) {

        // we are done with graph so we should free it
        struct AdjListNode* adjListNode = graph->array[j].head;
        struct AdjListNode* tempAdjListNode;

        while (adjListNode != NULL) {
            tempAdjListNode = adjListNode->next;
            free(adjListNode);
            adjListNode = tempAdjListNode;
        }

    }

    free(graph->array);
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
int main(int argc, char *argv[]) {


    // Anzahl der Knoten
    int V = 0;
    int startNode = -1;
    int endNode = -1;
    long maxWeight = 0;
    bool firstLine = true;

    // Input lesen
    int fromNode, toNode;
    long weight;
    char line[64];

    struct InputLinkedList* inputLinkedList = (struct InputLinkedList*) malloc(sizeof(struct InputLinkedList));
    inputLinkedList->count = 0;
    inputLinkedList->head = NULL;

    struct CampList* campList = (struct CampList*) malloc(sizeof(struct CampList));
    campList->counter = 0;
    campList->head = NULL;


    /* read at least 63 characters or unitl newline charater is encountered with */
    /*    fgets(line, sizeof line, stdin) */
    /* if the first character is a newline, then it's an empty line of input */
    while ((fgets(line, sizeof line, stdin) != NULL) && (line[0] != '\n')) {
        /* parse the read line with sscanf */
        if (sscanf(line, "%d%d%li", &fromNode, &toNode, &weight) == 3) {
            //printf("%d %d %li\n", fromNode, toNode, weight);
            if (firstLine) {
                startNode = fromNode;
                endNode = toNode;
                maxWeight = weight;
                firstLine = false;
            } else {
                insertLineIntoLinkedList(inputLinkedList, fromNode, toNode, weight);
                // Anzahl der Knoten zählen
                if (fromNode > V) {
                    V = fromNode;
                }
                if (toNode > V) {
                    V = toNode;
                }
            }
            fflush(stdout);
        } else {
            //printf("%d\n", fromNode);
            insertShIntoShList(campList, fromNode);
        }
    }


    V++;


    if (firstLine == true) {
        printf("No start node, end node or max weight specified\n");
        exitEarly(inputLinkedList, campList);
    }

    if (startNode == endNode) {
        printf("%d\n", startNode);
        exitEarly(inputLinkedList, campList);
    }

    if (endNode > V) {
        exitEarly(inputLinkedList, campList);
    }

    struct Graph* graph = createGraph(V);
    struct Graph* reversedGraph = createGraph(V);

    struct InputLine* currentLine = inputLinkedList->head;
    struct InputLine* tempLine;

    // Construct graph from saved input and free input line
    while (currentLine != NULL) {
        addEdge(graph, currentLine->from, currentLine->to, currentLine->weight);
        addEdge(reversedGraph, currentLine->to, currentLine->from, currentLine->weight);
        tempLine = currentLine->next;
        free(currentLine);
        currentLine = tempLine;
    }

    free(inputLinkedList);

    long *dist = dijkstra(graph, startNode, maxWeight);

    freeGraph(graph, V);

    long *reverseDist = dijkstra(reversedGraph, endNode, maxWeight);

    freeGraph(reversedGraph, V);


    for (int i = 0; i < V; ++i) {

        if (dist[i] <= maxWeight && reverseDist[i] <= maxWeight) {
            if (searchShInList(campList, i)) {
                printf("%d\n", i);
            }
        }
    }

    free(dist);
    free(reverseDist);

    freeCampList(campList);

    return 0;
}
