#include <cstdio>
#include <fstream>
#include <iostream>
#include <time.h>
//#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <math.h>
#include <float.h>
#include <map>
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>
#include<sstream>
#include <list>
#include <stack>


#define numberOfVertex  9205
#define Max_Iteration_Number 10000
#define Alpha 0.15
#define END_WEIGHT 1e-3
#define InitPageRankValue 6

using namespace std;
static int CPUiter = 0;

class Graph
{
    int V;    // No. of vertices
    list<int> *adj;    // An array of adjacency lists

    // Fills Stack with vertices (in increasing order of finishing
    // times). The top element of stack has the maximum finishing 
    // time
    void fillOrder(int v, bool visited[], stack<int> &Stack);

    // A recursive function to print DFS starting from v
    void DFSUtil(int v, bool visited[]);
public:
    Graph(int V);
    void addEdge(int v, int w);

    // The main function that finds and prints strongly connected
    // components
    void printSCCs();

    // Function that returns reverse (or transpose) of this graph
    Graph getTranspose();
};
