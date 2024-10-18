//Copyright: Will James, 2016. All rights reserved.
#include <algorithm>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <deque>
#include <assert.h>
#include <getopt.h>
#include <stdlib.h>
#include <queue>
#include <functional>
#include <unordered_map>
#include <iterator>
#include <utility>
#include <iomanip>
#include <ctime>
//#include "antcolony.h"

 #define DDEBUG
 #ifdef DEBUG
 #define _(args) args
 #else
 #define _(args)
 #endif
 
//#define _(args) args
/*
 #ifdef __APPLE__
 freopen("inputFilename.txt", "r", stdin);
 freopen("outputFilename.txt", "w", stdout);
 #endif
 */


using namespace std;

enum Mode {MST,OPTTSP,FASTTSP,ANTCOL};
struct MstNode{
    
    long int x,y;
    double edgeLength = numeric_limits<double>::max();
    bool visited = false;
    long int parent;
    char terrain;
    
};

struct TspNode{
    
    int x,y;
    int index;
    
};

double getMSTDist(const MstNode &n1, const MstNode &n2){
    
    return sqrt(pow((n1.x - n2.x), 2.0)+ pow((n1.y - n2.y), 2.0));
    
}

//Use the pythagorean formula to calculate the distance between the two nodes
double getTSPDist(const TspNode &n1, const TspNode &n2){
    
    return sqrt(pow((n1.x - n2.x), 2.0)+ pow((n1.y - n2.y), 2.0));
    
}

double getTotalTSPDist(const deque<TspNode> & v){
    
    double dist = 0.00;
    
    for (unsigned long i = 1; i < v.size(); ++i){
        
        dist += getTSPDist(v[i-1], v[i]);
        
        if (i == v.size() -1 ){
            
            dist += getTSPDist(v[i], v[0]);
            
        }
        
    }
    return dist;
    
}

double getTotalTSPDistVector(const vector<TspNode> & v){
    
    double dist = 0.00;
    
    for (unsigned long i = 1; i < v.size(); ++i){
        
        dist += getTSPDist(v[i-1], v[i]);
        
        if (i == v.size() -1 ){
            
            dist += getTSPDist(v[i], v[0]);
            
        }
        
    }
    return dist;
    
}

double getDistWithoutLoop(const deque<TspNode> & v){
    
    double dist = 0.00;
    
    for (unsigned long i = 1; i < v.size(); ++i){
        
        dist += getTSPDist(v[i-1], v[i]);
        
        
    }
    return dist;
    
}

void printTSPSol(const deque<TspNode> & tspDeq, int numNodes){
    
    
    //cout<< tspDeq.size()<< endl;
    
    double totalDist = getTotalTSPDist(tspDeq);
    ostringstream oss;
    oss.str() = "";
    oss << fixed << setprecision(2);
    oss<<totalDist << '\n';
    
    
    for (int i = 0; i < numNodes; ++i){
        
        oss<< tspDeq[i].index;
        if (i < numNodes-1){
            oss<< ' ';
        }
    }
    oss << '\n';
    cout << setprecision(2) << fixed;
    cout<<oss.str();
    //cout << endl;
    
}

void printTSPSolVector(const vector<TspNode> & tspDeq, int numNodes){
    
    
    //cout<< tspDeq.size()<< endl;
    
    double totalDist = getTotalTSPDistVector(tspDeq);
    ostringstream oss;
    oss.str() = "";
    oss << fixed << setprecision(2);
    oss<<totalDist << '\n';
    
    
    for (int i = 0; i < numNodes; ++i){
        
        oss<< tspDeq[i].index;
        if (i < numNodes-1){
            oss<< ' ';
        }
    }
    oss << '\n';
    cout << setprecision(2) << fixed;
    cout<<oss.str();
    //cout << endl;
    
}

bool promising(const deque<TspNode>& path, const deque<TspNode>& unvisited, double upperBound){
    
    double mstDist = 0.00;
    int numVisited = 1;
    int previous = 0;
    int numNodes = (int)unvisited.size();
    
    
    double tempDist = numeric_limits<double>::max();
    //int nearestNodeIndex = 0;
    vector<bool> tfv(numNodes);
    vector<double> dv(numNodes);
    fill(dv.begin(), dv.end(), numeric_limits<double>::max());
    
    if(!tfv.empty() && !dv.empty()){
        dv.front() = 0;
        tfv.front() = true;
    }
    while (numVisited < numNodes){
        
        for (int i = 1; i < numNodes; ++i){
            
            tempDist = getTSPDist(unvisited[previous], unvisited[i]);
            if(tempDist < dv[i] && tfv[i] == false){
                
                dv[i] = tempDist;
                
            }
            
        }
        tempDist = numeric_limits<double>::max();
        for (int i = 1; i < numNodes; ++i){
            
            
            if (tfv[i] == false && dv[i] < tempDist){
                
                previous = i;
                tempDist = dv[i];
                
            }
            
        }
        
        mstDist += tempDist;
        tfv[previous] = true;
        tempDist = numeric_limits<double>::max();
        
        ++numVisited;
        
    }
    
    double endToMST = 0;
    double startToMST = 0;
    
    if (!path.empty() && !unvisited.empty()){
        
        
        endToMST = numeric_limits<double>::max();
        startToMST = numeric_limits<double>::max();
        for(int i = 0; i < numNodes; ++i){
            double dist = getTSPDist(path.front(), unvisited[i]);
            double dist2 = getTSPDist(path.back(), unvisited[i]);
            if (dist < startToMST){
                startToMST = dist;
                
            }
            if (dist2 < endToMST){
                endToMST = dist2;
                
            }
            
        }
        
        
        if (path.size() == 1){
            
            endToMST = 0.0;
            
        }
    }
    
    
    double pathLength = getDistWithoutLoop(path);
    
    double lowerBound = mstDist + pathLength + endToMST + startToMST;
    
    if (lowerBound >= upperBound){
        
        
        //printTSPSol(path, (int)path.size());
        //printTSPSol(unvisited, (int)unvisited.size());
        //cout<<mstDist<<endl;
        return false;
        
    }
    
    else{
        
        
        return true;
    }
    
}

void genPerms(deque<TspNode>& path, deque<TspNode>& unvisited, deque<TspNode> &finalPath, double & upperBound) {
    
    if (!promising(path, unvisited, upperBound)){
        return;
        
    }
    
    if(unvisited.empty()){
        double newPathLength = getTotalTSPDist(path);
        if (newPathLength < upperBound){
            upperBound = newPathLength;
            finalPath = path;
            //printTSPSol(finalPath, (int)finalPath.size());
        }
        return;
    }
    for(unsigned i = 0; i < unvisited.size(); ++i) {
        path.push_back(unvisited.front());
        unvisited.pop_front();
        genPerms(path, unvisited, finalPath, upperBound);
        unvisited.push_back(path.back());
        path.pop_back();
    } // for
}
void doMST(vector<MstNode> & v, int numNodes){
    
    int previous = 0,bestNode = 0,numVisited = 1;
    double totalDist= 0.0;
    double tempLength= numeric_limits<double>::max();
    if (!v.empty()){
        
        v[0].visited =true;
        v[0].edgeLength = 0;
        
    }
    while (numVisited < numNodes){
        for (int i = 1; i < numNodes; ++i){
            if (!(v[i].terrain == 'L' && v[previous].terrain == 'S')&&
                !(v[i].terrain == 'S' && v[previous].terrain == 'L')&&
                !v[i].visited){
                tempLength = getMSTDist(v[previous], v[i]);
                if (tempLength < v[i].edgeLength){
                    v[i].edgeLength = tempLength;
                    v[i].parent = previous;
                }//if
            }//if
        }//for
        tempLength = numeric_limits<double>::max();
        for (int i = 1; i<numNodes; ++i){
            if (!v[i].visited && v[i].edgeLength < tempLength){
                bestNode = i;
                tempLength = v[i].edgeLength;
            }
        }
        
        previous = bestNode;
        totalDist += tempLength;
        v[bestNode].visited = true;
        tempLength = numeric_limits<double>::max();
        ++numVisited;
    }//while
    ostringstream oss;
    oss.str() = "";
    oss << setprecision(2) << fixed;
    oss<<totalDist << '\n';
    for (int i = 1; i < numNodes; ++i){
        if (i > v[i].parent){
            oss<< v[i].parent << ' '<< i <<'\n';
        }
        else {
            oss<< i << ' '<< v[i].parent <<'\n';
        }
    }
    cout << setprecision(2) << fixed;
    cout<<oss.str();
}
//this is a variation of farthest insertion, or just regular insertion
void doFASTTSP(deque<TspNode> &tspDeq, int numNodes){
    
    
    int nearest = numeric_limits<int>::max();
    double nearDist = numeric_limits<double>::max();
    
    vector<TspNode> path;
    
    path.push_back(tspDeq.front());
    tspDeq.pop_front();
    path.push_back(tspDeq.front());
    tspDeq.pop_front();
    path.push_back(tspDeq.front());
    tspDeq.pop_front();
    
    int n = (int)tspDeq.size();
    
    for(int i = 0; i < n; ++i){
        
        for (int j = 0; j <(int)path.size(); ++j){
            
            double dist = numeric_limits<double>::max();
            if(j == (int)path.size() - 1){
                dist = (getTSPDist(path[j], tspDeq[i]) + getTSPDist(path[0], tspDeq[i])) - getTSPDist(path[j],path[0]);
                
            }
            else{
                dist = (getTSPDist(path[j], tspDeq[i]) + getTSPDist(path[j+1], tspDeq[i])) - getTSPDist(path[j],path[j+1]);
            }
            if(dist < nearDist){
                
                nearDist = dist;
                nearest = j;
                
            }
            
        }
        path.insert(path.begin() + nearest+1, tspDeq[i]);
        nearDist = numeric_limits<double>::max();
        nearest = numeric_limits<int>::max();
    }
    
    
    printTSPSolVector(path, numNodes);
    
}

//this is an implementation of branch and bound for tsp
void doOPTTSP(deque<TspNode> &tspDeq, int numNodes){
    
    double tempDist = numeric_limits<double>::max();
    int nearestNodeIndex = 0;
    
    for(int i = 0; i< numNodes-1; ++i){
        for (int j = i+1; j < numNodes; ++j){
            
            if (getTSPDist(tspDeq[i], tspDeq[j]) < tempDist){
                
                tempDist = getTSPDist(tspDeq[i], tspDeq[j]);
                nearestNodeIndex = j;
            }
            
        }
        swap(tspDeq[i+1], tspDeq[nearestNodeIndex]);
        tempDist = numeric_limits<double>::max();
    }
    
    double upperBound = getTotalTSPDist(tspDeq);
    deque<TspNode> path;
    deque<TspNode> final = tspDeq;
    deque<TspNode> unvisited = tspDeq;
    //printTSPSol(tspDeq, numNodes);
    path.push_front(unvisited.front());
    unvisited.pop_front();
    
    genPerms(path, unvisited, final, upperBound);
    
    printTSPSol(final, numNodes);
    
    
}

class InnerVector {
public:
    vector<double> inner;
    
    //~InnerVector(){}
    InnerVector(int s){
        inner.resize(s);
        
    }
    
    
};

struct Antdata {
    
    deque<int> tour;
    //double tourLength;
    int vertex; // (ant's starting vertex in the graph)
    
};

deque<int> conversionHelper(deque<TspNode> const & t){
    
    deque<int> td;
    
    td.resize(t.size());
    
    for(int i = 0; i < t.size(); ++i){
        
        td[i] = t[i].index;
        
    }
    return td;
}

deque<TspNode> conversionInttoDeq(deque<int> const & t, deque<TspNode> const & tD){
    
    deque<TspNode> tp;
    
    tp.resize(tD.size());
    
    for(int i = 0; i < t.size(); ++i){
        
        tp[i] = tD[t[i]];
        
    }
    return tp;
}


int updateWeights(vector<InnerVector> & wM, const vector<Antdata> & tours, deque<TspNode> const & tpD, int numNodes, int numAnts){
    
    //Here are our parameters for pheremone level,
    //When updating the weights for the edges, we add weight to the edges that ants travel on the most, if those edges correspond to shorter tours
    
    double K = 100; //Pheremone Level
    double P = .2; //Pheremone evaporation constant
    
    double minTourDist = numeric_limits<double>::max();
    int optimalTourIndex = 0;
    
    //iterate through the ants and their tours
    
    for(int i = 0; i < numAnts; ++i){
        double edgeDist = 0;
        double additionalWeight = 0;
        double tourDist = 0;
        
        deque<int> ad1 = tours[i].tour;
        
        deque<TspNode> antITSPNodes = conversionInttoDeq(ad1, tpD);
        
        tourDist = getTotalTSPDist(antITSPNodes);
        
        if (tourDist < minTourDist){
            minTourDist = tourDist;
            optimalTourIndex = i; //ant i had the optimal tour
        }
        
        //iterate through the vertices in ant i's tour, getting the total distance of the tour
        // originall the number of vertices was the same as the number of ants, but I changed it
        
        for (int k = 0; k < numNodes; ++k){
            
            TspNode node1;
            TspNode node2;
            if (k == numNodes-1){
                node1 = antITSPNodes[k];
                node2 = antITSPNodes[0];
                
            }
            else{
                node1 = antITSPNodes[k];
                node2 = antITSPNodes[k+1];
                
            }
            edgeDist = getTSPDist(node1,node2);
            
            //edge distance should determine how much pheremone gets put on each edge, no- equal pheremone gets put on an edge if its chosen
            
            additionalWeight = K/tourDist;
            double initialWeight = wM[node1.index].inner[node2.index];
            
            wM[node1.index].inner[node2.index] = initialWeight * P + additionalWeight;
            
        }//end for (int k = 0; k < numNodes; ++k)
        
    }//end for(int i = 0; i < numAnts; ++i)
    
    return optimalTourIndex;
}


void doAntColony(deque<TspNode> tspDeq, int numNodes){
    
    //Here are our parameters that need to be set.
    //Future statistical analysis could determine how exactly they are interacting with each other
    
    int numIterations = 10;
    double P = .3; //pheremone Exponent parameter
    double D = .5; //distance Exponent parameter
    int ratioOfAnts = 20; //ratio of ants is 20 to 1
    
    
    //Create Weight Matrix, which is a representation of all the edges of a directed graph. Edges are stored as to represent a directed graph, in that lookup is not symmetrical.
    
    InnerVector iv1(numNodes);
    vector<InnerVector> weightMatrix(numNodes,iv1);
    
    //initialize variables to keep track of the optimal tour
    
    int tourWithMinLength = 0;
    //double lengthOfMinTour = 0;
    deque<TspNode> tspMinTour;
    
    
    
    for (int i=0; i < numIterations; ++i){ //iterate an arbitrary amount of times
        
        
        
        int numAnts = numNodes * ratioOfAnts;
        
        vector<int> antPositions(numAnts);
        
        //here, store the tours for all the ants (there are more ants than there are edges)
        //ants always start and end on the same edges
        //A possible implementation could be making it so that ants can start and end on different edges.
        
        vector<Antdata> antTours(numAnts);
        
        //We assign the ants' initial vertices relatively equally
        
        int vertex = 0;
        for (int i =0; i<numAnts; ++i){
            
            antTours[i].vertex = vertex;
            
            if (vertex >= numNodes-1){
                vertex = 0;
            }
            ++vertex;
        }
        
        // Iterate through the number of ants, and for each ant j, store its tour
        
        for (int j=0; j < numAnts; ++j){
            deque<TspNode> antJTour;
            
            //start at ant j's vertex
            
            int startingVertex = antTours[j].vertex;
            TspNode tourNode1 = tspDeq[startingVertex];
            //int tN1 = startingVertex;
            
            antJTour.push_back(tourNode1);
            
            //We Start individual Ant's Tour by making a copy of the deque so we can delete nodes.
            //We delete tourNode1 so ant j can start finding its edges
            
            deque<TspNode>cityDeque = tspDeq;
            auto i = cityDeque.begin() + startingVertex;
            cityDeque.erase(i);
            
            // create a tour for ant j, adding nodes until the cityDeque is empty and ant j's tour is full
            
            while (!cityDeque.empty()) {
                double bestProbability = 0;
                
                //calculate the denominator for probability, which is all the avaiable nodes (the nodes that haven't been added to the tour already
                
                double allowed = 0;
                for (int i = 0; i < cityDeque.size(); ++i){
                    TspNode node = cityDeque[i];
                    double distance = getTSPDist(node, antJTour.back());
                    int index1 = node.index;
                    int index2 = antJTour.back().index;
                    double addedWeight = weightMatrix[index1].inner[index2];
                    allowed += (pow(addedWeight,P) * pow(1/distance, D));
                    
                }//end for (int i = 0; i < cityDeque.size(); ++i)
                
                //choose the next node that is the best to add to ant J's tour
                
                auto nodePointer = cityDeque.cbegin();
                
                for(auto it = cityDeque.cbegin(); it != cityDeque.cend(); ++it){
                    
                    TspNode node = *it;
                    double distance = getTSPDist(node, antJTour.back());
                    int index1 = node.index;
                    int index2 = antJTour.back().index;
                    double newWeight = weightMatrix[index1].inner[index2];
                    double pEdge = (pow(newWeight, P) * pow(1/distance, D));
                    double probability = pEdge/allowed;//write a way to do probability (needs a denominator)
                    
                    if( probability > bestProbability ){
                        bestProbability = probability;
                        nodePointer = it;
                    }
                }//end for(auto it = cityDeque.cbegin(); it != cityDeque.cend(); ++it)
                
                //add the chosen node to the tour
                
                TspNode newNode = *nodePointer;
                antJTour.push_back(newNode);
                cityDeque.erase(nodePointer);
                
            }//end while (!cityDeque.empty())
            
            //save antJ's tour
            
            Antdata antJData;
            
            
            antJData.tour = conversionHelper(antJTour);
            //antJData.tourLength = getTotalTSPDist(antJTour);
            
            antTours[j] = antJData;
            
            
            
            //antTours[j].tour = antJTour;
            
            
        }//end for (int j=0; j < numAnts; ++j)
        
        //here, we update the edges based on the tours taken by every ant
        //more ants should go to places where the edges are pheremone-laden and this should make those edges' weights higher.
        
        tourWithMinLength = updateWeights(weightMatrix, antTours, tspDeq, numNodes, numAnts);
        
        tspMinTour = conversionInttoDeq(antTours[tourWithMinLength].tour,tspDeq);
        
        //lengthOfMinTour = getTotalTSPDist(tspMinTour);
        
        //lengthOfMinTour = getTotalTSPDist (antTours[tourWithMinLength].tour);
        
        
        //tspMinTour = tourWithMinLength;
        //cout << getTotalTSPDist(tspMinTour) << endl;
        
    }//end for (int i=0; i < numIterations; ++i)
    
    printTSPSol(tspMinTour, numNodes);
    
}

int main(int argc, char ** argv) {
    
    _(ifstream arq(getenv("MYARQ")); cin.rdbuf(arq.rdbuf());)
    
    cin.sync_with_stdio(false);
    
    struct option longOpts[] = {
        {"mode", required_argument, NULL, 'm'},
        {"help", no_argument, NULL, 'h'},
    };
    
    //Project Data Structures and Vars
    //------------------------------------------------------------
    Mode mode = MST;
    vector<MstNode> mstVector;
    deque<TspNode> tspVector;
    int numNodes;
    
    //------------------------------------------------------------
    opterr = false;
    
    int opt = 0, index = 0;
    string arg;
    //Main Processing Loop
    while ((opt = getopt_long(argc, argv,
                              "m:h", longOpts, &index)) != -1) {
        switch (opt) {
            case 'm':{
                arg = optarg;
                if (arg == "MST"){
                    mode = MST;
                    
                }
                else if (arg == "OPTTSP"){
                    mode = OPTTSP;
                    
                }
                else if (arg == "FASTTSP"){
                    mode = FASTTSP;
                }
                else if (arg == "ANTCOL"){
                    mode = ANTCOL;
                }
                else {
                    
                    cerr<< "invalid mode\n";
                    exit(1);
                    
                }
            }
                break;
                
            case 'h':{
                cout<< "Please refer to instructions.pdf\n";
                exit(0);
            }
                break;
                
            default:{
                cerr<< "I didn't recognize one of your flags";
                exit(1);
            }
        }
    }
    
    if (mode == MST){
        
        cin >> numNodes;
        
        MstNode tempN;
        
        bool isValidMap = false;
        
        while (cin>>tempN.x){
            
            cin>>tempN.y;
            
            if (tempN.x < 0 && tempN.y < 0){
                tempN.terrain = 'S';
                
            }
            
            else if ((tempN.x == 0 && tempN.y <=0) || (tempN.y == 0 && tempN.x <= 0)){
                
                tempN.terrain = 'B';
                
                isValidMap = true;
                
            }
            
            else {
                
                tempN.terrain = 'L';
                
            }
            
            mstVector.push_back(tempN);
            
        }
        
        if (!isValidMap){
            
            cerr<< "Cannot construct MST\n";
            exit(1);
            
        }
        else {
            
            doMST(mstVector, numNodes);
            
        }
    }//endif
    
    else {
        cin >> numNodes;
        
        tspVector.resize(numNodes);
        
        int i = 0;
        TspNode tempN;
        while (cin>>tempN.x){
            
            cin>>tempN.y;
            
            tempN.index = i;
            tspVector[i] = tempN;
            
            
            ++i;
            
            
        }
        
        if (mode == FASTTSP){
            
            doFASTTSP(tspVector, numNodes);
            
        }
        
        else if (mode == OPTTSP){
            
            doOPTTSP(tspVector, numNodes);
            
        }
        
        else if (mode == ANTCOL){
            
            //doOPTTSP(tspVector, numNodes);
            doAntColony(tspVector, numNodes);
            
        }
        
        
        
    }
    
    return 0;
}

