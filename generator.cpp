#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include<fstream>
using namespace std;

#define dbg(x) cerr << #x << "=" << x << endl

struct Edge{
  int nxtNode;
  int probability;
};
typedef vector<Edge> Edges;
typedef vector<Edges> Graph;
typedef vector<int> Path;

Graph createDAG(int numofNode, int maxgroupSize, int maxpathSize){
  random_device rnd;
  mt19937 mt(rnd());
  uniform_int_distribution<> rand100(1, 99);
  uniform_int_distribution<> randgroupSize(2, maxgroupSize);
  uniform_int_distribution<> randpathSize(1, maxpathSize);

  // make DAG topology
  Graph resultGraph(numofNode, Edges());
  for(int crtNode = 0; crtNode < numofNode; ){
    int crtgroupSize = randgroupSize(mt);
    crtgroupSize = min(crtgroupSize, numofNode - crtNode);
    int startNode = crtNode;
    int endNode = startNode + crtgroupSize - 1;

    for(int midNode = startNode + 1; midNode < endNode; ){
      int crtpathSize = randpathSize(mt);
      crtpathSize = min(crtpathSize, endNode - midNode);

      for(int i = midNode; i < midNode + crtpathSize; i++){
        if(i == midNode)resultGraph[startNode].push_back((Edge){i, -1});
        if(i == midNode + crtpathSize - 1)resultGraph[i].push_back((Edge){endNode, -1});
        if(i < midNode + crtpathSize - 1)resultGraph[i].push_back((Edge){i+1, -1});
      }
      
      midNode += crtpathSize;
    }

    if(startNode + 1 == endNode)resultGraph[startNode].push_back((Edge){endNode, -1});

    if(endNode + 1 < numofNode){
      resultGraph[endNode].push_back((Edge){endNode + 1, -1});
    }

    crtNode += crtgroupSize;
  }

  // decide probability
  for(int crtNode = 0; crtNode < numofNode; crtNode++){
    int branchSize = resultGraph[crtNode].size();
    vector<int> probnums;
    probnums.push_back(0);
    probnums.push_back(100);
    for(int i = 0; i < branchSize - 1; i++) probnums.push_back(rand100(mt));
    sort(probnums.begin(), probnums.end());
    for(int i = 0; i < branchSize; i++){
      Edge &e = resultGraph[crtNode][i];
      e.probability = (probnums[i+1] - probnums[i]);
    }
  }

  return resultGraph;
}

void graphOutput(Graph g){
  int numofNode = g.size();

  ofstream outputfile1("topology.txt");
  outputfile1 << "Number of nodes = " << numofNode << endl;
  for(int crtNode = 0; crtNode < numofNode; crtNode++){
    for(Edge e : g[crtNode]){
      outputfile1 << crtNode + 1 << " " << e.nxtNode + 1 << endl; // 1-indexed
    }
  }
  outputfile1.close();

  ofstream outputfile2("topology_prob.txt");
  outputfile2 << "Number of nodes = " << numofNode << endl;
  for(int crtNode = 0; crtNode < numofNode; crtNode++){
    for(Edge e : g[crtNode]){
      outputfile2 << crtNode + 1 << " " << e.nxtNode + 1 << " : prob = " << e.probability << endl; // 1-indexed
    }
  }
  outputfile2.close();
}

void subpathSearch(const Graph& givenGraph, int crtNode, vector<Path>& resultPaths, vector<double>& subpathProbs,
  vector<int> crtNodes, double crtProb, int pathlenLb, int pathlenUb){
  if((int) crtNodes.size() >= pathlenLb && (int) crtNodes.size() <= pathlenUb){
    resultPaths.push_back(crtNodes);
    subpathProbs.push_back(crtProb);
  }
  for(Edge e : givenGraph[crtNode]){
    crtNodes.push_back(e.nxtNode);
    subpathSearch(givenGraph, e.nxtNode, resultPaths, subpathProbs, crtNodes, crtProb * ((double)e.probability / 100.0), pathlenLb, pathlenUb);
    crtNodes.pop_back();
  }
}

void extractSubpath(const Graph& givenGraph, int numofNode, int pathlenLb, int pathlenUb){
  vector<Path> resultPaths;
  vector<double> subpathProbs;

  for(int crtNode = 0; crtNode < numofNode; crtNode++){
    subpathSearch(givenGraph, crtNode, resultPaths, subpathProbs, vector<int>(1,crtNode), 1.0, pathlenLb, pathlenUb);
  }

  vector<int> selectedPathidx;
  random_device rnd;
  mt19937 mt(rnd());
  uniform_int_distribution<> rand100(0, 99);
  for(int i = 0; i < (int)resultPaths.size(); i++){
    double p = (double)rand100(mt) / 100.0;
    dbg(p); dbg(subpathProbs[i]);
    if(p < subpathProbs[i]){
      selectedPathidx.push_back(i);
    }
  }

  dbg(selectedPathidx.size()); dbg(resultPaths.size());

  ofstream outputfile1("haploSet.csv");
  for(int i : selectedPathidx){
    Path p = resultPaths[i];
    double prob = subpathProbs[i];
    for(int node : p){
      outputfile1 << node + 1 << ","; // 1-indexed
    }
    outputfile1 << "0" <<endl; // ENDMARKER
  }
  outputfile1.close();

  ofstream outputfile2("haploSet_prob.txt");
  for(int i : selectedPathidx){
    Path p = resultPaths[i];
    double prob = subpathProbs[i];
    for(int node : p){
      outputfile2 << node + 1 << ","; // 1-indexed
    }
    outputfile2 << prob;
    outputfile2 << endl;
  }
  outputfile2.close();
}

int main(){

  Graph givenGraph = createDAG(25, 5, 2);
  graphOutput(givenGraph);
  extractSubpath(givenGraph,givenGraph.size(),1,100);

  return 0;
}
