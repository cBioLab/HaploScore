#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <fstream>
using namespace std;

#define dbg(x) cerr << #x << "=" << x << endl

struct Edge{
  int nxtNode;
  int probability;
};
typedef vector<Edge> Edges;
typedef vector<Edges> Graph;
typedef vector<int> Path;

string getfileName(string prf, string suf, int fileidx){
  return "./testdata/" + prf + to_string(fileidx) + suf;
}

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

void graphOutput(Graph g,int fileidx){
  int numofNode = g.size();

  ofstream outputfile1(getfileName("topology", ".txt", fileidx));
  outputfile1 << "Number of nodes = " << numofNode << endl;
  for(int crtNode = 0; crtNode < numofNode; crtNode++){
    for(Edge e : g[crtNode]){
      outputfile1 << crtNode + 1 << " " << e.nxtNode + 1 << endl; // 1-indexed
    }
  }
  outputfile1.close();

  ofstream outputfile2(getfileName("topology_prob", ".txt", fileidx));
  outputfile2 << "Number of nodes = " << numofNode << endl;
  for(int crtNode = 0; crtNode < numofNode; crtNode++){
    for(Edge e : g[crtNode]){
      outputfile2 << crtNode + 1 << " " << e.nxtNode + 1 << " : prob = " << e.probability << endl; // 1-indexed
    }
  }
  outputfile2.close();
}

void subpathEnum(const Graph& givenGraph, int crtNode, vector<Path>& resultPaths,
 vector<double>& subpathProbs, Path crtNodes, double crtProb){

  resultPaths.push_back(crtNodes);
  subpathProbs.push_back(crtProb);
  for(Edge e : givenGraph[crtNode]){
    crtNodes.push_back(e.nxtNode);
    subpathEnum(givenGraph, e.nxtNode, resultPaths, subpathProbs, crtNodes, crtProb * ((double)e.probability / 100.0));
    crtNodes.pop_back();
  }
}

void extractSubpath(const Graph& givenGraph, int numofNode, int fileidx){
  random_device rnd;
  mt19937 mt(rnd());
  uniform_int_distribution<> rand100(0, 99);

  vector<Path> allsubPaths;
  vector<double> subpathProbs;
  for(int crtNode = 0; crtNode < numofNode; crtNode++){
    subpathEnum(givenGraph, crtNode, allsubPaths, subpathProbs, vector<int>(1,crtNode), 1.0);
  }

  vector<Path> selectedPaths;
  vector<double> selectedPathprobs;
  for(int i = 0; i < (int)allsubPaths.size(); i++){
    double p = (double)rand100(mt) / 100.0;
    if(p < subpathProbs[i]){
      selectedPaths.push_back(allsubPaths[i]);
      selectedPathprobs.push_back(subpathProbs[i]);
    }
  }

  // output haploSet.csv
  ofstream outputfile1(getfileName("haploSet", ".csv", fileidx));
  for(int i = 0;i < (int)selectedPaths.size(); i++){
    Path p = selectedPaths[i];
    for(int node : p){
      outputfile1 << node + 1 << ","; // 1-indexed
    }
    outputfile1 << "0" <<endl; // ENDMARKER
  }
  outputfile1.close();

  ofstream outputfile2(getfileName("haploSet_prob", ".txt", fileidx));
  for(int i = 0;i < (int)selectedPaths.size(); i++){
    Path p = selectedPaths[i];
    double prob = selectedPathprobs[i];
    for(int node : p){
      outputfile2 << node + 1 << ","; // 1-indexed
    }
    outputfile2 << prob;
    outputfile2 << endl;
  }
  outputfile2.close();

  // cutting path
  vector<Path> cutPaths;
  double cuttingProb = 0.02;
  for(int i = 0; i < (int)selectedPaths.size(); i++){
    Path path;
    for(int node : selectedPaths[i]){
      double p = (double)rand100(mt) / 100.0;
      if(path.size() > 0 && p < cuttingProb){
        cutPaths.push_back(path);
        path.clear();
      }
      path.push_back(node);
    }
    if(path.size() > 0) cutPaths.push_back(path);
  }

  // output haploSet_cut.csv
  ofstream outputfile3(getfileName("haploSet_cut", ".csv", fileidx));
  for(int i = 0;i < (int)cutPaths.size(); i++){
    Path p = cutPaths[i];
    for(int node : p){
      outputfile3 << node + 1 << ","; // 1-indexed
    }
    outputfile3 << "0" <<endl; // ENDMARKER
  }
  outputfile3.close();
}

void searchHaplo(const Graph& givenGraph, int crtNode, Path crtNodes, Path& foundHaplo, int haploSize){
  if((int)crtNodes.size() >= haploSize){
    foundHaplo = crtNodes;
    return ;
  }
  for(Edge e : givenGraph[crtNode]){
    crtNodes.push_back(e.nxtNode);
    searchHaplo(givenGraph, e.nxtNode, crtNodes, foundHaplo, haploSize);
    if(foundHaplo.size() > 0)return ;
    crtNodes.pop_back();
  }
}

void createHaplo(const Graph& givenGraph, int fileidx, int haploSize){
  Path foundHaplo;
  searchHaplo(givenGraph, 0, Path(1,0), foundHaplo, haploSize);

  // output Haplo
  ofstream outputfile(getfileName("haplo", ".csv", fileidx));
  for(int node : foundHaplo){
    outputfile << node + 1 << ","; // 1-indexed
  }
  outputfile << "0" <<endl; // ENDMARKER
  outputfile.close();
}

int main(int argc, char *argv[]){
  if(argc < 2){
    cerr << "few arguments." << endl;
    return -1;
  }

  int fileidx = atoi(argv[1]);
  Graph givenGraph = createDAG(25, 5, 2);
  graphOutput(givenGraph, fileidx);
  extractSubpath(givenGraph, givenGraph.size(), fileidx);
  createHaplo(givenGraph, fileidx, 12);

  return 0;
}
