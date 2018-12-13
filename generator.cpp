#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <fstream>
#include <cassert>
using namespace std;

#define dbg(x) cerr << #x << "=" << x << endl

struct Edge{
  int nxtNode;
  int probability;
};
typedef vector<Edge> Edges;
typedef vector<Edges> Graph;
typedef vector<int> Haplo;
typedef vector<Haplo> HaploSet;

string getfileName(string prf, string suf, int fileidx){
  return "./testdata/" + prf + to_string(fileidx) + suf;
}

void outputHaplo(Haplo h, string filename){
  ofstream ofs(filename);
  for(int node : h){
    ofs << node + 1 << ","; // 1-indexed
  }
  ofs << "0" << endl; // END-MARKER
  ofs.close();
}

void outputHaploSet(HaploSet hs, string filename){
  ofstream ofs(filename);
  for(Haplo h : hs){
    for(int node : h){
      ofs << node + 1 << ","; // 1-indexed
    }
    ofs << "0" << endl; // END-MARKER
  }
  ofs.close();
}

void outputGraph(Graph g,int fileidx){
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

void subhaploEnum(const Graph& givenGraph, int crtNode, HaploSet& resultHaplos,
 vector<double>& subhaploProbs, Haplo crtHaplo, double crtProb){

  resultHaplos.push_back(crtHaplo);
  subhaploProbs.push_back(crtProb);
  for(Edge e : givenGraph[crtNode]){
    crtHaplo.push_back(e.nxtNode);
    subhaploEnum(givenGraph, e.nxtNode, resultHaplos, subhaploProbs, crtHaplo, crtProb * ((double)e.probability / 100.0));
    crtHaplo.pop_back();
  }
}

void extractSubhaplo(const Graph& givenGraph, int numofNode, int fileidx){
  random_device rnd;
  mt19937 mt(rnd());
  uniform_int_distribution<> rand100(0, 99);

  HaploSet allsubHaplos;
  vector<double> subhaploProbs;
  for(int crtNode = 0; crtNode < numofNode; crtNode++){
    subhaploEnum(givenGraph, crtNode, allsubHaplos, subhaploProbs, Haplo(1,crtNode), 1.0);
  }

  HaploSet selectedHaplos;
  vector<double> selectedHaploprobs;
  for(int i = 0; i < (int)allsubHaplos.size(); i++){
    double p = (double)rand100(mt) / 100.0;
    if(p < subhaploProbs[i]){
      selectedHaplos.push_back(allsubHaplos[i]);
      selectedHaploprobs.push_back(subhaploProbs[i]);
    }
  }

  // output haploSet.csv
  outputHaploSet(selectedHaplos, getfileName("haploSet", ".csv", fileidx));

  // FOR DEBUG OUTPUT
  {
    ofstream outputfile2(getfileName("haploSet_prob", ".txt", fileidx));
    for(int i = 0;i < (int)selectedHaplos.size(); i++){
      Haplo p = selectedHaplos[i];
      double prob = selectedHaploprobs[i];
      for(int node : p){
        outputfile2 << node + 1 << ","; // 1-indexed
      }
      outputfile2 << prob;
      outputfile2 << endl;
    }
    outputfile2.close();
  }

  // cutting haplo
  HaploSet cutHaplos;
  double cuttingProb = 0.02; // Cutting probability = 2%
  for(int i = 0; i < (int)selectedHaplos.size(); i++){
    Haplo h;
    for(int node : selectedHaplos[i]){
      double p = (double)rand100(mt) / 100.0;
      if((int)h.size() > 0 && p < cuttingProb){
        cutHaplos.push_back(h);
        h.clear();
      }
      h.push_back(node);
    }
    if((int)h.size() > 0) cutHaplos.push_back(h);
  }

  outputHaploSet(cutHaplos, getfileName("haploSet_cut", ".csv", fileidx));
}

void searchHaplo(const Graph& givenGraph, int crtNode, Haplo crtHaplo, Haplo& foundHaplo, int haploSize){
  if((int)crtHaplo.size() >= haploSize){
    foundHaplo = crtHaplo;
    return ;
  }
  for(Edge e : givenGraph[crtNode]){
    crtHaplo.push_back(e.nxtNode);
    searchHaplo(givenGraph, e.nxtNode, crtHaplo, foundHaplo, haploSize);
    if(foundHaplo.size() > 0)return ;
    crtHaplo.pop_back();
  }
}

void createHaplo(const Graph& givenGraph, int fileidx, int haploSize){
  Haplo foundHaplo;
  searchHaplo(givenGraph, 0, Haplo(1,0), foundHaplo, haploSize);
  if((int)foundHaplo.size() < haploSize){
    fprintf(stderr, "Haplo Not Found.\n");
    assert(false);
  }

  outputHaplo(foundHaplo, getfileName("haplo", ".csv", fileidx));
}

int main(int argc, char *argv[]){
  if(argc < 2){
    cerr << "few arguments." << endl;
    return -1;
  }

  int fileidx = atoi(argv[1]);
  Graph givenGraph = createDAG(30, 5, 2);
  outputGraph(givenGraph, fileidx);
  extractSubhaplo(givenGraph, givenGraph.size(), fileidx);
  createHaplo(givenGraph, fileidx, 15);

  return 0;
}
