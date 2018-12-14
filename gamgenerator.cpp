#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <fstream>
#include <cassert>
using namespace std;

#define dbg(x) cerr << #x << "=" << x << endl

struct GAMelem{
  int node;
  int basenum;

};
typedef vector<GAMelem> Haplo;
typedef vector<Haplo> GAM;

string getfileName(string prf, string suf, int fileidx){
  return "./testdata" + prf + to_string(fileidx) + suf;
}

void outputGAM(const &GAM hs, string filename){
  ofstream ofs(filename);
  for(Haplo h : hs){
    for(GAMelem g : h){
      ofs << g.node + 1 << ","; // 1-indexed
    }
    ofs << "0" << endl; // END-MARKER
  }
  ofs.close();
}

Haplo createHaplo(const GAM& gam){
  int haploNum = gam.size();
  int haploLen = gam[0].size();
  random_device rnd;
  mt19937 mt(rnd());
  uniform_int_distribution<> randh(0, haploNum - 1);
  uniform_int_distribution<> rand100(0, 99);
  int crtHaplo = randh(mt); // select random haplo
  int crtNode = gam[crtHaplo][0].node;
  double piC = 0.7; // current piC
  Haplo resultHaplo;

  for(int i = 0; i < haploLen; i++){
    //resultHaplo.push_back(crtNode);
    vector<int> canditate;
    for(int j = 0; j < haploNum; j++){
      if(j == crtHaplo || gam[j][i].node != crtNode) continue;
      canditate.push_back(j);
    }
    if(canditate.size() > 0) {
      double piR = (1.0 - piC) / (double) canditate.size();
      for(int cnt = 0; cnt < gam[crtHaplo][i].basenum; cnt++){
        double p = (double) rand100(mt) / 100.0;
        if(p >= piC){
          int nxtHaploidx = (int)floor((p - piC) / piR);
          nxtHaploidx = min(nxtHaploidx, (int)canditate.size() - 1);
          crtHaplo = canditate[nxtHaploidx];
        }
      }
    }
    if(i + 1 < haploLen)crtNode = gam[crtHaplo][i+1].node;
  }

  return resultHaplo;
}

GAM gamGenerate(const GAM& gam){
  GAM resultSet;
  for(int i = 0; i < (int) gam.size(); i++){
    resultSet.push_back(createHaplo(gam));
  }
  return resultSet;
}

GAM getInput(){
  GAM inputgam;

  string s;
  while(getline(cin,s)){
    int node, basenum;
    sscanf(s.c_str(), "%d,%d", &node, &basenum);
    //inputgam.push_back((Haplo){node, basenum});
  }
  return inputgam;
}

int main(int argc, char *argv[]){
  GAM inputgam = getInput();
  GAM resultgam = gamGenerate(inputgam);
  outputGAM(resultgam, getfileName());

  return 0;
}
