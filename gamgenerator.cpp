#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <fstream>
#include <cassert>
#include <map>
#include <set>
using namespace std;

#define dbg(x) cerr << #x << "=" << x << endl
//#define DEBUG

struct GAMelem{
  int node, basenum;
  string base;
};
bool operator==(const GAMelem& a, const GAMelem& b){
  if(a.node == b.node && a.base == b.base)return true;
  else return false;
}
typedef vector<GAMelem> Haplo;
typedef vector<Haplo> GAM;

vector<string> split(const string &s, char c1, char c2) {
  vector<string> elems;
  string item;
  int lev = 0;
  for (char ch: s) {
    if (ch == c1) {
      lev++;
    } else if(ch == c2) {
      elems.push_back(item);
      lev--;
      item.clear();
    } else if(lev > 0){
      item += ch;
    }
  }
  return elems;
}

string getfileName(string prf, string suf, int fileidx){
  return "./testdata" + prf + to_string(fileidx) + suf;
}

GAM inputGAM(){
  ifstream inputfile("input.txt");

  string pathinfo, haploinfo;
  Haplo inputhaplo;
  getline(inputfile, pathinfo);
  getline(inputfile, haploinfo);
  inputfile.close();
  map<int, int> node2idx;
  GAM inputgam;

  // get GAM data
  haploinfo = haploinfo.substr(1,haploinfo.size()-2);
  vector<string> gs = split(haploinfo, '(', ')');
  for(string str : gs){
    int node = 0, basenum;
    string base;
    int idx = 0;
    while(idx < (int)str.size() && str[idx] != ','){
      node *= 10;
      node += str[idx++] - '0';
    }
    idx += 3;
    while(idx < (int)str.size() && str[idx] != '\''){
      base += str[idx++];
    }
    basenum = base.size();
    inputhaplo.push_back((GAMelem){node, basenum, base});
    node2idx[node] = inputhaplo.size() - 1;
  }

  // get path info
  pathinfo = pathinfo.substr(1,pathinfo.size()-2);
  vector<string> hs = split(pathinfo, '[', ']');
  for(string str : hs){
    int idx = 0;
    Haplo path; 
    while(1){
      int node = 0;
      while(idx < (int)str.size() && isdigit(str[idx])){
        node *= 10;
        node += str[idx++] - '0';
      }
      path.push_back(inputhaplo[node2idx[node]]);
      if(idx == (int)str.size() )break;
      while(idx < (int)str.size() && !isdigit(str[idx])){
        idx++;
      }
    }
    inputgam.push_back(path);
  }

  return inputgam;
}

void outputGAM(const GAM& hs, string filename){
  ofstream ofs(filename);
  for(Haplo h : hs){
    ofs << h[0].node;
    for(int i = 1; i < (int)h.size(); i++){
      ofs << "," << h[i].node;
    }
    ofs << endl;
  }
  ofs << "0" << endl; // END-MARKER

  set<pair<int,string> > seqset;
  for(Haplo h : hs){
    for(GAMelem e : h){
      seqset.insert(pair<int,string>(e.node,e.base));
    }
  }
  for(pair<int,string> p : seqset){
    ofs << p.first << "," << p.second << endl;
  }
  ofs << "0" << endl;

  ofs.close();
}

Haplo createHaplo(const GAM& gam){
  int haploNum = gam.size();
  random_device rnd;
  mt19937 mt(rnd());
  uniform_int_distribution<> randh(0, haploNum - 1);
  uniform_int_distribution<> rand100(0, 99);
  int crtHaplo = randh(mt); // select random haplo
  int crtNode = gam[crtHaplo][0].node;
  double piC = 0.7; // current piC
  Haplo resultHaplo;

  for(int i = 0; ; ){
    resultHaplo.push_back((GAMelem){crtNode, gam[crtHaplo][i].basenum, gam[crtHaplo][i].base});
    vector<pair<int,int> > canditate;
    for(int j = 0; j < haploNum; j++){
      for(int k = 0; k < (int)gam[j].size(); k++){
        if(gam[j][k].node != crtNode) continue;
        canditate.push_back(pair<int,int>(j,k));
        break;
      }
    }
    double piR = (1.0 - piC) / (double) canditate.size();
    for(int cnt = 0; cnt < gam[crtHaplo][i].basenum; cnt++){
      double p = (double) rand100(mt) / 100.0;
      if(p >= piC){
        int nxtHaploidx = (int)floor((p - piC) / piR);
        nxtHaploidx = min(nxtHaploidx, (int)canditate.size() - 1);
        crtHaplo = canditate[nxtHaploidx].first;
        i = canditate[nxtHaploidx].second;
      }
    }
    if(i + 1 < (int)gam[crtHaplo].size()) {
      crtNode = gam[crtHaplo][i+1].node; i++;
    } else break;
  }

  return resultHaplo;
}

GAM gamGenerate(const GAM& gam){
  GAM resultSet;
  while(1){
    Haplo newh = createHaplo(gam);
    bool ok = true;
    for(Haplo h : gam){
      if(newh == h){
        ok = false;
        break;
      }
    }
    if(ok){
      resultSet.push_back(newh);
      break;
    }
  }
  return resultSet;
}

int main(int argc, char *argv[]){
  GAM inputgam = inputGAM();
#ifdef DEBUG
  // for DEBUG
  for(Haplo h : inputgam){
    cerr << " ===== " << endl;
    int len = 0;
    for(GAMelem g : h){
      cerr << g.node << " " ;
      len += g.basenum;
    }
    cerr << endl;
    dbg(len);
  }
  cerr<<endl;
#endif

  GAM resultgam = gamGenerate(inputgam);
#ifdef DEBUG
  // for DEBUG
  for(Haplo h : resultgam){
    cerr << " ===== " << endl;
    int len = 0;
    for(GAMelem g : h){
      cerr << g.node << " ";
      len += g.basenum;
    }
    cerr << endl;
  }
#endif
  outputGAM(resultgam, "output.csv");

  return 0;
}
