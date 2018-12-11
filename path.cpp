//パスの確率計算
#include <gbwt/dynamic_gbwt.h>

#define DEBUG

using namespace std;
using namespace gbwt;

typedef struct{
  double pi_c,pi_r;
}Param;

vector<int> subpath_count; //subpathの本数を記録
vector<int> memoIndex; //メモ化されたどの位置を参照するか記録
vector<vector<int>> memo_B; //Bのメモ化した値を記録 (|A_b| * R)
int haplo_len;
int n_of_haplos;
Param params;

double calcHaploScore();
int A(int ep,int r);
int B(int n,int ep,int r);

void setParams(double p_of_chain){
  params.pi_c = p_of_chain;
  params.pi_r = (1 - params.pi_c) / n_of_haplos;
}

void setHaplos(string& haploFile,vector<vector_type>& haplos){
    ifstream ifs(haploFile);
    string str;
    if (!ifs) exit(-1);
    while (getline(ifs, str)){
        string token;
        istringstream stream(str);
        vector_type v;
        while (getline(stream, token, ',')){
            v.push_back(stoi(token.c_str()));
        }
        if(v.size()) haplos.push_back(v);
    }
}

void setPaths(string& haploFile,string& haploSetFile,vector_type& haplo,DynamicGBWT& dynamic_index){
  vector<vector_type> h;
  vector<vector_type> haploSet;
  setHaplos(haploFile,h);
  setHaplos(haploSetFile,haploSet);
  haplo = h[0];
  n_of_haplos = haploSet.size();
  for(int i=0;i<n_of_haplos;i++){
    dynamic_index.insert(haploSet[i]);
  }
  haplo_len = haplo.size();
  subpath_count.reserve(haplo_len*(haplo_len+1)/2);
  memoIndex.reserve(haplo_len*(haplo_len+1)/2);

  #ifdef DEBUG
  cerr << "haplotype sequence" << endl;
  for(int i=0;i<haplo_len;i++) cerr << haplo[i] << " ";
  cerr << endl;
  cerr << "haplo set" << endl;
  for(auto haplo_itr = haploSet.begin();haplo_itr != haploSet.end();haplo_itr++){
    for(auto itr = (*haplo_itr).begin();itr != (*haplo_itr).end();itr++){
      cerr << *itr << " ";
    }
    cerr << endl;
  }
  #endif
}

void countPaths(vector_type& haplo,DynamicGBWT& dynamic_index){
  int count = 0;
  for(int i=0;i<haplo_len;i++){
    SearchState s = find(dynamic_index,haplo[i]);
    subpath_count.push_back(s.size());
    memoIndex.push_back(count);
    for(int j=i+1;j<haplo_len;j++){
      if(s.size() == 0){
        subpath_count.push_back(0);
        memoIndex.push_back(count);
        continue;
      }
      s = extend(dynamic_index,s,haplo[j]);
      int range = s.size();
      if(subpath_count.back() != range) count++;
      subpath_count.push_back(range);
      memoIndex.push_back(count);
    }
    count++;
  }
  memo_B.reserve(count);
  vector<int> v(haplo_len,-1);
  for(int i=0;i<count;i++){
    memo_B.push_back(v);
  }
}

int convertIndex(int i,int j){
  int n = haplo_len;
  return i*n - i*(i-1)/2 + (j-i);
}

//for debug subpath_count
void printMatrix(){
  #ifdef DEBUG
  cerr << "subpath_count" << endl;
  for(int i=0;i<haplo_len;i++){
    for(int j=i;j<haplo_len;j++){
      cerr << subpath_count[convertIndex(i,j)] << " ";
    }
    cerr << endl;
  }
  cerr << "memoIndex" << endl;
  for(int i=0;i<haplo_len;i++){
    for(int j=i;j<haplo_len;j++){
      cerr << memoIndex[convertIndex(i,j)] << " ";
    }
    cerr << endl;
  }
  #endif
}

int main(int argc,char* argv[]){
  string haploFile = "haplo.csv";
  string haploSetFile = "haploSet.csv";
  string haploScoreFile = "haploScore.csv";
  double p_of_chain = 0.7;

  int opt;
  while((opt = getopt(argc, argv, "h:p:H:s:")) != -1){
    switch(opt){
    case 'h':
      {
    	  haploFile = optarg;
	      break;
      }
    case 'H':
      {
	      haploSetFile = optarg;
	      break;
      }
    case 'p':
      {
	      p_of_chain = stod(optarg);
	      break;
      }
    case 's':
      {
        haploScoreFile = optarg;
        break;
      }
    default:
      {
	      fprintf(stderr,"error! \'%c\' \'%c\'\n", opt, optopt);
	      return -1;
	      break;
      }
    }
  }

  Verbosity::set(Verbosity::SILENT);

  DynamicGBWT dynamic_index;
  vector_type haplo;

  setPaths(haploFile,haploSetFile,haplo,dynamic_index);
  setParams(p_of_chain); 

  countPaths(haplo,dynamic_index);

  printMatrix();

  ofstream outputfile(haploScoreFile);
  double haplo_score = calcHaploScore();
  outputfile << haplo_score << endl;
  outputfile.close();

  return 0;
}

double calcHaploScore(){
  double score = 0.0;
  double mosaic_score = pow(params.pi_c,haplo_len-1);
  double rho = params.pi_r / params.pi_c;
  for(int i=0;i<haplo_len;i++){
    int mosaic_count = A(haplo_len-1,i);
    score += mosaic_count * mosaic_score;
    mosaic_score *= rho;
  }
  return score;
}

//h[0,ep],recombination数rのモザイクのパターン数
int A(int ep,int r){
  if(r > ep) return 0;
  if(r == 0) return subpath_count[convertIndex(0,ep)];
  int count = 0;
  for(int i=0;i<=ep;i++){
    if(i < r) continue;
    count += B(i,ep,r);
  }
  return count;
}

//h[0,ep],recombination数rかつ，最右のセグメントがh[n,ep]であるモザイクのパターン数
int B(int n,int ep,int r){
  if(memo_B[memoIndex[convertIndex(n,ep)]][r] != -1) return memo_B[memoIndex[convertIndex(n,ep)]][r];
  if(r == 0 && n == 0) return memo_B[memoIndex[convertIndex(n,ep)]][r] = A(ep,r);
  if(r == 0 || n == 0) return memo_B[memoIndex[convertIndex(n,ep)]][r] = 0;
  int count = subpath_count[convertIndex(n,ep)] * A(n-1,r-1);
  for(int i=0;i<n;i++){
    count -= B(i,ep,r-1);
  }
  return memo_B[memoIndex[convertIndex(n,ep)]][r] = count;
}
