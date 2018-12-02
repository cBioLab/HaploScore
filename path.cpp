//パスの確率計算
#include <gbwt/dynamic_gbwt.h>

using namespace std;
using namespace gbwt;

typedef struct{
  double pi_c,pi_r;
}Param;

vector<int> subpath_count; //subpathの本数を記録
vector<int> memoIndex; //メモ化されたどの位置を参照するか記録
vector<vector<int>> memo_B; //Bのメモ化した値を記録 (|A_b| * R)
int haplo_len;
Param params;

double calcHaploScore();
int A(int ep,int r);
int B(int n,int ep,int r);

void setParams(){
  params.pi_c = 0.9;
  params.pi_r = 0.01; //todo : |H|から計算するように変更
}

//todo : ファイル読み込み
void setPaths(vector_type& haplo,DynamicGBWT& dynamic_index){
  vector_type h{2,3,4};
  vector_type h1{1,2,3,0};
  vector_type h2{2,3,4,0};
  vector_type h3{3,4,5,0};
  dynamic_index.insert(h1);
  dynamic_index.insert(h2);
  dynamic_index.insert(h3);
  haplo = h;
  haplo_len = haplo.size();
  subpath_count.reserve(haplo_len*(haplo_len+1)/2);
  memoIndex.reserve(haplo_len*(haplo_len+1)/2);
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
  for(int i=0;i<count;i++){
    vector<int> v(haplo_len,-1);
    memo_B.push_back(v);
  }
}

int convertIndex(int i,int j){
  int n = haplo_len;
  return i*n - i*(i-1)/2 + (j-i);
}

//for debug subpath_count
void printMatrix(){
  for(int i=0;i<haplo_len;i++){
    for(int j=i;j<haplo_len;j++){
      cout << subpath_count[convertIndex(i,j)] << " ";
    }
    cout << endl;
  }
  for(int i=0;i<haplo_len;i++){
    for(int j=i;j<haplo_len;j++){
      cout << memoIndex[convertIndex(i,j)] << " ";
    }
    cout << endl;
  }
}

int main(){
  Verbosity::set(Verbosity::SILENT);

  DynamicGBWT dynamic_index;
  vector_type haplo;

  setParams();
  setPaths(haplo,dynamic_index);

  countPaths(haplo,dynamic_index);

  printMatrix();

  double haplo_score = calcHaploScore();

  cout << haplo_score << endl;

  return 0;
}

double calcHaploScore(){
  double score = 0.0;
  double mosaic_score = pow(params.pi_c,haplo_len-1);
  double rho = params.pi_r / params.pi_c;
  for(int i=0;i<haplo_len;i++){
    int mosaic_count = A(haplo_len-1,i);
    cout << mosaic_count << endl;
    score += mosaic_count * mosaic_score;
    mosaic_score *= rho;
  }
  return score;
}

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
