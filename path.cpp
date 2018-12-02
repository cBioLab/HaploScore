//パスの確率計算
#include <gbwt/dynamic_gbwt.h>

using namespace std;
using namespace gbwt;

typedef struct{
  double pi_c,pi_r;
}Param;

double calcHaploScore(vector<int>& count_memo,int haplo_len,Param& param);
int A(vector<int>& count_memo,int haplo_len,int ep,int r);
int B(vector<int>& count_memo,int haplo_len,int n,int ep,int r);

void setParams(Param& param){
  param.pi_c = 0.9;
  param.pi_r = 0.01; //todo : |H|から計算するように変更
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
}

void countPaths(vector<int>& count_memo,vector_type& haplo,DynamicGBWT& dynamic_index){
  int haplo_len = haplo.size();
  for(int i=0;i<haplo_len;i++){
    SearchState s = find(dynamic_index,haplo[i]);
    count_memo.push_back(s.size());
    for(int j=i+1;j<haplo_len;j++){
      if(s.size() == 0){
        count_memo.push_back(0);
        continue;
      }
      s = extend(dynamic_index,s,haplo[j]);
      count_memo.push_back(s.size());
    }
  }
}

int convertIndex(int i,int j,int n){
  return i*n - i*(i-1)/2 + (j-i);
}

//for debug count_memo
void printMatrix(vector<int>& matrix, int n){
  for(int i=0;i<n;i++){
    for(int j=i;j<n;j++){
      cout << matrix[convertIndex(i,j,n)] << " ";
    }
    cout << endl;
  }
}

int main(){
  Verbosity::set(Verbosity::SILENT);

  DynamicGBWT dynamic_index;
  vector_type haplo;
  vector<int> count_memo;
  Param params;

  setParams(params);
  setPaths(haplo,dynamic_index);

  int haplo_len = haplo.size();
  count_memo.reserve(haplo_len*(haplo_len+1)/2);
  countPaths(count_memo,haplo,dynamic_index);

  printMatrix(count_memo,haplo_len);

  double haplo_score = calcHaploScore(count_memo,haplo_len,params);

  cout << haplo_score << endl;

  return 0;
}

double calcHaploScore(vector<int>& count_memo,int haplo_len,Param& param){
  double score = 0.0;
  double mosaic_score = pow(param.pi_c,haplo_len-1);
  double rho = param.pi_r / param.pi_c;
  for(int i=0;i<haplo_len;i++){
    int mosaic_count = A(count_memo,haplo_len,haplo_len-1,i);
    cout << mosaic_count << endl;
    score += mosaic_count * mosaic_score;
    mosaic_score *= rho;
  }
  return score;
}

int A(vector<int>& count_memo,int haplo_len,int ep,int r){
  if(r > ep) return 0;
  if(r == 0) return count_memo[convertIndex(0,ep,haplo_len)];
  int count = 0;
  for(int i=0;i<=ep;i++){
    if(i < r) continue;
    count += B(count_memo,haplo_len,i,ep,r);
  }
  return count;
}

int B(vector<int>& count_memo,int haplo_len,int n,int ep,int r){
  if(r == 0 && n == 0) return A(count_memo,haplo_len,ep,r);
  if(r == 0 || n == 0) return 0;
  int count = count_memo[convertIndex(n,ep,haplo_len)] * A(count_memo,haplo_len,n-1,r-1);
  for(int i=0;i<n;i++){
    count -= B(count_memo,haplo_len,i,ep,r-1);
  }
  return count;
}
