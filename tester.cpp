#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <fstream>
using namespace std;

#define dbg(x) cerr << #x << "=" << x << endl

int main(){

  system("mkdir testdata");
  const int numofTest = 10;
  for(int i = 1; i <= numofTest; i++){
    string haploSetfile1 = "haploSet" + to_string(i) + ".csv";
    string haploSetfile2 = "haploSet_cut" + to_string(i) + ".csv";
    string haplofile = "haplo" + to_string(i) + ".csv";

    system(("./generator " + to_string(i)).c_str());
    system(("./path -H " + haploSetfile1 + " -h " + haplofile + " >> result.csv").c_str());
    system("echo , >> result.csv");
    system(("./path -H " + haploSetfile2 + " -h " + haplofile + " >> result.csv").c_str());
    system("echo \n >> result.csv");
  }

  return 0;
}
