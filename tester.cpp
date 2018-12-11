#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

#define dbg(x) cerr << #x << "=" << x << endl

int main(){

  const int numofTest = 1000;

  for(double p = 0.1; p < 1; p += 0.1){
    system("rm -rf testdata");
    system("rm -rf results");
    system("mkdir testdata");
    system("mkdir results");
    ofstream outputfile("result" + to_string(p) + ".csv");
    for(int i = 1; i <= numofTest; i++){
      string haploSetfile1 = "./testdata/haploSet" + to_string(i) + ".csv";
      string haploSetfile2 = "./testdata/haploSet_cut" + to_string(i) + ".csv";
      string haploScorefile1 = "./results/haploScore" + to_string(i) + ".csv";
      string haploScorefile2 = "./results/haploScore_cut" + to_string(i) + ".csv";
      string haplofile = "./testdata/haplo" + to_string(i) + ".csv";

      system(("./generator " + to_string(i)).c_str());
      system(("./path -H " + haploSetfile1 + " -h " + haplofile + " -s " + haploScorefile1 + " -p " + to_string(p)).c_str());
      system(("./path -H " + haploSetfile2 + " -h " + haplofile + " -s " + haploScorefile2 + " -p " + to_string(p)).c_str());
      dbg(to_string(p));
      ifstream inputfile1(haploScorefile1), inputfile2(haploScorefile2);
      double haploScore1, haploScore2;
      inputfile1 >> haploScore1; inputfile1.close();
      inputfile2 >> haploScore2; inputfile2.close();
      double arglog = log(haploScore2 / haploScore1);
      outputfile << arglog;
      if(i < numofTest) outputfile << ",";
    }
    outputfile << endl;
    outputfile.close();
  }

  return 0;
}
