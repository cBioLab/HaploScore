#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

#define dbg(x) cerr << #x << "=" << x << endl

int main(){

  const int numofTest = 1000;

  system("rm *.csv");
  system("rm -rf testdata");
  system("mkdir testdata");

  for(int i = 1; i <= numofTest; i++){
    system(("./generator " + to_string(i)).c_str());
  }
  fprintf(stderr, "[+] Generated. \n");

  for(int p = 1; p < 10; p += 1){
    system("rm -rf results");
    system("mkdir results");
    ofstream outputfile1("result_cut" + to_string(p) + ".csv");
    ofstream outputfile2("result_shf" + to_string(p) + ".csv");
    for(int i = 1; i <= numofTest; i++){
      string haploSetfile1 = "./testdata/haploSet" + to_string(i) + ".csv";
      string haploSetfile2 = "./testdata/haploSet_cut" + to_string(i) + ".csv";
      string haploSetfile3 = "./testdata/haploSet_shf" + to_string(i) + ".csv";
      string haploScorefile1 = "./results/haploScore" + to_string(i) + ".csv";
      string haploScorefile2 = "./results/haploScore_cut" + to_string(i) + ".csv";
      string haploScorefile3 = "./results/haploScore_shf" + to_string(i) + ".csv";
      string haplofile = "./testdata/haplo" + to_string(i) + ".csv";


      system(("./path -H " + haploSetfile1 + " -h " + haplofile + " -s " + haploScorefile1 + " -p 0." + to_string(p)).c_str());
      system(("./path -H " + haploSetfile2 + " -h " + haplofile + " -s " + haploScorefile2 + " -p 0." + to_string(p)).c_str());
      system(("./path -H " + haploSetfile3 + " -h " + haplofile + " -s " + haploScorefile3 + " -p 0." + to_string(p)).c_str());
      ifstream inputfile1(haploScorefile1), inputfile2(haploScorefile2), inputfile3(haploScorefile3);
      double haploScore1, haploScore2, haploScore3;
      inputfile1 >> haploScore1; inputfile1.close();
      inputfile2 >> haploScore2; inputfile2.close();
      inputfile3 >> haploScore3; inputfile3.close();

      outputfile1 << haploScore2 - haploScore1;
      outputfile2 << haploScore3 - haploScore1;
      if(i < numofTest){
        outputfile1 << ",";
        outputfile2 << ",";
      }
    }
    outputfile1 << endl;
    outputfile1.close();
    outputfile2 << endl;
    outputfile2.close();

    fprintf(stderr, "[+] Ended %d process. \n", p);
  }
  system("python3 plot.py");
  fprintf(stderr, "[+] Ended plotting. \n");

  return 0;
}
