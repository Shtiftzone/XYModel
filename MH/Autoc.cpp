#include "ranlxs.h"
#include "MultiHistRW.h"
#include "Autocorr.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

int main(int argc, char** argv) {
  rlxs_init(0,872348673);

  ifstream input(argv[1]);
  const int Ls = atoi(argv[2]);
  const int DIM=atoi(argv[3]);

  ofstream outf;
  
  int prova;
  double vol = pow(Ls,DIM);
  //Fill BETA, Action, O
  printf("#Files used:\n");

  string line;
  int nrun=0;
  while (getline(input, line)) nrun++;
  //cout << nrun << endl;
  input.clear();
  input.seekg (0, input.beg);
  
  vector<double> Action[nrun], O[nrun];
  double BETA[nrun];
  string data;
  char fname[256];
  cout << "BOH" << endl;
  for (int i=0; i<nrun; i++){
    getline(input, data);
    //cout << data << endl;
    
    sscanf(data.c_str(),"%lf %d\n",(BETA+i),&prova);
    if( DIM > 3 && Ls>72) sprintf(fname, "data/eul%dd/conf-%1.6f-%d.dat", DIM, BETA[i],prova);
    else if( DIM > 3 && Ls>48) sprintf(fname, "data/eul%dd/conf-%1.5f-%d.dat", DIM, BETA[i],prova);
    else sprintf(fname, "data/eul%dd/conf-%1.4f-%d.dat", DIM, BETA[i],prova);
    cout << fname << endl;
    ifstream in(fname);
    if(!in) {
      cerr<<"Cannot open file: "<<fname<<endl;
      exit(1);
    }
    double pE, pM;
    double Chip,Chin;
    while(in>>pM>>pE>>Chip>>Chin){
      //cout << Ls << " " << vol << " " << -pE << " " << abs(pM)/vol << endl;
      //cout << -pE << " " << abs(pM)/vol << " " << abs(Chip-Chin)/vol<< endl;
      //cout << -pE << " " << vol << " " << abs(pM)/vol << endl;
      
      Action[i].push_back(-pE);
      //O[i].push_back(abs(Chip-Chin)/vol);
      O[i].push_back(abs(pM)/vol);
    }

  }
  
  outf.open(argv[4],ios::out);

  for (int i=0; i<nrun; i++) {

     double autocorr, autocorr_err;
     AutoCorr(O[i], autocorr, autocorr_err);
     outf << DIM << " " 
	  << Ls << " " 
	  << BETA[i] << " " 
	  << autocorr << endl;
  }
  outf.close();

  return 0;

}
