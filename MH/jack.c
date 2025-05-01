#include "ranlxs.h"
#include "MultiHistRW.h"
//#include "MCData.h"
#include "Autocorr.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

// Setting default values
static int nrun = 4;
char file[255];
double BETA={1.37};
// Dimensioni  del reticolo
double Ls=8., Lt=8.;


void JackKnife(const vector<double> &dt, vector<double> &dt2, int blsize, double &avr, double &err){

 int size = dt.size();
 int nblocks = size/blsize; //Numero di blocchi

 double psum1[nblocks]; //Somme parziali del blocco
 double psum2[nblocks]; //Somme parziali del blocco
 double avp1[nblocks]; //Medie togliendo un blocco
 double avp2[nblocks]; //Medie togliendo un blocco
 double tsum1=0.0; //Somma totale
 double tsum2=0.0; //Somma totale

 for(int j=0; j<nblocks; j++) {
   psum1[j]=0.0;
   psum2[j]=0.0;
   for(int i=j*blsize; i<blsize*(j+1); i++) {
     psum1[j]+=dt[i];
     psum2[j]+=dt2[i];
   }
   tsum1+=psum1[j];
   tsum2+=psum2[j];
 }

 avr = err = 0.;
 for(int j=0; j<nblocks; j++) {
   avp1[j]=(tsum1-psum1[j])/(blsize*(nblocks-1));
   avp2[j]=(tsum2-psum2[j])/(blsize*(nblocks-1));
   avr += avp1[j]/avp2[j];
   err += avp1[j]*avp1[j] / (avp2[j]*avp2[j]);
 }

 avr /= static_cast<double>(nblocks);
 err -= avr * avr * nblocks;
 err *= static_cast<double>(nblocks-1) / static_cast<double>(nblocks);
 err = sqrt(err);

 /*for(int j=blsize*nblocks; j<size; j++) {
   tsum1+=dt1[j];
   tsum2+=dt2[j];

 }
 avr = tsum/static_cast<double>(size);*/

}


int main(int argc, char*argv[]) {
   int i,blsize;
   rlxs_init(0,52342373);

   if (argc<2) {
       printf ("Not enough parameters\n", argc);
       exit(1);
   }
   else {
       blsize = atoi(argv[1]);
       //printf ("%f %f %f %d\n",Ls, Lt, BETA, nrun);
   }

 vector<double> O[nrun],O2[nrun];

 //Fill BETA, Action, O
 for (int i=0; i<nrun; i++){
     sprintf(file, "mu_%1.0fx%1.0f_%g_%d_%d",Lt,Ls,BETA,nrun,i);
     ifstream in(file);
     if(!in) {
         cerr<<"Cannot open file: "<<file<<endl;
         exit(1);
     }
     double deltaS;
     while(in>>deltaS){
         O[i].push_back( exp(- deltaS));
         O2[i].push_back(exp(  deltaS));
     }

     cout<<BETA<<" "<<i<<" ";
     double val, err;
     AutoCorr(O[i], val, err);
     cerr<<val<<" "<<err<<" ";
     AutoCorr(O2[i], val, err);
     cerr<<val<<" "<<err<<" ";
     //int blsize=static_cast<int>(val*3.); //block size for jacknife 3*autocorr time

     if(blsize<1) blsize=1;
     JackKnife(O[i], O2[i], blsize, val, err);
     cout<<val<<" "<<err<<" "<<endl;

 }

}

