#include "Autocorr.h"

#include <cmath>
#include <iostream>

void AutoCorr(const vector<double> &dt, double &autocorr, double &err) {

  int len = dt.size();
  if(!len){  //Empty data set;
    autocorr=err=0.;
    return;
  }

  double C0, Ct, rho, tint;

  double avr=0.0; // Media f(x)
  double f2=0.0;  // Somma di f(x)f(x+t)
  double f=0.0;   // Somma di f(x)+f(x+t)

  //t=0
  for (int i=0; i<len; i++) {
    avr+=dt[i];
    f2+=dt[i]*dt[i];
  }
  f = 2.0 * avr;
  avr/=static_cast<double>(len);
  C0 = (f2/static_cast<double>(len))-avr*avr;
  tint = 0.5;

  bool valid=false; int M;
  for (M=1; M<len; M++) {
    f2=f=avr=0.0;
    for (int i=0; i<len-M; i++) {
      f2+=dt[i]*dt[i+M];
      f+=dt[i]+dt[i+M];
      avr+=dt[i];
    }
    avr/=static_cast<double>(len-M);
    Ct = (f2/static_cast<double>(len-M))+avr*(avr-(f/static_cast<double>(len-M)));
    rho = Ct/C0;
    tint += rho;
    //std::cerr<<M<<"\t"<<tint<<"\t"<<rho<<std::endl;
    //Check for end condition M/5>=tint
    if(M>4.0*tint) {valid=true; break;}
  }

  if(valid) {
    autocorr=tint;
    err=sqrt(2.0*(2*M+1)/len)*tint;
  } else {
    autocorr=err=0.0; //Not Valid value
  }

}
