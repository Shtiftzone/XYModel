//Implementation of class MCData

#include "MCData.h"
#include "ranlxs.h"

#include <iostream>
#include <cmath>

#include "commondef.h"

void MCData::print(vector<MCObs> O) {
  vector<MCObs>::iterator Ocur = O.begin();
  vector<MCObs>::iterator Olast = O.end();
  for(;Ocur!=Olast;Ocur++) {
    AutoCorr(*Ocur);
    JackKnifeAvr(data[*Ocur], rep[*Ocur].block, rep[*Ocur].avr, rep[*Ocur].avr_err);
    JackKnifeSusc(data[*Ocur], data[*Ocur], rep[*Ocur].block, rep[*Ocur].susc, rep[*Ocur].susc_err);
    std::cout << rep[*Ocur];
  }
}

void MCData::AutoCorr(vector<MCObs> O) {
  vector<MCObs>::iterator Ocur = O.begin();
  vector<MCObs>::iterator Olast = O.end();
  for(;Ocur!=Olast;Ocur++) 
    AutoCorr(*Ocur);
}

void MCData::printdata(MCObs O) {
  vector<double>::iterator cur = data[O].begin();
  vector<double>::iterator last = data[O].end();
  for(int i=0; cur!=last; cur++) {
    double d = (*cur);
    std::cerr << rep[O].therm+(i++)*rep[O].samplingrate << " " << d << std::endl;
  }
}

void MCData::clear() {
  data.clear();
  rep.clear();
}

void MCData::Avr(MCObs O, double &avr, double &avr_err) {
  AutoCorr(O);
  JackKnifeAvr(data[O], rep[O].block, rep[O].avr, rep[O].avr_err);
  avr = rep[O].avr;
  avr_err = rep[O].avr_err;
}

double MCData::SimpleAvr(MCObs O) {
  double avr=0.0;
  int len = data[O].size();
  for (int i=0; i<len; i++) 
    avr += data[O][i];
  return avr/len;
}

double MCData::SimpleErr(MCObs O) {
  double avr, avr2;
  avr=avr2=0.0;
  int len = data[O].size();
  for (int i=0; i<len; i++) {
    double d=data[O][i];
    avr += d;
    avr2 += d*d;
  }
  avr/=len;
  avr2/=len;
  return sqrt(avr2-avr*avr);
  
}

void MCData::Rw(MCData &M, MCObs O, double deltabeta, int k) {
  if (!M.data.count(Action)) return;
  vector<double> &Nrg = M.data[Action];

  int len = M.data[O].size();
  if (!len) return;

  int sr = M.rep[O].samplingrate;
  int tO = sr*(M.rep[O].therm/sr);
  int tE = M.rep[Action].therm;

  int dE, dO;
  dE=dO=0;
  if(tE>tO) dO=(tE-tO)/sr;
  else dE=(tO-tE);

  len-=dO;
  if(len*sr>(Nrg.size()-dE))
    len=(Nrg.size()-dE)/sr;
  
  double EAvr, EAvr_err;
  M.Avr(Action, EAvr, EAvr_err);

  double lognorm; //Logaritmo della normalizzazione
  lognorm = -deltabeta*(Nrg[dE]-EAvr);
  for (int i=1; i<len; i++) {
    double le = -deltabeta*(Nrg[dE+i*sr]-EAvr);
    if(le<lognorm)
      lognorm += log1p(exp(le-lognorm));
    else
      lognorm = le + log1p(exp(lognorm-le));
  }

  for (int i=0; i<len; i++) {
    //    std::cout<<"idx: " << dE+i*sr << " " << (dO+i)*sr << std::endl;
    double logO = static_cast<double>(k)*log(M.data[O][dO+i]) - deltabeta*(Nrg[dE+i*sr]-EAvr);
    logO -= lognorm;
    data[O].push_back(len*exp(logO));
  }

  rep[O].samplingrate = sr;
  rep[O].therm = M.rep[O].therm+dO*sr;

}

void MCData::Bootstrap(MCData &M, MCObs O) {
  int len = M.data[O].size();
  if (!len) return;
  
  clear();

  if(O!=Action) {  
    int sr = M.rep[O].samplingrate;
    int tO = sr*(M.rep[O].therm/sr);
    int tE = M.rep[Action].therm;

    int dE, dO;
    dE=dO=0;
    if(tE>tO) dO=(tE-tO)/sr;
    else dE=(tO-tE);

    len-=dO;
    if(len*sr>(M.data[Action].size()-dE))
      len=(M.data[Action].size()-dE)/sr;
    
    
    float dt[len];
    ranlxs(dt, len);

    double ac1=M.rep[Action].autocorr;    
    double ac2=M.rep[O].autocorr;
    int ac = static_cast<int>((ac1>ac2)?ac1:ac2);
    if (ac<1) ac=1;

    for(int i=0; ac*i<len; i++) {
      int idx = static_cast<int>(len*dt[i])+dO;
      data[O].push_back(M.data[O][idx]);
      data[Action].push_back(M.data[Action][tO-tE+idx*sr]);
    }

    rep[O].samplingrate = M.rep[O].samplingrate;
    rep[O].therm = M.rep[O].therm+sr*dO;
    rep[Action].samplingrate = M.rep[O].samplingrate;
    rep[Action].therm = M.rep[Action].therm+dE;
  } else {
    float dt[len];
    ranlxs(dt, len);
    int ac = static_cast<int>(M.rep[Action].autocorr);
    if (ac<1) ac=1;
    //Resize Action vector
    //data[Action].resize(M.data[Action].size());
    for(int i=0; ac*i<len; i++) {
      int idx = static_cast<int>(len*dt[i]);
      data[Action].push_back(M.data[Action][idx]);
    }
    rep[Action].samplingrate = M.rep[Action].samplingrate;
    rep[Action].therm = M.rep[Action].therm;    
  }
}

void MCData::AutoCorr(MCObs O) {
  if(!data.count(O)) return; //Empty data set;

  vector<double> &dt = data[O];
  int len = dt.size();
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

  //t=1..len/2-1
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
    //std::cerr<<M<<"\t"<<tint<<std::endl;
    //Check for end condition M/5>=tint
    if(M>4.0*tint) {valid=true; break;}
  }

  double min_tint = 1.0;
  if(valid) {
    if(tint<min_tint) tint = min_tint; 
    rep[O].autocorr=tint;
    rep[O].autocorr_err=sqrt(2.0*(2*M+1)/len)*tint;
    rep[O].block=static_cast<int>(3.0*tint);
  } else {
    rep[O].autocorr=rep[O].autocorr_err=0.0; //Not Valid value
    rep[O].block=static_cast<int>(3.0*min_tint); //Standard Block Size
  }
}

void MCData::JackKnifeAvr(const vector<double> &dt, int blsize, double &avr, double &err){

  int size = dt.size();
  int nblocks = size/blsize; //Numero di blocchi

  double psum[nblocks]; //Somme parziali del blocco
  double avp[nblocks]; //Medie togliendo un blocco
  double tsum=0.0; //Somma totale

  for(int j=0; j<nblocks; j++) {
    psum[j]=0.0;
    for(int i=j*blsize; i<blsize*(j+1); i++) {
      psum[j]+=dt[i];
    }
    tsum+=psum[j];
  }

  avr = err = 0.;
  for(int j=0; j<nblocks; j++) {
    avp[j]=(tsum-psum[j])/(blsize*(nblocks-1));
    avr += avp[j];
    err += avp[j]*avp[j];
  }

  avr /= static_cast<double>(nblocks);
  err -= avr * avr * nblocks;
  err *= static_cast<double>(nblocks-1) / static_cast<double>(nblocks);
  err = sqrt(err);

  for(int j=blsize*nblocks; j<size; j++)
    tsum+=dt[j];
  avr = tsum/static_cast<double>(size);

}

void MCData::JackKnifeAvr(MCObs O, double &avr, double &err){

  if(!data.count(O)) { //Empty data set;
    avr=err=0.0;
    return;
  }

  vector<double> &dt = data[O];
  if(rep[O].block==datarep::DR_NOTVALID) AutoCorr(O);
  int blsize=rep[O].block;
  int size = dt.size();
  int nblocks = size/blsize; //Numero di blocchi

  double psum[nblocks]; //Somme parziali del blocco
  double avp[nblocks]; //Medie togliendo un blocco
  double tsum=0.0; //Somma totale

  for(int j=0; j<nblocks; j++) {
    psum[j]=0.0;
    for(int i=j*blsize; i<blsize*(j+1); i++) {
      psum[j]+=dt[i];
    }
    tsum+=psum[j];
  }

  avr = err = 0.;
  for(int j=0; j<nblocks; j++) {
    avp[j]=(tsum-psum[j])/(blsize*(nblocks-1));
    avr += avp[j];
    err += avp[j]*avp[j];
  }

  avr /= static_cast<double>(nblocks);
  err -= avr * avr * nblocks;
  err *= static_cast<double>(nblocks-1) / static_cast<double>(nblocks);
  err = sqrt(err);

  for(int j=blsize*nblocks; j<size; j++)
    tsum+=dt[j];
  avr = tsum/static_cast<double>(size);

}

void MCData::JackKnifeSusc(const vector<double> &dt1, const vector<double> &dt2,
			   const int blsize, double &susc, double &err){

  int sz1 = dt1.size(); int sz2 = dt2.size();
  int sz = (sz1>sz2)?sz2:sz1;
  std::cout << "Size: " << sz << std::endl;
  int nblocks = sz/blsize;

  double xpsum[nblocks], ypsum[nblocks], zpsum[nblocks]; //Somme parziali del blocco
  double avp[nblocks]; //Medie togliendo un blocco
  double xtsum=0.0, ytsum=0.0, ztsum=0.0; //Somme totali
  
  for(int j=0; j<nblocks; j++) {
    xpsum[j]=ypsum[j]=zpsum[j]=0.0;
    for(int i=j*blsize; i<blsize*(j+1); i++) {
      xpsum[j]+=dt1[i]*dt2[i];
      ypsum[j]+=dt1[i];
      zpsum[j]+=dt2[i];
    }
    xtsum+=xpsum[j];
    ytsum+=ypsum[j];
    ztsum+=zpsum[j];
  }

  susc = err = 0.;
  for(int j=0; j<nblocks; j++) {
    avp[j]=((xtsum-xpsum[j])-(ytsum-ypsum[j])*(ztsum-zpsum[j])/(blsize*(nblocks-1)))/(blsize*(nblocks-1));
    susc += avp[j];
    err += avp[j]*avp[j];
  }

  susc /= static_cast<double>(nblocks);
  err -= susc * susc * nblocks;
  err *= static_cast<double>(nblocks-1) / static_cast<double>(nblocks);
  err = sqrt(err);

  for(int j=blsize*nblocks; j<sz; j++) {
    xtsum+=dt1[j]*dt2[j];
    ytsum+=dt1[j];
    ztsum+=dt2[j];
  }
  susc = (xtsum-(ytsum*ztsum)/sz)/sz;

}

void MCData::JackKnifeSusc(MCObs O1, MCObs O2,
			   double &susc, double &err){

  if(!data.count(O1) || !data.count(O2)) { //Empty data set;
    susc=err=0.0;
    return;
  }

  vector<double> &dt1=data[O1];
  vector<double> &dt2=data[O2];
  if(rep[O1].block==datarep::DR_NOTVALID) AutoCorr(O1);
  if(rep[O2].block==datarep::DR_NOTVALID) AutoCorr(O2);
  int sr1=rep[O1].samplingrate;
  int sr2=rep[O2].samplingrate;
  int bl1=rep[O1].block*sr1;
  int bl2=rep[O2].block*sr2;
  int blsize=(bl1>bl2)?bl1:bl2;
  int sz1 = dt1.size()*rep[O1].samplingrate; 
  int sz2 = dt2.size()*rep[O2].samplingrate;
  int sz = (sz1>sz2)?sz2:sz1;
  //  std::cout << "Size: " << sz << std::endl;
  int nblocks = sz/blsize;

  double xpsum[nblocks], ypsum[nblocks], zpsum[nblocks]; //Somme parziali del blocco
  double avp[nblocks]; //Medie togliendo un blocco
  double xtsum=0.0, ytsum=0.0, ztsum=0.0; //Somme totali
  
  for(int j=0; j<nblocks; j++) {
    xpsum[j]=ypsum[j]=zpsum[j]=0.0;
    for(int i=j*blsize; i<blsize*(j+1); i++) {
      xpsum[j]+=dt1[i/sr1]*dt2[i/sr2];
      ypsum[j]+=dt1[i/sr1];
      zpsum[j]+=dt2[i/sr2];
    }
    xtsum+=xpsum[j];
    ytsum+=ypsum[j];
    ztsum+=zpsum[j];
  }

  susc = err = 0.;
  for(int j=0; j<nblocks; j++) {
    avp[j]=((xtsum-xpsum[j])-(ytsum-ypsum[j])*(ztsum-zpsum[j])/(blsize*(nblocks-1)))/(blsize*(nblocks-1));
    susc += avp[j];
    err += avp[j]*avp[j];
  }

  susc /= static_cast<double>(nblocks);
  err -= susc * susc * nblocks;
  err *= static_cast<double>(nblocks-1) / static_cast<double>(nblocks);
  err = sqrt(err);

  for(int j=blsize*nblocks; j<sz; j++) {
    xtsum+=dt1[j/sr1]*dt2[j/sr2];
    ytsum+=dt1[j/sr1];
    ztsum+=dt2[j/sr2];
  }
  susc = (xtsum-(ytsum*ztsum)/sz)/sz;

}
