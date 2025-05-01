#ifndef MULTIHISTRW_H 
#define MULTIHISTRW_H 

#include <vector>

using std::vector;

double LnZ(double beta, int nrun, double *Beta, vector<double> *d, double *LogZ);
double LnO_n(vector<double> *O, double n, double beta, int nrun, double *Beta, vector<double> *d, double *LogZ);
void MultiHistRw(int nrun, double *beta, vector<double> *d, double *LogZ);

#endif
