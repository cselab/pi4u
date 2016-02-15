#ifdef __cplusplus
extern "C" {
#endif

#ifndef MAXDIM
#define MAXDIM 64
#endif

// netname: file to store NN
// Xin: evaluation points
// Tin: target values
// D: number of inputs
// F: number of outputs (1)
// N: size of training set
// Eout: error for each point
void build_network(char *netname, double *Xin, double *Tin, int D, int F, int N, double *Eout);

// netname: file with stored NN
// Xin: evaluation point
// D: number of inputs
double predict_network(char *netname, double *Xin, int D);

#ifdef __cplusplus
}
#endif
