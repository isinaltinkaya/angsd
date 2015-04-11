#pragma once
#include "abc.h"


typedef struct{
  int **dat;
}genoCalls;


class abcCallGenotypes:public abc{
private:
  int doGeno;
  float postCutoff;
  gzFile outfileZ;
  double geno_minMM;
  int geno_minDepth;
  int geno_maxDepth;
  int minInd;
public:
  abcCallGenotypes(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcCallGenotypes();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  void getGeno(funkyPars *pars);
  void printGeno(funkyPars *pars);
  
};
