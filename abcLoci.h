#pragma once
#include "abc.h"

class abcLoci:public abc{
private:
  kstring_t bufstr;
  int doLoci;
  FILE *outfile;
	int iPos;
  char *ref;
  char *rFasta;//contains the new fastasequence for the current loci
  int hasData;
	int currentChr;
  int NbasesPerLine;
	int rmTrans;
	double lphred[256];//these are log phread scores log(10^(1:255)/(-10))
	double iupacRatio;
  int seed;
  std::string locus;
  std::string cName;
  std::string lStart;
  std::string lStop;
  std::string gfn;
  std::size_t pos1;
  std::size_t pos2;
  std::string fn;
  const char* postfix;
	int rStart;
	int rStop;
	int rSites;
  char *ind;

  typedef struct{

    char*wFasta; //fasta waiter

    int start; //start of locus
    int stop; //stop of locus
    int refid; //chr refid locus belongs to
    int guests; //guest loci waiting in the waiting room

  }waitLoci;

  waitLoci wl;

public:
  int doLociFasta;
  int doIndFasta;
  int doCount;

  BGZF *outfileZ;

  int l;
  char *buffer;
  char *original;

  void run(funkyPars  *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  void openfile(const char *outfiles);
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  void changeChr(int refId);
  abcLoci(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcLoci();
};
