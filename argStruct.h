#pragma once
#include <cstring>
#include <vector>
#include <map>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <cstdio>
#include "sample.h"
//little struct for keeping information of regions to extract
typedef struct{
  int refID;
  int start;
  int stop;
}regs;

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

typedef std::map<const char *,int,ltstr> aMap;


typedef struct{
  int nInd;//number of inds inferred from filelists
  int argc;
  char **argv;
  int inputtype;//
  int *usedArgs; //array of ints telling if args been used
  FILE *argumentFile; //logfile
  char *outfiles; //prefix output
  bam_hdr_t *hd;
  // const aHead *hd;
  const aMap *revMap;
  char *infile;//contains, the -bam fname,-glf fname, -pileup fname
  std::vector<char *> nams;//contains either the above or the contents of -bam;
  std::vector<regs> regions;//regions to use -r/-rf when using seqdata bcf/vcf bam/cram
  int nLines;//nLines;
  int nReads;//number of reads to pop from each BAM/CRAM
  int show;
  char *fai;
  bam_sample_t *sm;//for dealing with readgroups
  char *ref;
  char *anc;
  char *cmdline;
  char *version;
}argStruct;

argStruct *setArgStruct(int argc,char **argv);
void destroy_argStruct(argStruct *args);
