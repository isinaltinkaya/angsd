/*
 *
 * fasta calling functions from abcWriteFasta are modified for loci
 * contact: isinaltinkaya@gmail.com
 *
 */

#include <ctype.h> //<-used for isupper/islower
#include <cmath> //in order to construct large array/objects
#include <stdlib.h>
#include <cstdlib>
#include <htslib/bgzf.h>
#include <assert.h>
#include <ctime>
#include <htslib/kstring.h>
#include <string>

#include "aio.h"
#include "shared.h" //<-contains the struct defintino for funkypars
#include "analysisFunction.h" //<-contains some utility functions, usefull for parsing args
#include "abcLoci.h"//contains the analysis class defition associated with this file
#include "abc.h"


void abcLoci::printArg(FILE *argFile){
	fprintf(argFile,"------------------------\n%s:\n",__FILE__);
	fprintf(argFile,"-doLoci\t\t%d (Which analysis should we perform?)\n",doLoci);
	fprintf(argFile,"\t\t1: Count and print basetypes in combination with strand\n");
	fprintf(argFile,"\t\t2: Call fasta for loci\n");
	fprintf(argFile,"\n");

	fprintf(argFile,"-doIndFasta (create fasta file for individual)\t%d\n",doIndFasta);
	fprintf(argFile,"-doLociFasta (create fasta file for each loci to be used in phylogenetic analyses)\t%d\n",doLociFasta);
	fprintf(argFile,"\t\t1: use a random (non N) base (needs -doCounts 1)\n");
	fprintf(argFile,"\t\t2: use the most common (non N) base (needs -doCounts 1)\n");
	fprintf(argFile,"\t\t3: use the base with highest ebd (under development) \n");
	fprintf(argFile,"\t\t4: output iupac codes (under development) \n");
	fprintf(argFile,"\t\t-basesPerLine\t%d\t(Number of bases perline in output file)\n",NbasesPerLine);
	fprintf(argFile,"\t\t-rmTrans\t%d\t remove transitions as different from -ref bases (0:no,1:yes)\n",rmTrans);
	fprintf(argFile,"\t\t-ref\t%s\t reference fasta, only used with -rmTrans 1\n",ref);
	fprintf(argFile,"\t\t-iupacRatio\t%.3f\t (Remove nucleotide below total depth ratio for IUPAC assignment)\n",iupacRatio);
	fprintf(argFile,"\t\t-seed\t%d\t use non random seed\n",seed);
	fprintf(argFile,"\n");
}

void abcLoci::getOptions(argStruct *arguments){

	doLoci=angsd::getArg("-doLoci",doLoci,arguments);
	if(doLoci==0){
		shouldRun[index]=0;
		return;
	}

	doLociFasta=angsd::getArg("-doLociFasta",doLociFasta,arguments);
	doIndFasta=angsd::getArg("-doIndFasta",doIndFasta,arguments);

	if ((doLociFasta) && (doIndFasta)){
		fprintf(stderr,"Error: -doLociFasta and -doIndFasta can not be used together\n");
		exit(0);
	}

	doCount=angsd::getArg("-doCounts",doCount,arguments);
	iupacRatio=angsd::getArg("-iupacRatio",iupacRatio,arguments);
	seed=angsd::getArg("-seed",seed,arguments);
	NbasesPerLine = angsd::getArg("-basesPerLine",NbasesPerLine,arguments);

	rmTrans=angsd::getArg("-rmTrans",rmTrans,arguments);
	ref=angsd::getArg("-ref",ref,arguments);
	if(rmTrans && ref==NULL){
		fprintf(stderr,"\t-> Must supply reference with -rmTrans 1\n");
		exit(0);
	}

	if ((doLoci == 1) && ((doLociFasta) || (doIndFasta))) {
		fprintf(stderr, "\tError: doLoci 1 cannot be used with doLociFasta or doIndFasta\n");
		printArg(stderr);
		exit(0);
	}

	if ((doLoci == 2) && ((doLociFasta == 0) && (doIndFasta == 0))) {
		fprintf(stderr, "\t-> Must supply -doLociFasta or -doIndFasta\n");
		printArg(stderr);
		exit(0);
	}

	if(doLociFasta){
		if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
			fprintf(stderr,"Error: bam or soap input needed for -doLociFasta \n");
			exit(0);
		}
		if((doLociFasta==2||doLociFasta==1) && doCount==0){
			fprintf(stderr,"Error: -doLociFasta 1 or 2 needs allele counts (use -doCounts 1)\n");
			exit(0);
		}
	}

	if(doIndFasta){
		if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
			fprintf(stderr,"Error: bam or soap input needed for -doIndFasta \n");
			exit(0);
		}
		if(doLociFasta==3 && arguments->nInd!=1){
			fprintf(stderr,"Error: -doIndFasta 3 only works for a single individual\n");
			exit(0);
		}
		if((doLociFasta==2||doLociFasta==1) && doCount==0){
			fprintf(stderr,"Error: -doIndFasta 1 or 2 needs allele counts (use -doCounts 1)\n");
			exit(0);
		}
	}
}


//constructor
abcLoci::abcLoci(const char *outfiles,argStruct *arguments,int inputtype){
	doLoci = 0; //defaults= dont do analysis
	doLociFasta = 0;
	doIndFasta = 0;
	outfile = NULL;
	rFasta = NULL;
	rSites=0;
	rmTrans = 0;
	currentChr=-1;
	iupacRatio = 0.0;
	NbasesPerLine=0;
	seed=0;
	ref = NULL;

	if(arguments->argc==2){
		if(!strcasecmp(arguments->argv[1],"-doLoci")){
			printArg(stdout);
			exit(0);
		} else
			return;
	}

	getOptions(arguments);

	if(doLoci==0){
		shouldRun[index] = 0;
		return;
	}

	printArg(arguments->argumentFile);

	if(seed) srand(seed);
	else srand(time(0));

	if(doLoci==0) return ;

	if(doLociFasta){
		//get filename from -i
		gfn = arguments->infile;
		//match filename in the filesystem
		pos1 = gfn.find_last_of("/\\");
		//get individual name
		gfn = gfn.substr(pos1+1);
		//match only with name without extension
		//this will exclude anything coming after a dot
		pos2 = gfn.find_last_of(".");
		fn = gfn.substr(0,pos2);
	}


	//initalize outputfiles
	if(doLoci==1) outfile = aio::openFile(outfiles,".results");

	if(doLoci==2){
		postfix=".fa.gz";
		for(int i=0;i<256;i++)
			lphred[i] =log(1.0-pow(i,-1.0*i/10.0));
		bufstr.s=NULL;
		bufstr.m=0;
		bufstr.l=0;

		if(doIndFasta){
			outfileZ = NULL;
			outfileZ = aio::openFileBG(outfiles,postfix);
		}

	}


}

//destructor
abcLoci::~abcLoci(){

	if(doLoci==0) return;
	if(ref) free(ref);

	if(outfileZ!=NULL) bgzf_close(outfileZ); 
	if(bufstr.s!=NULL) free(bufstr.s);
	if(outfile!=NULL) fclose(outfile);
}


void abcLoci::print(funkyPars *pars){

	if(doLoci==0) return;

	// count the number of A, C, G, T, N's for both strand
	if(doLoci==1){

		//rawseqdata is in chunkyT struct (bambi_interface.h)
		chunkyT *chk = pars->chk;
		// start of the region -0 indexed
		rStart = chk->regStart;
		// end of the region -0 indexed
		rStop = chk->regStop;
		// number of sites in rf
		rSites=rStop-rStart;

		//loop over sites
		for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){

			if(pars->keepSites[s]==0)continue;

			fprintf(stderr,"\r\t-> Printing at chr: %s pos:%d loci %d contains %d sites",header->target_name[pars->refId],pars->posi[s]+1,pars->chunkNumber,rSites);
			int bases[2][5] = {{0,0,0,0,0},{0,0,0,0,0}};      
			//loop over samples
			for(int i=0;i<chk->nSamples;i++){
				//all seqdata associated with single bamfile is in a tNode
				tNode *nd = chk->nd[s][i];
				//loop over the individual bases
				for(int l=0;l<nd->l;l++){
					char c = nd->seq[l]; //this is the base
					char q = nd->qs[l]; //this is the associated qscore, fancy shit
					int strand = isupper(nd->seq[l])==0; //strand is defined as either small/big letters
					bases[strand][refToInt[c]]++;

				}

				fprintf(outfile,"%s\t%d",header->target_name[pars->refId],pars->posi[s]+1);//position is zero index internally
				//print the basecount
				for(int i=0;i<2;i++)
					for(int j=0;j<5;j++)
						fprintf(outfile,"\t%d",bases[i][j]);
				fprintf(outfile,"\n");

			}
		}
	}
}


void abcLoci::changeChr(int refId) {
	return;
}

void abcLoci::run(funkyPars *pars){
	return;
}

void abcLoci::clean(funkyPars *pars){
	return;
}
