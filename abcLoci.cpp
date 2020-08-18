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
#include "abcFilter.h"
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

void writeIndFasta(kstring_t *bufstr,size_t len,char *nam,char *d, int start, int stop, int nbpl){

	fprintf(stderr,"\t-> [%s] writing loci: %s:%d-%d\n",__FUNCTION__,nam, start, stop);
	ksprintf(bufstr,">%s_%d-%d",nam, start, stop);
	if(nbpl != 0){
		for(size_t i=0;i<len;i++){
			if(i % nbpl == 0)
				aio::kputc('\n',bufstr);
			aio::kputc(d!=NULL?d[i]:'N',bufstr);
		}
	}else{
		aio::kputc('\n',bufstr);
		for(size_t i=0;i<len;i++){
			aio::kputc(d!=NULL?d[i]:'N',bufstr);
		}
	}
	aio::kputc('\n',bufstr);

}


void writeLociFasta(kstring_t *bufstr,size_t len,char *nam,char *d, int start, int stop, int nbpl, char *ind){

	fprintf(stderr,"\t-> [%s] writing loci: %s:%d-%d\n",__FUNCTION__,nam, start, stop);
	ksprintf(bufstr,">%s",ind);
	if(nbpl != 0){
		for(size_t i=0;i<len;i++){
			if(i % nbpl == 0)
				aio::kputc('\n',bufstr);
			aio::kputc(d!=NULL?d[i]:'N',bufstr);
		}
	}else{
		aio::kputc('\n',bufstr);
		for(size_t i=0;i<len;i++){
			aio::kputc(d!=NULL?d[i]:'N',bufstr);
		}
	}
	aio::kputc('\n',bufstr);

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

	if (doIndFasta){
		writeIndFasta(&bufstr, wl.stop-wl.start,header->target_name[wl.refid],wl.wFasta,wl.start+1,wl.stop,NbasesPerLine);
		aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
	}

	if (doLociFasta){
		fprintf(stderr,"IM WRITING %s %d %d", header->target_name[wl.refid], wl.start+1, wl.stop);
		outfileZ = aio::appFileBG(locus.c_str(),postfix);
		writeLociFasta(&bufstr,wl.stop-wl.start,header->target_name[wl.refid],wl.wFasta,wl.start+1,wl.stop,NbasesPerLine,ind);
		aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
	}

	changeChr(-1);

	if(outfileZ!=NULL) bgzf_close(outfileZ); 
	if(bufstr.s!=NULL) free(bufstr.s);
	if(outfile!=NULL) fclose(outfile);
}


void abcLoci::print(funkyPars *pars){


	extern abc **allMethods;
	filt *fl = ((abcFilter *) allMethods[0])->fl;

	extern size_t total_number_of_sites_filtered;

	//FIXME
	//fprintf(stderr,"offs size %s\n", fl->offs.find(header->target_name[pars->refId])->second);
	fprintf(stderr,"offs BEGIN %s\n", fl->keeps[]);
	//fprintf(stderr,"\nHERE %lu\n", total_number_of_sites_filtered);

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

		//FIXME
		//fprintf(stderr,"\nhere1 %d\n", fl->offs.find(header->target_name[pars->refId])->second);
		//fprintf(stderr,"\nhere1 %lu\n", fl->keeps);
		//fprintf(stderr,"LOVE %d" ,pars->posi[0]);

		fprintf(stderr, "HELP %d", pars->numSites);
		//loop over sites
		for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){

		//fprintf(stderr,"\nhere2 %d\n", fl->offs);
			//FIXME
			if(pars->keepSites[s]==0)continue;
			//if(fl->keeps[pars->posi[s]]==0)continue;
			//fprintf(stderr,"\nPOSI %d\n", pars->posi[s]);
			fprintf(stderr,"\nHERE1 %d\n", fl->keeps[pars->posi[s]]);
			//fprintf(stderr,"\nHERE2 %d\n", pars->keepSites[s]);

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
	// print fasta
	if(doLoci==2){

		if (doIndFasta){
			hasData=1;

			//rawseqdata is in chunkyT struct (bambi_interface.h)
			chunkyT *chk = pars->chk;
			currentChr=pars->refId;

			// if it's a different locus
			//FIXME
			//fprintf(stderr, "START AT %d\n", chk->regStart);
			if(iPos == 0){
				// init region fasta
				free(rFasta);
				rFasta=(char*)malloc(header->target_len[currentChr]);
				memset(rFasta,'N', header->target_len[currentChr]);
			}

			if(doIndFasta==1){//random number read
				for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[currentChr];s++){
					if(pars->keepSites[s]==0) continue;
					iPos=pars->posi[s];
					if(pars->nInd==1)
						rFasta[iPos] = intToRef[ angsd::getRandomCount(pars->counts[s],0) ];
					else
						rFasta[iPos] = intToRef[ angsd::getRandomCountTotal(pars->counts[s],pars->nInd) ];
				}     
			}

			else if(doIndFasta==2) {//most common
				for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
					if(pars->keepSites[s]==0) continue;
					iPos=pars->posi[s]-rStart;
					if(pars->nInd==1)
						rFasta[iPos] = intToRef[ angsd::getMaxCount(pars->counts[s],0) ];
					else
						rFasta[iPos] = intToRef[ angsd::getMaxCountTotal(pars->counts[s],pars->nInd) ];
				}
			}

			else if(doIndFasta==3){
				for(int i=0;i<pars->nInd;i++){
					for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
						if(pars->keepSites[s]==0)continue;
						iPos=pars->posi[s]-rStart;
						tNode *tn = pars->chk->nd[s][i];
						if(tn==NULL) continue;
						double ebds[]= {0.0,0.0,0.0,0.0};
						for(int b=0;b<tn->l;b++){
							int bof = refToInt[tn->seq[b]];
							if(bof==4) continue;
							ebds[bof] += exp(lphred[tn->qs[b]]+lphred[tn->mapQ[b]]);
						}

						for(int b=0;0&&b<4;b++)
							fprintf(stderr,"b:%d %f\n",b,ebds[b]);
						int wh = angsd::whichMax(ebds,4);
						if(wh==-1) wh=4;//catch no information
						rFasta[iPos] = intToRef[wh];
					}
				}
			}else if(doIndFasta==4){
				//supplied by kristian ullrich
				for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
					if(pars->keepSites[s]==0) continue;
					iPos=pars->posi[s]-rStart;
					if(pars->nInd==1){
						rFasta[iPos] = intToIupac[ angsd::getIupacCount(pars->counts[s],0,iupacRatio) ];
					}else{
						rFasta[iPos] = intToIupac[ angsd::getIupacCountTotal(pars->counts[s],pars->nInd,iupacRatio) ];
					}
				}
			}
			//Do transitions removal
			if(rmTrans){
				assert(pars->ref!=NULL);
				for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){

					int ob = refToInt[rFasta[pars->posi[s]]];
					int rb = refToInt[pars->ref[pars->posi[s]]];
					//A <-> G, C <-> T
					if(ob==0&&rb==2)
						rFasta[pars->posi[s]]=intToRef[4];
					else if(ob==2&&rb==0)
						rFasta[pars->posi[s]]=intToRef[4];
					else if(ob==1&&rb==3)
						rFasta[pars->posi[s]]=intToRef[4];
					else if(ob==3&&rb==1)
						rFasta[pars->posi[s]]=intToRef[4];
				}
			}

			if(wl.wFasta){
				//if not first run
				if (wl.start == rStart && wl.stop == rStop && wl.refid == pars->refId){
					//still looping through chunks
					free(wl.wFasta);
					wl.wFasta=(char*)malloc(rSites);
					//memset(wl.wFasta,'N',rSites);
					//snprintf(wl.wFasta, sizeof(wl.wFasta), "%s", rFasta);
					strncpy(wl.wFasta, rFasta, rSites);

				}else{
					//catch reg change 
					fprintf(stderr,"REGCHANGE");
					//write previous reg
					//if(rFasta!=NULL){//proper case we have data
					//writeIndFasta(&bufstr,rSites,header->target_name[currentChr],wl.wFasta,wl.start+1,wl.stop,NbasesPerLine);
					//FIXME
					fprintf(stderr,"IM WRITING %s %d %d", header->target_name[wl.refid], wl.start+1, wl.stop);
					writeIndFasta(&bufstr,wl.stop-wl.start,header->target_name[wl.refid],wl.wFasta,wl.start+1,wl.stop,NbasesPerLine);
					aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
					iPos=0;
					wl.start = rStart;
					wl.stop = rStop;
					wl.refid = pars->refId;
					free(wl.wFasta);
					wl.wFasta=(char*)malloc(rSites);
					//memset(wl.wFasta,'N',rSites);
					strncpy(wl.wFasta, rFasta, rSites);
					//snprintf(wl.wFasta, sizeof(wl.wFasta), "%s", rFasta);
					//}
				}

			}else{

				//if first run, init
				fprintf(stderr, "INIT");
				wl.start = rStart;
				wl.stop = rStop;
				wl.refid = pars->refId;
				wl.wFasta=(char*)malloc(rSites);
				//memset(wl.wFasta,'N',rSites);
				strncpy(wl.wFasta, rFasta, rSites);
				//snprintf(wl.wFasta, sizeof(wl.wFasta), "%s", rFasta);
			}
		}

		// give nice output to be directly used in phylogenetic analyses
		// create fasta file for each locus
		// same loci from different individuals will be in the same loci fasta file
		// sequences will be named by individual names (bam file names)
		if(doLociFasta){
			hasData=1;

			char ind[fn.size() + 1];
			strcpy(ind, fn.c_str());


			chunkyT *chk = pars->chk;
			// start of the region -0 indexed
			rStart = chk->regStart;
			// end of the region -0 indexed
			rStop = chk->regStop;
			// number of sites in rf
			rSites=rStop-rStart;
			currentChr=pars->refId;

			locus = "";
			cName = header->target_name[currentChr];
			lStart = std::to_string(rStart+1);
			lStop = std::to_string(rStop);
			locus += cName + std::string("_") + lStart + std::string("-") + lStop;

			postfix=".fa.gz";
			outfileZ = aio::appFileBG(locus.c_str(),postfix);


			// if it's a different locus
			if(iPos == 0){
				free(rFasta);
				rFasta=(char*)malloc(rSites);
				memset(rFasta,'N',rSites);
			}


			if(doLociFasta==1){//random number read
				for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
					if(pars->keepSites[s]==0) continue;
					iPos=pars->posi[s]-rStart;
					if(pars->nInd==1)
						rFasta[iPos] = intToRef[ angsd::getRandomCount(pars->counts[s],0) ];
					else
						rFasta[iPos] = intToRef[ angsd::getRandomCountTotal(pars->counts[s],pars->nInd) ];
				}     
			}

			else if(doLociFasta==2) {//most common
				for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
					if(pars->keepSites[s]==0) continue;
					iPos=pars->posi[s]-rStart;
					if(pars->nInd==1)
						rFasta[iPos] = intToRef[ angsd::getMaxCount(pars->counts[s],0) ];
					else
						rFasta[iPos] = intToRef[ angsd::getMaxCountTotal(pars->counts[s],pars->nInd) ];
				}
			}

			else if(doLociFasta==3){
				for(int i=0;i<pars->nInd;i++){
					for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
						if(pars->keepSites[s]==0) continue;
						iPos=pars->posi[s]-rStart;
						tNode *tn = pars->chk->nd[s][i];
						if(tn==NULL) continue;
						double ebds[]= {0.0,0.0,0.0,0.0};
						for(int b=0;b<tn->l;b++){
							int bof = refToInt[tn->seq[b]];
							if(bof==4) continue;
							ebds[bof] += exp(lphred[tn->qs[b]]+lphred[tn->mapQ[b]]);
						}

						for(int b=0;0&&b<4;b++)
							fprintf(stderr,"b:%d %f\n",b,ebds[b]);
						int wh = angsd::whichMax(ebds,4);
						if(wh==-1) wh=4;//catch no information
						rFasta[iPos] = intToRef[wh];
					}
				}
			}else if(doLociFasta==4){
				//supplied by kristian ullrich
				for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
					if(pars->keepSites[s]==0) continue;
					iPos=pars->posi[s]-rStart;
					if(pars->nInd==1){
						rFasta[iPos] = intToIupac[ angsd::getIupacCount(pars->counts[s],0,iupacRatio) ];
					}else{
						rFasta[iPos] = intToIupac[ angsd::getIupacCountTotal(pars->counts[s],pars->nInd,iupacRatio) ];
					}
				}
			}
			//Do transitions removal
			if(rmTrans){
				assert(pars->ref!=NULL);
				for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
					int ob = refToInt[rFasta[pars->posi[s]]];
					int rb = refToInt[pars->ref[pars->posi[s]]];
					//A <-> G, C <-> T
					if(ob==0&&rb==2)
						rFasta[pars->posi[s]]=intToRef[4];
					else if(ob==2&&rb==0)
						rFasta[pars->posi[s]]=intToRef[4];
					else if(ob==1&&rb==3)
						rFasta[pars->posi[s]]=intToRef[4];
					else if(ob==3&&rb==1)
						rFasta[pars->posi[s]]=intToRef[4];
				}
			}

			if(rFasta!=NULL)
				if(hasData)
					if(wl.wFasta){
						//if not first run
						if (wl.start == rStart && wl.stop == rStop && wl.refid == pars->refId){
							//still looping through chunks
							free(wl.wFasta);
							wl.wFasta=(char*)malloc(rSites);
							strncpy(wl.wFasta, rFasta, rSites);

						}else{
							//catch reg change 
							fprintf(stderr,"REGCHANGE");
							//write previous reg
							fprintf(stderr,"IM WRITING %s %d %d", header->target_name[wl.refid], wl.start+1, wl.stop);
							writeLociFasta(&bufstr,wl.stop-wl.start,header->target_name[wl.refid],wl.wFasta,wl.start+1,wl.stop,NbasesPerLine,ind);
							aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
							if(outfileZ!=NULL) {
								fprintf(stderr,"CLOSE");
								bgzf_close(outfileZ);
								outfileZ = NULL;
							}
							iPos=0;

							wl.start = rStart;
							wl.stop = rStop;
							wl.refid = pars->refId;
							free(wl.wFasta);
							wl.wFasta=(char*)malloc(rSites);
							//memset(wl.wFasta,'N',rSites);
							strncpy(wl.wFasta, rFasta, rSites);
							//snprintf(wl.wFasta, sizeof(wl.wFasta), "%s", rFasta);

						}

					}else{
						//if first run, init
						fprintf(stderr, "INIT");
						wl.start = rStart;
						wl.stop = rStop;
						wl.refid = pars->refId;
						wl.wFasta=(char*)malloc(rSites);
						//memset(wl.wFasta,'N',rSites);
						strncpy(wl.wFasta, rFasta, rSites);
					}
		}
		//if(outfileZ!=NULL) {
		//bgzf_close(outfileZ);
		//outfileZ = NULL;
		//}

	}

	if(doLoci==3){
		//one snp per tag
		return;
	}
}


void abcLoci::changeChr(int refId) {

	if(doIndFasta){
		if(refId!=-1){//-1 = destructor
			free(rFasta);
			rFasta=(char*)malloc(rSites);
			memset(rFasta,'N',rSites);
		}else{
			free(rFasta);
		}
	}

	if(doLociFasta){
		if(refId!=-1){//-1 = destructor
			free(rFasta);
			rFasta=(char*)malloc(rSites);
			memset(rFasta,'N',rSites);
		}else{
			free(rFasta);
		}
	}
}

void abcLoci::run(funkyPars *pars){
	return;
}

void abcLoci::clean(funkyPars *pars){
	return;
}
