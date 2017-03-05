//-#include <stdlib.h>
#include <string.h>
#include <math.h>
//-#include <time.h>
#include <sys/timeb.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <iostream>
#include "def.h"
#include "run.h"
#include "io.h"
using namespace std;
#ifndef	RAND_MAX
#define	RAND_MAX	2147483647
#endif
namespace Run {
	struct Option	option={0,0,0,0};
	void (*usage)();
	char
		programname[MAXLINLEN],
		configfile[MAXLINLEN],
		outputfile[MAXLINLEN],
		outputname[MAXLINLEN],
		inputfile[MAXLINLEN];
	struct Time	time;
	//! Test random numer generator
	void testrnd()	{
		for(int i=0;i<20;i++){
			long irnd=rand();
			double rnd=(double)irnd/(double)(RAND_MAX);
			REAL frnd=(REAL)rnd;
			cout<<i<<": irnd="<<irnd
			<<", RAND_MAX="<<RAND_MAX
			<<", rnd="<<rnd
			<<", frnd="<<frnd<<endl;
		}
	}
	//! Seeds random number
//-	#define	SEQUENCES	10000000
//-	void	timeseed() {
//-	//	REAL	r;
//-	//	randm_(&seed,&r);
//-		srandom(time(0) % SEQUENCES);
//-	//	srand(time(0) % SEQUENCES);
//-	}
	//! Generates random number between 0 and 1
	REAL rnd()
	{	return (REAL)((double)rand()/RAND_MAX);
	}
	/*! From Numerical Recepes in C, Ch.7.2 
	 * Zero mean and unit variance
	 */
	REAL	gauss() {
		static int iset = 0;
		static REAL gset;
		REAL fac,rsq,v1,v2;
		if (iset == 0)
		{	do {
				v1=2.*rnd()-1.;
				v2=2.*rnd()-1.;
				rsq=v1*v1+v2*v2;
			}	while (rsq >= 1.0 || rsq == 0);
			fac=(REAL)sqrt(-2.*log(rsq)/rsq);
			gset=v1*fac;
			iset=1;
			return v2*fac;
		}
		else
		{
			iset=0;
			return gset;
		}
	}
	//! \brief Initialize all variables, and execution options
	void	init(int argc, char *argv[])
	{	using namespace IO;
		struct timeb	worldtime;
		option.verbose=0;
		option.debug=0;
		option.mesh=0;
		option.xterm=0;
		option.restart=0;
		ftime(&worldtime);
		srandom(worldtime.time);//seed random numbers
		//srandom(777);//DEBUG
		readcmdline(argc, argv);
		getTime(configfile);
	}
	//! Parses command-line options
	void	readcmdline(int argc, char *argv[]) {
		char *p,*q;
		sprintf(configfile,"%s.xml",argv[0]);
		sprintf(outputfile,"%s.dat.gz",argv[0]);
		inputfile[0]='\0';
		strcpy(outputname,argv[0]);
		for (p=argv[0];!isalpha(*p)&&*p!='\0';p++);
		strcpy(programname,p);
		for (int i=1; i<argc; i++)
		{	if (*argv[i]=='-')
			switch (*(p=argv[i]+1))
			{	case 'f':
					if(i++>=argc)goto Usage;
					if(*argv[i]=='-')goto Usage;
					strcpy(configfile,argv[i]);
					break;
				case 'm':
					option.mesh=1;
					break;
				case 'v':
					option.verbose=1;
					cout<<"Verbose Mode\n";
					break;
				case 'V':
					option.debug=1;
					cout<<"Debug Mode\n";cout.flush();
					break;
				case 'o':
					if(i++>=argc)goto Usage;
					if(*argv[i]=='-')goto Usage;
					strcpy(outputfile,argv[i]);
					break;
				case 'h':
					goto Usage;
					break;
			} else {
				sprintf(inputfile,"%s",argv[i]);
			}
		}
		option.restart=inputfile[0]=='\0'?0:1;
		return;
		Usage:
		(*usage)();
		exit(1);
	}
}
