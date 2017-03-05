/* 
 * Splits domain in two by a plane with
 * normal direction idir and position at pos
 * Author: andrei.v.smirnov@gmail.com
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#define REAL	double
#define MAXLINLEN	510
#define WORDLENGTH	64
#define DIM	3
#define POS	0.0 // partition plane position
#define DIR	2 // partition plane direction

float rnd() { //Returns a random number between 0 and 1
	return (float)((double)rand()/RAND_MAX);
}

int main(int argc, char *argv[]) {
	uLong	nbuf=MAXLINLEN*sizeof(int);
	Byte *buf=(Byte*)calloc((uInt)nbuf,1);
	char *inpname,*outname,*oldname,*newname;
	gzFile	inp,out;
	int nmolecules=0;
	REAL pos=POS;
	//Check arguments
	if(argc<5||argc>6){
		printf("Usage: %s input.gz output.gz old new [pos]\n",argv[0]);
		printf("where:'old' is the old moleule name to be renamed into 'new' at the height above pos\n");
		return 1;
	}
	inpname=argv[1];
	//Open input file
	if((inp = gzopen(inpname, "rb"))==NULL)
	{	fprintf(stderr,"CAN'T OPEN %s\n",inpname);fflush(stderr);
		return 2;
	}
	//Open output file:
	outname=argv[2];
	if((out = gzopen(outname, "wb"))==NULL)
	{	fprintf(stderr,"CAN'T OPEN %s\n",outname);fflush(stdout);
		return 3;
	}
	oldname=argv[3];
	newname=argv[4];
	if(argc==6)pos=atof(argv[5]);
	printf("Replacing name '%s' with '%s' for all across plane at position=%g in direction=%d\n",oldname,newname,pos,DIR);
	//Perform I/O:
	while(!gzeof(inp)) {
		// READ LINE:
		char *s=(char*)buf,
			name[WORDLENGTH];
		REAL X[DIM],V[DIM],energy;
		int i=0;
		gzgets(inp, (char*)buf,nbuf);
		if((char)buf[0]=='#' || isdigit((char)buf[0])){
			gzputs(out, (char*)buf);
			continue;
		}
		sscanf(s,"%s",&name);// molecule name
		while(!isspace(*++s));
		while(isspace(*++s));
		for (i=0;i<DIM;i++) {//Read coordinates:
			REAL x;
			sscanf(s,"%lg",&x);
			X[i]=x;
			while(!isspace(*++s));
			while(isspace(*++s));
		}
		for (i=0;i<DIM;i++) {//Read velocity:
			REAL v;
			sscanf(s,"%lg",&v);
			V[i]=v;
			while(!isspace(*++s));
			while(isspace(*++s));
		}
		{	//Read energy
			sscanf(s,"%lg",&energy);
			while(!isspace(*++s));
			while(isspace(*++s));
		}
		// PERFORM CONVERSION
		if(X[DIR]>=pos&&strcmp(name,oldname)==0) strcpy(name,newname);
		// WRITE LINE:
		gzprintf(out,"%s",name);//write molecule name.
		for(i=0;i<DIM;i++)gzprintf(out,"\t%g",X[i]);
		for(i=0;i<DIM;i++)gzprintf(out,"\t%g",V[i]);
		gzprintf(out,"\t%g\n",energy);
		nmolecules++;
	}
	printf("Converted %d molecules from '%s' to '%s'\n",nmolecules,inpname,outname);
	gzclose(inp);
	gzclose(out);
}


