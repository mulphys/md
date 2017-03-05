#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <zlib.h>

#define MAXLINLEN	510

using namespace std;
const int dim=3;
const int nmolecules=10;
const int natoms=1; //number of atoms for each molecule
const int atype=2; // type of atoms
const float moleculesize=200.0;//size of the molecule in Angstroms
const float boxsize=1000.0; //size of the simulation domain in Angstroms
const float vmax=10.0;//maximum velocity in m/s

float rnd() { //Returns a random number between 0 and 1
	return (float)((double)rand()/RAND_MAX);
}

main(int argc, char *argv[]) {
	char filename[MAXLINLEN];
	if(argc!=2) {
		cout << "Usage: " << argv[0] <<" outputfile\n";
		exit(1);
	}
	sprintf(filename,"%s.dat.gz",argv[1]);
	gzFile out;
	if((out = gzopen(filename, "wb"))==NULL)
	{	cerr << "CAN'T OPEN "<<filename<<endl;cerr.flush();
		exit(1);
	}
	cout<<"Writing "<<nmolecules<<" molecules to "<<filename<<" ... ";
	gzprintf(out,"# Initial data\n%d\n",nmolecules);
	int nx=(int)pow(nmolecules,1.0/3.0)+1,//number of molecules in one direction
		nx2=nx*nx;
	float dx=boxsize/nx; //molecular separation
	cout
		<<"Number of molecules: "<<nmolecules<<endl
		<<"Number of molecules in each direction: "<<nx<<endl
		<<"Molecular separation: "<<dx<<endl;
	srand(1);//seed random generator (see man rand) 
cout<<"nmolecules="<<nmolecules<<", nx="<<nx<<endl;///DDD
	for(int imolecule=0;imolecule<nmolecules;imolecule++) {
		gzprintf(out,"%g",moleculesize);
		int n=imolecule%nx2, m=n/nx, ind[dim];
		ind[0]=imolecule/nx2;
		ind[1]=m;
		ind[2]=n%nx;
		for(int i=0;i<dim;i++) {
			float x=-0.5*boxsize + ind[i]*dx; //random position between 0 and boxsize
			gzprintf(out,"\t%g",x);
		}
		for(int i=0;i<dim;i++) {//velocity of the molecule in direction ix
			float 
				v=(-1.0 + 2.0*rnd())*vmax;//random velocity between -vmax and +vmax
			gzprintf(out,"\t%g",v);
		}
		gzprintf(out,"\t%d\n",natoms);
		//Output Atoms:
		if(natoms>0) {
			for(int iatom=0;iatom<natoms;iatom++) {
				gzprintf(out,"\t%d\n",(int)atype);
			}
		} else cerr<<"NO ATOMS SPECIFIED FOR MOLECULE "<<imolecule<<endl;
	}
	cout << nmolecules << " molecules written\n";
	gzclose(out);
}


