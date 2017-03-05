#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <zlib.h>

#define MAXLINLEN	510
#define PI 3.14159265358979323846
#define BoltzmannConstant	1.38066e-20 //g*(nm/ns)^2/K
#define AtomicMassUnit	1.66053886e-24 // in g

//#define IGNITION //define igniion region (if compbustion)

using namespace std;
const int dim=3;
const float 
	boxsize=100.0, //size of the simulation domain in nanometers
	temperature=300; //K: temperature
#ifdef IGNITION
	const float ignition_temperature=20000.0,//K
		xignmin[]={-0.5*boxsize,-0.5*boxsize,-0.5*boxsize},
		xignmax[]={-0.4*boxsize,0.5*boxsize,0.5*boxsize};//ignition box
#endif
struct Species {
	char *id; 
	int	number;
	float mass;// in atomic units
}	species[]={
//	{"H2O", 10, 18.0},
	{"H2" , 600, 2.0},
	{"O2" , 300, 32.0}
//		{"H2S", 10, 34.},
//		{"N2", 10, 28.}
};

float rnd() { //Returns a random number between 0 and 1
	return (float)((double)rand()/RAND_MAX);
}

main(int argc, char *argv[]) {
	int 
		nspecies=sizeof(species)/sizeof(*species),
		nmolecules=0;
	cout<<"Number of species: "<<nspecies<<endl;
	for(int i=0;i<nspecies;i++) {
		Species *specie=species+i;
		nmolecules+=specie->number;
		specie->mass*=AtomicMassUnit;// convert to grams
	}
	cout<<"Total number of molecules: "<<nmolecules<<endl;
 	float *conc=new float[nspecies];//concentration of each specie
	for(int i=0;i<nspecies;i++) conc[i]=(float)species[i].number/(float)nmolecules;
	for(int i=0;i<nspecies;i++) {
		Species *specie=species+i;
		cout<<i<<": id="<<specie->id<<", number of molecules="<<specie->number<<", concentration="<<conc[i]<<endl;
	}
	//Write output file:
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
	int *spcount=new int[nspecies];//count of each specie molecules
	for(int i=0;i<nspecies;i++) spcount[i]=0;
	float x[dim],v[dim];
	for(int imolecule=0;imolecule<nmolecules;imolecule++) {
		//Select specie type:
		bool selected=false;
		int isp=0;
		do {
			float r=rnd(),
				threshold=0.0;//random number between 0 and 1
			for(isp=0;isp<nspecies;isp++) {
				threshold+=conc[isp];
				if(r<=threshold) {
					if(spcount[isp]<species[isp].number) {
						spcount[isp]++;
						selected=true;
						break;
					}
				}
			}
		} while(!selected);
		//Set coordinates:
		int n=imolecule%nx2, m=n/nx, ind[dim];
		ind[0]=imolecule/nx2;
		ind[1]=m;
		ind[2]=n%nx;
#ifdef IGNITION
		int ignbox=0;//count of ignition box hits
#endif
		for(int i=0;i<dim;i++) {
			x[i]=-0.5*boxsize + (ind[i]+0.5)*dx; //position between 0 and boxsize
#ifdef IGNITION
			if(x[i]>=xignmin[i]&&x[i]<=xignmax[i]) ignbox++;
#endif
		}
		//Set velocity:
		float local_temperature=temperature;
#ifdef IGNITION
		if(ignbox==3) {
			local_temperature=ignition_temperature;
//isp=0;///DDD
		}
#endif
		float
			mass=species[isp].mass,
			velocity=(float)pow(8.0*BoltzmannConstant*local_temperature/(PI*mass),0.5);
		for(int i=0;i<dim;i++) {//velocity of the molecule in direction ix
			v[i]=(-1.0 + 2.0*rnd())*velocity;//random velocity between -vmax and +vmax
		}
		gzprintf(out,"%s",species[isp].id);//write specie index (type of molecule).
		for(int i=0;i<dim;i++)gzprintf(out,"\t%g",x[i]);
		for(int i=0;i<dim;i++)gzprintf(out,"\t%g",v[i]);
		gzprintf(out,"\n");
	}
	cout << nmolecules << " molecules written\n";
	gzclose(out);
}


