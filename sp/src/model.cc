#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
using namespace std;
#include "def.h"
#include "run.h"
#include "list.h"
#include "collection.h"
#include "model.h"
#include "species.h"
using namespace Species;
#include "domain.h"

//-ATOM::ATOM(){	
//-	type=undefined;
//-	for(int i=0;i<DIM;i++) {
//-		x[i]=0.0;
//-	}
//-}
Molecule::Molecule() {
	type=-1;
#ifdef LOCAL
	gridcell=-1;
	ptr=NULL;
#endif
	energy=0.0;
}
void Molecule::copy(Molecule *molecule) {
	type=molecule->Type();
#ifdef LOCAL
	int	gridcell=molecule->GridCell();//compressed 3D index to grid cell
#endif
	for(int i=0;i<DIM;i++) {
		x[i]=molecule->Coordinate(i);
		v[i]=molecule->Velocity(i);
	}	
}
REAL Molecule::KineticEnergy() {
	REAL
		mass=species[type].Mass(), 
		v2=0.0;
	for(int i=0;i<DIM;i++) {
		REAL vel=v[i];
		v2+=vel*vel;
	}
	return 0.5*mass*v2;
}
//! Set internal energy according to the temperature
void Molecule::Temperature(REAL temperature //!< temperature in K
) {
	int ispecie=Type();
	REAL // Add energy per internal degree of freedom:
		energy=temperature*BoltzmannConstant/AtomicMassUnit,// au*(nm/ns)^2
		cp=species[ispecie].Cp(),// specific heat
		dof=2.0*(cp*GasConstantInv-1.0),//=2*cv/R=degrees of freedom
		idof=MAX(dof-3.0,0.0); // internal degrees of freedom = dof>3.0?dof-3.0:0.0;
	InternalEnergy(idof*energy);
}
//! Moving a molecule by one time step
void Molecule::Move() {
	REAL dt=Run::time.step,
		mass=species[type].Mass();
	if (mass<SMALL) {cerr<<"MOLECULAR MASS "<<mass<<" IS TOO SMALL\n";return;}
	for (int i=0; i<DIM; i++)
	{	REAL
	//		drag=0.0, //Cdrag*(fvel[i]-v[i]), 
			vnew=v[i];//DDD  +drag/mass*dt;
		x[i]+=0.5*(vnew+v[i])*dt;
		v[i]=vnew;
	}
}
namespace Potential {
	REAL sigma=1.0, eta=1.0, cutoff=2.0;
//-	REAL potential(List<ATOM> *nodes) { //Compute potential for the current node
//-		ATOM *current=nodes->Current(),*next;
//-		REAL	*x=current->Coordinates();
//-		REAL potential=0.0;
//-		//Loop through other nodes:
//-		while((next=nodes->Next())!=current) {
//-			REAL *y=next->Coordinates(),
//-				d=0.0;
//-			for(int i=0;i<DIM;i++) {
//-				REAL r=x[i]-y[i];
//-				d+=r*r;
//-			}
//-			potential+=value((REAL)sqrt(d));
//-		}
//-		return potential;
//-	}
	REAL invdist(REAL x) {
		if(x<0.5) return log(2.0*x);
		else return -log(1.0-2.0*(x-0.5));
		
	}
}
REAL Geom::normalize(REAL a[]) {
	REAL d=length(a);
	for(int i=0;i<DIM;i++)a[i]/=d;//make unit length
	return d;
}
void Geom::distvec(REAL a[], REAL b[], REAL c[]) {
	for (int i=0;i<DIM;i++) c[i]=b[i]-a[i];
}
REAL Geom::distance(REAL a[], REAL b[]) {
	REAL d=0.0;
	for (int i=0;i<DIM;i++) {
		REAL r=b[i]-a[i];
		d+=r*r;
	}
	return sqrt(d);
}
REAL Geom::distance(int dim, REAL a[], REAL b[]) {
	REAL d=0.0;
	for (int i=0;i<dim;i++) {
		REAL r=b[i]-a[i];
		d+=r*r;
	}
	return sqrt(d);
}
REAL Geom::sclp(REAL a[], REAL b[]) {
	REAL c=0.0;
	for (int i=0;i<DIM;i++) c+=a[i]*b[i];
	return c;
}
int Geom::hash(int i, int j, int n) {
	return i>j?i*n+j:j*n+i;
}
REAL Geom::length(REAL a[]) {
	REAL d=0.0;
	for(int i=0;i<DIM;i++) {
		REAL r=a[i];
		d+=r*r;
	}
	return sqrt(d);
}
REAL Geom::area(REAL a[], REAL b[]) {
	REAL c[DIM];
	VECP(c,a,b);
	return LENGTH(c);
}
