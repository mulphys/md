/*
 * Particle interaction acceleration using a grid of rectangular cells
 */

#include <math.h>
#include <cstdlib>
#include <iostream>
using namespace std;
#include "def.h"
#include "run.h"
#include "list.h"
#include "collection.h"
#include "species.h"
#include "grid.h"
#include "container.h"
#include "model.h"

namespace GridContainer {// Keeps local node lists to avoid N^2 pair interaction loop
	REAL xmin[DIM],xmax[DIM];
	REAL cellsize;
	int ncells[DIM+1],// number of cells in each space direction: +1 to avoid if in the loop
		mcells=0,mcells1=0;//total sum of all ncells[]
	Container<Molecule> *nodes;//nodes[ncells[0]*ncells[1]*...*ncells[DIM-1]]
	Pool<Molecule> *pool;
	void init(
		REAL ymin[], REAL ymax[], 
		REAL radius, 
		Pool<Molecule> *newpool
	) {
		REAL margin=radius;
		pool=newpool;
		cellsize=radius;
		mcells=1;//total number of grid cells
		for(int i=0;i<DIM;i++) {
			REAL 
				x0=ymin[i]-margin,
				x1=ymax[i]+margin;
			int n=((int)((x1-x0)/cellsize))+1;
			ncells[i]=n;
			mcells*=n;
			xmin[i]=x0;
			xmax[i]=x1;
		}
		ncells[DIM]=1;//to avoid an if statement inside a loop
		nodes=new Container<Molecule>[mcells];
		if(nodes==NULL) {
			cout<<"Grid::init: Can't allocate nodes\n";cout.flush();
			exit(1);
		}
		mcells1=mcells/ncells[0];
		if(Run::option.debug||Run::option.verbose) {
			cout<<"Grid::init:cellsize:"<<cellsize<<"\n\txmin:";
			for(int i=0;i<DIM;i++) cout<<' '<<xmin[i];
			cout<<"\n\txmax:";
			for(int i=0;i<DIM;i++) cout<<' '<<xmax[i];
			cout<<"\n\tncells:";
			for(int i=0;i<DIM;i++) cout<<' '<<ncells[i];
			cout<<endl;cout.flush();
		}
	}
	int index(int ind[]) {//code index: ind->n
		int n=0;
		for(int i=0;i<DIM;i++) {
			n=n*ncells[i]+ind[i];
		}
		return n;
	}
	void index(int n, int ind[]) {//decode index: 
		if(n<0) {
			for(int i=0;i<DIM;i++) ind[i]=-1;
			return;
		}
		int m=mcells1;
		for(int i=0;i<DIM;i++) {
			ind[i]=n/m;
			n%=m;
			m/=ncells[i+1];
		}
	}
#ifdef LOCAL
	void put(Molecule *node) {
		REAL *x=node->Coordinates();//x[DIM]: node coordinates
//-		int *cellind=node->GridCell(),//cellind[DIM]: grid cell index
		int icell=node->GridCell(),//cellind[DIM]: grid cell index
			cellind[DIM],newind[DIM];//grid index
		bool moved=false;
		index(icell,cellind);
		for(int i=0;i<DIM;i++) {
			int j=(int)((x[i]-xmin[i])/cellsize);
			if(j>=ncells[i]||j<0) {
				cerr<<"Molecule ("<<x[0]<<','<<x[1]<<','<<x[2]<<") is outside of the grid: index: 0<"<<j<<"<="<<ncells[i]<<" erased\n";
				Ptr<Molecule> *pointer=node->getPtr();
				if(pointer!=NULL) {//Unlink the molecule from the grid
					nodes[icell].remove(pointer);
					node->putPtr(NULL);
				}
				node->Type(Species::VOIDSPECIE);//mark molecule to be erased
#ifdef DEBUG
				exit(1);///DDD
#else
				return; 
#endif
			}
			if(j!=cellind[i]) moved=true;
			newind[i]=j;
		}
		if(moved) {//reassign node to a different grid-cell
			Ptr<Molecule> *pointer=node->getPtr();//get list pointer 
			if(pointer!=NULL) {
				nodes[icell].unlink(pointer);
				icell=index(newind);
				nodes[icell].link(pointer);
			}	else {
				icell=index(newind);
				nodes[icell].insert(node);
				node->putPtr(nodes[icell].getPtr());
			}
			node->GridCell(icell);
		}
	}
#endif //END LOCAL INTERACTIONS
	Container<Molecule> *get(int ind[]) { 
		return nodes+index(ind);
	}
	int *dimensions() { return ncells; }
}

