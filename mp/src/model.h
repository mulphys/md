//MODELS:
#define HARDBALLS //!< hardball collision as opposed to interaction via a potential
/*!  Define the ADSORPTION flag below to enable 
 *   the time delay for the adsorbed species
 *   to be released from the surface 
 */
//#define ADSORPTION //!< turns on adsorbtion processing
enum AtomType {
	undefined=0,
	hydrogen=1,
	helium=2,
	carbon=12,
	oxygen=16,
	maxAtomTypes
};
union State {
//-	none=0,
//-	idle,
//-	busy,
	unsigned int locked: 1;
	unsigned int done: 1;
};
class Molecule {
	int type;//!< points to the species[type];
#ifdef OMP
	int id;//!< molecule id: context dependent:
	       //!< index (negative indicates locking)
	       //!< state: 0-free, 1-busy, 2-done
	State state;
#endif
#ifdef COLLISIONTIME
	REAL time,dt;// particle time and timestep
	Molecule *target;//collision partner
#endif
#ifdef LOCAL
	int	gridcell;//!< compressed 3D index to grid cell
	Ptr<Molecule> *ptr;//!< pointer to the grid cell slot 
#endif
	REAL 
		energy, //!<  internal energy 
#ifdef ADSORPTION
		residencetime,   //!<  release-time if adsorbed 
#endif
		x[DIM], //!<  position of center of mass
		v[DIM]; //!<  velocity of center of mass
//-	List<ATOM> atoms;
	public:
	Molecule();
	inline void Type(int t) {type=t;}
	inline int Type() {return type;}
#ifdef COLLISIONTIME
	inline REAL Time(){return time;}
	inline void Time(REAL t){time=t;}
	inline REAL TimeStep(){return dt;} 
	inline void TimeStep(REAL step){dt=step;}
	inline void Target(Molecule *a){target=a;}
	inline Molecule *Target(){return target;} 
#endif
#ifdef OMP
//-	inline bool islocked()   {return id<0;}
	inline bool islocked()   {return state.locked==1;}
//-	inline void lock()       {if(id>=0)id=-id-1;}
	inline void lock()       {state.locked=1;}
//-	inline void unlock()     {if(id<0)id=-id-1;}
	inline void unlock()     {state.locked=0;}
	inline bool isdone() {return state.done==1;}
	inline void done() {state.done==1;}
	inline void ready() {state.done=0;}
//-	inline State state()   {return id<0?(State)(-id):none; }
//-	inline void state(State state) {id=-((int)state); }
//-	inline int  index()      {return id<0?-id-1-nstates:id-nstates;}
	inline int  index()      {return id;}
//-	inline void index(int i) {id=i+nstates;}
	inline void index(int i) {id=i;}
#endif
#ifdef LOCAL
//-	inline int *GridCell() { return gridcell; }
	inline int GridCell() { return gridcell; }
	inline void GridCell(int icell) { gridcell=icell; }
	inline Ptr<Molecule> *getPtr() { return ptr; }
	inline void putPtr(Ptr<Molecule> *p) { ptr=p; }
#endif
//-	inline Pointer<Molecule> *getPointer() { return pointer; }
//-	inline void putPointer(Pointer<Molecule> *p) { pointer=p; }
	inline void InternalEnergy(REAL e) { energy=e; }
	inline REAL InternalEnergy() { return energy; }
#ifdef ADSORPTION
	/*! Release time after being adsorbed on the surface:
		NOTE: energy is overwrittent with 'time' when on the surface:
		The energy of reseased species will be determined by the surface temperature
	*/
	inline REAL ReleaseTime() {return residencetime;}
	inline void ReleaseTime(REAL t) {residencetime=t;}
#endif
	inline REAL Coordinate(int i) { 
		if(i>=DIM) {
			cerr<<"SubAtom coordinate index "<<i
			<<" exceeds "<<DIM<<endl;cerr.flush();
			return 0.0;
		}
		return x[i]; 
	}
	inline void Coordinate(int i, REAL y) { 
		if(i>=DIM) {
			cerr<<"SubAtom coordinate index "<<i
			<<" exceeds "<<DIM<<endl;cerr.flush();
		}
		x[i]=y; 
	}
	inline REAL* Coordinates() { return x; }
	inline void Coordinates(REAL y[]) { for (int i=0;i<DIM;i++) y[i]=x[i]; }
	inline void setCoordinates(REAL y[]) { for (int i=0;i<DIM;i++) x[i]=y[i]; }
	inline REAL Velocity(int i) { 
		if(i>=DIM) {
			cerr<<"SubAtom velocity index "<<i
			<<" exceeds "<<DIM<<endl;cerr.flush();
			return 0.0;
		}
		return v[i]; 
	}
	inline void Velocity(int i, REAL u) { 
		if(i>=DIM) {
			cerr<<"Molecule velocity index "<<i
			<<" exceeds "<<DIM<<endl;cerr.flush();
		}
		v[i]=u; 
	}
	inline REAL* Velocity() { return v; }
	inline void Velocity(REAL u[]) { for (int i=0;i<DIM;i++) u[i]=v[i]; }
	inline void setVelocity(REAL u[]) { for (int i=0;i<DIM;i++) v[i]=u[i]; }
	//! Retrieve kinetic energy
	REAL KineticEnergy();
	//! Set internal energy according to temperature:
	void Temperature(REAL temperature); 
	void Move();
	void copy(Molecule *molecule);
};
/*! \namespace Potential
 * \brief Interaction-potential functions, such as LJ 
 */
namespace Potential {
	const REAL small=1.0e-20, large=1.0e20;
	extern REAL cutoff;
	extern REAL sigma, eta;
	inline void strength(REAL strength) {
		eta=strength;
	}
	inline REAL strength() { return eta; }
	inline void lengthscale(REAL lengthscale) {
		sigma=lengthscale;
	}
	inline REAL lengthscale() { return sigma; }
	inline void Cutoff(REAL newcutoff) { cutoff=newcutoff; }
	inline REAL Cutoff() { return cutoff; }
	inline REAL value(REAL r) {
		if(r<small) r=small;
		// LJ with r-cutoff:
		REAL x=sigma/r,
		val=r>cutoff?0.0:eta*(pow(x,12.0) - pow(x,6.0));
//		val=r>cutoff?0.0:eta*(pow(x,16.0) - pow(x,8.0));
		return val<large?val:large;
//Octave check:
// X=[1.8:0.2:4.0];eta=1.0;sig=2.0;P=eta*((sig./X).^(16)-(sig./X).^8);plot(X,P);
//		return r<sigma?eta:0.0;
	}
	inline REAL force(REAL r) {
		return eta/sigma*(-12.0*pow(sigma/r,13.0) + 6.0*pow(sigma/r,7));
	}
	REAL value(REAL r);
	REAL derivative(REAL r);
//-	REAL potential(List<ATOM> *nodes);
};
/*! \namespace Geom
 * \brief Provides some trigeometry functions
 */
namespace Geom {
	void distvec(REAL a[], REAL b[], REAL c[]);
	REAL distance(REAL a[], REAL b[]);
	REAL distance(int dim, REAL a[], REAL b[]);
	REAL sclp(REAL a[], REAL b[]);
	REAL length(REAL a[]);
	REAL area(REAL a[], REAL b[]);
	REAL normalize(REAL a[]);
	int hash(int i, int j, int n);
}
