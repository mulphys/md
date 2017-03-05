//MODELS:
#define HARDBALLS //!< hardball collision as opposed to interaction via a potential
/*!  Define the ADSORPTION flag below to enable 
 *   the time delay for the adsorbed species
 *   to be released from the surface 
 */
//#define ADSORPTION //!< turns on adsorbtion processing
class Molecule;//to be declared later
enum AtomType {
	undefined=0,
	hydrogen=1,
	helium=2,
	carbon=12,
	oxygen=16,
	maxAtomTypes
};
//- DEPRECATED CLASS:
//-class ATOM {
//-	int id;// used for IO
//-	AtomType	type; 
//-	REAL x[DIM];
//-	public:
//-	ATOM();
//-	inline int Id(){ return id; }
//-	inline void setId(int newid){ id=newid; }
//-	inline AtomType Type(){ return type; }
//-	inline REAL Mass() { return (REAL) type; }  //Approximate mass of atom
//-	inline void setType(AtomType newtype){ type=newtype; } 
//-	inline void Coordinate(int i, REAL pos) {	
//-		if(i>=DIM) {
//-			cerr<<"SubAtom coordinate index "
//-			<<i<<" exceeds "<<DIM<<endl;cerr.flush();
//-			return;
//-		}
//-		x[i]=pos;
//-	}
//-	inline REAL Coordinate(int i)	{	
//-		if(i>=DIM) {
//-			cerr<<"SubAtom coordinate index "<<i
//-			<<" exceeds "<<DIM<<endl;cerr.flush();
//-			return 0.0;
//-		}
//-		return x[i];
//-	}
//-	inline REAL *Coordinates() { return x; }
//-	inline void Coordinates(REAL y[]) { for(int i=0;i<DIM;i++)y[i]=x[i]; }
//-	inline void setCoordinates(REAL y[]) { for(int i=0;i<DIM;i++)x[i]=y[i]; }
//-	inline REAL Time() { return x[DIM-1]; }
//-};
class Molecule {
	int type;//!< points to the species[type];
//-	int natoms;
#ifdef LOCAL
	int	gridcell;//!< compressed 3D index to grid cell
	Ptr<Molecule> *ptr;//!< pointer to the grid cell slot 
#endif
	REAL 
		energy, //!<  internal energy 
#ifdef ADSORPTION
		time,   //!<  release-time if adsorbed 
#endif
		x[DIM], //!<  position of center of mass
		v[DIM]; //!<  velocity of center of mass
//-	List<ATOM> atoms;
	public:
	Molecule();
	inline void Type(int t) {type=t;}
	inline int Type() {return type;}
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
	inline REAL ReleaseTime() {return time;}
	inline void ReleaseTime(REAL t) {time=t;}
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
