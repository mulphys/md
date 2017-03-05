/*! \file domain.h
 *	Main definitions of the Domain class and its
 *	associate class Boundary
 */
enum BoundaryTypes {
	insideBoundary=0,
	periodicBoundary,
	elasticBoundary,
	openBoundary,
	maxBoundaryTypes
};
/*! \class Boundary
 * \brief Domain Boundary Definition
 *
 */
struct Boundary {
	int idir,//!< boundary orientation in 3D space: 0,1,2
		iside;//!< side of the boundary: 0,1
	REAL *xmin,*xmax;//!< min/max bounds
	BoundaryTypes	type; //!< defines boundary type
	bool adiabatic;//!< flag for adiabatic boundary
	REAL	
		area, //!< boundary surface area
		temperature;//!< boundary temperature
	List<Gas> gases;//!< list of gases across the boundary
	Reaction *reactions;//!< boundary reactions with all species
	//! Boundary constructor: sets defaults
	Boundary() {
		idir=iside=-1;
		xmin=xmax=NULL;
		type=elasticBoundary;
		adiabatic=true;
		temperature=0.0;
		area=0.0;
		reactions=NULL;
	}
	//! Boundary destructor
	~Boundary() {
//+		delete reactions;
	}
	inline REAL Area() {return area;}
	inline void Area(REAL a){ area=a; }
	inline REAL Temperature() {return temperature;}
	inline void Temperature(REAL t){ temperature=t; }
	void init(int i, int j, REAL ymin[], REAL ymax[], int nspecies);
	void Inject(
		int specietype, //!< injected specie type
		REAL vavx,//!< average velocity in on direction
		REAL velx,//!< actual velocity normal to the boundary
		REAL temperature, //!< temperature
		Collection<Molecule> *molecules //!< collection of molecules
	);
};
/*! \class Bulk
 * \brief Domain Bulk Properties
 *
 * Holds thermodynamics properties and the list of species
 * in the bulk of the domain. Initialized only when the new
 * run. If Run::option.restart=1, it is ingored and the bulk
 * properties and species are retrieved from the input file.
 */
class Bulk {
	REAL temperature,volume,
		xmin[DIM],xmax[DIM];//!< box-dimensions
public:
	List<Gas> gases;//!< list of gases initially present in the bulk
	inline REAL Temperature() {return temperature; }
	inline void Temperature(REAL temp) {temperature=temp;}
	inline REAL Volume() {return volume; }
	inline void Volume(REAL vol) {volume=vol;}
	void init(REAL xmin[], REAL xmax[]);
	//! Adding a molecule to the pool with not registering on the grid
	Molecule *inject(int specie, REAL velocity, Collection<Molecule> *molecules);
	//! Additn a molecule both to the pool and to the grid
};
/*! \class Domain
 * \brief Computational Domain 
 *
 * Implements main data types inside a computational domain
 *
 */
class Domain {
	int *distribution;//!< molecule counts for all species
	REAL xmin[DIM],xmax[DIM];//!< domain bounds
	char boundaryName[maxBoundaryTypes][WORDLENGTH];//!< to be depricated
	Bulk bulk;//!< bulk of the domain: initialized only if Run::option.restart=0
	Boundary boundaries[DIM][2]; //!< six boundaries of a hexahedral domain
	Collection<Molecule> molecules; //!< Collection of molecules
	Pool<Molecule> *pool;//!< Molecule pointer pool 
	public:
	Domain(char *filename);
	~Domain();
	inline void BoundaryType(enum BoundaryTypes b, int idir, int inside) {
		if(idir<0||idir>=DIM||inside<0||inside>2) {
			fprintf(stderr,"ERROR: Invalid boundary index: %d, %d in BoundaryType\n",idir,inside);
			exit(1);
		} 
		boundaries[idir][inside].type=b; 
	};
	inline enum BoundaryTypes BoundaryType(int idir, int inside) {
		if(idir<0||idir>=DIM||inside<0||inside>2) {
			fprintf(stderr,"ERROR: Invelid boundary index: %d, %d in BoundaryType\n",idir,inside);
			exit(1);
		} 
		return boundaries[idir][inside].type; 
	};
	inline void setMinBound(int i, REAL x){ xmin[i]=x;};
	inline void setMaxBound(int i, REAL x){ xmax[i]=x;};
	inline REAL minBound(int i){ return xmin[i];};
	inline REAL maxBound(int i){ return xmax[i];};
//-	inline int getBounds(REAL *x0, REAL *x1) { x0=xmin; x1=xmax;};
	int computeBounds(REAL x0[], REAL x1[]);
//?	inline List<Molecule> *Molecules(){ return &molecules; };
	inline Collection<Molecule> *Molecules(){ return &molecules; };
//-	void load(char *filename, List<Molecule> *molecules, REAL xmin[], REAL xmax[]);
	int run(int niter);
	enum BoundaryTypes boundary(Molecule *molecule);
	/*! Inject cross-boundary species.
    * The algorithm uses the density and temperature to caclulate the
    * frequency of injection of molecules of specie, \f$s\f$ at the boundary.
    * The injection frequency, \f$f_s\f$ per unit area, \f$A\f$ is computed as
    * \f[
    * f_s = \frac{N}{\Delta\,t}
    * \f]
    * where \f$\Delta\,t_s\f$ is the time interval at which 
    * a molecule hits the boundary area \f$A\f$, and
    * \f$N\f$ is the number of molecules in a volume with the base \f$A\f$
    * and length \f$\Delta\,x\f$ as shown in the figure:
    * \image html 'injectionfrequency.png'
    * The time inteval, \f$\Delta\,t\f$, between the collisions
    * can be related to the component of velocity of the molecule of species \f$s\f$,
    * in direction \f$x\f$, \f$v_{sx}\f$, as follows:
    * \f[
    * \Delta\,t = \frac{2\Delta\,x}{v_{sx}}
    * \f]
    * The number of molecules, \f$N\f$ can be related to density, \f$\rho\f$, as:
    * \f[
    * N = \frac{\rho}{\mu}\Delta\,x\,A
    * \f]
    * where \f$\mu\f$ is the mass of one molecule and \f$A = \Delta\,y\,\Delta\,z\f$. 
    * Thus, the frequency is:
    * \f[
    * f_s = \frac{\rho}{\mu}\Delta\,x\,A\frac{v_{sx}}{2\Delta\,x}
    *     = \frac{\rho\,v_{sx}}{2\mu}A
    * \f]
    * 
    */
	void injection();
	void interaction();
	//! Computes interaction between two molecules with the exchange of 
	// energy and momentum according to hard-ball collisions, and
	// the possibility of chemical reactions. 
	Interaction interact(Molecule *a, Molecule *b);
/*!< Initializing molecules according to bulk properties (restart=0)
 * The velocity of each molecule is set from given temperature according to:
 * \f[ \frac{3}{2}kT = \frac{mV_{av}^2}{2} 
 * \f]
 * Since the average velocity relates to the average velocity in one direction, x, as:
 * \f[ V_{av}^2 = 3V_{x,av}^2 
 * \f]
 * we obtain for \f$V_{av,x} = kT\f$
 */
	void init();
	void load(char *filename);//!< Initializing molecules according to bulk properties (restart=0)
	void save(char *taskname);
};

