#define ABOUT "MOLECULAR DYNAMICS SIMULATOR BY A.V.SMIRNOV\nandrei.v.smirnov@gmail.com\n"
//#define DEBUG //uncomment this for debugging checks
#define DOCTYPE	"remody"
#define REAL	double 
#define TINY 1e-10
#define SMALL 1e-30
#define LARGE 1e30
#define MAXLINLEN	510
#define WORDLENGTH	64
#define ESC	0x1B
#define PI 3.14159265358979323846
#define SQRTPIo8	0.62666  //sqrt(pi/8)
#define ERROR(message) \
		{	fprintf(stderr,"\nERROR: %s\n",message);exit(1);}
#define TRUE	1
#define FALSE	0

// OpenMP 
#define OMP // comment this to use sequential mode (saves memory)
#define WAIT	for(int i=0;i<999999&&locked;i++)
const int maxcounters=4;// used to run multiple threads

// Numerical constants
const REAL 
	onethird=1.0/3.0,
	twothirds=2.0/3.0,
	Sqrt2=sqrt(2.0),
	Sqrt2Over2=0.5*Sqrt2,
	sqrt1p5=sqrt(1.5),
	sqrt1p25=sqrt(1.25),
	Sqrt3=sqrt(3.0),
	OneOverSqrt2=1.0/Sqrt2,
	OneOverSqrt3=1.0/Sqrt3,
	Sqrt3Over2=0.5*Sqrt3;

// Numerical Scheme Flags
#define LOCAL  //LOCAL INTERACTION SCHEME
#define COLLISIONTIME  //ESTIMATE THE TIME STEP AS MINIMUM TIME TO COLLISION 

#define RND	(REAL)((double)rand()/RAND_MAX)
// MACROS:
#define MAX(a,b)	((a)>(b)?a:b)
#ifndef ZERO
#define	ZERO(n,x)	for (int i=0; i<(n); i++)x[i]=0.0
#endif
#define ADD(n,x,y)	for (int i=0; i<(n); i++) (x)[i] += (y)[i]
#define SUB(n,x,y)	for (int i=0; i<(n); i++) (x)[i] -= (y)[i]
#define ADDC(n,x,y,d)	for (int i=0; i<(n); i++) (x)[i] += (y)[i]*(d)
#define SUBC(n,x,y,d)	for (int i=0; i<(n); i++) (x)[i] -= (y)[i]*(d)
#define MUL(n,x,y)	for (int i=0; i<(n); i++) (x)[i] *= (y)[i]
#define MULC(n,x,c)	for (int i=0; i<(n); i++) (x)[i] *= (c)
#define DIV(n,x,y)	for (int i=0; i<(n); i++) (x)[i] /= (y)[i]
#define MULVEC(x,s)	for (int i=0; i<DIM; i++) (x)[i] *= (s)
#define ADDVEC(x,y)	for (int i=0; i<DIM; i++) (x)[i] += (y)[i]
#define COPYVEC(x,y)	for (int i=0; i<DIM; i++) (x)[i] = (y)[i]
#define ZEROVEC(x)	for (int i=0; i<DIM; i++) (x)[i] = 0.
#define SCLP(A,B)	(A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
#define VECP(A,B,C)	A[0]=B[1]*C[2]-B[2]*C[1];\
							A[1]=B[2]*C[0]-B[0]*C[2];\
							A[2]=B[0]*C[1]-B[1]*C[0];
#define LENGTH(A)	sqrt(SCLP(A,A))
#define VOIDSPECIE	nspecies  //void specie index is the max spicie
#define GAUSS Run::gauss()

// Physical constants:
#define DIM	3
#define BoltzmannConstant	1.38066e-23 //!< \def Boltzman Constant in kg*(nm/ns)^2/K
#define AtomicMassUnit	1.66053886e-27 //!< \def atomic mass unit in kg/au
#define AvogadroNumber	6.02214179e23 //!< Avogadro Number in 1/mol
//#define GasConstant	8.31457070 //[70] J/(K*mol): ideal gas constant = BoltzmanConstant*AvogadroNumber
const REAL 
	GasConstant=BoltzmannConstant*AvogadroNumber,
	GasConstantInv=1.0/GasConstant;
// Compute degrees of freedom from specific heat Cp:
#define DOFK	3 // Kinetic degrees of freedom = 3
//#define DOF(cp)	MAX((cp)*GasConstantInv+0.5,DOFK) //= 3 + (2(Cp/R-1)-3)/2 = 3 + Cp/R -1 - 1.5 = Cp/R + 1.5 - 1 = Cp/R + 0.5
#define DOF(cp)	MAX(2.0*((cp)*GasConstantInv-1.0),DOFK) //= 3 + (2(Cp/R-1)-3)/2 = 3 + Cp/R -1 - 1.5 = Cp/R + 1.5 - 1 = Cp/R + 0.5
#define DOFI(dof)	(dof-DOFK) //Internal DOF from total DOF where 3.0 is kinetic dof

