#define ABOUT "MOLECULAR DYNAMICS SIMULATOR BY A.V.SMIRNOV\nandrei.v.smirnov@gmail.com\n"
#define DEBUG //uncomment this for debugging checks
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
#define LOCAL  //LOCAL INTERACTION SCHEME
#define DIM	3
#define BoltzmannConstant	1.38066e-23 //!< \def Boltzman Constant in kg*(nm/ns)^2/K
#define AtomicMassUnit	1.66053886e-27 //!< \def atomic mass unit in kg/au
#define AvogadroNumber	6.02214179e23 //!< Avogadro Number in 1/mol
//#define GasConstant	8.31457070 //[70] J/(K*mol): ideal gas constant = BoltzmanConstant*AvogadroNumber
const REAL 
	GasConstant=BoltzmannConstant*AvogadroNumber,
	GasConstantInv=1.0/GasConstant;
#define RND	(REAL)((double)rand()/RAND_MAX)
// MACROS:
#define MAX(a,b)	a>b?a:b;
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
