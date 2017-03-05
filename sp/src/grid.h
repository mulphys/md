class Molecule;
/*! \namespace GridContainer
 * \brief Domain Segmentation for Interaction Acceleration 
 *
 * Splits space into cubic cells to consider local interactions
 * between ajacent cells only
 */
namespace GridContainer {// Keeps local node lists to avoid N^2 pair interaction loop
	extern REAL cellsize;
	extern int ncells[DIM+1];// number of cells in each space direction
	extern Container<Molecule> *nodes;//nodes[ncells[0]*ncells[1]*...*ncells[DIM-1]]
	extern Pool<Molecule> *pool;
	extern int *dimensions();
	extern int index(int gridcell[]);
	extern void index(int icell, int gridcell[]);
	extern void init(REAL ymin[], REAL ymax[], REAL radius, Pool<Molecule> *pool);
	extern void put(Molecule *node);
	extern Container<Molecule> *get(int index[]); 
	extern bool checkPool(int icell, char *msg);
};

