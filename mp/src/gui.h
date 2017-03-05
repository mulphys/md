/*
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 */


#define	WINDOW_SIZE	800
#define	MAX_WINDOW_WIDTH	800
#define	MAX_WINDOW_HEIGHT	600

#define	LEFT_BUTTON	1
#define	MIDDLE_BUTTON	2
#define	RIGHT_BUTTON	3
/*! \namespace Gui
 * \brief OpenGL window output routines
 */
namespace	Gui {
	enum	Elements {
		points=0,
		nodes,
		edges,
		faces,
		cells,
		frame,
		mesh,
		boundary_nodes,
		boundary_edges,
		boundary_faces,
		boundary_cells,
		maxelements
	};
	enum	Color {
		red=0,
		green,
		blue,
		skyblue,
		brown,
		magenta,
		yellow,
		color6,
		color7,
		color8,
		color9,
		color10,
		color11,
		color12,
		color13,
		color14,
		color15,
		maxcolor
	};
	enum ColorSchemes {
		colorByType=0,
		colorByBoundary,
		colorByTime,
		colorByMass,
		maxColorSchemes
	};
	enum	PointParameters
	{
		variable=maxcolor,
		vmin,vmax,size,
		sizevar,massvar,
		maxpointprm
	};
	enum	LineParameters
	{
		thickness=maxcolor,
		length,
		maxlineprm
	};
	enum	AxesParameters
	{
		axeswidth=0,
		xaxislength,
		yaxislength,
		zaxislength,
		arrowheight,
		arrowidth,
		maxaxesprm
	};
	enum	SurfDispType
	{
		gridlines=0,
		solidsurface,
		maxsurfdisptypes
	};
	extern char *configfile;
	extern REAL
		xmin[],xmax[],
		lx,ly,lz,lmin,lmax;
	extern REAL
		step,
		vecval[maxlineprm],
		axes[maxaxesprm],
		rgbcolor[maxcolor][3];

	struct	ElementDisp
	{	int
			lighting;
		REAL
			rgbcolor[3];//0-red, 1-green, 2-blue
	};	
	extern ElementDisp disp[maxelements];
	struct	ColorScale {
		int	relative	:1;
	};
	struct	Scene {
		struct Color {
			ColorSchemes scheme; // 0: by boundary, 1: by time
			REAL minvalue,maxvalue;
		}	color;
		struct Mesh {
			REAL
				node[maxpointprm],
				line[maxlineprm];
		}	mesh;
		struct Frame {
			REAL
				line[maxlineprm];
		}	frame;
	};
//	public:
	enum Movement
	{	stay=0,
		rotate,
		moveuv,
		movew
	};
	enum ShowPars
	{	showRun=0,
		showAxes,
		showSpheres,
		showBonds,
		showGrid,
		showNodes,
		showVariables,
		showBoundaryVertexes,
		showBoundaryVectors,
		showToolVertexes,
		showBoundaryFaceCenters,
		showBoundaryFaces,
		showBoundaryGrid,
		showToolGrid,
		showFrame,
		showCellCenters,
		showFaceCenters,
		showIsoSurfaces,
		dumpWindow,
		maxshowpars
	};
	extern int
		showpar[maxshowpars],
		finished,
		animation;
	extern REAL
		xo[], //origin of coordinate system
		dx,dy,dz,
		lastx, lasty,	lastz;
	
	// Added
	extern int mouseButtons[3];
	extern REAL zoom;
	extern REAL rotx;
	extern REAL roty;
	extern REAL tx;
	extern REAL ty;
	
	
	struct	WindowGeom
	{
		int
			width,height;
	};
	void	Exit();
	void	Quit();
	void	readconf();
	void	helpDisplay();
	int	readparam (
		char	*s,
		char	*param[],
		int	maxparam,
		char	*filename,
		REAL	val[]
	);
	void	initdisp();
	void	Materials(int argc, char *argv[]);
	void	display();
	void	displayAxes();
	void	writeData();
	void	helpCommand();
	void	printGridXLimits();
	void	printGridVecLimits();
	void	getXLimits(REAL	*xmin, REAL *xmax);
	void	getScaling(REAL	&lmin, REAL &lmax);
	void	printPartVelLimits();
	void	printParticleVariables();
	void	showGridElements
	(	
		int ne, 
		double	*X,//coordinates
		REAL color[]
	);
	void	showVectorComponent
	(
		int	ivar,
		int	icomp,
		double	*X,//coordinates
		double	vmn,
		double	vmx
	);
	void	showVector
	(
		int	ivar,
		double	*X//coordinates
	);
	void	commandMode();
	void	printVariables
	(
		int	n
//		struct Variables	*V
	);

void	init(int argc, char *argv[], Domain *domain);
void	finish();
void	helpDisplay();
void	readconf();
void	initdisp();
void	InitMaterials(void);
/*  Initialize material property and depth buffer.
 */
void	showGridElements
(	int	ne,//number of elements to show
	double	*X,//coordinates
	REAL	color[]
);
void	showVector
(
	int	ivar,
	double	*X//coordinates
);
void	showVectorComponent
(
	int	ivar,
	int	icomp,
	double	*X,//coordinates
	double	vmn,
	double	vmx
);
void	displayAxes();
/* ARGSUSED1 */

void	getScaling
(
	double	&l0,
	double	&l1
);
void	getXLimits
(	double	*x0,
	double	*x1
);
void	printPartVelLimits();
void	printVariables
(
	int	n
);
void	Exit();
void	Quit();
void	selectVariable();
void	printParticleVariables();
void	setBackgroundRun();
void	setForegroundRun();
void	commandMode();

/* Main */

/******** GLOBAL GUI-FUNCTIONS **********/

void	display();
void	animate();
void	writeData();
void	printGridXLimits();
void	printGridVecLimits();
void	helpCommand();
void	reshape(int w, int h);
void	mouse(int button, int state, int x, int y);
void	motion(int x, int y);
/* ARGSUSED3 */
void	keyboard(unsigned char key, int x, int y);
void	menu(int value);
void	run();
void	drawSegment
(
	REAL	*x,
	REAL	*y 
);
void	displaymessage(char *);
void	setElementColor(int element_type, REAL rgbcolor[]);
void	getElementColor(int element_type, REAL rgbcolor[]);
} //END NAMESPACE
/** @namespace Palette
 * \brief Defines color palette
 *
 */
namespace Palette {
//-			alow =PI/(2.*(vav-vmn)),
//-			ahigh=PI/(2.*(vmx-vav));
	extern void init(REAL vmin, REAL vmax);
	extern void pickcolor(REAL var, REAL color[]);
}

