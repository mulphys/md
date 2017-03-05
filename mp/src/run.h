/*! \brief Execution options.
 */
struct	Option {
	unsigned int	restart	:1;//< new run (0) or restrart from inputfile (1)
	unsigned int	verbose	:1;//< show as much output as possible
	unsigned int	debug	:1;//!< run in debugging mode
	unsigned int	mesh	:1;//!< used to start the mesh generator
	unsigned int	tags	:1;//!< used in the trim module to trigger tag assignement
	unsigned int	xterm	:1;//!< used in the trim module to trigger tag assignement
};
/*! \namespace Run
 * \brief Processing command-line and execution control parameters
 */
namespace Run {
	extern struct Option	option;//!< command-line options
	extern char programname[];
	extern char configfile[];
	extern char inputfile[];
	extern char outputfile[];
	extern char outputname[];
#ifdef OMP
	extern int nthreads;
#endif
	//! Execution time-control structure
	struct Time
	{	double
			start,	//!< start-time
			end,	//!< end-time
			prev,	//!< previous time
			current,	//!< current time
			step,	//!< time-step
#ifdef COLLISIONTIME
			next, //!< next time
#else
			nextstep,	//!< next time-step
#endif
			step0,	//!< previous time-step
			output;	//!< data dump time inteval
	};
	extern struct Time time;
	void	init(int argc, char *argv[]);
	void	readcmdline(int argc, char *argv[]);
	int elapsed();//!< time (sec) elapsed since the start of the run
	void seed();//!< seed random number generator
	REAL rnd();//!< random number generator
	REAL gauss();//!< normal random number generator
	extern void (*usage)();
}
