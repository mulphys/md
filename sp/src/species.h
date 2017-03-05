/*! \namespace Species
 * \brief Description of Specie and Reaction classes
 */
namespace Species {
	extern int nspecies;
	class Specie {
		char id[WORDLENGTH];
		REAL mass,size,//molecular characteristics
			cp,//heat capacity
//-			dof,// =2*c_v/R=2*(c_p/R-1) http://en.wikipedia.org/wiki/Heat_capacity#Heat_capacity
			inertia[DIM];//tensor of inertia (geometrical param).
		public:
		Specie() {
			id[0]='\0';
			mass=size=cp=0.0;
			for(int i=0;i<DIM;i++) inertia[i]=0.0;
		}
//-		inline void DOF(REAL f) {dof=f;}
//-		inline REAL DOF() {return dof;}
		inline void Cp(REAL c) {cp=c;}
		inline REAL Cp() {return cp;}
		inline char* Id() {return id;}
		inline void Mass(REAL m) { mass=m; }
		inline REAL Mass() { return mass; }
		inline REAL Size() { return size; }
		inline void Size(REAL s) { size=s; }
	};
	//! Gas is a specie with some additional parameters
	struct Gas {
		REAL density; //!< density in kg/m^3
		Specie *specie;
		Gas(Specie *s) {specie=s;}
	};
	//! Reaction determines the two products only
	struct Reaction {
		//! Possible outcomes of a reaction
		struct Outcome {
			int product[2];//!< two indexes of the products of a reaction
			REAL activationEnergy,probability,enthalpy,
				time; //!< time of reaction: used primarily for surface reactions,
		           //!< where there is a time delay between species adsorption
		           //!< on the surface and product formation.
			Outcome *next;//!< linked list of outcomes
			Outcome() {
				for(int j=0;j<2;j++) product[j]=VOIDSPECIE;//VOID reaction
				probability=enthalpy=time=0.0;
				next=NULL;
			}
			Outcome(
				int ip0,//!< index of the first reaction product  
				int ip1,//!< index of the second reaction product 
				REAL a, //!< activation energy
				REAL p, //!< outcome probability
				REAL h  //!< enthalpy of reaction
			) {
				if(ip0<0||ip0>nspecies) {
					cout<<"ERROR: Invalid first product specie index "<<ip0<<" of maximum "<<nspecies<<endl;
					exit(1);
				}
				product[0]=ip0;
				if(ip1<0||ip1>nspecies) {
					cout<<"ERROR: Invalid second product specie index "<<ip1<<" of maximum "<<nspecies<<endl;
					exit(1);
				}
				product[1]=ip1;
				activationEnergy=a;
				probability=p;
				enthalpy=h;
				time=0.0;
				next=NULL;
			}
#ifdef ADSORPTION
			Outcome(
				int ip0, //!< index of the first reaction product  
				int ip1, //!< index of the second reaction product 
				REAL a, //!< activation energy
				REAL p, //!< outcome probability
				REAL h, //!< enthalpy of reaction
				REAL t  //!< delay time
			) {
				if(ip0<0||ip0>nspecies) {
					cout<<"ERROR: Invalid first product specie index "<<ip0<<" of maximum "<<nspecies<<endl;
					exit(1);
				}
				product[0]=ip0;
				if(ip1<0||ip1>nspecies) {
					cout<<"ERROR: Invalid second product specie index "<<ip1<<" of maximum "<<nspecies<<endl;
					exit(1);
				}
				product[1]=ip1;
				activationEnergy=a;
				probability=p;
				enthalpy=h;
				time=t;
				next=NULL;
			}
#endif
			void Products(int ip0, int ip1) {
				if(ip0<0||ip0>=nspecies) {
					cout<<"Product index "<<ip0<<" outside of species index\n";
					exit(1);
				}
				product[0]=ip0;
				if(ip1<0||ip1>=nspecies) {
					cout<<"Product index "<<ip1<<" outside of species index\n";
					exit(1);
				}
				product[1]=ip1;
			}
			int Product(int i) {
				if(i<0||i>1) {
					cout<<"ERROR: product index "<<i<<" outside of bounds\n";
					exit(1);
				}
				return product[i];
			}
			inline void Time(REAL t) {time=t;}
			inline REAL Time() {return time;}
			inline void ActivationEnergy(REAL a) { activationEnergy=a; }
			inline REAL ActivationEnergy() { return activationEnergy; }
			inline void Probability(REAL r) { probability=r; }
			inline REAL Probability() { return probability; }
			inline void Enthalpy(REAL h) { enthalpy=h; }
			inline REAL Enthalpy() { return enthalpy; }
		}	*outcomes,*current;//!< all outcomes, and the current outcome
		Reaction() {current=outcomes=NULL;}
		~Reaction() {
			Erase(outcomes);
			outcomes=NULL;
		}
		void Erase(Outcome *outcome) {
			if(outcome->next!=NULL) Erase(outcome->next);
			delete outcome;
		}
		void Add(int ip0, int ip1, REAL a, REAL p, REAL h) {
			//! Add reaction products, probability, and enthalpy:
			if(outcomes==NULL) {
				outcomes=new Outcome(ip0,ip1,a,p,h);
				return;
			}
			Outcome *outcome=outcomes;
			for(;outcome->next!=NULL;outcome=outcome->next);
			outcome->next=new Outcome(ip0,ip1,a,p,h);
		}
#ifdef ADSORPTION
		void Add(int ip0, int ip1, REAL a, REAL p, REAL h, REAL t) {
			//! Add reaction products, probability, enthalpy, and time:
			if(outcomes==NULL) {
				outcomes=new Outcome(ip0,ip1,a,p,h,t);
				return;
			}
			Outcome *outcome=outcomes;
			for(;outcome->next!=NULL;outcome=outcome->next);
			outcome->next=new Outcome(ip0,ip1,a,p,h,t);
		}
#endif
		Outcome *First() { current=outcomes; return current; } 
		Outcome *Next() { 
			current=current->next; 
			if(current!=NULL) return current; 
			current=outcomes; return NULL;
		} 
	};
	enum Interaction {
		missed=0,
		collided,
		reacted, //!< two reactants -> two products
		annihalated  //!< two reactants -> one product
	};
	extern Specie *species;
	extern Reaction *reactions;
}

