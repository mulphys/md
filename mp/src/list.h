/*
 * List is a loop-linked list of objects.
 * When inserting/deleting objects list uses
 * system memory allocation deallocation
 * functions.
 * Advantages:
 * The size of list is not specified at the beginning and 
 * can grow to fill up the available system memory.
 * Drawbacks:
 * System calls to new/delete functions are slow.
 * For more efficient allocation/deallocation use
 * Collection and Container classes which use their 
 * own allocation/deallocation routines from a fixed i
 * size memory pool (see collection.*).
*/

// CLASSES:
template <class Element>
struct Pointer
{	Element	*element;
	Pointer	*next,*prev;
};
template <class Element>
class List {
	int	nelements;//number of elements in the list
	struct Pointer<Element> *first,*last,*current;
	public:
	List();
	List(int n);
	~List();
	inline int number(){ return nelements; }
	inline void setFirst(){ first=last->next; }
	inline void setFirstLast(){ first=last; }
	inline void goFirst(){ current=first; }
	inline void goFirstNext(){ current=first->next; }
	inline void FirstNext(){ first=first->next; }
	inline void FirstPrev(){ first=first->prev; }
	inline bool isFirstLast(){ if(first==last) return true; else return false; }
	inline void goLastNext(){ current=last->next; }
	inline void goLast(){ current=last; }
	inline void goNext(){ current=current->next; }
	inline void goPrev(){ current=current->prev; }
	inline bool isFirst(){ return current==first?true:false; }
	inline bool isLast(){ if(current==last) return true; else return false; }
	inline void	LastNext(){ last=last->next; }
	inline void	LastPrev(){ last=last->prev; }
	inline void go(Pointer<Element> *p) { current=p; }//not safe
	inline Element *First(){ return first->element; }
	inline Element *Last(){ return last->element; }
	inline Element *Next(){ current=current->next; return current->element; }
	inline Element *getNext(){ return current->next->element; }
	inline Element *Prev(){ current=current->prev; return current->element; }
	inline Element *getPrev(){ return current->prev->element; }
	inline Element *Current(){ return current->element; }
	inline Pointer<Element> *getPointer() { return current; }
	inline Pointer<Element> *getFirstPointer() { return first; }
	inline Pointer<Element> *getLastPointer() { return last; }
//-	void init(); //DDD: the same as List();
	int init(int n); //initialize with n elements
	int insert(Element *element);
	int insert();
	int append(Element *element);
	int append();
	int link(Pointer<Element> *p);
	void unlink(Pointer<Element> *p);
	int prepend(Element *element);
	bool locate(Element *element);
	bool locate(Pointer<Element> *pointer);
	void moveAfterFirst();
	void swapAfterFirst();
	void moveBehindFirst();
	int erase();
	bool erase(Element *element);
	void eraseall();
	bool remove();
	void cleanup();
	void remove(Pointer<Element> *pointer);
};
template <class Element>
List<Element>::List()
{	nelements=0;
	current=first=last=NULL;
}
//-template <class Element>
//-void List<Element>::init()///DDD
//-{	nelements=0;
//-	current=first=last=NULL;
//-}
template <class Element>
List<Element>::~List() {
	cleanup();
}
template <class Element>
void List<Element>::eraseall() {
	while (erase());
	current=first=last=NULL;
	nelements=0;
}
template <class Element>
void List<Element>::cleanup() {
	while (remove());
	current=first=last=NULL;
	nelements=0;
}
template <class Element>
int	List<Element>::init(int n)
{	//initializes list with n elements
	for(int i=0;i<n;i++) if(append()==0) return 0;
	nelements=n;
	return n;
}
template <class Element>
int	List<Element>::prepend(Element *element)
{	//appends before last
	if(last==NULL) {
		last=new Pointer<Element>;
		last->next=last->prev=last;
		last->element=element;
		current=first=last;
		nelements=1;
		return 1;
	}
	Pointer<Element> *tmp=last->prev;
	last->prev=new Pointer<Element>;
	if(last->prev==NULL) return 0;//failed to append: the memory is full
	tmp->next=last->prev;
	last->prev->next=last;
	last=last->prev;
	last->prev=tmp;
	last->element=element;
	nelements++;
	first=last->next;
	return 1;
}
template <class Element>
int	List<Element>::append(Element *element)
{	//appends after last
	if(last==NULL) {
		last=new Pointer<Element>;
		last->next=last->prev=last;
		last->element=element;
		current=first=last;
		nelements=1;
		return 1;
	}
	Pointer<Element> *lastnext=last->next;
	last->next=new Pointer<Element>;
	if(last->next==NULL) return 0;//failed to append: the memory is full
	lastnext->prev=last->next;
	last->next->prev=last;
	last=last->next;
	last->next=lastnext;
	last->element=element;
	nelements++;
	first=lastnext;
	return 1;
}
template <class Element>
int	List<Element>::append()
{	//appends before last
	Element *element=new Element;
	return append(element);
}
template <class Element>
int	List<Element>::link(Pointer<Element> *pointer)
{	//appends after last
	if(pointer==NULL) return 0;
	if(last==NULL) {
		last=pointer;
		last->next=last->prev=last;
		current=first=last;
		nelements=1;
		return 1;
	}
	Pointer<Element> *lastnext=last->next;
	last->next=pointer;
	lastnext->prev=last->next;
	last->next->prev=last;
	last=last->next;
	last->next=lastnext;
	nelements++;
//-	first=lastnext;
	return 1;
}
template <class Element>
int	List<Element>::insert(Element *element)
{	//inserts after current
	if(last==NULL) {
		last=new Pointer<Element>;
		last->next=last->prev=last;
		last->element=element;
		current=first=last;
		nelements=1;
		return 1;
	}
	Pointer<Element> *tmp=current->next;
	current->next=new Pointer<Element>;
	if(current->next==NULL) return 0;//failed to insert: the memory is full
	tmp->prev=current->next;
	current->next->prev=current;
	current=current->next;
	current->next=tmp;
	current->element=element;
	last=current;
	nelements++;
	return 1;
}
template <class Element>
int	List<Element>::insert()
{	//inserts after current
	Element *element=new Element;
	return insert(element);
}
template <class Element>
bool	List<Element>::remove() {	// deletes current, but does not delete the element
	if(current==NULL) return false;//failed to remove: the list is empty
	if(current->next==current) {
		first=last=NULL;
		delete current;
		current=NULL;
		nelements=0;
		return true;
	}
	current->next->prev=current->prev;
	current->prev->next=current->next;
	if(current==last)last=last->prev;
	if(current==first)first=first->next;
	Pointer<Element> *old=current;
	current=current->prev;
	delete old;
	nelements--;
	return true;
}
template <class Element>
int	List<Element>::erase() {	// deletes current
	if(current==NULL) return 0;//failed to erase: the list is empty
	if(current->next==current)
	{	delete current->element;
		first=last=NULL;
		delete current;
		current=NULL;
		nelements=0;
		return 1;
	}
	current->next->prev=current->prev;
	current->prev->next=current->next;
	if(last==current)last=last->prev;
	Pointer<Element> *old=current;
	current=current->prev;
	delete old->element;
	delete old;
	nelements--;
	first=last->next;
	return 1;
}
template <class Element>
void	List<Element>::unlink(Pointer<Element> *pointer) {	// deletes current
	if(last==NULL||pointer==NULL) {
//-		cerr<<"Failed to remove from the list: the list is empty\n";
		return;
	}
//-	if(locate(pointer)) {
//-		remove();///DDD: slow but safe
	if(pointer==current) current=current->next;
	if(pointer==last)last=last->prev;
	if(pointer==first)first=first->next;
	pointer->next->prev=pointer->prev;
	pointer->prev->next=pointer->next;
	if(--nelements==0){ first=NULL; last=NULL;}
}
template <class Element>
void	List<Element>::remove(Pointer<Element> *pointer) {	// deletes current
	if(current==NULL||pointer==NULL) {
//-		cerr<<"Failed to remove from the list: the list is empty\n";
		return;
	}
	if(locate(pointer)) remove();///DDD: slow but safe
	else {
		cerr<<"Can't locate pointer in the list\n";
		exit(1);
	}///DDD

//+	current=pointer; //fast but unsafe
//+	remove();
}
template <class Element>
bool	List<Element>::erase(Element *element) {
	if(!locate(element))return false;
	erase();
	return true;
}
template <class Element>
bool List<Element>::locate(Element *element) {
	if(current==NULL) return false;
	current=first;
	do {
		if(current->element==element) return true;
		current=current->next;
	}	while(current!=first);
	return false;
}
template <class Element>
bool List<Element>::locate(Pointer<Element> *pointer) {
	if(current==NULL) return false;
	current=first;
	do {
		if(current==pointer) return true;
		current=current->next;
	}	while(current!=first);
	return false;
}
template <class Element>
void List<Element>::swapAfterFirst() {
	Element *tmp=current->element;
	current->element=first->next->element;
	first->next->element=tmp;
	current=first->next;
}
template <class Element>
void List<Element>::moveAfterFirst() {
	Pointer<Element> *firstnext=first->next;
	if(current==firstnext)return;
	first->next=current;
	firstnext->prev=current;
	current->next->prev=current->prev;
	current->prev->next=current->next;
	current->prev=first;
	current->next=firstnext;
}
template <class Element>
void List<Element>::moveBehindFirst() {
	Pointer<Element> *firstprev=first->prev;
	if(current==firstprev)return;
	first->prev=current;
	firstprev->next=current;
	current->next->prev=current->prev;
	current->prev->next=current->next;
	current->next=first;
	current->prev=firstprev;
}
//-DEPRECATED:
//-class Molecule;
//-namespace GridList {// Keeps local node lists to avoid N^2 pair interaction loop
//-	extern REAL cellsize;
//-	extern int ncells[DIM+1];// number of cells in each space direction
//-	extern List<Molecule> *nodes;//nodes[ncells[0]*ncells[1]*...*ncells[DIM-1]]
//-	extern int *dimensions();
//-	extern int index(int gridcell[]);
//-	extern void index(int icell, int gridcell[]);
//-	extern void init(REAL ymin[], REAL ymax[], REAL radius);
//-	extern void put(Molecule *node);
//-	extern List<Molecule> *get(int index[]); 
//-};
