/*
 * This file implements classes Collection and Container
 * Collection is used to keep the dynamic particle list
 * of variable size where the particles (molecules) can be
 * added/deleted without system memory allocation calls.
 * Collection is a fixed size dynamic list.
 * Container is used to arrange the dynamic grid-cell-to-particle
 * pointer arrays (grid.*) used to accelerate particle interaction 
 * scheme.
 * As in collection the Container allocates a fixed size pointer
 * array (pool) and then retrieves/discards pointers to the pool
 * without making expensive system memory allocation/deallocation
 * calls.
 *
*/


// CLASSES:
template <class Element>
struct Item {//The Element is imbedded	
	Element	element;
	Item	*next,*prev;
};
template <class Element>
class Collection {//Uses Items
	int mitems,nitems;//maximum and current number of items
	bool locked;
	struct Item<Element>
		*items,//origin of the items array
		*dead,// dead items.
		*first,*last,*current,
		*current1[maxcounters];
	public:
	Collection();
	~Collection();
	inline int maxCounters(){return maxcounters;}
	inline int number(){return nitems;}
	inline int maxnumber(){return mitems;}
	inline void setFirst(){ first=last->next; }
	inline void setFirstLast(){ first=last; }
	inline void goFirst(){ 
		current=first; 
#ifdef OMP
		for(int i=0;i<maxcounters;i++)current1[i]=first;
#endif
	}
	inline void goFirst(int i){ current1[i]=first; }
	inline void goFirstNext(){ current=first->next; }
	inline void goFirstNext(int i){ current1[i]=first->next; }
	inline void FirstNext(){ first=first->next; }
	inline void FirstPrev(){ first=first->prev; }
	inline bool isFirstLast(){ if(first==last) return true; else return false; }
	inline void goLastNext(){ current=last->next; }
	inline void goLastNext(int i){ current1[i]=last->next; }
	inline void goLast(){ current=last; }
	inline void goLast(int i){ current1[i]=last; }
	inline void goNext(){ current=current->next; }
	inline void goNext(int i){ current1[i]=current1[i]->next; }
	inline void goPrev(){ current=current->prev; }
	inline void goPrev(int i){ current1[i]=current1[i]->prev; }
	inline bool isFirst(){ return current==first?true:false; }
	inline bool isFirst(int i){ return current1[i]==first?true:false; }
	inline bool isLast(){ if(current==last) return true; else return false; }
	inline bool isLast(int i){ if(current1[i]==last) return true; else return false; }
	inline void	LastNext(){ last=last->next; }
	inline void	LastPrev(){ last=last->prev; }
	inline void go(Item<Element> *p) { current=p; }//not safe
	inline Element *First(){ return &first->element; }
	inline Element *Last(){ return &last->element; }
	inline Element *Next(){ current=current->next; return &current->element; }
	inline Element *Next(int i){ current1[i]=current1[i]->next; return &current1[i]->element; }
	inline Element *getNext(){ return &current->next->element; }
	inline Element *getNext(int i){ return &current1[i]->next->element; }
	inline Element *Prev(){ current=current->prev; return &current->element; }
	inline Element *Prev(int i){ current1[i]=current1[i]->prev; return &current1[i]->element; }
	inline Element *getPrev(){ return &current->prev->element; }
	inline Element *getPrev(int i){ return &current1[i]->prev->element; }
	inline Element *Current(){ return &current->element; }
	inline Element *Current(int i){ return &current1[i]->element; }
	inline Item<Element> *getItem() { return current; }
	inline Item<Element> *getItem(int i) { return current1[i]; }
	inline Item<Element> *getFirstItem() { return first; }
	inline Item<Element> *getLastItem() { return last; }
	bool init(int n);
	Element *insertAfter();
	Element *insert0();
	Element *insert();
	Element *insert(int icounter);
	bool insert0(Element *element);
	bool insert(Element *element);
	bool insert(Element *element, int icounter);
	Element *append0();
	Element *append();
	Element *append(int icounter);
	bool append0(Element *element);
	bool append(Element *element);
	bool append(Element *element, int icounter);
	bool remove0();
	bool remove();
	bool remove(int icounter);
};
template <class Element>
Collection<Element>::Collection() {
	mitems=nitems=0;
	locked=false;
	dead=first=last=current=NULL;
	for(int i=0;i<maxcounters;i++)current1[i]=NULL;
}
template <class Element>
Collection<Element>::~Collection() {
//+	delete items;
	mitems=nitems=0;
	dead=first=last=current=NULL;
	for(int i=0;i<maxcounters;i++)current1[i]=NULL;
}
template <class Element>
bool Collection<Element>::init(int n) {
	dead=new Item<Element>[n];
	if(dead==NULL) {
		cout<<"CAN'T ALLOCATE "<<n<<" ITEMS\n";cout.flush();
		exit(1);
	}
	for(int i=1;i<n-1;i++) {
		current=dead+i;
		current->next=dead+i+1;
		current->prev=dead+i-1;
		for(int i=0;i<maxcounters;i++) {
			current1[i]=dead+i;
			current1[i]->next=dead+i+1;
			current1[i]->prev=dead+i-1;
		}
	}
	dead->next=dead+1;
	dead->prev=dead+n-1;
	dead[n-1].next=dead;
	dead[n-1].prev=dead+n-2;
	nitems=0;mitems=n;
	first=last=current=dead;
#ifdef OMP
	for(int i=0;i<maxcounters;i++) current1[i]=dead;
#endif
}
template <class Element>
Element	*Collection<Element>::insertAfter() {
	//Insert after current:
	if(dead->next==dead) {//keep the last dummy dead item to save on if's
		cout<<"CAN'T INSERT AFTER THE ELEMENT: Collection OF "<<nitems<<" IS FULL\n";
		cout.flush();
		return NULL;
	}
	Item<Element> *newitem=dead;
	if(first!=dead) {
		//Resurrect one:
		dead=dead->next;
		dead->prev=newitem->prev;
		dead->prev->next=dead;
		//Insert behind of current:
		Item<Element> *tmp=current->prev;
		current->prev=newitem;
		tmp->next=newitem;
		newitem->next=current;
		newitem->prev=tmp;
//?		current=newitem;
	}	else {//empty collection: initialize first:
		dead=dead->next;
		dead->prev=first->prev;
		dead->prev->next=dead;
		first->next=first->prev=first;
		newitem=last=current=first;
#ifdef OMP
		for(int i=0;i<maxcounters;i++)current1[i]=current;
#endif
	}
//? last=current; 
	nitems++;
	return &newitem->element;
}
template <class Element>
Element	*Collection<Element>::insert0() {
	//Insert before current:
	if(dead->next==dead) {//keep the last dummy dead item to save on if's
		cout<<"CAN'T INSERT ELEMENT: Collection OF "<<nitems<<" IS FULL\n";
		cout.flush();
		return NULL;
	}
	Item<Element> *newitem=dead;
	if(first!=dead) {
		//Resurrect one:
		dead=dead->next;
		dead->prev=newitem->prev;
		dead->prev->next=dead;
		//Insert ahead of current:
		Item<Element> *tmp=current->next;
		current->next=newitem;
		tmp->prev=newitem;
		newitem->prev=current;
		newitem->next=tmp;
//?		current=newitem;
	}	else {//empty collection: initialize first:
		dead=dead->next;
		dead->prev=first->prev;
		dead->prev->next=dead;
		first->next=first->prev=first;
		newitem=last=current=first;
#ifdef OMP
		for(int i=0;i<maxcounters;i++)current1[i]=current;
#endif
	}
//? last=current; 
	nitems++;
	return &newitem->element;
}
template <class Element>
Element	*Collection<Element>::insert() {
#ifdef OMP
	WAIT; locked=true;
#endif
	Element *newelement=insert0();
#ifdef OMP
	newelement->index(number());
	locked=false;
#endif
	return newelement;
}
template <class Element>
Element	*Collection<Element>::insert(int icounter) {
#ifdef OMP
	WAIT; locked=true;
	if(current1[icounter]!=NULL)current=current1[icounter];
#endif
	Element *element=insert0();
#ifdef OMP
	for(int i=0;i<maxcounters;i++) if(current1[i]==NULL)current1[i]=current;
	locked=false;
#endif
	if(element==NULL) return NULL;
	return element;
}
template <class Element>
bool	Collection<Element>::insert(Element *element) {
	Element *newelement=insert();
	if(newelement==NULL) return false;
	newelement->copy(element);
	return true;
}
template <class Element>
bool	Collection<Element>::insert(Element *element,int icounter) {
#ifdef OMP
	WAIT;
	locked=true;
	if(current1[icounter]!=NULL)current=current1[icounter];
#endif
	Element *newelement=insert0();
#ifdef OMP
	for(int i=0;i<maxcounters;i++) if(current1[i]==NULL)current1[i]=current;
	locked=false;
#endif
	if(newelement==NULL) return false;
	newelement->copy(element);
	return true;
}
template <class Element>
Element *Collection<Element>::append0() { //Append after last:
	if(dead->next==dead) {//keep the last dummy dead item to save on if's
		cout<<"CAN'T APPEND ELEMENT: Collection OF "<<nitems<<" IS FULL\n";
		cout.flush();
		return NULL;
	}
	if(first!=dead) {
		Item<Element> *old=dead;
		//Resurrect one:
		dead=dead->next;
		dead->prev=old->prev;
		dead->prev->next=dead;
		//Insert behind last:
		Item<Element> *tmp=last->next;
		last->next=old;
		tmp->prev=old;
		old->prev=last;
		old->next=tmp;
		last=old;
	}	else {//empty collection: initailize first:
		dead=dead->next;
		dead->prev=first->prev;
		dead->prev->next=dead;
		first->next=first->prev=first;
		last=current=first;
#ifdef OMP
		for(int i=0;i<maxcounters;i++)current1[i]=current;
#endif
	}
//? last=current; 
	nitems++;
	return &last->element;
}
template <class Element>
Element *Collection<Element>::append() { //Append after last:
#ifdef OMP
	WAIT; locked=true;
#endif
	Element *element=append0();
#ifdef OMP
	locked=false;
#endif
	return element;
}
template <class Element>
Element *Collection<Element>::append(int icounter) { //Append after last:
#ifdef OMP
	WAIT; locked=true;
	if(current1[icounter]!=NULL)current=current1[icounter];
#endif
	Element *element=append0();
#ifdef OMP
	for(int i=0;i<maxcounters;i++) if(current1[i]==NULL)current1[i]=current;
	locked=false;
#endif
	if(element==NULL) return NULL;
	return element;
}
template <class Element>
bool	Collection<Element>::append(Element *element) {
	Element *newelement=append();
	if(newelement==NULL) return false;
	newelement->copy(element);
	return true;
}
template <class Element>
bool	Collection<Element>::append(Element *element,int icounter) {
#ifdef OMP
	WAIT; locked=true;
	if(current1[icounter]!=NULL)current=current1[icounter];
#endif
	Element *newelement=append0();
#ifdef OMP
	for(int i=0;i<maxcounters;i++) if(current1[i]==NULL)current1[i]=current;
	locked=false;
#endif
	if(newelement==NULL) return false;
	newelement->copy(element);
	return true;
}
template <class Element>
bool	Collection<Element>::remove0() {
	if(current==dead||nitems<=0) {
		cout<<"CAN'T REMOVE ITEM FROM AN EMPTY COLLECTION\n";
		cout.flush();
		return false;
	}
	if(current->next==current) {
		current->next=dead;
		current->prev=dead->prev;
		dead->prev->next=current;
		dead->prev=current;
		first=last=dead=current;
#ifdef OMP
		for(int i=0;i<maxcounters;i++)current1[i]=current;
#endif
	} else {
		current->next->prev=current->prev;
		current->prev->next=current->next;
		Item<Element> *deadprev=dead->prev;
		deadprev->next=current;
		dead->prev=current;
		if(current==first)first=first->next;
		if(current!=last) current=current->next;
		else {
			last=last->prev;
			current=last;
		}
		dead->prev->next=dead;
		dead=dead->prev;
		dead->prev=deadprev;
	}
	nitems--;
	return true;
}
template <class Element>
bool	Collection<Element>::remove() {
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
	for(int i=0;i<maxcounters;i++) if(current1[i]==current)current1[i]=NULL;
#endif
	returned=remove0();
#ifdef OMP
	for(int i=0;i<maxcounters;i++) if(current1[i]==NULL)current1[i]=current;
	locked=false;
#endif
	return returned;
}
template <class Element>
bool	Collection<Element>::remove(int icounter) {
	bool returned=true;
#ifdef OMP
	WAIT; locked=true;
	if(icounter>=0&&icounter<maxcounters&&current1[icounter]!=NULL)current=current1[icounter];
	for(int i=0;i<maxcounters;i++)if(current1[i]==current)current1[i]=NULL;
#endif
	returned=remove0();
#ifdef OMP
	for(int i=0;i<maxcounters;i++)if(current1[i]==NULL)current1[i]=current;
	locked=false;
#endif
	return returned;
}
/* ---------------------------------------- */
template <class Element>
struct Ptr {//The Element is deferenced
	Element	*element;
	Ptr	*next,*prev;
};
template <class Element>
struct Pool {// Memory pool of pointers
	int mptrs,nptrs;
	Ptr<Element> *hook;
	Pool(int nptrs);
	~Pool();
	inline int size() {return nptrs;}
	Ptr<Element> *get();
	bool put(Ptr<Element> *ptr);
	bool check(char *msg);
};
template <class Element>
Pool<Element>::Pool(int n) {
	hook=new Ptr<Element>[n];
	if(hook==NULL) {
		cout<<"CAN'T ALLOCATE "<<n<<" ITEMS FOR THE POOL\n";
		cout.flush();
		exit(1);
	}
	//Link all items in a loop:
	for(int i=1;i<n-1;i++) {
		Ptr<Element> *ptr=hook+i;
		ptr->next=hook+i+1;
		ptr->prev=hook+i-1;
		ptr->element=NULL;
	}
	Ptr<Element> *first=hook;
	first->next=hook+1;
	first->prev=hook+n-1;
	first->element=NULL;
	Ptr<Element> *last=hook+n-1;
	last->next=hook;
	last->prev=hook+n-2;
	last->element=NULL;
	mptrs=nptrs=n;
	if(Run::option.verbose||Run::option.debug) cout<<"Pool size: "<<nptrs<<endl;
}
template <class Element>
Pool<Element>::~Pool() {
//+	delete hook;
	mptrs=nptrs=0;
	hook=NULL;
}
template <class Element>
Ptr<Element>	*Pool<Element>::get() {
	if(hook->next==hook) {//leave one dummy element to save on if's
		cout<<"\nPool::get:POOL-ALLOCATION FAILED: POOL SIZE: "<<nptrs<<endl;
		cout.flush();
		return NULL;//alternatively, initialize data dump at this point and abort
	}
	Ptr<Element> *out=hook->next;
	hook->next=out->next;
	out->next->prev=hook;
	out->next=out->prev=NULL;///DDD: probably not needed, but just to be on the safe side
	nptrs--;
	return out;
}
template <class Element>
bool	Pool<Element>::put(Ptr<Element> *ptr) {
	if(nptrs==mptrs) {//leave one dummy element to save on if's
		cout<<"POOL IS FILLED UP\n";
		cout.flush();
		return false;// or dump data and abort
	}
	if(ptr==NULL) {
		cout<<"\nNULL POINTER RETURNED TO POOL\n";
		cout.flush();
		return false;
	}
//-if(ptr==hook->next){cout<<"Pool::put:hook->next=ptr="<<ptr<<endl;cout.flush();exit(1);}///DDD
	//Disconnect ptr from the old loop:
	ptr->prev->next=ptr->next;
	ptr->next->prev=ptr->prev;
	//Insert ptr into the pool-loop before the hook:
	ptr->prev=hook;
	ptr->next=hook->next;
	ptr->next->prev=ptr;
	hook->next=ptr;
	nptrs++;
	return true;
}
template <class Element>
	bool Pool<Element>::check(char *msg) {
	Ptr<Element> *tmp=hook;///DDD
	int n=0;///DDD
	do {///DDD
		n++;tmp=tmp->next;///DDD
		if(n>mptrs+100) break;
	} while(tmp!=hook);///DDD
	cout<<"Pool::check("<<msg<<"):nptrs="<<nptrs<<", n="<<n<<", hook="<<hook<<", hook->next="<<hook->next<<"     ";cout.flush();///DDD
	if(n!=nptrs){cout<<"broken pool!\n";cout.flush();exit(1);}///DDD
	cout<<endl;///DDD
}
/*
 * Container is a collection of lists of pointers to any element list
 * The pointers are allocated/deallocated from the fixed size memory pool.
 *
 */
template <class Element>
class Container {//Uses the Pool of Ptrs
	int	nptrs;// number of elements in the list
	bool locked;
//-	Pool<Element> *pool;
	Ptr<Element>
		*first,*last,*current,
		*first1[maxcounters],*current1[maxcounters];
	public:
	Container();
	~Container();
	inline bool islocked(){ return locked; }
	inline void lock() { locked=true; }
	inline void unlock() { locked=false; }
	inline int number(){ return nptrs; }
	inline void setFirst(){ first=last->next; }
	inline void setFirst(int i){ first1[i]=last->next; }
	inline void setFirstLast(){ first=last; }
	inline void setFirstLast(int i){ first1[i]=last; }
	inline void goFirst(){ 
		current=first; 
#ifdef OMP
		for(int i=0;i<maxcounters;i++)current1[i]=first1[i]=first;
#endif
	}
	inline void goFirst(int i){ current1[i]=first1[i]; }
	inline void goFirstNext(){ current=first->next; }
	inline void goFirstNext(int i){ current1[i]=first1[i]->next; }
	inline void FirstNext(){ first=first->next; }
	inline void FirstNext(int i){ first1[i]=first1[i]->next; }
	inline void FirstPrev(){ first=first->prev; }
	inline void FirstPrev(int i){ first1[i]=first1[i]->prev; }
	inline bool isFirstLast(){ if(first==last) return true; else return false; }
	inline bool isFirstLast(int i){ if(first1[i]==last) return true; else return false; }
	inline void goLastNext(){ current=last->next; }
	inline void goLastNext(int i){ current1[i]=last->next; }
	inline void goLast(){ current=last; }
	inline void goLast(int i){ current1[i]=last; }
	inline void goNext(){ current=current->next; }
	inline void goNext(int i){ current1[i]=current1[i]->next; }
	inline void goPrev(){ current=current->prev; }
	inline void goPrev(int i){ current1[i]=current1[i]->prev; }
	inline bool isFirst(){ return current==first?true:false; }
	inline bool isFirst(int i){ return current1[i]==first1[i]?true:false; }
	inline bool isLast(){ if(current==last) return true; else return false; }
	inline bool isLast(int i){ if(current1[i]==last) return true; else return false; }
	inline void	LastNext(){ last=last->next; }
	inline void	LastPrev(){ last=last->prev; }
	inline void go(Ptr<Element> *p) { current=p; }//not safe
	inline void go(Ptr<Element> *p, int i) { current1[i]=p; }//not safe
	inline Element *First(){ return first->element; }
	inline Element *First(int i){ return first1[i]->element; }
	inline Element *Last(){ return last->element; }
	inline Element *Next(){ current=current->next; return current->element; }
	inline Element *Next(int i){ current1[i]=current1[i]->next; return current1[i]->element; }
	inline Element *getNext(){ return current->next->element; }
	inline Element *getNext(int i){ return current1[i]->next->element; }
	inline Element *Prev(){ current=current->prev; return current->element; }
	inline Element *Prev(int i){ current1[i]=current1[i]->prev; return current1[i]->element; }
	inline Element *getPrev(){ return current->prev->element; }
	inline Element *getPrev(int i){ return current1[i]->prev->element; }
	inline Element *Current(){ return current->element; }
	inline Element *Current(int i){ return current1[i]->element; }
	inline Ptr<Element> *getPtr() { return current; }
	inline Ptr<Element> *getPtr(int i) { return current1[i]; }
	inline Ptr<Element> *getFirstPtr() { return first; }
	inline Ptr<Element> *getFirstPtr(int i) { return first1[i]; }
	inline Ptr<Element> *getLastPtr() { return last; }
	bool insert0(Element *element);
	bool insert(Element *element);
	bool insert(Element *element,int i);
	Ptr<Element> *insert1(Element *element,int i);
	bool insert0();
	bool insert();
	bool insert(int i);
	bool append0(Element *element);
	bool append(Element *element);
	bool append(Element *element, int i);
	bool append0();
	bool append();
	bool append(int i);
	Ptr<Element> *append1(Element *element,int i);
	bool link0(Ptr<Element> *&ptr);
	bool link(Ptr<Element> *&ptr);
//-	bool link0(Ptr<Element> *&ptr, int i);
//-	bool link(Ptr<Element> *&ptr, int i);
	bool unlink0(Ptr<Element> *&ptr);
	bool unlink(Ptr<Element> *&ptr);
//-	bool unlink(Ptr<Element> *&ptr, int i);
	bool remove0();
	bool remove();
	bool remove(int i);
	bool remove0(Ptr<Element> *&ptr);
	bool remove(Ptr<Element> *ptr);
	bool remove(Ptr<Element> *ptr, int i);
	bool checkPool(char *msg, Pool<Element> *pool);
};
