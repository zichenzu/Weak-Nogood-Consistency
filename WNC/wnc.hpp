#include <gecode/int/branch.hh>
#include <gecode/int.hh> 
#include <gecode/WNC/head.hh>
namespace Gecode{  
 
 
	extern IntArgs PartA; 
	extern int ex_nSym, ex_psize; 
	extern IntArgs ex_VarM,ex_ValM, ex_PosM,ex_nSymR, checkA;
	extern int ex_dosize;
	extern	IntArgs XP,VP;
	extern	int re;
	extern int conN;
	template<class View>
	 class WNC : public Propagator {
	 protected: 
		ViewArray<View> x;
		int id;
		int y0,v0,p0;
		int y1,v1,p1;
	 public: 
		//post
		WNC(Space& home,  ViewArray<View>& x0,  int _id, int _y0,int _v0,int _p0,int _y1,int _v1,int _p1);
		//copy
		WNC(Space& home, bool share, WNC<View>& p);
		virtual Propagator* copy(Space& home, bool share);
		//cost
		virtual PropCost cost(const Space&, const ModEventDelta&) const;
		//propagation
		virtual ExecStatus propagate(Space& home, const ModEventDelta&);
		//post
		static ExecStatus post(Space& home, ViewArray<View>& x0,  int _id, int _y0,int _v0,int _p0,int _y1,int _v1,int _p1);
		//dispose
		virtual size_t dispose(Space& home);
		void symGoal(Space& home, int index, int val);
	};

	template<class View>
	void
	WNC<View>::symGoal(Space& home, int index, int val) {
		 
					 
		// special for IntVars but could be generalised
		if (ex_nSym == 0) { return;} // i.e. no goal 	
				   // better way to return Truth?
		else { 
 
			
			for (int i = 0; i < ex_nSym; i++)
			{
				 
				Index_class k=(*_symmetries)(ex_nSymR[i],index,val);
				if(!x[k.index].in(k.val) || (k.index==index&&k.val==val)) continue;
				
				if(ex_VarM[i]==-1&&ex_ValM[i]==-1 &&  checkA[k.index*ex_dosize+k.val]==-1) 
				{
					checkA[k.index*ex_dosize+k.val]=1;
					XP[re]=k.index;
					VP[re++]=k.val;
				}
				else  if(ex_VarM[i]!=-1&&ex_ValM[i]!=-1 && x[ex_VarM[i]].in(ex_ValM[i]))
				{
					conN++;
					WNC(home, x,ex_nSymR[i], ex_VarM[i], ex_ValM[i], ex_PosM[i], k.index,k.val,ex_psize);
				}
			}
			
		 
			return; 
		}  
			 
	};
		
	// posting
	template<class View>
	inline
	WNC<View>::WNC(Space& home,   ViewArray<View>& x0, int _id, int _y0,int _v0,int _p0,int _y1,int _v1,int _p1) : 
			Propagator(home), x(x0), y0(_y0), v0(_v0), p0(_p0), y1(_y1), v1(_v1), p1(_p1) , id(_id)
	{
		 
		x[y0].subscribe(home,*this,Int::PC_INT_VAL); 
		//x.subscribe(home,*this,Int::PC_INT_VAL);  
	} 
	template<class View>
	ExecStatus
	WNC<View>::post(Space& home, ViewArray<View>& x0, int id, int y0,int v0,int p0,int y1,int v1,int p1) {
	
	
	
		 
		(void) new (home) WNC(home,x0, id, y0, v0, p0, y1,v1,p1); 
		return ES_OK; 
	} 
	// disposal 
	template<class View>
	forceinline size_t
	WNC<View>::dispose(Space& home) {
		 
 
		(void) Propagator::dispose(home);
		return sizeof(*this); 
	} 
	// copying
	template<class View>
	forceinline
	WNC<View>::WNC(Space& home, bool share, WNC<View>& p) : Propagator(home,share,p),
						 id(p.id),y0(p.y0),v0(p.v0),p0(p.p0),y1(p.y1),v1(p.v1),p1(p.p1){
												 
		 
		x.update(home,share,p.x); 
		 
	} 
	template<class View> Propagator*  WNC<View>::copy(Space& home, bool share) {
		return new (home) WNC(home,share,*this); 
	} 
	// cost computation 
	template<class View>
	PropCost
	WNC<View>::cost(const Space&, const ModEventDelta&) const {
		return PropCost::linear(PropCost::LO,int(100000000000));
	} 
	
	 
	// propagation 
	template<class View>
	ExecStatus
	WNC<View>::propagate(Space& home, const ModEventDelta&) { 
		
		if(y0==-1 || v0==-1) std::cout<<"error now\n";
		if(!x[y0].in(v0))
		{
			x[y0].cancel(home,*this,Int::PC_INT_VAL);
			return home.ES_SUBSUMED(*this);
		}
		if(x[y0].assigned())
		{		
			x[y0].cancel(home,*this,Int::PC_INT_VAL);
			int s; 
			for(s=p0+1;s<p1;s++) 
			{
				 
				Index_class g=(*_symmetries)(id,PartA[2*s],PartA[1+2*s]);
				if(!x[g.index].in(g.val))//now the symmetry is broken
				{	
					return home.ES_SUBSUMED(*this);
				}
				if(!x[g.index].assigned())//now find the first pointer
				{
					y0=g.index;
					v0=g.val;
					p0=s;  
					x[g.index].subscribe(home,*this,Int::PC_INT_VAL);
					break;
				}
			 
			}
			if(s==p1)//can prune
			{
				if(x[y1].in(v1)){
					GECODE_ME_CHECK(x[y1].nq(home,v1));
					symGoal(home,y1,v1);
				}
				return home.ES_SUBSUMED(*this);
			}
		}
	  
		return ES_NOFIX;

	} 


	void wnc(Space& home, ViewArray<Int::IntView>& xv, int id, int y0,int v0,int p0,int y1,int v1,int p1) { 
		// constraint post function  
		if (WNC<Int::IntView>::post(home,xv,id, y0, v0, p0, y1,v1,p1) != ES_OK)
			home.fail(); 
	} 
}