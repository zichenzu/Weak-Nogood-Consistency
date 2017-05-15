
//***********************************************************************//

#include <gecode/int/branch.hh>
#include <gecode/int.hh>
#include <gecode/WNC/wnc.hpp>
namespace Gecode { 
	int nNum;
    IntArgs Aa,Bb,Cc,Dd,reF;
	IntArgs PartA; 
	int ex_nSym;
	IntArgs ex_VarM,ex_ValM, ex_PosM,ex_nSymR;
	IntArgs checkA;
	int ex_psize;
	int ex_dosize;
	IntArgs XP,VP;
	int re;
	int conN;
	Index_class _symmetries (int id, int index, int val){
	   
	 if(reF[id]==1)
	 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
		else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;}
	 else if(reF[id]==2)
	 {	if(index%nNum==Aa[id]-1) index=index+Bb[id]-Aa[id];
		else if(index%nNum==Bb[id]-1) index=index-Bb[id]+Aa[id];}
	 else if (reF[id]==3)
	 {	if(index>=nNum*(Aa[id]-1)&&index<nNum*Aa[id]) index=index+(Bb[id]-Aa[id])*nNum;
		else if(index>=nNum*(Bb[id]-1)&&index<nNum*Bb[id]) index=index-(Bb[id]-Aa[id])*nNum;
		if(index%nNum==Cc[id]-1) index=index+Dd[id]-Cc[id];
		else if(index%nNum==Dd[id]-1) index=index-Dd[id]+Cc[id];
	 }
	 return Index_class(index,val); 
	}
	//symmetry class
	class SBDS_Sym  {
		private:
			Space& _home;
			int _nSym;
			

		public:
			IntArgs _nSymR;

			int dosize;
			IntArgs VarM;
			IntArgs ValM;
			IntArgs PosM; 
			~SBDS_Sym() {} ;
			Space& getManager() {return _home;};
		  
			int getNSym() {return _nSym;};
			 

		   

			void symGoal(ViewArray<Int::IntView>& vars, int index, int val ,int p_size );//-----------------------------this is the normal SBDS adds conditional contrains
		 
			SBDS_Sym(Space& home, int nSym,int var_n,int dosize) ;
			void SymAdjust(ViewArray<Int::IntView>& vars, int index, int val, int p_size) ;
			SBDS_Sym(Space& home, SBDS_Sym OldSym) ;
	};



	void SBDS_Sym::symGoal(ViewArray<Int::IntView>& vars, int index, int val, int p_size) {
				 
		// special for IntVars but could be generalised
		if (_nSym == 0) { return;} // i.e. no goal 	
				   // better way to return Truth?
		else { 
			int i = 0; 

			for (i = 0; i < _nSym; i++)
			{
				Index_class k=(*_symmetries)(_nSymR[i],index,val);
				if(!vars[k.index].in(k.val) || (k.index==index&&k.val==val) ) continue; 
				if(VarM[i]!=-1&&ValM[i]!=-1 && vars[VarM[i]].in(ValM[i]))
				{	
					conN++;
					wnc(_home, vars,_nSymR[i], VarM[i], ValM[i], PosM[i], k.index,k.val,p_size);
				}
				else if(VarM[i]==-1&&ValM[i]==-1 && checkA[k.index*dosize+k.val]==-1) 
				{	
					 checkA[k.index*dosize+k.val]=1;

					XP[re]=k.index;
					VP[re++]=k.val;
					symGoal(vars,k.index, k.val,p_size);
					
				}
				 
			}
	 
			return; 
		}  
	};
				 
	SBDS_Sym::SBDS_Sym(Space& home, int nSym,int var_n,int dosize0) :
					_home(home),
					_nSym(nSym),
					
					dosize(dosize0)
	{
		 
		_nSymR=IntArgs(_nSym);
		VarM=IntArgs(_nSym);
		ValM=IntArgs(_nSym);
		PosM=IntArgs(_nSym);
		
		ex_VarM=IntArgs(_nSym);
		ex_ValM=IntArgs(_nSym);
		ex_PosM=IntArgs(_nSym);
		ex_nSymR=IntArgs(_nSym);
		checkA=IntArgs(var_n*dosize0); 
		XP=IntArgs(var_n*dosize0);
		VP=IntArgs(var_n*dosize0);
		ex_dosize=dosize;
		for(int i=0;i<_nSym;i++)
		{
			_nSymR[i]=i;
			VarM[i]=-1;
			ValM[i]=-1;
			PosM[i]=-1;
		}
	};
	SBDS_Sym::SBDS_Sym(Space& home,SBDS_Sym OldSym) :
					 _nSym(OldSym.getNSym()),
					_home(home),
					dosize(OldSym.dosize),
					_nSymR(OldSym.getNSym()),
					VarM(OldSym.getNSym()),
					ValM(OldSym.getNSym()),
					PosM(OldSym.getNSym())
	{
		for(int i=0;i<_nSym;i++)
		{
			_nSymR[i]=OldSym._nSymR[i];
			VarM[i]=OldSym.VarM[i];
			ValM[i]=OldSym.ValM[i];
			PosM[i]=OldSym.PosM[i];
		}




	};
	void SBDS_Sym::SymAdjust(ViewArray<Int::IntView>& vars, 
				   int index, int val, int p_size) 
	{
		int oldNSym = _nSym; 

	  
		_nSym = 0;
	 


		for (int i = 0 ; i < oldNSym ; i++) {
	 
			
			if(VarM[i]!=-1&& ValM[i]!=-1 && !vars[VarM[i]].in(ValM[i]))	continue;
			if(VarM[i]==-1&& ValM[i]==-1)
			{
				Index_class k=(*_symmetries)(_nSymR[i],index,val);
				if(!vars[k.index].in(k.val)) continue;
				else
					if(k.index!=index||k.val!=val)
					{
						VarM[_nSym]=k.index;
						ValM[_nSym]=k.val;
						PosM[_nSym]=p_size-1;
						_nSymR[_nSym++]=_nSymR[i];
					}
					else
					{
						VarM[_nSym]=-1;
						ValM[_nSym]=-1;
						PosM[_nSym]=-1;
						_nSymR[_nSym++]=_nSymR[i];
					}
			}
			else if(vars[VarM[i]].in(ValM[i])&& vars[VarM[i]].assigned())
			{
				int s; 
				for(s=PosM[i]+1;s<p_size;s++)
				{
				 
					Index_class g=(*_symmetries)(_nSymR[i],PartA[2*s],PartA[1+2*s]);
					if(!vars[g.index].in(g.val))//now the symmetry is broken
					{	
						break;
					}
					if(!vars[g.index].assigned())//now find the first pointer
					{
						VarM[_nSym]=g.index;
						ValM[_nSym]=g.val;
						PosM[_nSym]=s;
						_nSymR[_nSym++]=_nSymR[i];
						break;
					}
			 
				}
				if(s==p_size)
				{
					VarM[_nSym]=-1;
					ValM[_nSym]=-1;
					PosM[_nSym]=s-1;
					_nSymR[_nSym++]=_nSymR[i];
				}
			}					
			else //if(_nSym!=i)
			{
				VarM[_nSym]=VarM[i];
				ValM[_nSym]=ValM[i];
				PosM[_nSym]=PosM[i];
				_nSymR[_nSym++]=_nSymR[i];
			} 
 
		}
	 
	 
	 
	};

	//**************************************************//
	//**************************************************//
	//**************************************************//
	//**************************************************//
	//**************************************************//
	//*************   branch ***************************//

	 template<class Val>
	  class GECODE_VTABLE_EXPORT SBDSChoice : public PosValChoice<Val> {
	  
	  public:
		/// Initialize choice for brancher \a b, position \a p, value \a
		/// n, and set of literals \a literals (of size \a nliterals)
		SBDSChoice(const Brancher& b, unsigned int a, const Pos& p, const Val& n);
		/// Destructor
		~SBDSChoice(void);

		virtual size_t size(void) const;

	  };

	  template<class Val>
	  inline
	  SBDSChoice<Val>::SBDSChoice(const Brancher& b, unsigned int a, const Pos& p, 
								  const Val& n)
		: PosValChoice<Val>(b,a,p,n)
		{}

	  template<class Val>
	  SBDSChoice<Val>::~SBDSChoice(void) {
		
	  }
	  
	  template<class Val>
	  size_t
	  SBDSChoice<Val>::size(void) const {
		return sizeof(SBDSChoice<Val>);
	  }

	  


	 template<class View, int n, class Val, unsigned int a>
	class LReSBDSBrancher : public ViewValBrancher<View,n,Val,a> {
		typedef typename ViewBrancher<View,n>::BranchFilter BranchFilter;
	  public:
		/// Array of symmetry implementations
		SBDS_Sym SymObject;
		int dosize;
		int start;
		int p_size;
	  protected:
		/// Constructor for cloning \a b
		LReSBDSBrancher(Space& home, bool share, LReSBDSBrancher& b);
		/// Constructor for creation
		
		
		LReSBDSBrancher(Home home, 
					 ViewArray<View>& x,
					 ViewSel<View>* vs[n], 
					 ValSelCommitBase<View,Val>* vsc,
					 int nSym, int dosize0,
					 BranchFilter bf,IntVarValPrint vvp);
	  public:
		/// Return choice
		virtual const Choice* choice(Space& home);
		/// Return choice
		virtual const Choice* choice(const Space& home, Archive& e);
		/// Perform commit for choice \a c and alternative \a b
		virtual ExecStatus commit(Space& home, const Choice& c, unsigned int b);
		/// Perform cloning
		virtual Actor* copy(Space& home, bool share);
		/// Perform dispose
		virtual size_t dispose(Space& home);
		/// Delete brancher and return its size
		static BrancherHandle post(Home home,
								   ViewArray<View>& x,
								   ViewSel<View>* vs[n],
								   ValSelCommitBase<View,Val>* vsc,
								   int nSym,int dosize0,
								   BranchFilter bf,IntVarValPrint vvp);
	  };

	  template<class View, int n, class Val, unsigned int a>
	  LReSBDSBrancher<View,n,Val,a>
	  ::LReSBDSBrancher(Home home, ViewArray<View>& x,
					 ViewSel<View>* vs[n],
					 ValSelCommitBase<View,Val>* vsc,
					 int nSym, int dosize0,
					 BranchFilter bf,IntVarValPrint vvp)
		: ViewValBrancher<View,n,Val,a>(home, x, vs, vsc, bf,vvp),
		  SymObject(home,nSym,x.size(),dosize0),
	dosize(dosize0)
	  {	start=0;
	  PartA=IntArgs(2*x.size());
	  p_size=0;
		 home.notice(*this, AP_DISPOSE);

	  }

	  template<class View, int n, class Val, unsigned int a>
	  forceinline BrancherHandle
	  LReSBDSBrancher<View,n,Val,a>::
	  post(Home home, ViewArray<View>& x,
		   ViewSel<View>* vs[n], ValSelCommitBase<View,Val>* vsc,
		   int nSym,int dosize0,
		   BranchFilter bf,IntVarValPrint vvp) {
		return *new (home) LReSBDSBrancher<View,n,Val,a>(home,x,vs,vsc,nSym,dosize0,bf,vvp);
	  }

	  template<class View, int n, class Val, unsigned int a>
	  forceinline
	  LReSBDSBrancher<View,n,Val,a>::
	  LReSBDSBrancher(Space& home, bool shared, LReSBDSBrancher<View,n,Val,a>& b)
		: ViewValBrancher<View,n,Val,a>(home,shared,b), 
		  SymObject(home,b.SymObject),
		  dosize(b.dosize) {
		  p_size=b.p_size;
			start=b.start;
		  

	  }
	  
	  template<class View, int n, class Val, unsigned int a>
	  Actor*
	  LReSBDSBrancher<View,n,Val,a>::copy(Space& home, bool shared) {
		return new (home) LReSBDSBrancher<View,n,Val,a>(home,shared,*this);
	  }


	  // Compute choice
	  template<class View, int n, class Val, unsigned int a>
	  const Choice*
	  LReSBDSBrancher<View,n,Val,a>::choice(Space& home) {

		 
		const Choice* c = ViewValBrancher<View,n,Val,a>::choice(home);
		const PosValChoice<Val>* pvc = static_cast<const PosValChoice<Val>* >(c);
		
		// Compute symmetries.

		int choicePos = pvc->pos().pos;
		int choiceVal = pvc->val();
		return new SBDSChoice<Val>(*this,a,choicePos,choiceVal);
	  }

	 template<class View, int n, class Val, unsigned int a>
	  const Choice*
	  LReSBDSBrancher<View,n,Val,a>::choice(const Space& home, Archive& e) {

		(void) home;
		int p; e >> p;
		Val v; e >> v;
		return new SBDSChoice<Val>(*this,a,p,v);
	  }
	template<class View, int n, class Val, unsigned int a>
	  size_t
	  LReSBDSBrancher<View,n,Val,a>::dispose(Space& home) {
		home.ignore(*this,AP_DISPOSE);
		(void) ViewValBrancher<View,n,Val,a>::dispose(home);
		return sizeof(LReSBDSBrancher<View,n,Val,a>);
	  }
	  template<class View, int n, class Val, unsigned int a>
	  ExecStatus
	  LReSBDSBrancher<View,n,Val,a>
	  ::commit(Space& home, const Choice& c, unsigned int b) {
	  

		 const SBDSChoice<Val>& pvc
		  = static_cast<const SBDSChoice<Val>&>(c);
		int pos = pvc.pos().pos;
		int val = pvc.val();
 
		for(int i=0;i<checkA.size();i++)
			checkA[i]=-1;
		re=0;
		 
		conN=0;
		if(b==0)
		{			
		 
			PartA[p_size*2]=pos;
			PartA[p_size*2+1]=val;
			p_size++;
			ex_psize=p_size;
			SymObject.SymAdjust(this->x,pos,val,p_size);//ajust symmetries
			 
			ex_nSym=SymObject.getNSym();
			for(int i=0;i<ex_nSym;i++)
			{
				ex_VarM[i]=SymObject.VarM[i];
				ex_ValM[i]=SymObject.ValM[i];
				ex_PosM[i]=SymObject.PosM[i];
				ex_nSymR[i]=SymObject._nSymR[i];
			}
			 GECODE_ME_CHECK(this->x[pos].eq(home,val));
			 StatusStatistics c;
			 home.status(c); 
			 
			 //ExecStatus fromBase = ViewValBrancher<View,n,Val,a>::commit(home, c, b);
			 //GECODE_ES_CHECK(fromBase);
		}
		else
		{
			ex_psize=p_size;
			ex_nSym=SymObject.getNSym();
			for(int i=0;i<ex_nSym;i++)
			{
				ex_VarM[i]=SymObject.VarM[i];
				ex_ValM[i]=SymObject.ValM[i];
				ex_PosM[i]=SymObject.PosM[i];
				ex_nSymR[i]=SymObject._nSymR[i];
			}
			GECODE_ME_CHECK(this->x[pos].nq(home,val));
			   
			SymObject.symGoal(this->x,pos, val,p_size);
			while(re!=0)
			{
				int re1=re;
				re=0;
				for(int i=0;i<re1;i++)
				{
					GECODE_ME_CHECK(this->x[XP[i]].nq(home,VP[i]));
					SymObject.symGoal(this->x,XP[i], VP[i],p_size);
				}
			  
			} 
			
			
			StatusStatistics c;
			home.status(c); 
			//ExecStatus fromBase = ViewValBrancher<View,n,Val,a>::commit(home, c, b);
			//GECODE_ES_CHECK(fromBase);
		}


		return ES_OK;
	  }
	  
	//**************************************************//
	//**************************************************//
	//**************************************************//
	//**************************************************//
	//**************************************************//
	//*************   post ***************************//

	 BrancherHandle
	  branch(Home home, const IntVarArgs& x,
			 IntVarBranch vars, IntValBranch vals,
			 int n,int m,int dosize0,IntBranchFilter bf=NULL) {
		using namespace Int;
		if (home.failed()) return BrancherHandle();
		ViewArray<IntView> xv(home,x);
		ViewSel<IntView>* vs[1] = { 
		  Branch::viewselint(home,vars) 
		};
		int nSym=n*m*n*m;
	 
		int a=n*(n-1)/2;
		int b=m*(m-1)/2;

		nSym=n*m*n*m;
		std::cout<<nSym<<"\t"<<n<<","<<m<<"\n";
		nNum=m;
		Aa=IntArgs(nSym);
		Bb=IntArgs(nSym);
		Cc=IntArgs(nSym);
		Dd=IntArgs(nSym);
		reF=IntArgs(nSym);
		 
 
		IntArgs _A(a+b);
		IntArgs _B(a+b);
		int r=0;
		int r1=0;
		 
		 
		for(int i=0;i<n-1;i++)
		{
			int j=i+1;
			for(int j=i+1;j<n;j++)
			{
				_A[r1]=i+1;
				_B[r1++]=j+1;

			}
			Aa[r]=i+1;
			Bb[r]=i+2;
			 
			reF[r++]=1;
		}
		
		int re0=r;

		for(int i=0;i<m-1;i++)
		{  
		   int j=i+1;
		   for(int j=i+1;j<m;j++)
		   {
			 _A[r1]=i+1;
			 _B[r1++]=j+1;
			 
		   } 
			Aa[r]=i+1;
			Bb[r]=i+2;
			reF[r++]=2;
		}
		 
		for(int i=0;i<a;i++)
			for(int j=0;j<b;j++)
			{ 	Aa[r]=_A[i];
				Bb[r]=_B[i];
				Cc[r]=_A[a+j];
				Dd[r]=_B[a+j];
			 
				reF[r++]=3;
			 
			}
		 
		 

		nSym=r;
		std::cout<<"---------------------nSym is "<<nSym<<std::endl;


		return LReSBDSBrancher<IntView,1,int,2>::post
			(home,xv,vs,Branch::valselcommitint(home,x.size(),vals),nSym,dosize0,bf,NULL);
	}
	   
	   
   
   
}


  

