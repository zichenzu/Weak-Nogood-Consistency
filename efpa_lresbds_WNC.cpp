                                             
                                                                              
#include <gecode/driver.hh>                                                                                        
#include <gecode/int.hh>                                       
               
                                                                                           
#include <gecode/WNC/LReSBDSBrancher_WNC.hpp> 
                                                                    
                                                                                              
using namespace Gecode;                                                                                 
IntArgs VH;     
 
int _symmetries (int a,int b,int c, int d, int index){
  
	if(index>=nNum*(a-1)&&index<nNum*a) index=index+(b-a)*nNum;
	else if(index>=nNum*(b-1)&&index<nNum*b) index=index-(b-a)*nNum;
	if(index%nNum==c-1) index=index+d-c;
	else if(index%nNum==d-1) index=index-d+c;

	return index; }                                                                                      
/**                                                               
 * \brief %Options for %EFPA problems                              
 *                             
 */        
class EFPAOptions : public Options {                
public: 
  int q,l,d,v,a1,a2,a3,a4;        ///< Derived parameters           
  /// Derive additional parameters
  
  /// Initialize options for example with name \a s  
  EFPAOptions(const char* s,
              int q0, int l0, int d0,int v0) 
    : Options(s), q(q0), l(l0), d(d0),v(v0) {                     
    
  }          
  /// Parse options from arguments \a argv (number is \a argc)
  void parse(int& argc, char* argv[]) {     
    Options::parse(argc,argv);        
    if (argc < 5) 
      return; 
    q = atoi(argv[1]);    
    l = atoi(argv[2]);   
    d = atoi(argv[3]);     
    v = atoi(argv[4]);
  }    
  /// Print help message    
  virtual void help(void) {         
    Options::help();
    std::cerr << "\t(unsigned int) default: " << q << std::endl
              << "\t\tparameter v" << std::endl
              << "\t(unsigned int) default: " << l << std::endl
              << "\t\tparameter k" << std::endl
              << "\t(unsigned int) default: " << d << std::endl
              << "\t\tparameter lambda" << std::endl;
  }  
}; 
   
/**
 * \brief %Example: Balanced incomplete block design (%EFPA)
 *
 * See problem 28 at http://www.csplib.org/. 
 *
 * \ingroup Example
 *
 */
class EFPA : public Script {
protected:
  /// Options providing access to parameters
  const EFPAOptions& opt;
  /// Matrix of Boolean variables
  int v; ///< Number of sequences
  int q; ///< Number of symbols
  int l; ///< Number of sets of symbols for a sequence (\f$\lambda\f$)
  int d; ///< Hamming distance between any pair of sequences
  int n; ///< Length of sequence (\f$q\cdot\lambda\f$)
  int nseqpair;  ///< Number of sequence pairs (\f$\frac{v(v-1)}{2}\f$)
  IntVarArray  c; ///< Variables for sequences
  BoolVarArray diff; ///< Differences between sequences

public:
  /// Symmetry breaking variants
  enum { 
    SYMMETRY_NONE,      ///< No symmetry breaking
    SYMMETRY_LEX,       ///< DoubleLex-constraints on rows/columns 
    SYMMETRY_SNAKE,       ///< SnakeLex-constraints on rows/columns
    SYMMETRY_LDSB,       ///< LDSB on rows/columns
    SYMMETRY_LRESBDS      ///< LReSBDS on rows/columns
  }; 
 
  /// Actual model 
  EFPA(const EFPAOptions& o)                            
    : opt(o),  v(opt.v),
      q(opt.q),
      l(opt.l),
      d(opt.d),   
      n(q*l),
      nseqpair((v*(v-1))/2),
      c(*this, n*v, 1,q),
      diff(*this, n*nseqpair, 0, 1)
  { 
   
	 
    // Matrix access
    // q*lambda=n columns, and v rows
    Matrix<IntVarArray> cm(c, n, v);
    // q*lambda=n columns, and nseqpair rows
    Matrix<BoolVarArray> diffm(diff, n, nseqpair);       

    // Counting symbols in rows 
    {     
      IntArgs values(q);
      for (int i = q; i--; ) values[i] = i+1;
      IntSet cardinality(l, l);
      for (int i = v; i--; ) 
        count(*this, cm.row(i), cardinality, values, opt.icl());
    }     
           
    // Difference variables      
    {
      int nseqi = 0;  
      for (int a = 0; a < v; ++a) {  
        for (int b = a+1; b < v; ++b) { 
          for (int i = n; i--; ) {
            rel(*this, cm(i, a), IRT_NQ, cm(i, b), diffm(i, nseqi));
          }
          ++nseqi;        
        }
      }
      assert(nseqi == nseqpair);
    }  

    // Counting the Hamming difference
    {
      for (int i = nseqpair; i--; ) {   
        linear(*this, diffm.row(i), IRT_EQ, d);  
      }
    } 
 
	
	//build a rowwise heuristic
	VH=IntArgs (c.size());         
	for(int i=0;i<v;i++)
		if(i%2==0)
			for(int j=0;j<n;j++) 
				VH[i*n+j]=i*n+j;
		else
			for(int j=0;j<n;j++)
				VH[i*n+j]=i*n+n-1-j;   
   
    if (opt.symmetry() == SYMMETRY_LDSB) {           
		Symmetries s;    
		s << rows_interchange(cm);     
		s << columns_interchange(cm);
		branch(*this, c, INT_VAR_NONE(), INT_VAL_MIN(), s); 
    } 
	else if (opt.symmetry() == SYMMETRY_LRESBDS)      
    {    
		      
		branch(*this,c,INT_VAR_NONE(), INT_VAL_MIN(),v,n,q+1); 

    }         
    else 
	  if (opt.symmetry() == SYMMETRY_SNAKE)         
      {      
 	  
		//rowwise snakelex symmetry breaking   
		IntVarArgs rows;    
		for(int i=0;i<v;i++)
			for(int j=0;j<n;j++)
				rows<<cm(n-1-j,i);
		Matrix<IntVarArgs> s(rows,n,v);    
		
		//row symmetry
		for (int i=0; i<v-1; i++)   
        {
			if(i%2==0) 
			{
				rel(*this, cm.row(i), IRT_LQ, cm.row(i+1));
				if(i+2<v)
				rel(*this, cm.row(i), IRT_LQ, cm.row(i+2));
			}
			else
			{
				rel(*this, s.row(i), IRT_LQ, s.row(i+1));
				
				if(i+2<v)
				rel(*this, s.row(i), IRT_LQ, s.row(i+2));
			
			}
 
		} 
		  
		//column symmetry  
        for (int j=0; j<n-1; j++)
        {
			IntVarArgs t1;    
			IntVarArgs t2;    
			for(int i=0;i<v;i++)
				if(i%2==0)
				{	t1<<cm(j,i); 
					t2<<cm(j+1,i);
				}
				else
				{	t2<<cm(j,i); 
					t1<<cm(j+1,i);
				}
			rel(*this, t1, IRT_LQ, t2);  

			
		}
 
		branch(*this,c,INT_VAR_NONE(), INT_VAL_MIN(),0);  
     
		  
		  

      }
     	  
      else { 
		  if (opt.symmetry() == SYMMETRY_LEX) {  
			 
			for (int i=1; i<v; i++)  
			 rel(*this, cm.row(i-1), IRT_LQ, cm.row(i));
		  
			 
			for (int j=1; j<n; j++)
			  rel(*this, cm.col(j-1), IRT_LQ, cm.col(j));    
	  
		  } 

          branch(*this,c,INT_VAR_NONE(), INT_VAL_MIN());   
      }
   
  }   
 
  /// Print solution 
  virtual void
  print(std::ostream& os) const {
 
  }

  /// Constructor for cloning \a s 
  EFPA(bool share, EFPA& s)
    : Script(share,s), opt(s.opt) {      
    c.update(*this,share,s.c);
  }

  /// Copy during cloning  
  virtual Space* 
  copy(bool share) { 
    return new EFPA(share,*this);
  }
 
};
  
/** \brief Main-function    
 *  \relates EFPA       
 */        
int
main(int argc, char* argv[]) {            
  EFPAOptions opt("EFPA",5,3,3,4);  
       
  opt.symmetry(EFPA::SYMMETRY_LRESBDS);  
  opt.symmetry(EFPA::SYMMETRY_NONE,"none");  
  opt.symmetry(EFPA::SYMMETRY_LEX,"double"); 
  opt.symmetry(EFPA::SYMMETRY_SNAKE,"snake");
  opt.symmetry(EFPA::SYMMETRY_LDSB,"ldsb");
  opt.symmetry(EFPA::SYMMETRY_LRESBDS,"lresbds");
  opt.solutions(0);
  opt.c_d(1);
  opt.parse(argc,argv);

  /*
   * Other interesting instances:
   * EFPA(7,3,1), EFPA(6,3,2), EFPA(7,3,20), ...
   */

  Script::run<EFPA,DFS,EFPAOptions>(opt);
  return 0;
}
 
