
namespace Gecode {
 
	//index_class
	class Index_class{
		public: 
			int index;
			int val;
			~Index_class(){};
			Index_class(int _index,int _val){index=_index;val=_val;};
				
	};
  
	typedef Index_class (SymmetryFunction)(int ,int,int);
	int solr;
	Index_class _symmetries (int id, int index, int val);


}


  

