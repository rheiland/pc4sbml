
int read_sbml_file(const char* file);

#undef foobar
#ifdef foobar
class Molecular
{
	private:
	public: 
		Microenvironment* pMicroenvironment; 
	
		// model much of this from Secretion 
		Molecular(); 
 	
		// we'll set this to replace BioFVM's version		
		std::vector<double> internalized_total_substrates; 

		// for each substrate, a fraction 0 <= f <= 1 of the 
		// total internalized substrate is released back inot
		// the environment at death 
		std::vector<double> fraction_released_at_death; 

		// for each substrate, a fraction 0 <= f <= 1 of the 
		// total internalized substrate is transferred to the  
		// predatory cell when ingested 
		std::vector<double> fraction_transferred_when_ingested; 

		// Model_t *molecular_model;			// for LibSBMLSim (/usr/local/include/sbml/common/sbmlfwd.h)
		
		/* prototyping / beta in 1.5.0 */ 
		// Boolean, Integer, and Double parameters
/*		
		std::vector<bool> bools; 
		std::unordered_map<std::string,int> bool_name_map; 
		std::string& bool_name( int i ); 
		std::vector<std::string> bool_units; 
		void resize_bools( int n ); 
		int add_bool( std::string name , std::string units , bool value ); 
		bool& access_bool( std::string name ); 
		
		std::vector<int> ints; 
		std::unordered_map<std::string,int> int_name_map; 
		std::string& int_name( int i ); 
		std::vector<std::string> int_units; 
		int& access_int( std::string name ); 
		
		std::vector<int> doubles; 
		std::unordered_map<std::string,int> double_name_map; 
		std::string& double_name( int i ); 
		std::vector<std::string> double_units; 
		double& access_double( std::string name ); 
*/
	
		// use this to properly size the secretion parameters to the 
		// microenvironment in molecular.pMicroenvironment. 
		void sync_to_current_microenvironment( void ); // done 
		
//		void advance( Basic_Agent* pCell, Phenotype& phenotype , double dt ); 
		
		// use this to properly size the secretion parameters to the microenvironment in 
		// pMicroenvironment
		void sync_to_microenvironment( Microenvironment* pNew_Microenvironment ); // done 
		
		// use this 
		void sync_to_cell( Basic_Agent* pCell ); 
		
};
#endif  // foobar