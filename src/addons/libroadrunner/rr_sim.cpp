
#include "../../core/PhysiCell_cell.h"

// #include "libsbmlsim/libsbmlsim.h"

// Model_t* my_model;
SBMLDocument_t* sbml_doc;
// Model my_sbml_model;

//extern "C" {
//}

//SBMLSIM_EXPORT myResult* simulateSBMLFromFile(const char* file, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method) {


int read_sbml_file(const char* file)
{
  // SBMLDocument_t* sbml_doc;
  // my_sbml_model = Model();

  // Model_t* m;
  // myResult *rtn;

  unsigned int err_num;
  // double atol = 0.0;
  // double rtol = 0.0;
  //  rf. ~/dev/libsbmlsim-1.4.0/src/lib_main.c
  double atol = ABSOLUTE_ERROR_TOLERANCE;
  double rtol = RELATIVE_ERROR_TOLERANCE;
  // double facmax = 0.0;
    double facmax = DEFAULT_FACMAX;

  std::cout << "------------>  readSBMLFromFile\n";
  sbml_doc = readSBMLFromFile(file);  //rwh2
  if (sbml_doc == NULL) {
//    return create_myResult_with_errorCode(Unknown);
    std::cout << "CRAP: got NULL in read_sbml_file\n";
    return 1;
  }
  // SBMLDocument_free(sbml_doc);

  return 0;
}

//------------------------------------------------------------
void solve_sbml(PhysiCell::Cell* pCell)
{
  extern int oxygen_i, glucose_i;  //, energy_i;   // substrate index
  extern int energy_vi;   // custom var index

  // static double sim_time = 6.0;  // 10.0
  double sim_time = 6.0;  // 10.0
  // static double sbml_dt = 0.001;   // diffusion_dt = 0.01 (default; rf. PhysiCell_constants.h)

  // 9-26-19: from test driver:  steps=300 --> dt=0.000005
  // static double sbml_dt = 0.000005;
  // static double sbml_dt = 0.0001;
  double sbml_dt = 0.0001;

  // static double sbml_dt = 0.01; //0.001;  // diffusion_dt = 0.01 (default; rf. PhysiCell_constants.h)
  // static double sbml_dt = 0.1; //0.001;  // diffusion_dt = 0.01 (default; rf. PhysiCell_constants.h)
  /*
  static int print_interval = 10;
  static int print_amount = 1;
  static int method = MTHD_RUNGE_KUTTA;
  static boolean use_lazy_method = false;
  static double atol = 0.0;
  static double rtol = 0.0;
  static double facmax = 0.0;
  */

  int print_interval = 10;
  int print_amount = 1;
  int method = MTHD_RUNGE_KUTTA;
  boolean use_lazy_method = false;
  double atol = 0.0;
  double rtol = 0.0;
  double facmax = 0.0;

  std::vector<double> density[3];

  std::cout << "   --- sbml_sim.cpp:  ...->getNumSpecies() =" << pCell->phenotype.molecular.molecular_model->getNumSpecies() << std::endl;

  Model_t *mm = pCell->phenotype.molecular.molecular_model;  // arrgh. NOOOOoooooo, don't do this!

  // std::cout << " >>>> sbml_sim.cpp: mm oxygen ic="<<mm->getAttribute.getElementBySId("oxygen") << std::endl;
  //initialConcentration
  double oxy_val;
  double new_oxy_val;
  // int oxygen_i = microenvironment.find_density_index( "oxygen" ); 
  // int glucose_i = microenvironment.find_density_index( "glucose" ); 
  // int energy_i = microenvironment.find_density_index( "energy" ); 

  int vi = microenvironment.nearest_voxel_index(pCell->position);
  // dv = microenv->nearest_density_vector(vi);

#define debug_sbmlsim
#ifdef debug_sbmlsim
  std::cout << " >>>> sbml_sim.cpp: t="<<PhysiCell::PhysiCell_globals.current_time << ", solve_sbml for Cell ID= " << pCell->ID << std::endl;
  std::cout << " >>>> sbml_sim.cpp: Cell x,y pos ="<< pCell->position[0]<<", " << pCell->position[1]<<" " << std::endl;
  std::cout << " >>>> sbml_sim.cpp: vi="<< vi << std::endl;
#endif


  // mm->getElementBySId("Oxygen")->getAttribute("initialConcentration", oxy_val);
  // new_oxy_val = oxy_val + 2.42;
  // mm->getElementBySId("Oxygen")->setAttribute("initialConcentration", new_oxy_val);
  // mm->getElementBySId("Oxygen")->getAttribute("initialConcentration", oxy_val);
  // std::cout << " >>>> sbml_sim.cpp: mm oxygen ic="<< oxy_val << std::endl;

#define debug_sbmlsim2
#ifdef debug_sbmlsim2
  std::cerr << "       pCell->ID = " << pCell->ID << std::endl;
  // std::cerr << "       mm = " << mm << std::endl;
  std::cerr << "       vi = " << vi << std::endl;
  std::cerr << "       oxygen_i = " << oxygen_i << std::endl;
  std::cerr << "       setting oxygen IC = " << microenvironment(vi)[oxygen_i] << std::endl;
  // std::cout << "       setting glucose IC = " << microenvironment(vi)[glucose_i] << std::endl;
#endif

// http://sbml.org/Software/libSBML/5.18.0/docs/cpp-api/class_model.html#a58f9cf402bf4f5e090799270aad1c1cb
// hmm, can we safely access the microenv thread-safe?
//#pragma omp critical
// {
  // if (pCell->phenotype.molecular.molecular_model->getElementBySId("Oxygen") == NULL)
  // if (mm->getElementBySId("Oxygen") == NULL)
  if (false)
  {
    std::cerr << "      ...->getElementBySId('Oxygen') == NULL  !!" << std::endl;
  }
  else 
  {
    mm->getElementBySId("Oxygen")->setAttribute("initialConcentration", 42. );
    std::cerr << "            ...->getElementBySId('Oxygen') != NULL  !!" << std::endl;
    // pCell->phenotype.molecular.molecular_model->getElementBySId("Oxygen")->setAttribute("initialConcentration", 42. );
  }
// }
  // mm->getElementBySId("Oxygen")->setAttribute("initialConcentration", microenvironment(vi)[oxygen_i] );
  // mm->getElementBySId("Glucose")->setAttribute("initialConcentration", microenvironment(vi)[glucose_i] );

  // std::cout << "       setting glucose IC = " << microenvironment(vi)[glucose_i] << std::endl;
  // mm->getElementBySId("Glucose")->setAttribute("initialConcentration", microenvironment(vi)[glucose_i] );

  // std::cout << " >>>> sbml_sim.cpp: t="<<PhysiCell::PhysiCell_globals.current_time << ", solve_sbml for Cell ID= " << pCell->ID << std::endl;

  //                   \/ \/ \/ \/
  // myResult *mm_res = simulateSBMLModel( pCell->phenotype.molecular.molecular_model, sim_time, sbml_dt, print_interval, print_amount, method, use_lazy_method, atol, rtol, facmax);
  myResult *mm_res = simulateSBMLModel( mm, sim_time, sbml_dt, print_interval, print_amount, method, use_lazy_method, atol, rtol, facmax);
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#ifdef debug_sbmlsim
  //  std::cout << "solve_sbml for Cell ID= " << pCell->ID << std::endl;
  std::cout << "--- num_of_rows = " << mm_res->num_of_rows << std::endl;
  std::cout << "--- num_of_columns_sp = " << mm_res->num_of_columns_sp << std::endl;
#endif

  // Energy_0,Glucose_0,Hydrogen_0,Oxygen_0,  Energy_1, etc…   
  // So, to extract *just* the values of Energy, one would get every 4th value (modulo (# of species)), 
  // starting at the 0th (Energy) in the list.
  // std::cout << "values_sp = " <<mm_res->values_sp[0]<<","<<mm_res->values_sp[1]<<","<<mm_res->values_sp[2]<<","<<mm_res->values_sp[3]<<","<<mm_res->values_sp[4]<<","<<mm_res->values_sp[5]<<"," <<std::endl;
  int energy_species_idx  = 0;  // Energy
  int energy_species_idx_first  = 0;  // Energy
  // int g_idx_species = 1;  // Glucose
  // int idx_species = 2;  // Hydrogen
  // int o_idx_species = 3;  // Oxygen

  double o_last, g_last, e_last;
  float tmax = 598.;
  tmax = 29.;
  tmax = 0.;
  // std::cout << "PhysiCell::PhysiCell_globals.current_time = " <<PhysiCell::PhysiCell_globals.current_time<<std::endl;

  //rwh: let's only write results at the end
  // if (PhysiCell::PhysiCell_globals.current_time > tmax)    //rwh: let's only write results at the end
#undef dump_sbmlsim
  int glucose_species_idx = 1;
  int oxygen_species_idx = 3;

// BEWARE: thread-safe file i/o?
#ifdef dump_sbmlsim
  // if (PhysiCell::PhysiCell_globals.current_time >= 0.0 && PhysiCell::PhysiCell_globals.current_time < 1.0)
  if (PhysiCell::PhysiCell_globals.current_time >= 0.0 && PhysiCell::PhysiCell_globals.current_time < 1.0 && pCell->ID==0)
  {  
    std::cout << "solve_sbml for Cell ID= " << pCell->ID << std::endl;
    std::cout << "PhysiCell::PhysiCell_globals.current_time = " <<PhysiCell::PhysiCell_globals.current_time<<std::endl;

    std::ofstream o_file;
    o_file.open ("oxygen.dat");
    std::ofstream g_file;
    g_file.open ("glucose.dat");
    std::ofstream e_file;
    e_file.open ("energy.dat");
    // std::string fname = "species.csv";
    // myfile.open ("species.csv");
    // std::cout << "----------- " << __FILE__ << ":" <<__FUNCTION__ <<":  selected species in " << fname << std::endl;
    for (int idx=0; idx<=mm_res->num_of_rows-1; idx++ ) {
  //    std::cout << mm_res->values_sp[idx]<<", ";
      // myfile << "Writing this to a file.\n";
      // std::cout << mm_res->values_time[idx]<<", "<<mm_res->values_sp[idx_species]<<std::endl;
      o_file << mm_res->values_time[idx]<<", "<<mm_res->values_sp[oxygen_species_idx]<<std::endl;
      g_file << mm_res->values_time[idx]<<", "<<mm_res->values_sp[glucose_species_idx]<<std::endl;
      e_file << mm_res->values_time[idx]<<", "<<mm_res->values_sp[energy_species_idx]<<std::endl;

      o_last = mm_res->values_sp[oxygen_species_idx];
      g_last = mm_res->values_sp[glucose_species_idx];
      e_last = mm_res->values_sp[energy_species_idx];

      oxygen_species_idx += mm_res->num_of_columns_sp;  // increment to get next value for this species
      glucose_species_idx += mm_res->num_of_columns_sp;  // increment to get next value for this species
      energy_species_idx += mm_res->num_of_columns_sp;  // increment to get next value for this species
    }
//    std::cout << std::endl;
    o_file.close ();
    g_file.close ();
    e_last = mm_res->values_sp[energy_species_idx_first + (mm_res->num_of_rows-1) * mm_res->num_of_columns_sp];
    // e_file << "6.01, "<<e_last<<std::endl;
    e_file.close ();
  }
#endif

  // microenvironment(vi)[oxygen_i] = o_last;
  // microenvironment(vi)[glucose_i] = g_last;
  // microenvironment(vi)[energy_i] = e_last;
    // microenvironment(vi)[energy_i] = mm_res->values_sp[energy_species_idx ];  // update energy substrate with value from ODE
    // microenvironment(vi)[energy_i] = mm_res->values_sp[energy_species_idx_first + (mm_res->num_of_rows-1) * mm_res->num_of_columns_sp];   // update energy substrate with value from ODE
    pCell->custom_data[energy_vi] =mm_res->values_sp[energy_species_idx_first + (mm_res->num_of_rows-1) * mm_res->num_of_columns_sp]; 

  //  include/libsbmlsim/myResult.h 

  //   typedef struct myResult{
  //   LibsbmlsimErrorCode error_code;
  //   const char *error_message;
  //   int num_of_rows;
  //   int num_of_columns_sp;
  //   int num_of_columns_param;
  //   int num_of_columns_comp;
  //   const char *column_name_time;
  //   const char **column_name_sp;
  //   const char **column_name_param;
  //   const char **column_name_comp;
  //   double *values_time;
  //   double *values_sp;
  //   double *values_param;
  //   double *values_comp;
  // 	/* new code*/
  //   double* values_time_fordelay;
  //   int num_of_delay_rows;
  // } myResult;
}