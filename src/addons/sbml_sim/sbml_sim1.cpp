
#include "../../core/PhysiCell_cell.h"

#include "libsbmlsim/libsbmlsim.h"
// Model_t* my_model;
SBMLDocument_t* my_doc;
// Model my_sbml_model;

//extern "C" {
//}

//SBMLSIM_EXPORT myResult* simulateSBMLFromFile(const char* file, double sim_time, double dt, int print_interval, int print_amount, int method, int use_lazy_method) {


int read_sbml_file(const char* file)
{
  // SBMLDocument_t* my_doc;
  // my_sbml_model = Model();

  // Model_t* m;
  myResult *rtn;

  unsigned int err_num;
  // double atol = 0.0;
  // double rtol = 0.0;
  //  rf. ~/dev/libsbmlsim-1.4.0/src/lib_main.c
  double atol = ABSOLUTE_ERROR_TOLERANCE;
  double rtol = RELATIVE_ERROR_TOLERANCE;
  // double facmax = 0.0;
    double facmax = DEFAULT_FACMAX;

  std::cout << "------------>  readSBMLFromFile\n";
  my_doc = readSBMLFromFile(file);  //rwh2
  if (my_doc == NULL) {
//    return create_myResult_with_errorCode(Unknown);
    std::cout << "CRAP: got NULL in read_sbml_file\n";
    return 1;
  }

#ifdef get_model
  my_model = SBMLDocument_getModel(my_doc);  //rwh2 - needs to go into each cell
  // my_sbml_model = *my_model;
  std::cout << "---> sizeof(my_model) = " << sizeof(my_model) << std::endl;
  std::cout << "---> sizeof(*my_model) = " << sizeof(*my_model) << std::endl<< std::endl;
  std::cout << "---> my_model->getNumCompartments() = " << my_model->getNumCompartments() << std::endl;
  std::cout << "---> my_model->getNumParameters() = " << my_model->getNumParameters() << std::endl;
  std::cout << "---> my_model->getNumSpecies() = " << my_model->getNumSpecies() << std::endl;

// time:10 step:100 dt:0000024.
  double sim_time = 10.0;
  sim_time = 5.0;
  double dt = 0.000024;
  dt = 0.001;
  int print_interval = 10;
  int print_amount = 1;
  int method = MTHD_RUNGE_KUTTA;
  boolean use_lazy_method = false;

   // ---- do the actual simulation -------------
  std::cout << "--------->  Yehaw: calling simulateSBMLModel\n";
  rtn = simulateSBMLModel(my_model, sim_time, dt, print_interval, print_amount, method, use_lazy_method, atol, rtol, facmax);
  if (rtn == NULL) {
    rtn = create_myResult_with_errorCode(SimulationFailed);
    std::cout << "-------> CRAP: got NULL in sbml_sim.cpp: read_sbml_file(), simulateSBMLModel\n";
//    return 1;
  }
  else {
    std::cout << "---> rtn->num_of_rows = " << rtn->num_of_rows << std::endl;
    std::cout << "---> rtn->num_of_columns_sp = " << rtn->num_of_columns_sp << std::endl;
    std::cout << "---> rtn->num_of_columns_param = " << rtn->num_of_columns_param << std::endl;
    std::cout << "---> rtn->num_of_columns_comp = " << rtn->num_of_columns_comp << std::endl;
    for (int idx=0; idx<10; idx++ )
      std::cout << "---> rtn->column_name_time[ " << idx << " = " << rtn->column_name_time[idx] << std::endl;
    for (int idx=0; idx<4; idx++ )  // energy, glucose, hydrogen, oxygen
      std::cout << "---> rtn->column_name_sp[ " << idx << " = " << rtn->column_name_sp[idx] << std::endl;
//    for (int idx=0; idx<10; idx++ )
    for (int idx=0; idx<9; idx++ )
      std::cout << "---> rtn->column_name_param[ " << idx << " = " << rtn->column_name_param[idx] << std::endl;
//    for (int idx=0; idx<500; idx+=10 )
    for (int idx=0; idx<9; idx+=1 )
      std::cout << idx << " --> " << rtn->values_time[idx] << ", " << rtn->values_sp[idx] << std::endl;
//      std::cout << idx << " --> " << rtn->values_time[idx] << ", " << rtn->values_sp[idx] << std::endl;

  double *value_time_p  = rtn->values_time;
  double *value_sp_p    = rtn->values_sp;
  double *value_param_p = rtn->values_param;
  double *value_comp_p  = rtn->values_comp;

/* 
  if ((fp_s = my_fopen(fp_s, file_s, "w")) == NULL) {
    return;
  }
  if ((fp_p = my_fopen(fp_p, file_p, "w")) == NULL) {
    return;
  }
  if ((fp_c = my_fopen(fp_c, file_c, "w")) == NULL ) {
    return;
  }
  */

#ifdef print_species
  char delimiter = ',';
  /*  Species */
  for (int i = 0; i < rtn->num_of_rows; i++) {
    // fprintf(fp_s, "%.16g", *(value_time_p));
    printf( "%.16g", *(value_time_p));
    value_time_p++;
    for (int j = 0; j < rtn->num_of_columns_sp; j++) {
      // printf(fp_s, "%c%.16g", delimiter, *(value_sp_p));
      printf("%c%.16g", delimiter, *(value_sp_p));
      value_sp_p++;
    }
    printf("\n");
  }
#endif //print_species


  /*  Parameters */
  /* 
  value_time_p  = result->values_time;
  for (i = 0; i < result->num_of_rows; i++) {
    fprintf(fp_p, "%.16g", *(value_time_p));
    value_time_p++;
    for (j = 0; j < result->num_of_columns_param; j++) {
      fprintf(fp_p, "%c%.16g", delimiter, *(value_param_p));
      value_param_p++;
    }
    fprintf(fp_p, "\n");
  }
  */

  /*  Compartments */
  /* 
  value_time_p  = result->values_time;
  for (i = 0; i < result->num_of_rows; i++) {
    fprintf(fp_c, "%.16g", *(value_time_p));
    value_time_p++;
    for (j = 0; j < result->num_of_columns_comp; j++) {
      fprintf(fp_c, "%c%.16g", delimiter, *(value_comp_p));
      value_comp_p++;
    }
    fprintf(fp_c, "\n");
  }
  */

  }
#endif // get_model
  // SBMLDocument_free(my_doc);

  return 0;
}

void solve_sbml(PhysiCell::Cell* pCell)
{
  static double sim_time = 10.0;
  static double sbml_dt = 0.001;   // diffusion_dt = 0.01 (default; rf. PhysiCell_constants.h)
  // static double sbml_dt = 0.01; //0.001;  // diffusion_dt = 0.01 (default; rf. PhysiCell_constants.h)
  // static double sbml_dt = 0.1; //0.001;  // diffusion_dt = 0.01 (default; rf. PhysiCell_constants.h)
  static int print_interval = 10;
  static int print_amount = 1;
  static int method = MTHD_RUNGE_KUTTA;
  static boolean use_lazy_method = false;
  static double atol = 0.0;
  static double rtol = 0.0;
  static double facmax = 0.0;

  std::vector<double> density[3];
  
  Model_t *mm = pCell->phenotype.molecular.molecular_model;
  // std::cout << " >>>> sbml_sim.cpp: mm oxygen ic="<<mm->getAttribute.getElementBySId("oxygen") << std::endl;
  //initialConcentration
  double oxy_val;
  double new_oxy_val;
  int oxygen_i = microenvironment.find_density_index( "oxygen" ); 
	int glucose_i = microenvironment.find_density_index( "glucose" ); 

  int vi = microenvironment.nearest_voxel_index(pCell->position);
  // dv = microenv->nearest_density_vector(vi);
  std::cout << " >>>> sbml_sim.cpp: vi="<< vi << std::endl;
  microenvironment(vi)[oxygen_i] = 
  microenvironment(vi)[glucose_i] = 


  mm->getElementBySId("Oxygen")->getAttribute("initialConcentration", oxy_val);
  new_oxy_val = oxy_val + 2.42;
  mm->getElementBySId("Oxygen")->setAttribute("initialConcentration", new_oxy_val);
  mm->getElementBySId("Oxygen")->getAttribute("initialConcentration", oxy_val);
  std::cout << " >>>> sbml_sim.cpp: mm oxygen ic="<< oxy_val << std::endl;

  std::cout << " >>>> sbml_sim.cpp: t="<<PhysiCell::PhysiCell_globals.current_time << ", solve_sbml for Cell ID= " << pCell->ID << std::endl;
  myResult *mm_res = simulateSBMLModel( pCell->phenotype.molecular.molecular_model, sim_time, sbml_dt, print_interval, print_amount, method, use_lazy_method, atol, rtol, facmax);
//  std::cout << "solve_sbml for Cell ID= " << pCell->ID << std::endl;
  std::cout << "--- num_of_rows = " << mm_res->num_of_rows << std::endl;
  std::cout << "--- num_of_columns_sp = " << mm_res->num_of_columns_sp << std::endl;

  // Energy_0,Glucose_0,Hydrogen_0,Oxygen_0,  Energy_1, etc…   
  // So, to extract *just* the values of Energy, one would get every 4th value (modulo (# of species)), 
  // starting at the 0th (Energy) in the list.
  std::cout << "values_sp = " <<mm_res->values_sp[0]<<","<<mm_res->values_sp[1]<<","<<mm_res->values_sp[2]<<","<<mm_res->values_sp[3]<<","<<mm_res->values_sp[4]<<","<<mm_res->values_sp[5]<<"," <<std::endl;
  // int idx_species = 0;  // Energy
  // int idx_species = 1;  // Glucose
  // int idx_species = 2;  // Hydrogen
  int idx_species = 3;  // Oxygen

  float tmax = 598.;
  tmax = 0.;
  std::cout << "PhysiCell::PhysiCell_globals.current_time = " <<PhysiCell::PhysiCell_globals.current_time<<std::endl;
  if (PhysiCell::PhysiCell_globals.current_time > tmax)    //rwh: let's only write results at the end
  {  
    std::ofstream myfile;
    std::string fname = "species.csv";
    // myfile.open ("species.csv");
    myfile.open (fname);
    std::cout << "----------- " << __FILE__ << ":" <<__FUNCTION__ <<":  selected species in " << fname << std::endl;
    for (int idx=0; idx<=mm_res->num_of_rows-1; idx++ ) {
  //    std::cout << mm_res->values_sp[idx]<<", ";
      // myfile << "Writing this to a file.\n";
      // std::cout << mm_res->values_time[idx]<<", "<<mm_res->values_sp[idx_species]<<std::endl;
      myfile << mm_res->values_time[idx]<<", "<<mm_res->values_sp[idx_species]<<std::endl;
      idx_species += mm_res->num_of_columns_sp;  // increment to get next value for this species
    }
    std::cout << std::endl;
    myfile.close ();
  }

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