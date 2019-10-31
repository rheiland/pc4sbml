/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./heterogeneity.h"
#include "../modules/PhysiCell_settings.h"

// These are for C
// #define STATIC_RRC
#include "rrc_api.h"
#include "rrc_types.h"
// #include "rrc_utilities.h"
extern "C" rrc::RRHandle createRRInstance();


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
// #include <vector>
#include <string>

int oxygen_ID, glucose_ID, energy_ID; 
int energy_vi; 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints( "random_seed" ) ); 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// turn the default cycle model to live, 
	// so it's easier to turn off proliferation
	
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	// Make sure we're ready for 2D
	
	cell_defaults.functions.set_orientation = up_orientation;  
	
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// use default proliferation and death 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	
	cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	cell_defaults.parameters.o2_reference = 38.0; 
	
	// set default uptake and secretion 
	
	// oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0
	
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_ID] = 0; 
	// cell_defaults.phenotype.secretion.uptake_rates[oxygen_ID] = 10; 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_ID] = 0.01; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_ID] = 38; 

	// set the default cell type to no phenotype updates 
	
	cell_defaults.functions.update_phenotype = energy_based_cell_phenotype;   // rwh 

	cell_defaults.functions.update_migration_bias = cell_chemotaxis; 
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 
	
	cell_defaults.name = "cancer cell"; 
	cell_defaults.type = 0; 
	
	// add custom data 
	
	// cell_defaults.custom_data.add_variable( "alpha" , "dimensionless", 1.0 ); 
	// cell_defaults.custom_data.add_variable( "beta" , "dimensionless", 1.0 ); 
	// cell_defaults.custom_data.add_variable( "resistance" , "dimensionless", 1.0 ); 
	// cell_defaults.custom_data.add_variable( "use_rate" , "dimensionless", 1.0 ); 
	cell_defaults.custom_data.add_variable( "energy" , "dimensionless", 1.0 ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// domain parameters read from XML config file

	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	// initialize BioFVM 
	initialize_microenvironment(); 	

	oxygen_ID = microenvironment.find_density_index( "oxygen" ); 
	glucose_ID = microenvironment.find_density_index( "glucose" ); 
  	energy_ID = microenvironment.find_density_index( "energy" ); 
	std::cout << "---------- setup_microenv\n";
	std::cout << "    oxygen_ID = " << oxygen_ID << std::endl;
	std::cout << "    glucose_ID = " << glucose_ID << std::endl;
	std::cout << "    energy_ID = " << energy_ID << std::endl;

	double oxy = 38.0;  // IC
	double oxy_del = 9.0;
	double glu = 32.0;  // IC
	double glu_del = 7.5;
	double x;
	double xmin = -750.;
	int nregions = 5;
	double xdel = 1500./nregions;
	// 5 regions across x: [-750:-450:-150:150:450:750]
	std::cout << "setup_microenvironment: num voxels= " << microenvironment.number_of_voxels() << std::endl;
	double vmin = 100.;
	double vmax = -100.;
	double v;
	// ifstream infile( "oxy_irreg.dat" );

	std::ifstream in( "../data/oxy_irreg.dat");

    std::string line;
	int ivox = 0;
    while ( getline( in, line ) )                   // read a whole line of the file
    {
      std::stringstream ss( line );                     // put it in a stringstream (internal stream)
      std::string data;
      while ( getline( ss, data, ' ' ) )           // read (string) items up to a comma
      {
         v = stod( data );            // use stod() to convert to double
		 if (v < vmin) vmin = v;
		 if (v > vmax) vmax = v;
		 microenvironment(ivox)[oxygen_ID] = v;
		 ivox++;
      }
    }
	std::cout << "setup_microenvironment: oxy vmin,vmax= " << vmin << ", " << vmax << std::endl;


	// for( int n=0; n < microenvironment.number_of_voxels(); n++ )
	// {
	// 	// microenvironment(n)[oxygen_i] = oxy - (iregion-1)*oxy_del;
	// 	microenvironment(n)[oxygen_i] = n;
	// 	// microenvironment(n)[glucose_i] = glu - (iregion-1)*glu_del;
	// }
	
	return; 
}

void setup_tissue( void )
{
	// extern SBMLDocument_t *sbml_doc;
	rrc::RRVectorPtr vptr;
    rrc::RRCDataPtr result;  // start time, end time, and number of points

	// place a cluster of tumor cells at the center 
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = parameters.doubles( "tumor_radius" ); // 250.0; 
	// tumor_radius = 28.; // 250.0; 
	
	// Parameter<double> temp; 
	
	// int i = parameters.doubles.find_index( "tumor_radius" ); 
	
	Cell* pCell = NULL; 
	
	double x = 0.0; 
	double x_outer = tumor_radius; 
	double y = 0.0; 
	
	// some bookkeeping 
	energy_vi = cell_defaults.custom_data.find_variable_index( "energy" ); 

	// static int energy_i = cell_defaults.custom_data.find_variable_index( "energy" ); 
	// static int alpha_i = cell_defaults.custom_data.find_variable_index( "alpha" ); 
	// static int beta_i = cell_defaults.custom_data.find_variable_index( "beta" ); 
	// static int resistance_i = cell_defaults.custom_data.find_variable_index("resistance"); 
	// static int use_rate_i = cell_defaults.custom_data.find_variable_index( "use_rate" ); 
	
	int nrows = 0; 
	bool first_time = true;

	while( y < tumor_radius )
	{
		// if (nrows > 0) break;    //rwh
		x = 0.0; 
		if( nrows % 2 == 1 )
		{ x = 0.5*cell_spacing; }
		x_outer = sqrt( tumor_radius*tumor_radius - y*y ); 
		
		while( x < x_outer )
		{
			pCell = create_cell(); // tumor cell 
			pCell->assign_position( x , y , 0.0 );

			// std::cerr << "------------->>>>>  Creating rrHandle, loadSBML file\n\n";
			rrc::RRHandle rrHandle = createRRInstance();
			if (!rrc::loadSBML (rrHandle, "../Toy_Model_for_PhysiCell.xml")) {
				std::cerr << "------------->>>>>  Error while loading SBML file  <-------------\n\n";
			// 	printf ("Error message: %s\n", getLastError());
			// 	getchar ();
			// 	exit (0);
			}
			pCell->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell
			if (first_time)
			{
				first_time = false;
				int r = rrc::getNumberOfReactions(rrHandle);
				int m = rrc::getNumberOfFloatingSpecies(rrHandle);
				int b = rrc::getNumberOfBoundarySpecies(rrHandle);
				int p = rrc::getNumberOfGlobalParameters(rrHandle);
				int c = rrc::getNumberOfCompartments(rrHandle);

				std::cerr << "Number of reactions = " << r << std::endl;
				std::cerr << "Number of floating species = " << m << std::endl;  // 4
				std::cerr << "Number of boundary species = " << b << std::endl;  // 0
				std::cerr << "Number of compartments = " << c << std::endl;  // 1

				std::cerr << "Floating species names:\n";
				std::cerr << "-----------------------\n";
				std::cerr << stringArrayToString(rrc::getFloatingSpeciesIds(rrHandle)) <<"\n"<< std::endl;

				vptr = rrc::getFloatingSpeciesConcentrations(rrHandle);
				std::cerr << vptr->Count << std::endl;
				for (int kdx=0; kdx<vptr->Count; kdx++)
					std::cerr << kdx << ") " << vptr->Data[kdx] << std::endl;
			}
			
			//------------------------
			// pCell->custom_data[alpha_i] = UniformRandom(); 
			// pCell->custom_data[beta_i] = UniformRandom(); 
			// pCell->custom_data[resistance_i] = UniformRandom(); 

			
			if( fabs( y ) > 0.01 )
			// if( fabs( y ) > 0.01 )
			{
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( x , -y , 0.0 );
	
				// std::cerr << "------------->>>>>  Creating rrHandle, loadSBML file\n\n";
				rrc::RRHandle rrHandle = createRRInstance();
				if (!rrc::loadSBML (rrHandle, "../Toy_Model_for_PhysiCell.xml")) {
					std::cerr << "------------->>>>>  Error while loading SBML file  <-------------\n\n";
				// 	printf ("Error message: %s\n", getLastError());
				// 	getchar ();
				// 	exit (0);
				}
				pCell->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell
				
				// pCell->custom_data[alpha_i] = UniformRandom(); 
				// pCell->custom_data[beta_i] = UniformRandom(); 
				// pCell->custom_data[resistance_i] = UniformRandom(); 
			}
			
			if( fabs( x ) > 0.01 )
			{ 
				pCell = create_cell(); // tumor cell 
				pCell->assign_position( -x , y , 0.0 );

				// std::cerr << "------------->>>>>  Creating rrHandle, loadSBML file\n\n";
				rrc::RRHandle rrHandle = createRRInstance();
				if (!rrc::loadSBML (rrHandle, "../Toy_Model_for_PhysiCell.xml")) {
					std::cerr << "------------->>>>>  Error while loading SBML file  <-------------\n\n";
				// 	printf ("Error message: %s\n", getLastError());
				// 	getchar ();
				// 	exit (0);
				}
				pCell->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell

				// std::cout << ">>>>>>>>>>>>>>>>>>>>>>\n";
				// std::cout << "mm->getNumCompartments() =" << mm->getNumCompartments() << std::endl;
				// std::cout << "mm->getNumParameters() =" << mm->getNumParameters() << std::endl;
				// std::cout << "mm->getNumSpecies() =" << mm->getNumSpecies() << std::endl;
				// std::cout << "mm->getNumReactions() =" << mm->getNumReactions() << std::endl;
				// std::cout << ">>>>>>>>>>>>>>>>>>>>>>\n\n";

				// pCell->custom_data[alpha_i] = UniformRandom(); 
				// pCell->custom_data[beta_i] = UniformRandom(); 
				// pCell->custom_data[resistance_i] = UniformRandom(); 
		
				// if( fabs( y ) > 0.01 )
				if( fabs( y ) > 0.01 )
				{
					pCell = create_cell(); // tumor cell 
					pCell->assign_position( -x , -y , 0.0 );

					// std::cerr << "------------->>>>>  Creating rrHandle, loadSBML file\n\n";
					rrc::RRHandle rrHandle = createRRInstance();
					if (!rrc::loadSBML (rrHandle, "../Toy_Model_for_PhysiCell.xml")) {
						std::cerr << "------------->>>>>  Error while loading SBML file  <-------------\n\n";
					// 	printf ("Error message: %s\n", getLastError());
					// 	getchar ();
					// 	exit (0);
					}
					pCell->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell

					// int r = rrc::getNumberOfReactions(rrHandle);
					// int m = rrc::getNumberOfFloatingSpecies(rrHandle);
					// int b = rrc::getNumberOfBoundarySpecies(rrHandle);
					// int p = rrc::getNumberOfGlobalParameters(rrHandle);
					// int c = rrc::getNumberOfCompartments(rrHandle);

					// std::cerr << "Number of reactions = " << r << std::endl;
					// std::cerr << "Number of floating species = " << m << std::endl;  // 4
					// std::cerr << "Number of boundary species = " << b << std::endl;  // 0
					// std::cerr << "Number of compartments = " << c << std::endl;  // 1

					// std::cerr << "Floating species names:\n";
					// std::cerr << "-----------------------\n";
					// std::cerr << stringArrayToString(rrc::getFloatingSpeciesIds(rrHandle)) <<"\n"<< std::endl;

					// vptr = rrc::getFloatingSpeciesConcentrations(rrHandle);
					// std::cerr << vptr->Count << std::endl;
					// for (int kdx=0; kdx<vptr->Count; kdx++)
					// 	std::cerr << kdx << ") " << vptr->Data[kdx] << std::endl;


					// pCell->custom_data[alpha_i] = UniformRandom(); 
					// pCell->custom_data[beta_i] = UniformRandom(); 
					// pCell->custom_data[resistance_i] = UniformRandom(); 
				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		nrows++; 
	}
	
	return; 
}

void cell_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
{
	// static double bias = parameters.doubles("macrophage_migration_bias");
	// static double bias = 0.7;
	double bias = 0.7;   // avoid OpenMP "static" declaration race condition?
	// static int oxygen_i = microenvironment.find_density_index( "oxygen" ); 
	
	phenotype.motility.migration_bias = bias; 
	
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient( oxygen_ID ); 
	double denominator =  norm( phenotype.motility.migration_bias_direction ) + 1e-17; 
	
	phenotype.motility.migration_bias_direction /= denominator; 
	
	return; 
}

std::vector<std::string> energy_coloring_function( Cell* pCell )
{
	// color 0: cytoplasm fill 
	// color 1: outer outline 
	// color 2: nuclear fill 
	// color 3: nuclear outline 
	
	// some bookkeeping 
	// static int energy_i = pCell->custom_data.find_variable_index( "energy" ); 
	
	// start black 
	std::vector< std::string > output( 4, "white" ); 

	// std::cout << "--- energy_coloring_function: cell ID, energy = " << pCell->ID <<", "<< pCell->custom_data[energy_vi] << std::endl; 

	double interp = pCell->custom_data[energy_vi];
	if (interp > 1.0)
		interp = 1.0;

	int Red   = (int) floor( 255.0*interp ); 
	int Green = (int) floor( 255.0*interp ); 
	// int Blue  = (int) floor( 255.0 *(1-interp) ); 
	int Blue  = (int) floor( 255.0 *(1-interp) ); 
	
	char szTempString [128];
	sprintf( szTempString , "rgb(%u,%u,%u)", Red, Green, Blue );
	output[0].assign( szTempString );
	output[2].assign( szTempString );
	

	// if (pCell->custom_data[energy_vi] > 1.0)
	// 	output[0] = "rgb(0,255,0)";
	// else if (pCell->custom_data[energy_vi] > 0.8)
	// 	output[0] = "rgb(255,0,0)";
	// else if (pCell->custom_data[energy_vi] > 0.6)
	// 	output[0] = "rgb(255,255,255)";
	// else if (pCell->custom_data[energy_vi] > 0.4)
	// 	output[0] = "rgb(255,255,0)";
	// else if (pCell->custom_data[energy_vi] > 0.2)
	// 	output[0] = "rgb(0,255,255)";
	// else 
	// 	output[0] = "rgb(0,0,0)";
	
	return output; 
}

// cell_defaults.functions.update_phenotype = energy_based_cell_phenotype; 
void energy_based_cell_phenotype(Cell* pCell, Phenotype& phenotype , double dt)
{
	#pragma omp critical
	{
	// static int idx_glucose = 1;
	// static int idx_oxygen = 3;
	// #pragma omp threadprivate(idx_glucose)   // alternatively, try this??
	int idx_glucose = 1;
	int idx_oxygen = 3;
	rrc::RRVectorPtr vptr;
	rrc::RRCDataPtr result;  // start time, end time, and number of points

	// pC->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell
	vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	// std::cout << "energy_based_cell_phenotype: --- before updating:" << std::endl;
	// for (int idx=0; idx<vptr->Count; idx++)
	// 	std::cout << idx << ", " << vptr->Data[idx] << std::endl;

	// vptr->Data[idx_oxygen] += 0.1;
	// rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);

	// vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	// // std::cout << vptr->Count << std::endl;
	// std::cout << "--- after updating oxygen:" << std::endl;
	// for (int idx=0; idx<vptr->Count; idx++)
	// 	std::cout << idx << ", " << vptr->Data[idx] << std::endl;

//	int oxygen_i = microenvironment.find_density_index( "oxygen" ); 
//	int glucose_i = microenvironment.find_density_index( "glucose" ); 
	// int energy_i = microenvironment.find_density_index( "energy" ); 

	int vi = microenvironment.nearest_voxel_index(pCell->position);
	double oxy_val = microenvironment(vi)[oxygen_ID];
	double glucose_val = microenvironment(vi)[glucose_ID];
	// std::cout << "oxy_val at voxel of cell = " << oxy_val << std::endl;
	// std::cout << "glucose_val at voxel of cell = " << glucose_val << std::endl;

	vptr->Data[idx_oxygen] = oxy_val;
	vptr->Data[idx_glucose] = glucose_val;
	rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);

	result = rrc::simulateEx (pCell->phenotype.molecular.model_rr, 0, 10, 10);  // start time, end time, and number of points
	int index = 0;
	// Print out column headers... typically time and species.
	// for (int col = 0; col < result->CSize; col++)
	// {
	// 	std::cout << result->ColumnHeaders[index++];
	// 	if (col < result->CSize - 1)
	// 	{
	// 		std::cout << "\t";
	// 	}
	// }
	// std::cout << "\n";

	// index = 0;
	// // // Print out the data
	// for (int row = 0; row < result->RSize; row++)
	// {
	// 	for (int col = 0; col < result->CSize; col++)
	// 	{
	// 		std::cout << result->Data[index++];
	// 		if (col < result->CSize -1)
	// 		{
	// 			std::cout << "\t";
	// 		}
	// 	}
	// 	std::cout << "\n";
	// }

	int last_row_idx = result->CSize * (result->RSize - 1);
	// std::cout << "\n-----> Final: t=" << result->Data[last_row_idx] <<", Energy= " <<result->Data[last_row_idx+1] <<  std::endl;

	pCell->custom_data[energy_vi] = result->Data[last_row_idx+1];

	microenvironment(vi)[energy_ID] = result->Data[last_row_idx+1];
	}
}
