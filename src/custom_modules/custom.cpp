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


#include "./custom.h"

int oxygen_i, glucose_i, energy_i; 
int energy_vi; 

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

	oxygen_i = microenvironment.find_density_index( "oxygen" ); 
	glucose_i = microenvironment.find_density_index( "glucose" ); 
  	energy_i = microenvironment.find_density_index( "energy" ); 
	std::cout << "---------- setup_microenv\n";
	std::cout << "    oxygen_i = " << oxygen_i << std::endl;
	std::cout << "    glucose_i = " << glucose_i << std::endl;
	std::cout << "    energy_i = " << energy_i << std::endl;

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

	std::ifstream in( "oxy_irreg.dat");

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
		 microenvironment(ivox)[oxygen_i] = v;
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

void setup_microenvironment1( void )
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

	oxygen_i = microenvironment.find_density_index( "oxygen" ); 
	glucose_i = microenvironment.find_density_index( "glucose" ); 
  	energy_i = microenvironment.find_density_index( "energy" ); 
	std::cout << "---------- setup_microenv\n";
	std::cout << "    oxygen_i = " << oxygen_i << std::endl;
	std::cout << "    glucose_i = " << glucose_i << std::endl;
	std::cout << "    energy_i = " << energy_i << std::endl;

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
	for( int n=0; n < microenvironment.number_of_voxels(); n++ )
	{
		// x coordinate of the nth voxel's center
		x = microenvironment.mesh.voxels[n].center[0];
		for( int iregion=1; iregion <= nregions; iregion++ )
		{
			if (x < (xmin + iregion*xdel))
			{
				microenvironment(n)[oxygen_i] = oxy - (iregion-1)*oxy_del;
				microenvironment(n)[glucose_i] = glu - (iregion-1)*glu_del;
				break;
			}
			// oxy -= iregion-5;
			// glu -= iregion-2;
		}
	}
	
	return; 
}

void create_cell_types( void )
{
	// housekeeping 
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	cell_defaults.type = 0; 
	cell_defaults.name = "cell"; 

	cell_defaults.custom_data.add_variable( "energy" , "dimensionless", 1.0 ); 
	energy_vi = cell_defaults.custom_data.find_variable_index( "energy" ); 
	cell_defaults.functions.update_phenotype = energy_based_cell_phenotype; 
	
	// set default cell cycle model 
	// cell_defaults.functions.cycle_model = live; 
	
	// // set default_cell_functions; 
	// cell_defaults.functions.update_phenotype = NULL; 
	
	// // needed for a 2-D simulation: 
	// cell_defaults.functions.set_orientation = up_orientation; 
	// cell_defaults.phenotype.geometry.polarity = 1.0;
	// cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// // make sure the defaults are self-consistent (do this before following settings!)
	// cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	// cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 


	// // Set phenotype as you want

	// // turn off birth
	int start_index = live.find_phase_index( PhysiCell_constants::live );
	int end_index = live.find_phase_index( PhysiCell_constants::live );
	cell_defaults.phenotype.cycle.data.transition_rate(start_index,end_index) = 0.0; 

	// // turn off death 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	cell_defaults.phenotype.death.rates[ apoptosis_model_index ] = 0.0; 
	cell_defaults.phenotype.death.rates[ necrosis_model_index ] = 0.0; 

	// // turn off adhesion and repulsion 
	cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = 0.0; 
	cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength = 0.0; 
	
	// // set motilty parameters 
	cell_defaults.phenotype.motility.is_motile = true;
	cell_defaults.phenotype.motility.persistence_time = parameters.doubles("persistence_time"); 
	cell_defaults.phenotype.motility.migration_speed = parameters.doubles("migration_speed"); 
	// cell_defaults.phenotype.motility.migration_speed = 1.0; 
	cell_defaults.phenotype.motility.migration_bias_direction = { 1.0, 0.0, 0.0 };  
	cell_defaults.phenotype.motility.migration_bias = parameters.doubles("migration_bias"); 
	// cell_defaults.phenotype.motility.migration_bias = 0.0; 

	return; 
}

void setup_tissue( void )
{
	// extern SBMLDocument_t *sbml_doc;
	rrc::RRVectorPtr vptr;
    rrc::RRCDataPtr result;  // start time, end time, and number of points

	Cell* pC;
	float xval = -600.0;
	float yval = 1000.0;
	// create just 3 cells, equally spaced in y; they'll migrate left-to-right
	for (int idx=0; idx<3; idx++) {
		pC = create_cell(); 
		yval -= 500.;
		pC->assign_position( xval, yval, 0.0 ); 
		pC->set_total_volume( pC->get_total_volume() * 3.0); 

		// Model_t *mm = SBMLDocument_getModel(sbml_doc);
		// std::cout << "mm =" << mm << std::endl;
		// pC->phenotype.molecular.molecular_model = mm;  // assign the intracellular model to each cell

		std::cerr << "------------->>>>>  Creating rrHandle, loadSBML file\n\n";
		rrc::RRHandle rrHandle = createRRInstance();
		if (!rrc::loadSBML (rrHandle, "../Toy_Model_for_PhysiCell.xml")) {
			std::cerr << "------------->>>>>  Error while loading SBML file  <-------------\n\n";
		// 	printf ("Error message: %s\n", getLastError());
		// 	getchar ();
		// 	exit (0);
		}
		pC->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell
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
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector<std::string> output( 4 , "white" );   // "gray"
	// output[0] = "cyan";  // cytoplasm
	// output[2] = "red";  // nucleus
	
	return output; 
}

// cell_defaults.functions.update_phenotype = energy_based_cell_phenotype; 
void energy_based_cell_phenotype(Cell* pCell, Phenotype& phenotype , double dt)
{
	static int idx_oxygen = 3;
	rrc::RRVectorPtr vptr;
	rrc::RRCDataPtr result;  // start time, end time, and number of points

	// pC->phenotype.molecular.model_rr = rrHandle;  // assign the intracellular model to each cell
	vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	std::cout << "--- before updating:" << std::endl;
	for (int idx=0; idx<vptr->Count; idx++)
		std::cout << idx << ", " << vptr->Data[idx] << std::endl;

	// vptr->Data[idx_oxygen] += 0.1;
	// rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);

	vptr = rrc::getFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr);
	// std::cout << vptr->Count << std::endl;
	std::cout << "--- after updating oxygen:" << std::endl;
	for (int idx=0; idx<vptr->Count; idx++)
		std::cout << idx << ", " << vptr->Data[idx] << std::endl;

	int oxygen_i = microenvironment.find_density_index( "oxygen" ); 
	int glucose_i = microenvironment.find_density_index( "glucose" ); 
	int energy_i = microenvironment.find_density_index( "energy" ); 
	int vi = microenvironment.nearest_voxel_index(pCell->position);
	double oxy_val = microenvironment(vi)[oxygen_i];
	std::cout << "oxy_val at voxel of cell = " << oxy_val << std::endl;

	vptr->Data[idx_oxygen] = oxy_val;
	rrc::setFloatingSpeciesConcentrations(pCell->phenotype.molecular.model_rr, vptr);

	result = rrc::simulateEx (pCell->phenotype.molecular.model_rr, 0, 10, 10);  // start time, end time, and number of points
	int index = 0;
	// Print out column headers... typically time and species.
	for (int col = 0; col < result->CSize; col++)
	{
		std::cout << result->ColumnHeaders[index++];
		if (col < result->CSize - 1)
		{
			std::cout << "\t";
		}
	}
	std::cout << "\n";

	index = 0;
	// Print out the data
	for (int row = 0; row < result->RSize; row++)
	{
		for (int col = 0; col < result->CSize; col++)
		{
			std::cout << result->Data[index++];
			if (col < result->CSize -1)
			{
				std::cout << "\t";
			}
		}
		std::cout << "\n";
	}
}

std::vector<std::string> energy_coloring_function( Cell* pCell )
{
	// color 0: cytoplasm fill 
	// color 1: outer outline 
	// color 2: nuclear fill 
	// color 3: nuclear outline 
	
	// some bookkeeping 
	// static int energy_i = pCell->custom_data.find_variable_index( "energy" ); 
	
	std::vector< std::string > output( 4, "white" ); 

	std::cout << "--- coloring fn: cell ID, energy = " << pCell->ID <<", "<< pCell->custom_data[energy_vi] << std::endl; 
	if (pCell->custom_data[energy_vi] > 6.0)
		output[0] = "rgb(0,255,0)";
	else if (pCell->custom_data[energy_vi] > 4.0)
		output[0] = "rgb(255,0,0)";
	else if (pCell->custom_data[energy_vi] > 3.0)
		output[0] = "rgb(255,255,255)";
	else if (pCell->custom_data[energy_vi] > 2.0)
		output[0] = "rgb(255,255,0)";
	else if (pCell->custom_data[energy_vi] > 1.0)
		output[0] = "rgb(0,255,255)";
	else 
		output[0] = "rgb(0,0,0)";
/*
	if( pCell->phenotype.death.dead == false )
	{
		int red = (int) round( 255.0 * pCell->custom_data[alpha_i] ); 
		int green = (int) round( 255.0 * pCell->custom_data[beta_i] ); 
		int blue = (int) round( 255.0 * pCell->custom_data[resistance_i] ); 
		int grey = (int) round( 255.0 * pCell->custom_data[energy_i] ); 
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", red, green, blue );
		output[2].assign( szTempString ); // nucleus by alpha, beta, resistance "genes"
		sprintf( szTempString , "rgb(%u,%u,%u)", grey, grey, grey );
		output[0].assign( szTempString ); // cyto by energy 
		
		return output; 
	}

	// if not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic ) 
		// Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown // switch to black???
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	
	*/
	
	return output; 
}