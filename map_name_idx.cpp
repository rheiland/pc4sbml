// clang++ -std=c++11 map_name_idx.cpp -o map_name_idx
#include <ctime>
#include <algorithm>
#include <limits>
#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <random>

// using namespace std;

int main()
{
	std::map<std::string, int> species_map;
	std::vector<std::string> species_names;
	// std::stringstream sstm;

	// Build the <string, int> map: celltype+sbml_species --> libroadrunner index
	int idx_rr_species = 0;
	for(int idx_celltype = 0; idx_celltype < 3; ++idx_celltype)
	{
		std::string celltype = "celltype" + std::to_string(idx_celltype);
		for(int idx_species = 0; idx_species < 6; ++idx_species)
		{
			std::string sbml_species = "species" + std::to_string(idx_species);
			std::string key = celltype + "_" + sbml_species;
			species_names.push_back(key);
			species_map[key] = idx_rr_species;
			std::cout << key << " = " << idx_rr_species << std::endl;
			idx_rr_species++;
		}
	}
	// for (auto i: species_names)
	// 	std::cout << i << ' ';

	// Timer for a zillion finds
	std::clock_t start_time = clock();

	std::random_device random_device;
	std::mt19937 engine{random_device()};
	std::uniform_int_distribution<int> dist(0, species_names.size() - 1);

	int max_finds = 1000000;
	std::cout << "----- test time to find " << max_finds << " items in map..." << std::endl;
	for (int idx = 0; idx < max_finds; idx++)
	{
		std::string species = species_names[dist(engine)];
		// idx_rr = species_map.find("celltype1_species2")->second;
		int idx_rr = species_map.find(species)->second;
		// std::cout << idx << std::endl;
		if (idx < 5)
		{
			std::cout << species << " --> " << idx_rr << std::endl;
		}
	}
	std::clock_t stop_time = clock();
	std::cout << "The time it took to build map<string, int> is "
		<< double(stop_time - start_time) / CLOCKS_PER_SEC << " seconds" << '\n';

	// std::cout << "Please enter a character to exit:" << "\n";
	// char ch = 0;
	// std::cin >> ch;

	return 0;
}
