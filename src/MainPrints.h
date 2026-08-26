#include "MainStructs.h"

void PrintMainFlags(const MainFlags* flags, std::ostream& out = std::cout);

bool PrintStarGridCoords(std::vector<std::vector<SpotData>> *StarGrid, std::ostream &out = std::cout);

bool PrintSpotData(std::vector<std::vector<SpotData>> *StarGrid, 
                   std::vector<std::string> print_header, 
                   std::string separator = ",", 
                   std::ostream &out = std::cout);
                   
bool PrintOutput(std::map<std::string, std::string> outputs, 
                 std::vector<std::string> print_header, 
                 std::string separator = ",", 
                 std::ostream &out = std::cout);