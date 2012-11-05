#include "utils.hpp"
#include "physicalParams.hpp"
#include "values.hpp"

void Values::ReadFromFile( std::ifstream &ifs, int n ) {
  
  if (ifs) {
  
    RealType value = 0.;
    
    // string holding the current line
    std::string line = "";
    
    int cnt = 0;
    while ( getline(ifs, line ) ) {
      
      if ( line.substr(0, 2) == std::string("//") ) // ignore comments
	{}

      else {
	std::istringstream isstream(line);
	std::string token;
	getline(isstream, token);
	
	if ( !parse<RealType>(value, token) )
	  std::cerr << "Bad input: " << token << "\n";
	
	switch (n) {
	
	case 0: 
	  (*this).H( cnt ) += value;
	  break;
	  
	case 1:
	  (*this).U( cnt ) += value;
	  break;

	case 2:
	  (*this).V( cnt ) += value;
	  break;

	case 3:
	  (*this).Zb( cnt ) += value;
	  break;
	  
	default:
	  std::cerr << "Invalid value for variable number: "
		    << n << "\n";
	  break;
	}
	
      }

      ++cnt;

    }
    
  }

}
