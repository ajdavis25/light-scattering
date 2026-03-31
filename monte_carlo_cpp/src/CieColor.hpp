// CieColor.hpp
#ifndef CIECOLOR_HPP
#define CIECOLOR_HPP

#include <map>

extern std::map<int,double> cieX; 
extern std::map<int,double> cieY; 
extern std::map<int,double> cieZ;

void loadCIEData(); // fill from 380..780

#endif
