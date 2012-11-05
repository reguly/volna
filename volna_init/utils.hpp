#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <sstream>

// Template to convert a string into any type, returning 0 for
// invalid input.
//
// N.B. -- as an alternative, one could use boost lexical_cast.
template <class T>
extern bool parse(T &t, std::string& s)
{
  std::string garbage;
  std::istringstream iss(s);
  bool begins_good = !(iss >> t).fail() ;
  iss.clear();
  iss.seekg(0, std::ios::beg );
  iss >> t >> garbage;
  return begins_good && garbage.empty();
}


#endif
