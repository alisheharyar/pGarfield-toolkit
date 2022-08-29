#ifndef G_UTILITIES_H
#define G_UTILITIES_H

#include <algorithm>
#include <functional>
#include <string>
#include <sstream>
#include <vector>

namespace Garfield {

inline void ltrim(std::string& line) {
  line.erase(line.begin(), 
    std::find_if(line.begin(), line.end(),
                 [](int ch) {return !std::isspace(ch);}));
}

inline void rtrim(std::string& line) {
  line.erase(
    std::find_if(line.rbegin(), line.rend(), 
                 [](int ch) {return !std::isspace(ch);}).base(), 
    line.end());
}

inline std::vector<std::string> tokenize(const std::string& line) {
  std::vector<std::string> words;
  std::istringstream ss(line);
  for (std::string word; ss >> word;) {
    words.push_back(word);
  }
  return words;
}

}

#endif
