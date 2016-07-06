#include "slmm_util.hpp"

#include <cctype>
#include <sys/time.h>
#include <sstream>

namespace slmm {

double wall_time () {
  static const double us = 1.0e6;
  timeval t;
  gettimeofday(&t, 0);
  return (t.tv_sec*us + t.tv_usec)/us;
}

std::string& tolower (std::string& s) {
  for (auto& c: s)
    c = std::tolower(c);
  return s;
}

std::string format_strings_as_list (const char** strings, const Size n) {
  std::stringstream ss;
  ss << "{";
  for (Size i = 0; i < n-1; ++i) ss << strings[i] << ", ";
  ss << strings[n-1] << "}";
  return ss.str();
}

} // namespace slmm
