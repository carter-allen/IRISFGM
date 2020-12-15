#include "main.h"

#include <csignal>
#include <algorithm>
#include <iostream>
#include <list>
#include <string>
#include <vector>

// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>

using namespace Rcpp;

extern "C" void my_function_to_handle_aborts(int signal_number) {
  /*Your code goes here. You can output debugging info.
   If you return from this function, and it was called
   because abort() was called, your program will exit or crash anyway
   (with a dialog box on Windows).
   */
  // stop("abort()");
}

//' @backref src/rcpp_irisfgm.cpp
// [[Rcpp::export(.main)]]
int qubic(const CharacterVector& str) {
  // may treat abort() more friendly, see http://stackoverflow.com/a/3911102
  signal(SIGABRT, &my_function_to_handle_aborts);
  try {
    int argc = str.size();
    char **argv = new char*[str.size()];

    for (int i = 0; i < argc; ++i) {
      argv[i] = strdup(as<std::string>(str[i]).c_str());
    }

    int rtn = do_qubic(str.size(), argv);
    for (int i = 0; i < argc; ++i)
      free(argv[i]);
    delete[] argv;
    return rtn;
  }
  catch (double) {
    stop("Something wrong near r_main_d function, maybe out of memory");
  }
  return -1; // avoid warning
}
