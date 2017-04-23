#include <cstdio>
#include <cstdlib>

#include <iostream>

#include "simulate.hh"

/**
   Finds the argument for a command line option in a stream of characters.
   Simple command line processing functions taken from:
   http://stackoverflow.com/questions/865668/parse-command-line-arguments

   @param begin[in] Start of character stream being read from
   @param end[in] End of character stream being read from
   @param option[in] Command line option to search for

   @return Pointer to argument of command line option, or null if not found
*/

char* getCmdOption(char ** begin, char ** end, const std::string& option)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) return *itr;
  return 0;
}

/**
   Looks for a command line option in a stream of characters.
   Simple command line processing functions taken from:
   http://stackoverflow.com/questions/865668/parse-command-line-arguments

   @param begin[in] Start of character stream being read from
   @param end[in] End of character stream being read from
   @param option[in] Command line option to search for

   @return True if found else false
*/


bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

int main(int argc, char *argv[])
{
  if (cmdOptionExists(argv, argv + argc, "-h")) {
    std::cout << argv[0]
              << " options, where options are:\n"
      "-h: help - Print this message.\n"
      "-f filename: Use filename as input.\n"
      "-s integer: Random seed (0 - use time).\n"
      "-t: run tests.\n";
    ParameterMap parameterMap;
    parameterMap.print(0, std::string(""));
    exit(1);
  }

  // In some implementations executing this may dramatically speed up
  // writing to standard output.
  // std::ios::sync_with_stdio(false);

  const char *seed_str = getCmdOption(argv, argv + argc, "-s");
  std::vector<ParameterMap> parameterMaps;
  const char *input_file_str = getCmdOption(argv, argv + argc, "-f");
  if (input_file_str) {
    std::ifstream infile;
    infile.open (input_file_str, std::ifstream::in);
    if (infile.fail()) {
      infile.close();
      std::cerr << "Error opening " << input_file_str << std::endl;
      exit(1);
    }
    readParameters(infile, parameterMaps);
    infile.close();
  }

  if (cmdOptionExists(argv, argv + argc, "-t")) {
    ParameterMap parameterMap;
    runTests(parameterMap);
  }

  execSimulationSet(parameterMaps, seed_str ? atoi(seed_str) : 1);
}
