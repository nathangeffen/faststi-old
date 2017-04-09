#include <cstdio>
#include <cstdlib>

#include <iostream>

#include "simulate.hh"

/*
  Simple command line processing functions taken from:
  http://stackoverflow.com/questions/865668/parse-command-line-arguments
*/

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
    {
      return *itr;
    }
  return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

int main(int argc, char *argv[])
{
  if (cmdOptionExists(argv, argv + argc, "-h")) {
    printf("%s options, where options are:\n"
           "-h: help - Print this message.\n"
           "-f filename: Use filename as input.\n"
           "-s integer: Random seed (0 - use time).\n"
           "-t: run tests.\n",
           argv[0]);
    ParameterMap parameterMap;
    parameterMap.print(std::string(""));
    exit(1);
  }

  const char *seed_str = getCmdOption(argv, argv + argc, "-s");
  std::vector<ParameterMap> parameterMaps;
  const char *input_file_str = getCmdOption(argv, argv + argc, "-f");
  if (input_file_str) {
    std::ifstream infile;
    infile.open (input_file_str, std::ifstream::in);
    if (infile.fail()) {
      fprintf(stderr, "Error opening %s\n", input_file_str);
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
