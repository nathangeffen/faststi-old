#include "parameters.hh"

/**
   Populates a vector of ParameterMaps by reading parameters from a stream.

   @param input is the stream read from
   @param parameterMaps is the vector to populate
*/

void readParameters(std::istream& input,
                    std::vector<ParameterMap>& parameterMaps)
{
  std::string lineString;
  unsigned lineNumber = 0;

  parameterMaps.push_back(ParameterMap());
  while (std::getline(input, lineString)) {
    RangeType rangeType = NONE;
    bool varyWithPrevious = false;
    ++lineNumber;
    boost::trim(lineString);
    if (lineString.size() == 0) continue;
    if (lineString[0] == '#') continue;
    if (lineString[0] == '-') {
      parameterMaps.push_back(ParameterMap());
      continue;
    }
    if (lineString.back() == '@') rangeType = THREE_PARM;
    if (lineString.back() == '!') rangeType = LIST;
    if (rangeType > NONE) {
      lineString.pop_back(); // Remove Range indicator
      boost::trim(lineString);
      if (lineString.back() == '+') varyWithPrevious = true;
    }
    std::istringstream line(lineString);
    std::string parameterName;
    line >> parameterName;
    try {
      parameterMaps.back().replaceParameter(parameterName.c_str(),
                                            line, rangeType, varyWithPrevious);

    } catch (const std::exception& e) {
      throw ConfigurationError(e.what(), lineNumber);
    }
  }
}

/**
   Creates a vector of doubles with values populated from *from* to *to*
   stepping by *step*.

   @param from first value in range
   @param to boundary of range, not included
   @param step added to each value in range to generate next value
*/

std::vector<double> setRange(const double from,
                             const double to,
                             const double step)
{
  assert(to >= from || step < 0.0);
  assert(step != 0 || from == to);
  std::vector<double> range;
  for (double d = from; d < to; d += step) range.push_back(d);
  return range;
}

/**
   Creates a vector of doubles with values populated using uniform random reals.

   @param elements number of elements in range
   @param from beginning of uniform random range
   @param to end of uniform random range
*/


std::vector<double> setRangeUniRand(const size_t elements,
                                    const double from,
                                    const double to)
{
  std::uniform_real_distribution<double> uni(from, to);
  std::vector<double> range;
  for (size_t i = 0; i < elements; ++i) range.push_back(uni(rng));
  return range;
}
