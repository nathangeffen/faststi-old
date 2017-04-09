#include <cstdlib>

#include "linear.hh"

/**
   Gets a column in a matrix.

   @param matrix the matrix to get the column from
   @param col the zero-indexed column in the matrix to return
   @return matrix column as a vector of doubles
*/

std::vector<double> getCol(const DblMatrix& matrix, unsigned col)
{
  std::vector<double> output;
  for (unsigned i = 0; i < matrix.size(); ++i)
    output.push_back(matrix[i][col]);
  return output;
}

/**
   Calculates a vector subtracted from a scalar.

   @param d scalar from which to subtract
   @param v vector to subtract
   @return vector = d - v
*/

std::vector<double> subVector(double d,
                              const std::vector<double>& v)
{
  std::vector<double> output;
  for (size_t i = 0; i < v.size(); ++i)
    output.push_back(d - v[i]);
  return output;
}

/**
   Calculates a vector multiplied by a scalar.

   @param d scalar to multiply
   @param v vector to multiply
   @return vector = d * v
*/

std::vector<double> multVector(double d,
                               const std::vector<double>& v)
{
  std::vector<double> output;
  for (size_t i = 0; i < v.size(); ++i)
    output.push_back(d * v[i]);
  return output;
}

/**
   Calculates product of two vectors.

   @param v1 vector to multiply
   @param v2 vector to multiply
   @return vector = v1 X v2
*/

std::vector<double> multVectors(const std::vector<double>& v1,
                                const std::vector<double>& v2)
{
  std::vector<double> output;
  for (size_t i = 0; i < v1.size(); ++i)
    output.push_back(v1[i] * v2[i]);
  return output;
}

/**
   Calculates sum of two vectors.

   @param v1 vector to add
   @param v2 vector to add
   @return vector = v1 + v2
*/

std::vector<double> addVector(const std::vector<double>& v1,
                              const std::vector<double>& v2)
{
  std::vector<double> output;
  for (size_t i = 0; i < v1.size(); ++i)
    output.push_back(v1[i] + v2[i]);
  return output;
}

/**
   Calculates sum of elements of vector.

   @param v vector to sum
   @return sum of elements of
*/

double sumVector(const std::vector<double>& v)
{
  double total = 0.0;
  for (auto d: v)
    total += d;
  return total;
}
