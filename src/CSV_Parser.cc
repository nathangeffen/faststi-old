#include "CSV_Parser.hh"

/**
   Prints the values of a matrix.

   @param dblMatrix[in] Matrix of doubles to print
 */

void printDblMatrix(const DblMatrix& dblMatrix)
{
  for (auto& r: dblMatrix) {
    for (auto &d: r)
      std::cout << d << " ";
    std::cout << std::endl;
  }
}
