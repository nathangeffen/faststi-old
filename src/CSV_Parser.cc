#include "CSV_Parser.hh"

void printDblMatrix(const DblMatrix& dblMatrix)
{
  for (auto& r: dblMatrix) {
    for (auto &d: r)
      std::cout << d << " ";
    std::cout << std::endl;
  }
}
