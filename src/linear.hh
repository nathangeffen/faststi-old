#ifndef LINEAR_HH
#define LINEAR_HH

#include <vector>

typedef std::vector< std::vector<double > > DblMatrix;

std::vector<double> getCol(const DblMatrix&, unsigned);
std::vector<double> subVector(double, const std::vector<double>&);
std::vector<double> multVector(double, const std::vector<double>&);
std::vector<double> multVectors(const std::vector<double>&,
                                const std::vector<double>&);
std::vector<double> addVector(const std::vector<double>&,
                              const std::vector<double>&);
double sumVector(const std::vector<double>&);

#endif
