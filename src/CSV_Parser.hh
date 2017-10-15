#ifndef CSVPARSER_HH
#define CSVPARSER_HH

/**
   C++ interface to C code of csvparser.
*/

#include <iostream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

#include "linear.hh"

extern "C" {
#include "csvparser.h"
}

class CSVParser {
public:
  /**
     Reads a CSV file into memory.

     @param filename[in] Name of CSV file
     @param delim[in] CSV delimiter, defaults to comma
     @param hasHeader[in] Whether the file has a header, defaults to true
   */
  explicit CSVParser(const char* filename,
                     const char* delim = ",",
                     const bool hasHeader = true)
  {
    csvparser = CsvParser_new(filename, delim, hasHeader);
    assert(csvparser);
    if (hasHeader) {
      header_ = CsvParser_getHeader(csvparser);
      if (!header_) {
	std::cerr << "Error getting csv header in " << filename << std::endl;
	exit(1);
      }
    }
    CsvRow *row;
    while ((row = CsvParser_getRow(csvparser)) ) {
      rows_.push_back(row);
      std::vector< std::string > str_row;
      const char **rowFields = CsvParser_getFields(row);
      for (int i = 0 ; i < CsvParser_getNumFields(row) ; i++)
        str_row.push_back(rowFields[i]);
      string_rows.push_back(str_row);
    }
  };


  /**
     Converts the CSV cells into doubles.
   */
  DblMatrix toDoubles() {
    DblMatrix double_rows;
    for (auto& r: string_rows) {
      std::vector<double> double_row;
      for (auto& s: r) {
	boost::replace_all(s, ",", ".");
	boost::replace_all(s, "NA", "0");
	double_row.push_back(stod(s));
      }
      double_rows.push_back(double_row);
    }
    return double_rows;
  };

  ~CSVParser() {
    for (auto& row: rows_) CsvParser_destroy_row(row);
    CsvParser_destroy(csvparser);
  };
  std::vector< std::vector<std::string> > string_rows;
private:
  CsvParser *csvparser;
  CsvRow *header_ = NULL;
  std::vector<CsvRow *> rows_;
};

void printDblMatrix(const DblMatrix&);

#endif
