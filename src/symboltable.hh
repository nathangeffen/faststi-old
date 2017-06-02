#ifndef SYMBOLTABLE_HH
#define SYMBOLTABLE_HH

/**
   Manages a symbol table that holds the range values of varying parameter
   values. When multiple simulations are executed, on each new simulation, the
   varying parameters for the ParameterMap it uses are populated by the next
   entry in the symbol table.
 */

#include <vector>

class SymbolTable {
public:
  /** Specifes the number of entries int the symbol table which
      usually will correspond to the number of simulations.

      @param[in] numEntries Number of entries in the symbol table.
  */
  explicit SymbolTable(const unsigned numEntries) {
    table.resize(numEntries);
  }

  /**
     Populates all the entries in a symbol table with the entries in a range,
     cycling from the beginning of the range if the symbol table has more entries
     than the range.

     @param[in] name Name of the symbol table entry

     @param[in] parent The name of the parent in the symbol table such that this
     entry will only vary once the parent has cycled through its range. Empty
     string if no parent.

     @param[in] range List of numbers comprising the range
   */

  void addSymbol(const std::string& name, const std::string& parent,
                 const vector<double>& range)
  {
    assert( parent == "" || initialValues.find(parent) !=  initialValues.end() );
    assert(range.size() > 0);

    if (parent == "") {
      initialValues[name] = range.size();
    } else {
      initialValues[name] = range.size() * initialValues[parent];
    }
    size_t j = range.size();
    for (size_t i = 0; i < table.size(); ++i) {
      if (parent == std::string("") ||
          (i % initialValues[parent]) == 0) {
        ++j;
        if (j >= range.size()) j = 0;
      }
      table[i][name] = range[j];
    }
  }


  /**
     Gets reference to entry symbol table.

     @param[in] index The index of the entry to be retrieved from the symbol
     table.
   */

  std::unordered_map<std::string, double>& operator[](const size_t index)
  {
    assert(index < table.size());
    return table[index];
  }

  /** Gets the number of entries in the symbol table.
   */

  size_t size() const
  {
    return table.size();
  }


  /**
     Prints the symbol table to standard output. Mainly for debugging.
   */

  void output()
  {
    for (size_t i = 0; i < table.size(); ++i) {
      std::cout << "Iteration: " << i << std::endl;
      for (auto& s: table[i]) {
        std::cout << s.first << " " << s.second << std::endl;
      }
    }
  }

  std::vector< std::unordered_map<std::string, double> > table;
  std::unordered_map<std::string, size_t> initialValues;
};

#endif
