#ifndef SYMBOLTABLE_HH
#define SYMBOLTABLE_HH

#include <vector>

class SymbolTable {
public:
  SymbolTable(unsigned numEntries) {
    table.resize(numEntries);
  }

  void addSymbol(const std::string& name, const std::string& parent,
                 const vector<double>& range)
  {
    assert( parent == "" || initialValues.find(parent) !=  initialValues.end()  );
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


  std::unordered_map<std::string, double>& operator[](const size_t index)
  {
    assert(index < table.size());
    return table[index];
  }

  const std::unordered_map<std::string, double>&
  operator[](const size_t index) const
  {
    assert(index < table.size());
    return table[index];
  }

  size_t size() const
  {
    return table.size();
  }

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
