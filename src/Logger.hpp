#ifndef LOGGER_H_
#define LOGGER_H_

#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>

class Logger {
 public:
  std::ofstream cao;
  bool is_screen = true;

  Logger() {}

  ~Logger() {}

  template <class S>
  Logger &operator<<(const S &val) {
    cao << val;
    if (is_screen) std::cout << val;
    return *this;
  }

  Logger &operator<<(std::ostream &(*pfun)(std::ostream &)) {
    pfun(cao);
    if (is_screen) pfun(std::cout);
    return *this;
  }

  template <class S>
  void printSpace(std::ostream &os, const S &val) {
    os << val << " ";
  }

  template <typename... Args>
  void print(const Args &...args) {
    std::initializer_list<int>{(printSpace(cao, args), 0)...};
    cao << std::endl;
    if (is_screen) {
      std::cout.precision(6);
      std::cout.flags(std::ios::fixed | std::ios::right);
      std::initializer_list<int>{(printSpace(std::cout, args), 0)...};
      std::cout << std::endl;
    }
  }

  template <class T>
  void warn(const T &s) {
    if (is_screen)
      std::cout << std::endl
                << "\x1B[33m"
                << "WARNING: "
                << "\033[0m" << s << std::endl;
    cao << std::endl << "WARNING: " << s << std::endl;
  }

  template <class T>
  void error(const T &s) {
    if (is_screen)
      std::cout << std::endl
                << "\x1B[31m"
                << "ERROR: "
                << "\033[0m" << s << std::endl;
    cao << std::endl << "ERROR: " << s << std::endl;
    exit(EXIT_FAILURE);
  }

  template <class T>
  void done(const T &s) {
    if (is_screen)
      std::cout << std::endl
                << "\x1B[32m"
                << "DONE: "
                << "\033[0m" << s << std::endl;
    cao << std::endl << "DONE: " << s << std::endl;
    exit(EXIT_SUCCESS);
  }
};

#endif  // LOGGER_H_
