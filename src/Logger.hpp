/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/log.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef LOG_H_
#define LOG_H_

#include <cstring>
#include <fstream>
#include <iomanip>  // setw
#include <iostream>
#include <sstream>
#include <stdexcept>

class Logger {
 public:
  std::ofstream cao;
  bool is_screen = false;

  Logger() {}

  ~Logger() {}

  template <class S>
  void printSpace(std::ostream& os, const S& val) {
    os << val << " ";
  }

  template <typename... Args>
  void print(const Args&... args) {
    (..., printSpace(cao, args));
    cao << std::endl;
    if (is_screen) {
      std::cout.precision(6);
      std::cout.flags(std::ios::fixed | std::ios::right);
      (..., printSpace(std::cout, args));
      std::cout << std::endl;
    }
  }

  template <typename... Args>
  void warn(const Args&... args) {
    (..., printSpace(cao, args));
    cao << std::endl;
    if (is_screen) {
      std::cout << "\x1B[33m";
      (..., printSpace(std::cout, args));
      std::cout << "\033[0m" << std::endl;
    }
  }

  template <typename... Args>
  void error(const Args&... args) {
    (..., printSpace(cao, args));
    cao << std::endl;
    std::ostringstream oss;
    oss << "\x1B[31m";
    (..., printSpace(oss, args));
    oss << "\033[0m" << std::endl;
    throw std::runtime_error(oss.str());
  }
};

#endif  // LOG_H_
