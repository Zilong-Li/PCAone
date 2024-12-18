/*******************************************************************************
 * @file        https://github.com/Zilong-Li/phaseless/src/log.hpp
 * @author      Zilong Li
 * Copyright (C) 2023. The use of this code is governed by the LICENSE file.
 ******************************************************************************/
#ifndef LOG_H_
#define LOG_H_

#include <cstring>
#include <fstream>
#include <iomanip> // setw
#include <iostream>
#include <sstream>

class Logger
{
  public:
    std::ofstream cao;
    bool is_screen;

    Logger() {}

    Logger(std::string filename, bool screen = true)
    {
        is_screen = screen;
        cao.open(filename.c_str());
        if(!cao) throw std::runtime_error(filename + " : " + std::strerror(errno));
        cao.precision(3);
        cao.flags(std::ios::fixed | std::ios::right);
    }

    ~Logger() {}

    template<class S>
    Logger & operator<<(const S & val)
    {
        cao << val;
        if(is_screen) std::cout << val;
        return *this;
    }

    Logger & operator<<(std::ostream & (*pfun)(std::ostream &))
    {
        pfun(cao);
        if(is_screen) pfun(std::cout);
        return *this;
    };

    template<class S>
    void printSpace(std::ostream & os, const S & val)
    {
        if(std::is_integral_v<std::decay_t<decltype(val)>>)
            os << std::setw(2) << val;
        else if(std::is_floating_point_v<std::decay_t<decltype(val)>>)
            os << val;
        else
            os << val << " ";
    }

    template<typename... Args>
    void print(const Args &... args)
    {
        (..., printSpace(cao, args));
        cao << std::endl;
        if(is_screen)
        {
            std::cout.precision(3);
            std::cout.flags(std::ios::fixed | std::ios::right);
            (..., printSpace(std::cout, args));
            std::cout << std::endl;
        }
    }

    // only print to stderr
    template<typename... Args>
    void cerr(const Args &... args)
    {
        (..., printSpace(std::cerr, args));
        std::cerr << std::endl;
    }

    template<typename... Args>
    void warn(const Args &... args)
    {
        (..., printSpace(cao, args));
        cao << std::endl;
        if(is_screen)
        {
            std::cout << "\x1B[33m";
            (..., printSpace(std::cout, args));
            std::cout << "\033[0m" << std::endl;
        }
    }

    template<typename... Args>
    void error(const Args &... args)
    {
        (..., printSpace(cao, args));
        cao << std::endl;
        std::ostringstream oss;
        oss << "\x1B[31m";
        (..., printSpace(oss, args));
        oss << "\033[0m" << std::endl;
        throw std::runtime_error(oss.str());
    }

    template<typename... Args>
    void done(const Args &... args)
    {
        (..., printSpace(cao, args));
        cao << std::endl;
        if(is_screen)
        {
            std::cout << "\x1B[32m";
            (..., printSpace(std::cout, args));
            std::cout << "\033[0m" << std::endl;
        }
    }
};

#endif // LOG_H_
