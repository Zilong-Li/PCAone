#ifndef LOGGER_H_
#define LOGGER_H_

#include <fstream>
#include <iomanip>
#include <iostream>

class Logger
{
  public:
    std::ofstream cao;
    bool is_screen = true;

    Logger() {}

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
    }

    void print(std::string s)
    {
        if(is_screen) std::cout << std::setprecision(16) << s << std::endl;
        cao << s << std::endl;
    }

    void warning(std::string s)
    {
        if(is_screen)
            std::cout << std::endl
                      << "\x1B[33m"
                      << "WARNING: "
                      << "\033[0m" << s << std::endl;
        cao << std::endl << "WARNING: " << s << std::endl;
    }

    void error(std::string s)
    {
        if(is_screen)
            std::cout << std::endl
                      << "\x1B[31m"
                      << "ERROR: "
                      << "\033[0m" << s << std::endl;
        cao << std::endl << "ERROR: " << s << std::endl;
        exit(EXIT_FAILURE);
    }

    void done(std::string s)
    {
        if(is_screen)
            std::cout << std::endl
                      << "\x1B[32m"
                      << "DONE: "
                      << "\033[0m" << s << std::endl;
        cao << std::endl << "DONE: " << s << std::endl;
        exit(EXIT_SUCCESS);
    }
};

#endif // LOGGER_H_
