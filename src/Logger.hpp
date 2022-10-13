#ifndef LOGGER_H_
#define LOGGER_H_

#include <fstream>
#include <iostream>


class Logger
{
public:
    std::ofstream clog;
    bool isprint = true;

    template <class S>
    Logger& operator<<(const S& val)
    {
        clog << val;
        if (isprint)
            std::cout << val;
        return *this;
    }

    Logger& operator<<(std::ostream& (*pfun)(std::ostream&))
    {
        pfun(clog);
        if (isprint)
            pfun(std::cout);
        return *this;
    };

    Logger()
    {
    }

    ~Logger()
    {
    }
};


#endif // LOGGER_H_
