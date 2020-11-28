#include "Utils.hpp"

size_t count_lines(string fpath)
{
    std::ifstream in(fpath);
    size_t count = 0;
    string line;
    while (getline(in, line)) {
        count++;
    }
    return count++;
}

string timestamp()
{
    time_t t = time(NULL);
    char *s = asctime(localtime(&t));
    s[strlen(s) - 1] = '\0';
    string str(s);
    str = string("[") + str + string("] ");
    return str;
}