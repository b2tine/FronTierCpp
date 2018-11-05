#include "filecheck.h"
#include <iostream>

static bool getString(std::ifstream& fin, const char* mesg)
{
    std::streampos position;
    std::string line;

    position = fin.tellg();
    while (getline(fin, line))
    {
        if (line.find(mesg) != std::string::npos)
        {
            char ch;
            fin.seekg(-1, std::ios_base::cur);
            while (fin.get(ch) && ch != ':')
                fin.seekg(-2, std::ios_base::cur);
            fin.seekg(1, std::ios_base::cur);
            return true;
        }
    }
    fin.clear();
    fin.seekg(position);
    return false;
}

bool findAndLocate(std::ifstream& fin, const char* mesg)
{
    std::cout << mesg << " ";
    if (!getString(fin, mesg))
    {
        fin.seekg(0);
        if (!getString(fin, mesg))
        {
            std::cout << "Can't find string " << mesg << " in file!\n";
            return false;
        }
    }
    return true;
}

