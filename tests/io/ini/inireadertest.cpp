// #include "freelb.h"
// #include "freelb.hh"

#include <iostream> 
#include <vector>
#include "io/ini_reader.h"

using T = float;

int main()
{
    iniReader reader("test.ini");

    std::cout << reader.getValue<int>("int", "intx") << std::endl;

    std::cout << reader.getValue<std::string>("string", "stringx") << std::endl;

    std::cout << reader.getValue<T>("real", "realx") << std::endl;

    std::cout << reader.getValue<T>("real", "realy") << std::endl;

    std::vector<int> vec;
    reader.getVector<int>("vector", vec);
    for (auto &i : vec)
    {
        std::cout << i << " ";
    }
    std::cout << std::endl;

    return 0;
}