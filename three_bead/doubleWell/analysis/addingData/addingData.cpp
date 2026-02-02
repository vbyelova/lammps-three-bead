#include "addingData.h"

int main()
{
    std::set<int> totalParticles;
    totalParticles = loadPar("something.txt");

    for (int p : totalParticles)
    {
        std::cout << p << " ";
    }
    std::cout << std::endl;
    return 0;
}
