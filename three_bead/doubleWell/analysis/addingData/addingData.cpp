#include "addingData.h"
#include "percPerBC.h"

int main()
{
    /*
    percPerBC<3> test3D;
    std::set<int> totalParticles;
    totalParticles = loadPar("systemData.txt");

    std::vector<Frame> percVals;
    percVals = loadPercVals("percinfo.txt");

    std::ofstream file;
    file.open("percVals.txt");

    for (auto& frame : percVals)
    {
        std::map< std::array<int, 2>, std::array<int, 3>> frameData = frame.data;
        int percDim = test3D.percolationDimension(totalParticles, frame.data);
        //std::cout << percDim << std::endl;
        file << percDim << std::endl;
    }

    file.close();
    */
    std::set<int> percDimSet;
    percDimSet = testAddingData3D();
/*
    for (auto dim: percDimSet)
    {
        std::cout << dim << std::endl;
    }
*/
    return 0;


}
