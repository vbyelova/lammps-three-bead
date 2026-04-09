#include "addingData.h"
#include "percPerBC.h"

int main(int argc, char* argv[])
{

    percPerBC<3> test3D;

    std::string systemDataFile = argv[1];
    std::string percInfoFile = argv[2];
    std::string outputFile = argv[3];

    std::set<int> totalParticles;
    totalParticles = loadPar(systemDataFile);

    std::vector<Frame> percVals;
    percVals = loadPercVals(percInfoFile);

    std::ofstream file(outputFile);

    for (auto& frame : percVals)
    {
        std::map< std::array<int, 2>, std::array<int, 3>> frameData = frame.data;
        std::cout << "frame number " << frame.frameNumber << std::endl;

        int percDim = test3D.percolationDimension(totalParticles, frame.data);
        if (percDim != 3)
        {
            file << percDim << std::endl;
        }
        else
        {
            file << percDim << std::endl;
            break;
        }
    }


    file.close();
/*
    std::set<int> percDimSet;
    percDimSet = testAddingData3D();

    for (auto dim: percDimSet)
    {
        std::cout << dim << std::endl;
    }
*/
    return 0;


}
