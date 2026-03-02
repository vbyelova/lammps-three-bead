#include <iostream>
#include <list>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <array>
#include "percPerBC.h"

std::set<int> loadPar(const std::string& fileName)
{
    // load in number of particles
    int numPar {};

    std::string line;
    std::ifstream fileInput(fileName);

    // read input file
    // find number of particles
    while (getline (fileInput, line))
    {
        if (line.find("Npar") != std::string::npos)
        {
            std::size_t pos = line.find("=");
            if (pos != std::string::npos)
            {
                numPar = std::stoi(line.substr(pos + 1));
                break;
            }
        }
    }

    // make a set of particles
    std::set<int> totalParticles {};
    for(int i = 0; i <= numPar; i++)
    {
        totalParticles.insert(i);
    }
    std::cout << std::endl;
    fileInput.close();
    std::cout << "Loaded " << totalParticles.size() << " particles, max=" << *totalParticles.rbegin() << std::endl;
    return totalParticles;

}


struct Frame
{
    int frameNumber;
    std::map<std::array<int, 2>, std::array<int, 3>> data;
};
    std::vector<Frame> loadPercVals(const std::string& fileName)
    {
        int numBonds;
        std::vector<Frame> frames;
        std::ifstream fileInput(fileName);
        std::string line;


        while (getline(fileInput, line))
        {
            Frame currentFrame;
            currentFrame.frameNumber = std::stoi(line);
            if (!std::getline(fileInput, line)) break;
            numBonds = std::stoi(line);

            std::cout << "frame number " << currentFrame.frameNumber << " number of bonds " << numBonds << std::endl;


            for (int i = 0; i < numBonds; i++ )
            {
                std::array<int, 2> bondedAtoms;
                std::array<int, 3> percAcross;

                if (!std::getline(fileInput, line)) break;
                bondedAtoms[0] = std::stoi(line);

                if (!std::getline(fileInput, line)) break;
                bondedAtoms[1] = std::stoi(line);

                if (!std::getline(fileInput, line)) break;
                percAcross[0] = std::stoi(line);

                if (!std::getline(fileInput, line)) break;
                percAcross[1] = std::stoi(line);

                if (!std::getline(fileInput, line)) break;
                percAcross[2] = std::stoi(line);

                currentFrame.data[bondedAtoms] = percAcross;
            }
            frames.push_back(currentFrame);
        }
        fileInput.close();
        return frames;
    }

/*
std::set<int> testAddingData3D()
{
    percPerBC<3> test3D;
    std::set<int> vectors {};
    for (int i = 0; i <= 13; i++ )
    {
        vectors.insert(i);
    }

    std::set<int> totalParticles;

    std::vector<Frame> percVals;
    percVals = loadPercVals("testfile.txt");

    std::ofstream file;
    std::set<int> percDimSet;
    file.open("percVals.txt");

    for (auto& frame : percVals)
    {
        std::map< std::array<int, 2>, std::array<int, 3>> frameData = frame.data;
        int percDim = test3D.percolationDimension(vectors, frame.data);
        std::cout << percDim << std::endl;
        file << percDim << std::endl;
        percDimSet.insert(percDim);
    }

    file.close();
    return percDimSet;
}
*/
