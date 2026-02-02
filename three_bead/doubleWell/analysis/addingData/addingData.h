#include <iostream>
#include <list>
#include <fstream>
#include <string>
#include <map>
#include <vector>

std::list<int> loadPar(const std::string& fileName)
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
    for(int i = 1; i <= numPar; i++)
    {
        totalParticles.push_back(i);
    }

    std::cout << std::endl;
    return totalParticles;

}



struct Frame
{
    int frameNumber;
    std::map<std::array<int,2>, std::array<int,3>> data;
};
    std::vector<Frame> loadPercVals(const std::string& fileName)
    {
        std::vector<Frame> frames;
        std::ifstream fileInput(fileName);
        std::string line;
        Frame* currentFrame = nullptr;

        while getline((fileInput, line))
        {
            if (line.find("frame") != std::string::npos)
            {
                if (currentFrame != nullptr)
                {
                    frames.push_back(*currentFrame);
                    delete currentFrame;
                }

                currentFrame = new Frame();
                size_t pos = line.find("frame");

                frameNumber = std::stoi(line.substr(pos + 1));
                frames.pushback(frameNumber);
            }
        }

    }
