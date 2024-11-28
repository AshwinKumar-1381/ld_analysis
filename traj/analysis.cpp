// analysis.cpp

#include "analysis.h"

using namespace analysis;

analysis::atomsXYZ::atomsXYZ(){}
analysis::atomsXYZ::~atomsXYZ(){}

analysis::particleBin::particleBin(float L, float binWidth)
{
	this->L = L;
	this->binWidth = binWidth;
	nBins = int(L/binWidth);
	bin = new float [nBins];

	zero();
}

analysis::particleBin::~particleBin(){}

void analysis::particleBin::addParticle(float rx)
{
	bin[int(rx/binWidth)] += 1;	
}

void analysis::particleBin::normalize()
{
	for(int i = 0; i < nBins; i++)
		bin[i] /= (L*binWidth);
}

void analysis::particleBin::zero()
{
	for(int i = 0; i < nBins; i++)
		bin[i] = 0.0;
}

analysis::Trajectory::Trajectory(float timeStep, float frameWidth)
{
	fpathI = new char [500];
	fpathO = new char [500];
	pipeString = new char [500];
	pipeChar = new char [500];
	fileI = nullptr;
	fileO = nullptr;

	frame_nr = -1;
	this->timeStep = timeStep;
	this->frameWidth = frameWidth;
}

analysis::Trajectory::~Trajectory(){}

void analysis::Trajectory::openTrajectory()
{
	fileI = fopen(fpathI, "r");
	if(fileI == NULL)
	{
		printf("Could not open %s. Exiting ...\n", fpathI);
		exit(-1);
	}
	else
	{
		fgets(pipeString, 500, fileI);
		sscanf(pipeString, "%d", &nAtoms);
		rewind(fileI);
		printf("Trajectory file opened and ready to be read...\n");
	}
}

void analysis::Trajectory::closeTrajectory()
{
	printf("Closing trajectory file.\n");
	fclose(fileI);
}

void analysis::Trajectory::readThisFrame(atom_style *ATOMS)
{
	fgets(pipeString, 500, fileI);
	
	fgets(pipeString, 500, fileI);
	sscanf(pipeString, "%*s %d", &step);

	for(int i = 0; i < nAtoms; i++)
	{
		fgets(pipeString, 500, fileI);
		sscanf(pipeString, "%c %f %f %f", &ATOMS[i].id, &ATOMS[i].rxt1, &ATOMS[i].ryt1, &ATOMS[i].rzt1);
	}

	frame_nr++;
}