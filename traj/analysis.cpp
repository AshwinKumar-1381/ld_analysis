// analysis.cpp

#include "analysis.h"

using namespace analysis;

analysis::atomsXYZ::atomsXYZ(){}
analysis::atomsXYZ::~atomsXYZ(){}

analysis::Trajectory::Trajectory(float timeStep, float frameWidth)
{
	fpath = new char [500];
	pipeString = new char [500];
	pipeChar = new char [500];
	fileI = nullptr;
	fileO = nullptr;

	frame_nr = 0;
	this->timeStep = timeStep;
	this->frameWidth = frameWidth;
}

analysis::Trajectory::~Trajectory(){}

void analysis::Trajectory::openTrajectory()
{
	fileI = fopen(fpath, "r");
	if(fileI == NULL)
	{
		printf("Could not open %s. Exiting ...\n", fpath);
		exit(-1);
	}
	else
	{
		fgets(pipeString, 500, fileI);
		sscanf(pipeString, "%d", &nAtoms);
		rewind(fileI);
	}
}

void analysis::Trajectory::closeTrajectory()
{
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