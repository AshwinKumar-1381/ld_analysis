// analysis.h

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "library.h"

namespace analysis {

	typedef class atomsXYZ {

	public:

	char id;
	float rxt1, ryt1, rzt1;
	float rxt2, ryt2, rzt2; 
	float dx, dy, dz;

	atomsXYZ();
	~atomsXYZ();
	
	} atom_style;

	class particleBin {

	public:

	float L, binWidth;
	int nBins;
	float *bin;

	particleBin(float L, float binWidth);
	~particleBin();

	void addParticle(float rx = 0.0);
	void normalize();
	void zero();

	};

	class Trajectory {

	public:

	int nAtoms, frame_nr, step;
	float timeStep, time, frameWidth;

	char *fpathI, *fpathO, *pipeString, *pipeChar;
	FILE *fileI, *fileO;

	Trajectory(float timeStep = 1.0, float frameWidth = 1.0);
	~Trajectory();

	void openTrajectory();
	void closeTrajectory();
	void readThisFrame(atom_style *ATOMS);
	void readNextFrame(atom_style *ATOMS);

	void write2file();
	void write2file(particleBin *binA, particleBin *binB, int ctr);

	};
}

#endif /*ANALYSIS_H*/