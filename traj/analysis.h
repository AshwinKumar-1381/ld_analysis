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

	class Trajectory {

	public:

	int nAtoms, frame_nr, step;
	float time, timeStep, frameWidth;

	char *fpath, *pipeString, *pipeChar;
	FILE *fileI, *fileO;

	Trajectory(float timeStep = 1.0, float frameWidth = 1.0);
	~Trajectory();

	void openTrajectory();
	void closeTrajectory();
	void readThisFrame(atom_style *ATOMS);
	void readNextFrame(atom_style *ATOMS);

	};
}

#endif /*ANALYSIS_H*/