// analysis.h

#ifndef ANALYSIS_H
#define ANALYSIS_H

#define MAXCELL 9000000

#include "library.h"

namespace analysis {

	typedef class atomsXYZ {

	public:

	char id;
	float rxt1, ryt1, rzt1;
	float rxt2, ryt2, rzt2; 
	float px, py;
	float dx, dy, dz;

	atomsXYZ();
	~atomsXYZ();
	
	} atom_style;

	class System {

	public:

	float Lx, Ly, rcellx, rcelly;
	int Ncellx, Ncelly, ncells;
	int nAtoms;

	int MAPS[MAXCELL], HEAD[MAXCELL], LIST[MAXCELL];

	System(float Lx, float Ly, float nAtoms);
	~System();

	int cellindex(int ix, int iy);
	void buildCellMaps();
	void buildCellList(atom_style *ATOMS);
	void checkMinImage(float *dx = NULL, float *dy = NULL);

	};

	class particleBin {

	public:

	float Lx, Ly, binWidth;
	int nBins;
	float *bin;
	float *binavg;

	particleBin(float Lx, float Ly, float binWidth);
	particleBin();
	~particleBin();

	void addToBin(float rx = 0.0, float val = 1.0);
	void normalize(float *binToavg, float norm = 0);
	void zero(float *binTozero);
	void addBins(float *bin1, float *bin2);

	};

	class Trajectory {

	public:

	int nAtoms, frame_nr, step;
	float timeStep, time;
	float xCom, yCom, zCom; 

	char *fpathI, *fpathO, *pipeString, *pipeChar;
	FILE *fileI, *fileO;

	Trajectory(float timeStep = 1.0);
	~Trajectory();

	void openTrajectory();
	void closeTrajectory();
	void readThisFrame(atom_style *ATOMS);
	void readNextFrame(atom_style *ATOMS);
	void computeCom(atom_style *ATOMS);

	void write2file();
	void write2file(particleBin *binA, particleBin *binB, int ctr = 0);
	void write2file(particleBin **Pin, int ctr = 0);

	};
}

#endif /*ANALYSIS_H*/