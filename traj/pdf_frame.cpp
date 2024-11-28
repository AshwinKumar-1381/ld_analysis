// pdf_frame.cpp

#include "analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	Trajectory *TRAJ = new Trajectory(5e-4, 0.25);
	sprintf(TRAJ -> fpath, "../../LD/LD-cpp/Data15/traj3.xyz");

	float sampleTimes[5] = {1e3, 2e3, 3e3, 4e3, 5e3};

	TRAJ -> openTrajectory();
	atom_style *ATOMS = new atom_style[TRAJ -> nAtoms];

	while( !feof(TRAJ -> fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);
	}

	TRAJ -> closeTrajectory();
}