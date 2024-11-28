// pdf_tavg.cpp

#include "analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	// Trajectory params - dt frame_w

	Trajectory *TRAJ = new Trajectory(5e-4, 0.25);
	sprintf(TRAJ->fpath, "../../LD/LD-cpp/Data15/traj3.xyz");

	TRAJ -> openTrajectory();
	atom_style *ATOMS = new atom_style [TRAJ->nAtoms]; 

	while(!feof(TRAJ->fileI))
	{
		TRAJ -> readThisFrame(ATOMS);
		printf("Reading frame %d at step %d\n", TRAJ->frame_nr, TRAJ->step);
		
		
	}

	printf("%d, %s\n", TRAJ->nAtoms, TRAJ->fpath);

	return(0);
}