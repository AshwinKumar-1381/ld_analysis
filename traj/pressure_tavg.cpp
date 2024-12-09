// pressure_tavg.cpp

#include "irvingKirkwood.h"

int main(int argc, char *argv[])
{
	float Lx = 200.0, Ly = 100.0, binWidth = 6.0;
	float timeStart = 10e3, timeEnd = 15e3;

	Trajectory *TRAJ = new Trajectory(5e-4);
	sprintf(TRAJ->fpathI, "../../LD/LD-cpp/Data32/traj3.xyz");
	sprintf(TRAJ->fpathO, "../../LD/LD-cpp/Data32/pressure_tavg_6.dat");
	TRAJ -> openTrajectory();

	System *BOX = new System(Lx, Ly, TRAJ->nAtoms);
	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	pair_style *INTERACTION = new pair_style(1.0, 1.0, 1.122462048);

	int nProps = 4;
	particleBin **Pin = new particleBin*[nProps];

	for(int i = 0; i < nProps; i++)
		Pin[i] = new particleBin(Lx, Ly, binWidth);

	int nFrames = 0;
	while (!feof(TRAJ->fileI))
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->step*TRAJ->timeStep >= timeStart and TRAJ->step*TRAJ->timeStep <= timeEnd)
		{
			nFrames++;
			printf("Reading frame %d\n", nFrames);

			TRAJ -> computeCom(ATOMS);
			analysis::computeInteractionPressure(TRAJ, ATOMS, BOX, INTERACTION, Pin);

			for(int i = 0; i < nProps; i++)
			{
				Pin[i] -> normalize(Pin[i] -> bin, Ly);
				Pin[i] -> addBins(Pin[i] -> binavg, Pin[i] -> bin);
				Pin[i] -> zero(Pin[i] -> bin);
			}
		}
	}

	for(int i = 0; i < nProps; i++)
		Pin[i] -> normalize(Pin[i] -> binavg, float(nFrames));

	TRAJ -> write2file(Pin);
	TRAJ -> closeTrajectory();
}

void analysis::Trajectory::write2file(particleBin **Pin, int ctr)
{
	remove(fpathO);

	fileO = fopen(fpathO, "w");
	if( fileO == NULL)
	{
		printf("Cannot create file %s. Exiting...\n", fpathO);
		exit(-1);
	}

	fprintf(fileO, "x Pin_xx Pin_xy Pin_yx Pin_yy\n");

	for(int i = 0; i < Pin[0] -> nBins; i++)
		fprintf(fileO, "%f %f %f %f %f\n", (i+1)*Pin[0]->binWidth, Pin[0]->binavg[i], Pin[1]->binavg[i], Pin[2]->binavg[i], Pin[3]->binavg[i]);

	fclose(fileO);
}