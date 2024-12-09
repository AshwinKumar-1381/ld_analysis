// pdf_frame.cpp

#include "analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	float sampleTimes[5] = {11e3, 12e3, 13e3, 14e3, 15e3};
	float Lx = 200.0, Ly = 100.0;
	float binWidth = 4.0;

	// Create a trajectory object
	Trajectory *TRAJ = new Trajectory(5e-4);
	sprintf(TRAJ -> fpathI, "../../LD/LD-cpp/Data32/traj3.xyz");
	sprintf(TRAJ -> fpathO, "../../LD/LD-cpp/Data32/pdf_frames.dat");
	TRAJ -> openTrajectory();

	// Create an atoms object
	atom_style *ATOMS = new atom_style[TRAJ -> nAtoms];

	// Create bins
	particleBin *binA = new particleBin(Lx, Ly, binWidth);
	particleBin *binB = new particleBin(Lx, Ly, binWidth);

	int ctr = 0;
	while( !feof(TRAJ -> fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(int(TRAJ->step*TRAJ->timeStep) == int(sampleTimes[ctr]))
		{
			printf("Computing PDF at time = %d\n", int(TRAJ->step*TRAJ->timeStep));

			TRAJ -> computeCom(ATOMS);
			float delxCom = TRAJ -> xCom - 0.5*Lx;

			for(int i = 0; i < TRAJ -> nAtoms; i++)
			{
				if(ATOMS[i].id == 'N')
					binA -> addToBin(ATOMS[i].rxt1 - delxCom);

				else if(ATOMS[i].id == 'O')
					binB -> addToBin(ATOMS[i].rxt1 - delxCom);
			}

			binA -> normalize();
			binB -> normalize();
			TRAJ -> write2file(binA, binB, ctr);
			binA -> zero();
			binB -> zero();

			ctr++;
		}

		if(ctr == int(sizeof(sampleTimes)/sizeof(float)))
		{
			TRAJ -> closeTrajectory();
			return(0);
		}
	}
}

void analysis::Trajectory::write2file(particleBin *binA, particleBin *binB, int ctr)
{
	if(ctr == 0)
	{
		remove(fpathO);

		fileO = fopen(fpathO, "w");
		if(fileO == NULL)
		{
			printf("Cannot create %s. Exiting...\n", fpathO);
			exit(-1);
		}

		fprintf(fileO, "x pdfA pdfB\n");
	}
	else
	{
		fileO = fopen(fpathO, "a+");
		if(fileO == NULL)
		{
			printf("Cannot open %s. Exiting...\n", fpathO);
			exit(-1);
		}
	}

	fprintf(fileO, "Time %d\n", int(step*timeStep));

	for(int i = 0; i < binA->nBins; i++)
		fprintf(fileO, "%d %f %f\n", int((i+1)*binA->binWidth), binA->bin[i], binB->bin[i]);

	fclose(fileO);
}