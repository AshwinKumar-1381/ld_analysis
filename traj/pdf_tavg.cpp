// pdf_tavg.cpp

#include "analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	float Lx = 200.0, Ly = 100.0, binWidth = 4.0;
	float timeStart = 10e3, timeEnd = 15e3;

	Trajectory *TRAJ = new Trajectory(5e-4);
	sprintf(TRAJ -> fpathI, "../../LD/LD-cpp/Data55/traj3.xyz");
	sprintf(TRAJ -> fpathO, "../../LD/LD-cpp/Data55/pdf_tavg.dat");
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style [TRAJ->nAtoms];

	particleBin *binA = new particleBin(Lx, Ly, binWidth);
	particleBin *binB = new particleBin(Lx, Ly, binWidth);

	int nFrames = 0;
	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->step*TRAJ->timeStep >= timeStart and TRAJ->step*TRAJ->timeStep <= timeEnd)
		{
			TRAJ -> computeCom(ATOMS);
			float delxCom = TRAJ -> xCom - 0.5*Lx;

			for(int i = 0; i < TRAJ -> nAtoms; i++)
			{
				if(ATOMS[i].id == 'N')
					binA -> addToBin(ATOMS[i].rxt1 - delxCom);
				else
					binB -> addToBin(ATOMS[i].rxt1 - delxCom);
			}

			binA -> normalize(binA->bin);
			binB -> normalize(binB->bin);
			binA -> addBins(binA->binavg, binA->bin);
			binB -> addBins(binB->binavg, binB->bin);
			binA -> zero(binA->bin);
			binB -> zero(binB->bin);
			nFrames++;
		}
	}

	binA -> normalize(binA->binavg, nFrames);
	binB -> normalize(binB->binavg, nFrames);

	printf("PDF averaged from %d to %d over %d frames.\n", int(timeStart), int(timeEnd), nFrames);

	TRAJ -> write2file(binA, binB);
	TRAJ -> closeTrajectory();

	return(0);
}

void analysis::Trajectory::write2file(particleBin *binA, particleBin *binB, int ctr)
{
	remove(fpathO);

	fileO = fopen(fpathO, "w");
	if(fileO == NULL)
	{
		printf("Cannot create %s file. Exiting...\n", fpathO);
		exit(-1);
	}

	fprintf(fileO, "x pdfA pdfB\n");

	for(int i = 0; i < binA->nBins; i++)
		fprintf(fileO, "%d %f %f\n", int((i+1)*binA->binWidth), binA->binavg[i], binB->binavg[i]);

	fclose(fileO);
}