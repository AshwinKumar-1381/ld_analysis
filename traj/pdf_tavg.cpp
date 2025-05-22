// pdf_tavg.cpp

#include "analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	float Lx = 200.0, Ly = 50.0, binWidth = 4.0;
	long eq_steps = long(1e7);
	long startStep = long(7e8), endStep = long(10e8); 
	float dt = 5e-4;
	int frameW = 100000;

	Trajectory *TRAJ = new Trajectory(dt, frameW);
	sprintf(TRAJ -> fpathI, "/media/ashwin/One Touch/ashwin_md/biphasic/May2025/lmp/Data12/traj2.xyz");
	sprintf(TRAJ -> fpathO, "/media/ashwin/One Touch/ashwin_md/biphasic/May2025/lmp/Data12/pdf_tavg.dat");
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style [TRAJ->nAtoms];

	Bin1D *binA = new Bin1D(Lx, Ly, binWidth);
	Bin1D *binB = new Bin1D(Lx, Ly, binWidth);

	int nFrames = 0;
	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if((TRAJ->step >= eq_steps + startStep) and (TRAJ->step <= eq_steps + endStep))
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

	binA -> normalize(binA->binavg, NULL, nFrames);
	binB -> normalize(binB->binavg, NULL, nFrames);

	printf("PDF averaged from %d to %d over %d frames.\n", int(startStep*dt), int(endStep*dt), nFrames);

	TRAJ -> write2file(binA, binB);
	TRAJ -> closeTrajectory();

	return(0);
}

void analysis::Trajectory::write2file(Bin1D *binA, Bin1D *binB, int ctr, float timeToAvg)
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