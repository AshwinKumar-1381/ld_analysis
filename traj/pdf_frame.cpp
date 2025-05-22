// pdf_frame.cpp

#include "analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	float Lx = 200.0, Ly = 50.0, binWidth = 4.0;
	float dt = 5e-4, timeToAvg = 500.0; 		// (in time units)
	int frameWidth = 10000;
	float sampleTimes[2] = {5e3, 30e3};
	
	// Create a trajectory object
	Trajectory *TRAJ = new Trajectory(dt, frameWidth);
	sprintf(TRAJ -> fpathI, "/media/ashwin/One Touch/ashwin_md/Feb2025/Pe_1e0/Data64/traj2.xyz");
	sprintf(TRAJ -> fpathO, "/media/ashwin/One Touch/ashwin_md/Feb2025/Pe_1e0/Data64/pdf_frame.dat");
	TRAJ -> openTrajectory();

	// Create an atoms object
	atom_style *ATOMS = new atom_style[TRAJ -> nAtoms];

	// Create bins
	Bin1D *binA = new Bin1D(Lx, Ly, binWidth);
	Bin1D *binB = new Bin1D(Lx, Ly, binWidth);

	int ctr = 0;
	int frame_ctr = 0;
	while( !feof(TRAJ -> fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if((int(TRAJ->step*TRAJ->timeStep) >= int(sampleTimes[ctr]) - timeToAvg/2) and (int(TRAJ->step*TRAJ->timeStep) < int(sampleTimes[ctr]) + timeToAvg/2))
		{
			frame_ctr++;
			TRAJ -> computeCom(ATOMS);
			float delxCom = TRAJ -> xCom - 0.5*Lx;

			for(int i = 0; i < TRAJ -> nAtoms; i++)
			{
				if(ATOMS[i].id == 'N')
					binA -> addToBin(ATOMS[i].rxt1 - delxCom);

				else if(ATOMS[i].id == 'O')
					binB -> addToBin(ATOMS[i].rxt1 - delxCom);
			}

			binA -> normalize(binA->bin);
			binB -> normalize(binB->bin);
			binA -> addBins(binA->binavg, binA->bin);
			binB -> addBins(binB->binavg, binB->bin);
			binA -> zero(binA->bin);
			binB -> zero(binB->bin);
		}

		if(((int(TRAJ->step*TRAJ->timeStep) == int(sampleTimes[ctr]) + timeToAvg/2)) or (feof(TRAJ -> fileI)))
		{
			binA -> normalize(binA->binavg, NULL, frame_ctr);
			binB -> normalize(binB->binavg, NULL, frame_ctr);

			TRAJ -> write2file(binA, binB, ctr, timeToAvg);

			binA -> zero(binA->binavg);
			binB -> zero(binB->binavg);

			printf("Computed PDF at time = %d by averaging over %d frames\n", int(sampleTimes[ctr]), frame_ctr);
			frame_ctr = 0;
			ctr++;
		}

		if(ctr == int(sizeof(sampleTimes)/sizeof(float)))
		{
			TRAJ -> closeTrajectory();
			return(0);
		}
	}
}

void analysis::Trajectory::write2file(Bin1D *binA, Bin1D *binB, int ctr, float timeToAvg)
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

	fprintf(fileO, "Time %d\n", int(step*timeStep - timeToAvg/2));

	for(int i = 0; i < binA->nBins; i++)
		fprintf(fileO, "%d %f %f\n", int((i+1)*binA->binWidth), binA->binavg[i], binB->binavg[i]);

	fclose(fileO);
}