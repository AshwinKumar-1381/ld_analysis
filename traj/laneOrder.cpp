// laneOrder.cpp

#include "analysis.h"

using namespace analysis;

float computeLaneOrder(atom_style *ATOMS, System *BOX, float slabW_y);

int main(int argc, char*argv[])
{
	char *option = new char[20];
	sprintf(option, "time_avg");

	// Trajectory params
	float dt = 5e-4;
	long eq_steps = long(1e7); 
	long startStep = long(8e7), endStep = long(1e8);
	int frameW = int(1e3);

	// System params
	float Lx = 150.0, Ly = 30.0; 
	float slabW = 1.0;

	// Keywords and values for option = time_evolve 
	int Nfreq = 500; 	
	int Nsamples = 100; 
	int Nevery = 1;		

	Trajectory *TRAJ = new Trajectory(dt, frameW);
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data51/traj2.xyz");
	sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data51/laneOrder.dat");

	TRAJ -> openTrajectory();
	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms);

	int frameStart = int(startStep/frameW);
	int frameEnd = int(endStep/frameW); 

	if(strcmp(option, "time_evolve") == 0)
	{
		TRAJ -> write2file(0.0, 0.0, 1);

		int ctr = 0;
		float avg = 0.0;

		int nextFrame = frameStart + ctr*Nfreq;
		int frameStartAvg = nextFrame;

		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);
			TRAJ->time = (TRAJ->step - eq_steps) * TRAJ->timeStep; 

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				if(TRAJ->frame_nr == frameStart)
				{
					ctr++;
					nextFrame = frameStart + ctr*Nfreq;
					frameStartAvg = nextFrame - Nsamples*Nevery;

					float order = computeLaneOrder(ATOMS, BOX, slabW);
					TRAJ -> write2file(TRAJ -> time, order);
				}

				if(TRAJ->frame_nr == frameStartAvg) avg = 0.0;

				if((TRAJ->frame_nr > frameStartAvg) and (TRAJ->frame_nr <= nextFrame))
				{
					if((TRAJ->frame_nr - frameStartAvg)%Nevery == 0)
						avg += computeLaneOrder(ATOMS, BOX, slabW);
				}

				if(TRAJ->frame_nr == nextFrame)
				{
					ctr++;
					nextFrame = frameStart + ctr*Nfreq;
					frameStartAvg = nextFrame - Nsamples*Nevery;

					float order = float(avg/Nsamples);
					TRAJ -> write2file(TRAJ -> time, order);
				}
			}
		}

		printf("Computed order parameter for %d frames.\n", ctr);

		TRAJ -> closeTrajectory();
		TRAJ -> write2file(0.0, 0.0, 0);
	}

	else if(strcmp(option, "time_avg") == 0)
	{
		int ctr = 0;
		float avg = 0.0;

		while( !feof(TRAJ -> fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				avg += computeLaneOrder(ATOMS, BOX, slabW);
				ctr++;
			}

			if(TRAJ->frame_nr == frameEnd)
			{
				avg /= ctr;
				break;
			}
		}

		printf("Lane order parameter averaged over %d frames = %f\n", ctr, avg);
	}
}

float computeLaneOrder(atom_style *ATOMS, System *BOX, float slabW_y)
{
	int nBins = int(BOX->Ly/slabW_y);

	float phi[nBins];
	int bin_ctr[nBins];
	for(int i = 0; i < nBins; i++) 
	{
		phi[i] = 0.0;
		bin_ctr[i] = 0;
	}

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		int bin_id = int(ATOMS[i].ryt1 / slabW_y);

		if(ATOMS[i].id == 'N')
			phi[bin_id] += -1.0;

		else if(ATOMS[i].id == 'O')
			phi[bin_id] += 1.0;

		bin_ctr[bin_id] += 1;
	}

	float order = 0.0;
	for(int i = 0; i < nBins; i++)
		order += abs(phi[i] / bin_ctr[i]);

	return(order/nBins);
}

void analysis::Trajectory::write2file(float time, float order, int option)
{
	if(option == 1)
	{
		remove(fpathO);

		fileO = fopen(fpathO, "a+");
		if(fileO == NULL)
		{
			printf("Data file cannot be created. Exiting...\n");
			exit(-1);
		}

		fprintf(fileO, "time order\n");
	}

	else if(option == 0)
		fclose(fileO);

	else
		fprintf(fileO, "%f %f\n", time, order);
}