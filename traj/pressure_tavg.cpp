// pressure_tavg.cpp

#include "irvingKirkwood.h"

int main(int argc, char *argv[])
{
	float Lx = 200.0, Ly = 100.0, binWidth = 4.0;
	float timeStart = 10e3, timeEnd = 15e3;
	float PeA = 0.0, PeB = 50.0;

	Trajectory *TRAJ = new Trajectory(5e-4);
	sprintf(TRAJ -> fpathI, "../../LD/LD-cpp/Data55/traj3.xyz");
	sprintf(TRAJ -> fpathO, "../../LD/LD-cpp/Data55/pressure_tavg.dat");
	TRAJ -> openTrajectory();

	System *BOX = new System(Lx, Ly, TRAJ->nAtoms);
	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	pair_style *INTERACTION = new pair_style(1.0, 1.0, 1.122462048);

	int nProps = 4;
	Bin1D **Pin = new Bin1D*[nProps];
	Bin1D **Pkin = new Bin1D*[nProps];
	Bin1D *Pswim = new Bin1D(Lx, Ly, binWidth);

	for(int i = 0; i < nProps; i++)
	{
		Pin[i] = new Bin1D(Lx, Ly, binWidth);
		Pkin[i] = new Bin1D(Lx, Ly, binWidth);
	}

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
			analysis::computeKineticPressure(TRAJ, ATOMS, BOX, Pkin);
			analysis::computeSwimPressure(TRAJ, ATOMS, BOX, Pswim, PeA, PeB);

			for(int i = 0; i < nProps; i++)
			{
				Pin[i] -> normalize(Pin[i] -> bin, NULL, Ly);
				Pin[i] -> addBins(Pin[i] -> binavg, Pin[i] -> bin);
				Pin[i] -> zero(Pin[i] -> bin);

				Pkin[i] -> normalize(Pkin[i] -> bin, NULL);
				Pkin[i] -> addBins(Pkin[i] -> binavg, Pkin[i] -> bin);
				Pkin[i] -> zero(Pkin[i] -> bin);
			}

			Pswim -> normalize(Pswim -> bin, NULL);
			Pswim -> addBins(Pswim -> binavg, Pswim -> bin);
			Pswim -> zero(Pswim ->bin);
		}
	}

	for(int i = 0; i < nProps; i++)
	{
		Pin[i] -> normalize(Pin[i] -> binavg, NULL, float(nFrames));
		Pkin[i] -> normalize(Pkin[i] -> binavg, NULL, float(nFrames));
	}
	Pswim -> normalize(Pswim -> binavg, NULL, float(nFrames));

	printf("PDF averaged from %d to %d over %d frames.\n", int(timeStart), int(timeEnd), nFrames);

	TRAJ -> write2file(Pin, Pkin, Pswim);
	TRAJ -> closeTrajectory();
}

void analysis::Trajectory::write2file(Bin1D **Pin, Bin1D **Pkin, Bin1D *Pswim, int ctr)
{
	remove(fpathO);

	fileO = fopen(fpathO, "w");
	if( fileO == NULL)
	{
		printf("Cannot create file %s. Exiting...\n", fpathO);
		exit(-1);
	}

	fprintf(fileO, "x Pin_xx Pin_xy Pin_yx Pin_yy Pkin_xx Pkin_xy Pkin_yx Pkin_yy Pswim\n");

	for(int i = 0; i < Pin[0] -> nBins; i++)
		fprintf(fileO, "%f %f %f %f %f %f %f %f %f %f\n", (i+1)*Pin[0]->binWidth, Pin[0]->binavg[i], Pin[1]->binavg[i], Pin[2]->binavg[i], Pin[3]->binavg[i], Pkin[0]->binavg[i], Pkin[1]->binavg[i], Pkin[2]->binavg[i], Pkin[3]->binavg[i], Pswim->binavg[i]);

	fclose(fileO);
}