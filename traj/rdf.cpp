// computeRDF.cpp
// Calculates full and partial radial distribution functions - g(r)
// Bins particle pairs according to their radial separations

#define PI 3.14159265359 
#include "analysis.h"

using namespace analysis;

void computeRDF(float **RDF, atom_style *ATOMS, System *BOX, int nRDFtypes, float Rcut, float binW);

int main(int argc, char* argv[])
{
	// Trajectory params
	long eq_steps = long(1e7);
	long startStep = long(8e9), endStep = long(10e9);
	float dt = 5e-4;
	int frameW = int(1e5);

	// System params
	int nAtomTypes = 2;
	float xA = 0.5, xB = 0.5;
	float Lx = 150.0, Ly = 30.0, rho = 0.45;
	float Rcut = 10.0;
	int Nbins = 200;

	// Opening trajectory 
	Trajectory *TRAJ = new Trajectory(dt, frameW);
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data38/traj2.xyz");
	sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/lmp/Data38/laneRDF.dat");
	TRAJ -> openTrajectory();

	// RDF initialization
	float binW = float(Rcut/Nbins);
	int nRDFtypes = 1;
	if(nAtomTypes > 1) nRDFtypes = nAtomTypes*nAtomTypes + 1;

	float *rn = new float[Nbins];
	float **RDF = new float*[nRDFtypes];
	for(int i = 0; i < nRDFtypes; i++)
		RDF[i] = new float[Nbins];
	
	for(int i = 0; i < Nbins; i++)
	{
		rn[i] = binW*(i + 0.5);
		for(int j = 0; j < nRDFtypes; j++)
			RDF[j][i] = 0.0;
	}

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms);

	int ctr = 0;
	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->step >= (startStep + eq_steps) and TRAJ->step <= (endStep + eq_steps))
		{
			ctr++;
			computeRDF(RDF, ATOMS, BOX, nRDFtypes, Rcut, binW);
			printf("Processing step %ld, frame %d\n", TRAJ->step - eq_steps, ctr);
		}
	}

	// RDF normalization
	float FAC[nRDFtypes] = {1.0, xA*xA, xA*xB, xB*xA, xB*xB};
	float norm = 2 * PI * rho * TRAJ->nAtoms * binW * ctr;
	for(int i = 0; i < nRDFtypes; i++)
	{
		for(int j = 0; j < Nbins; j++)
			RDF[i][j] /= (FAC[i] * rn[j] * norm); 
	}

	printf("RDF computed and averaged over %d frames.\n", ctr);

	TRAJ -> write2file(rn, RDF, nRDFtypes, Nbins);
	TRAJ -> closeTrajectory();
}

void computeRDF(float **RDF, atom_style *ATOMS, System *BOX, int nRDFtypes, float Rcut, float binW)
{
	for(int i = 0; i < BOX->nAtoms; i++)
	{
		float rxi = ATOMS[i].rxt1;
		float ryi = ATOMS[i].ryt1;

		for(int j = i + 1; j < BOX->nAtoms; j++)
		{
			float rxj = ATOMS[j].rxt1;
			float ryj = ATOMS[j].ryt1;

			float dxij = rxj - rxi;
			float dyij = ryj - ryi;

			BOX -> checkMinImage(&dxij, &dyij);

			float drij2 = dxij*dxij + dyij*dyij; 
			if(drij2 < Rcut*Rcut) 					// Bin atom pairs
			{
				float drij = sqrt(drij2);
				int bin_id = int(drij/binW);

				RDF[0][bin_id] += 2.0;

				if(nRDFtypes > 1)
				{
					if(ATOMS[i].id == ATOMS[j].id)
					{
						if(ATOMS[i].id == 'N')
							RDF[1][bin_id] += 2.0;

						else if(ATOMS[i].id == 'O')
							RDF[4][bin_id] += 2.0;
					}

					else
					{
						RDF[2][bin_id] += 1.0;
						RDF[3][bin_id] += 1.0;	
					}	
				}
			}
		}
	}
}

void analysis::Trajectory::write2file(float *rn, float **RDF, int nRDFtypes, int Nbins)
{
	remove(fpathO);

	fileO = fopen(fpathO, "w");

	if(nRDFtypes == 1)
	{
		fprintf(fileO, "bin rn g_all\n");

		for(int i = 0; i < Nbins; i++)
			fprintf(fileO, "%d %f %f\n", i + 1, rn[i], RDF[0][i]);
	}

	else
	{
		fprintf(fileO, "bin rn g_all gAA gAB gBA gBB\n");

		for(int i = 0; i < Nbins; i++)
			fprintf(fileO, "%d %f %f %f %f %f %f\n", i + 1, rn[i], RDF[0][i], RDF[1][i], RDF[2][i], RDF[3][i], RDF[4][i]);
	}

	fclose(fileO);
}