// computeRDF_1D.cpp
// Calculates 1D full and partial pair distribution functions - g(x) and g(y)
// Bins particles according to their longitudinal (x-) and lateral (y-) distances
//
// Date created  : 25.04.25
// Last modified : 27.04.25

#define PI 3.14159265359 
#include "analysis.h"

using namespace analysis;

float*** buildRDF_1D(int nRDFtypes, int Nbins[]);
void computeRDF_1D(float ***RDF_x_y, atom_style *ATOMS, System *BOX, int nRDFtypes, float Rcut[], float binW[], char option[]);

int main(int argc, char* argv[])
{
	// ---------- Trajectory params ----------
	long eq_steps = long(0e7);
	long startStep = long(8e7), endStep = long(10e7);
	float dt = 5e-4;
	int frameW = int(1e4);

	// ---------- System params ----------
	int nAtomTypes = 2;
	float xA = 0.5, xB = 0.5;
	float Lx = 150.0, Ly = 30.0, rho = 0.45;

	// ---------- RDF params ----------
	float Rcut[2] = {10.0, 10.0};
	char option[10] = "asymm";
	int Nbins[2] = {200, 200};

	// ---------- Opening trajectory ----------
	Trajectory *TRAJ = new Trajectory(dt, frameW);
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/cpp/Data6/traj2.xyz");
	sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/lane/Apr2025/cpp/Data6/laneRDF_1D_%s.dat", option);
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(Lx, Ly, TRAJ->nAtoms);

	// ---------- RDF initialization ----------
	float binW[2] = {float(Rcut[0]/Nbins[0]), float(Rcut[1]/Nbins[1])};

	int nRDFtypes = 1;
	if(nAtomTypes > 1) nRDFtypes = nAtomTypes*nAtomTypes + 1;

	int fac = 0;
	if(strcmp(option, "asymm") == 0) fac = 1;
	for(int i = 0; i < 2; i++) Nbins[i] *= (1 + fac);

	auto RDF_x_y = buildRDF_1D(nRDFtypes, Nbins);

	int ctr = 0;
	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->step >= (startStep + eq_steps) and TRAJ->step <= (endStep + eq_steps))
		{
			ctr++;
			printf("Processing step %ld, frame %d\n", TRAJ->step - eq_steps, ctr);
			computeRDF_1D(RDF_x_y, ATOMS, BOX, nRDFtypes, Rcut, binW, option);
		}
	}

	// ---------- RDF normalization ----------
	float FAC[nRDFtypes] = {1.0, xA*xA, xA*xB, xB*xA, xB*xB};
	float normx = 1.0, normy = 1.0;

	if(strcmp(option, "symm") == 0)
	{
		normx = rho * TRAJ->nAtoms * (2.0*binW[0]) * (2.0*binW[1]) * ctr;
		normy = normx;
	}
	else if(strcmp(option, "asymm") == 0)
	{
		normx = rho * TRAJ->nAtoms * (1.0*binW[0]) * (2.0*binW[1]) * ctr;
		normy = rho * TRAJ->nAtoms * (2.0*binW[0]) * (1.0*binW[1]) * ctr;
	}

	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < nRDFtypes; j++)
		{
			for(int k = 0; k < Nbins[i]; k++)
			{
				if(i == 0) RDF_x_y[i][j][k] /= (FAC[j] * normx);
				else if(i == 1) RDF_x_y[i][j][k] /= (FAC[j] * normy);
			}
		}	
	}

	printf("RDF computed and averaged over %d frames.\n", ctr);

	TRAJ -> write2file(RDF_x_y, nRDFtypes, Nbins, binW, option);
	TRAJ -> closeTrajectory();
}

float*** buildRDF_1D(int nRDFtypes, int Nbins[])
{
	float ***RDF_x_y = new float**[2];
	for(int i = 0; i < 2; i++)
	{
		RDF_x_y[i] = new float*[nRDFtypes];
		for(int j = 0; j < nRDFtypes; j++)
		{
			RDF_x_y[i][j] = new float[Nbins[i]];
			for(int k = 0; k < Nbins[i]; k++)
				RDF_x_y[i][j][k] = 0.0;
		}
	}

	return(RDF_x_y);
}

void computeRDF_1D(float ***RDF_x_y, atom_style *ATOMS, System *BOX, int nRDFtypes, float Rcut[], float binW[], char option[])
{
	float **RDFx = RDF_x_y[0];
	float **RDFy = RDF_x_y[1];

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

			if(abs(dxij) < Rcut[0] and abs(dyij) < Rcut[1])
			{

				// printf("Atom %d(%c) and %d(%c) are separated by dx=%f and dy=%f\n", i, ATOMS[i].id, j, ATOMS[j].id, dxij, dyij);
				// ---------- Bins symmetric about origin ----------
				if(strcmp(option, "symm") == 0)
				{
					if(abs(dyij) < binW[1])
					{
						int bin_x = int(abs(dxij) / binW[0]);

						RDFx[0][bin_x] += 2.0;

						if(nRDFtypes > 1)
						{
							if(ATOMS[i].id == ATOMS[j].id)
							{
								if(ATOMS[i].id == 'N')
									RDFx[1][bin_x] += 2.0;

								else if(ATOMS[i].id == 'O')
									RDFx[4][bin_x] += 2.0;
							}

							else
							{
								RDFx[2][bin_x] += 1.0;
								RDFx[3][bin_x] += 1.0;
							}
						}
					}

					if(abs(dxij) < binW[0])
					{
						int bin_y = int(abs(dyij) / binW[1]);

						RDFy[0][bin_y] += 2.0;

						if(nRDFtypes > 1)
						{
							if(ATOMS[i].id == ATOMS[j].id)
							{
								if(ATOMS[i].id == 'N')
									RDFy[1][bin_y] += 2.0;

								else if(ATOMS[i].id == 'O')
									RDFy[4][bin_y] += 2.0;
							}

							else
							{
								RDFy[2][bin_y] += 1.0;
								RDFy[3][bin_y] += 1.0;
							}
						}
					}
				}

				// ---------- Bin asymmetric about origin ---------- 
				else if(strcmp(option, "asymm") == 0)
				{
					if(abs(dyij) < binW[1])
					{
						int bin_xi = int(( dxij + Rcut[0]) / binW[0]);
						int bin_xj = int((-dxij + Rcut[0]) / binW[0]);

						RDFx[0][bin_xi] += 1.0;
						RDFx[0][bin_xj] += 1.0;

						if(nRDFtypes > 1)
						{
							if(ATOMS[i].id == ATOMS[j].id)
							{
								if(ATOMS[i].id == 'N')
								{
									RDFx[1][bin_xi] += 1.0;
									RDFx[1][bin_xj] += 1.0;
								}

								else if(ATOMS[i].id == 'O')
								{
									RDFx[4][bin_xi] += 1.0;
									RDFx[4][bin_xj] += 1.0;	
								}
							}

							else
							{
								if(ATOMS[i].id == 'N')
								{
									RDFx[2][bin_xi] += 1.0;
									RDFx[3][bin_xj] += 1.0;	
								}

								else if(ATOMS[i].id == 'O')
								{
									RDFx[3][bin_xi] += 1.0;
									RDFx[2][bin_xj] += 1.0;
								}
							}
						}
					}

					if(abs(dxij) < binW[0])
					{
						int bin_yi = int(( dyij + Rcut[1]) / binW[1]);
						int bin_yj = int((-dyij + Rcut[1]) / binW[1]);

						RDFy[0][bin_yi] += 1.0;
						RDFy[0][bin_yj] += 1.0;

						if(nRDFtypes > 1)
						{
							if(ATOMS[i].id == ATOMS[j].id)
							{
								if(ATOMS[i].id == 'N')
								{
									RDFy[1][bin_yi] += 1.0;
									RDFy[1][bin_yj] += 1.0;
								}

								else if(ATOMS[i].id == 'O')
								{
									RDFy[4][bin_yi] += 1.0;
									RDFy[4][bin_yj] += 1.0;	
								}
							}

							else
							{
								if(ATOMS[i].id == 'N')
								{
									RDFx[2][bin_yi] += 1.0;
									RDFx[3][bin_yj] += 1.0;	
								}

								else if(ATOMS[i].id == 'O')
								{
									RDFx[3][bin_yi] += 1.0;
									RDFx[2][bin_yj] += 1.0;
								}
							}
						}						
					}
				}	
			}
		}
	}
}

void analysis::Trajectory::write2file(float ***RDF_x_y, int nRDFtypes, int Nbins[], float binW[], char option[])
{
	float *xn = new float[Nbins[0]];
	float *yn = new float[Nbins[1]];
	float **RDFx = RDF_x_y[0];
	float **RDFy = RDF_x_y[1];

	int fac = 0;
	if(strcmp(option, "asymm") == 0) fac = 1;

	for(int i = 0; i < Nbins[0]; i++) xn[i] = (i - fac*Nbins[0]/2 + 0.5)*binW[0];
	for(int i = 0; i < Nbins[1]; i++) yn[i] = (i - fac*Nbins[1]/2 + 0.5)*binW[1];

	remove(fpathO);

	fileO = fopen(fpathO, "w");

	if(nRDFtypes == 1)
	{
		fprintf(fileO, "bin xn g_all\n");

		for(int i = 0; i < Nbins[0]; i++)
			fprintf(fileO, "%d %f %f\n", i + 1, xn[i], RDFx[0][i]);

		fprintf(fileO, "bin yn, g_all\n");

		for(int i = 0; i < Nbins[1]; i++)
			fprintf(fileO, "%d %f %f\n", i + 1, yn[i], RDFy[0][i]);
	}

	else
	{
		fprintf(fileO, "bin xn g_all gAA gAB gBA gBB\n");

		for(int i = 0; i < Nbins[0]; i++)
			fprintf(fileO, "%d %f %f %f %f %f %f\n", i + 1, xn[i], RDFx[0][i], RDFx[1][i], RDFx[2][i], RDFx[3][i], RDFx[4][i]);

		fprintf(fileO, "bin yn g_all gAA gAB gBA gBB\n");

		for(int i = 0; i < Nbins[1]; i++)
			fprintf(fileO, "%d %f %f %f %f %f %f\n", i + 1, yn[i], RDFy[0][i], RDFy[1][i], RDFy[2][i], RDFy[3][i], RDFy[4][i]);
	}

	fclose(fileO);

	delete[] xn;
	delete[] yn;
}