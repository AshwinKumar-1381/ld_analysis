// makeRestart.cpp
// Generates a data file for LAMMPS restart simulations by selecting a particular frame from an .xyz trajectory file
// 
// Date created  : 05.05.25
// Last modified : 06.05.25

#include "analysis.h"

using namespace analysis;

int main(int argc, char* argv[])
{
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;

	char typeLabels[nAtomTypes] = {'N', 'O'};

	float dt = 5e-4;
	long eqStep = long(1e7);
	long chooseStep = long(8e9);

	Trajectory *TRAJ = new Trajectory(dt);
	sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/Apr2025/lane/lmp/Data22/traj2.xyz");
	sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/Apr2025/lane/lmp/Data22/lane.res");
	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ -> nAtoms]; 

	while( !feof(TRAJ -> fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ -> step == (eqStep + chooseStep))
		{
			printf("Step %ld extracted for restart simulations\n", TRAJ->step);
			break;
		}
	}

	TRAJ -> closeTrajectory();

	remove(TRAJ->fpathO);
	TRAJ->fileO = fopen(TRAJ->fpathO, "w");

	fprintf(TRAJ->fileO, "# makeRestart.cpp : LAMMPS data file for restart simulation generated using %ld step from trajectory file %s\n\n", TRAJ->step, TRAJ->fpathI);

	fprintf(TRAJ->fileO, "%d atoms\n", TRAJ->nAtoms);
	fprintf(TRAJ->fileO, "%f %f xlo xhi\n", 0.0, Lx);
	fprintf(TRAJ->fileO, "%f %f ylo yhi\n", 0.0, Ly);
	fprintf(TRAJ->fileO, "%d atom types\n\n", nAtomTypes);

	fprintf(TRAJ->fileO, "Atom Type Labels\n\n");

	for(int i = 0; i < nAtomTypes; i++)
		fprintf(TRAJ->fileO, "%d %c\n", i + 1, typeLabels[i]);

	fprintf(TRAJ->fileO, "\nMasses\n\n");

	for(int i = 0; i < nAtomTypes; i++)
		fprintf(TRAJ->fileO, "%d %f\n", i + 1, 1.0);

	fprintf(TRAJ->fileO, "\nAtoms # atomic\n\n");

	for(int i = 0; i < TRAJ->nAtoms; i++)
		fprintf(TRAJ->fileO, "%d %c %f %f %f\n", i + 1, ATOMS[i].id, ATOMS[i].rxt1, ATOMS[i].ryt1, 0.0);

	fclose(TRAJ->fileO);
}