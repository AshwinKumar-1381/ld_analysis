// irvingKirkwood.cpp

#include "irvingKirkwood.h"

void insertSnippet1(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, pair_style *INTERACTION, Bin1D **Pin, float rxi, float ryi, float rxj, float ryj)
{
	float delxCom = TRAJ->xCom - 0.5*BOX->Lx;
	float binW = Pin[0]->binWidth;

	float dx = rxi - rxj;
	float dy = ryi - ryj;

	BOX->checkMinImage(&dx, &dy);
	float r2ij = dx*dx + dy*dy;

	float F = INTERACTION->force(r2ij);
	float Fx = F*dx;
	float Fy = F*dy;

	int binIDi = int((rxi - delxCom)/binW);
	int binIDj = int((rxj - delxCom)/binW);

	if(binIDi == binIDj)
	{
		Pin[0] -> bin[binIDi] += (dx * Fx) / binW;
		Pin[1] -> bin[binIDi] += (dx * Fy) / binW;
		Pin[2] -> bin[binIDi] += (dy * Fx) / binW;
		Pin[3] -> bin[binIDi] += (dy * Fy) / binW;
	}

	else
	{
		int bin1, bin2;
		float rx1, rx2, a, b;

		if(binIDi < binIDj)
		{
			bin1 = binIDi;
			bin2 = binIDj;
			rx1 = rxi;
			rx2 = rxj;
		}
		else
		{
			bin1 = binIDj;
			bin2 = binIDi;
			rx1 = rxj;
			rx2 = rxi;
		}

		if(bin2 - bin1 == int(BOX->Lx/binW) - 1)
		{
			a = abs(rx1);
			b = abs(BOX->Lx - rx2);
			a /= (abs(dx) * binW);
			b /= (abs(dx) * binW);
		}
		else
		{
			a = abs((bin1+1)*binW - rx1);
			b = abs(bin2*binW - rx2);
			a /= (abs(dx) * binW);
			b /= (abs(dx) * binW);

			for(int k = bin1 + 1; k < bin2; k++)
			{
				Pin[0] -> bin[k] += (dx * Fx) / abs(dx);
				Pin[1] -> bin[k] += (dx * Fy) / abs(dx);
				Pin[2] -> bin[k] += (dy * Fx) / abs(dx);
				Pin[3] -> bin[k] += (dy * Fy) / abs(dx);
			}
		}

		Pin[0] -> bin[bin1] += (dx * Fx * a);
		Pin[1] -> bin[bin1] += (dx * Fy * a);
		Pin[2] -> bin[bin1] += (dy * Fx * a);
		Pin[3] -> bin[bin1] += (dy * Fy * a);

		Pin[0] -> bin[bin2] += (dx * Fx * b);
		Pin[1] -> bin[bin2] += (dx * Fy * b);
		Pin[2] -> bin[bin2] += (dy * Fx * b);
		Pin[3] -> bin[bin2] += (dy * Fy * b);
	}
}

void analysis::computeInteractionPressure(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, pair_style *INTERACTION, Bin1D **Pin)
{
	BOX->buildCellList(ATOMS); 

	for(int icell = 1; icell <= BOX->ncells; icell++)
	{
		int i = BOX->HEAD[icell];
		while (i != 0)
		{
			int ii = i - 1;
			float rxi = ATOMS[ii].rxt1;
			float ryi = ATOMS[ii].ryt1;

			int j = BOX->LIST[i];
			while (j != 0)
			{
				int jj = j - 1;
				float rxj = ATOMS[jj].rxt1;
				float ryj = ATOMS[jj].ryt1;

				float dx = rxi - rxj;
				float dy = ryi - ryj;
				float r2ij = dx*dx + dy*dy;

				if(r2ij <= INTERACTION->rcut*INTERACTION->rcut)
					insertSnippet1(TRAJ, ATOMS, BOX, INTERACTION, Pin, rxi, ryi, rxj, ryj);
				 
				j = BOX->LIST[j];
			}
			i = BOX->LIST[i];
		}
	}

	int nNbors = 4;
	for(int icell = 1; icell <= BOX->ncells; icell++)
	{
		int i = BOX->HEAD[icell];
		while (i != 0)
		{
			int ii = i - 1;
			float rxi = ATOMS[ii].rxt1;
			float ryi = ATOMS[ii].ryt1;

			int icell_index = nNbors*(icell - 1);
			for(int nbor = 1; nbor <= nNbors; nbor++)
			{
				int jcell = BOX->MAPS[icell_index + nbor];

				int j = BOX->HEAD[jcell];
				while(j != 0)
				{
					int jj = j - 1;
					float rxj = ATOMS[jj].rxt1;
					float ryj = ATOMS[jj].ryt1;

					float dx = rxi - ATOMS[jj].rxt1;
					float dy = ryi - ATOMS[jj].ryt1;

					BOX->checkMinImage(&dx, &dy);
					float r2ij = dx*dx + dy*dy;

					if(r2ij <= INTERACTION->rcut*INTERACTION->rcut)
						insertSnippet1(TRAJ, ATOMS, BOX, INTERACTION, Pin, rxi, ryi, rxj, ryj);

					j = BOX->LIST[j];
				}
			}
			i = BOX->LIST[i];
		}
	}
}	

void analysis::computeKineticPressure(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, Bin1D **Pkin)
{
	float delxCom = TRAJ->xCom - 0.5*BOX->Lx;

	for(int i = 0; i < TRAJ -> nAtoms; i++)
	{
		float vx = ATOMS[i].px;
		float vy = ATOMS[i].py;

		Pkin[0] -> addToBin(ATOMS[i].rxt1 - delxCom, vx * vx);
		Pkin[1] -> addToBin(ATOMS[i].rxt1 - delxCom, vx * vy);
		Pkin[2] -> addToBin(ATOMS[i].rxt1 - delxCom, vy * vx);
		Pkin[3] -> addToBin(ATOMS[i].rxt1 - delxCom, vy * vy);
	}
}

void analysis::computeSwimPressure(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, Bin1D *Pswim, float PeA, float PeB)
{
	float delxCom = TRAJ->xCom - 0.5*BOX->Lx;

	for(int i = 0; i < TRAJ -> nAtoms; i++)
	{
		float rx = abs(ATOMS[i].rxt1 - 0.5*BOX->Lx);

		if(ATOMS[i].id == 'N')
			Pswim -> addToBin(ATOMS[i].rxt1 - delxCom, rx*PeA);			
		else
			Pswim -> addToBin(ATOMS[i].rxt1 - delxCom, rx*PeB);
	}
}