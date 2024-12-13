// analysis.cpp

#include "analysis.h"

using namespace analysis;

/* ----------------- atomsXYZ members -----------------*/

analysis::atomsXYZ::atomsXYZ(){}
analysis::atomsXYZ::~atomsXYZ(){}

/* ----------------- System members -----------------*/

analysis::System::System(float Lx, float Ly, float nAtoms)
{
	this -> Lx = Lx;
	this -> Ly = Ly;
	this -> nAtoms = nAtoms;

	rcellx = 1.122462048;
	rcelly = rcellx;

	Ncellx = int(this->Lx/rcellx);
	Ncelly = int(this->Ly/rcelly);
	ncells = int(Ncellx*Ncelly);

	rcellx = this->Lx/Ncellx;
	rcelly = this->Ly/Ncelly; 

	for(int i = 0; i < MAXCELL; i++)
	{
		MAPS[i] = 0;
		LIST[i] = 0;
		HEAD[i] = 0;
	}

	buildCellMaps();
}

analysis::System::~System(){}

int analysis::System::cellindex(int ix, int iy)
{   
    if (ix == Ncellx) ix = 0; 
	else if (ix==-1) ix = Ncellx - 1; 
	if (iy == Ncelly) iy = 0;
	else if (iy == -1) iy = Ncelly - 1;

	if(Ncellx >= Ncelly)
		return (1 + ix + iy*Ncellx);
	else
		return (1 + ix + iy*Ncelly);
}

void analysis::System::buildCellMaps()
{
	for (int ix = 0; ix < Ncellx; ix++)
	{
		for (int iy = 0; iy < Ncelly; iy++)
		{
			int imap = 4*(cellindex(ix,iy)-1);
			MAPS[imap+1] = cellindex(ix-1,iy);
			MAPS[imap+2] = cellindex(ix-1,iy+1);
			MAPS[imap+3] = cellindex(ix,iy+1);
			MAPS[imap+4] = cellindex(ix+1,iy+1);
		}
	}
    
    printf("Successfully constructed MAPS array with %d cells.\n", int(Ncellx * Ncelly));
}

void analysis::System::buildCellList(atom_style *ATOMS)
{
	for(int i = 0; i < MAXCELL; i++)
	{
		LIST[i] = 0;
		HEAD[i] = 0;
	}

	for (int i = 1; i <= nAtoms; i++)
	{
		int ii = i - 1;
		int ix = int(ATOMS[ii].rxt1/rcellx);
		int iy = int(ATOMS[ii].ryt1/rcelly);
		
		int icell = cellindex(ix, iy);	        
		LIST[i] = HEAD[icell];
		HEAD[icell] = i;		    
	}	
}

void analysis::System::checkMinImage(float *dx, float *dy)
{
	if(dx != NULL)
	{
		if(*dx >= 0.5*Lx) *dx -= Lx;
		else if(*dx <= -0.5*Lx) *dx += Lx;
	}
	if(dy != NULL)
	{
		if(*dy >= 0.5*Ly) *dy -= Ly;
		else if(*dy <= -0.5*Ly) *dy += Ly; 
	}
}

/* ----------------- particleBin members -----------------*/

analysis::particleBin::particleBin(float Lx, float Ly, float binWidth)
{
	this -> Lx = Lx;
	this -> Ly = Ly;
	this -> binWidth = binWidth;
	nBins = int(Lx/binWidth);
	bin = new float[nBins];
	binavg = new float[nBins];

	zero(bin);
	zero(binavg);
}

analysis::particleBin::~particleBin(){}

void analysis::particleBin::addToBin(float rx, float val)
{
	if(val == 1.0)
		bin[int(rx/binWidth)] += 1;
	else
		bin[int(rx/binWidth)] += val;	
}

void analysis::particleBin::normalize(float *binTonorm, float norm)
{
	if(norm == 0.0) norm = Ly*binWidth;

	for(int i = 0; i < nBins; i++)
		binTonorm[i] /= norm;
}

void analysis::particleBin::zero(float *binTozero)
{
	for(int i = 0; i < nBins; i++)
		binTozero[i] = 0.0;
}

void analysis::particleBin::addBins(float *bin1, float *bin2)
{
	for(int i = 0; i < nBins; i++)
		bin1[i] += bin2[i];
}

/* ----------------- Trajectory members -----------------*/

analysis::Trajectory::Trajectory(float timeStep)
{
	fpathI = new char [500];
	fpathO = new char [500];
	pipeString = new char [500];
	pipeChar = new char [500];
	fileI = nullptr;
	fileO = nullptr;

	frame_nr = -1;
	this->timeStep = timeStep;
	xCom = yCom = zCom = 0.0;
}

analysis::Trajectory::~Trajectory(){}

void analysis::Trajectory::openTrajectory()
{
	fileI = fopen(fpathI, "r");
	if(fileI == NULL)
	{
		printf("Could not open %s. Exiting ...\n", fpathI);
		exit(-1);
	}
	else
	{
		fgets(pipeString, 500, fileI);
		sscanf(pipeString, "%d", &nAtoms);
		rewind(fileI);
		printf("\nTrajectory %s file opened and ready to be read...\n", fpathI);
	}
}

void analysis::Trajectory::closeTrajectory()
{
	printf("Closing trajectory file.\n");
	fclose(fileI);
}

void analysis::Trajectory::readThisFrame(atom_style *ATOMS)
{
	fgets(pipeString, 500, fileI);
	
	fgets(pipeString, 500, fileI);
	sscanf(pipeString, "%*s %d", &step);

	for(int i = 0; i < nAtoms; i++)
	{
		fgets(pipeString, 500, fileI);
		sscanf(pipeString, "%c %f %f %f %f %f", &ATOMS[i].id, &ATOMS[i].rxt1, &ATOMS[i].ryt1, &ATOMS[i].rzt1, &ATOMS[i].px, &ATOMS[i].py);
	}

	frame_nr++;
}

void analysis::Trajectory::computeCom(atom_style *ATOMS)
{
	xCom = yCom = zCom = 0.0;

	for(int i = 0; i < nAtoms; i++)
	{
		xCom += ATOMS[i].rxt1;
		yCom += ATOMS[i].ryt1;
		zCom += ATOMS[i].rzt1;
	}	

	xCom /= nAtoms;
	yCom /= nAtoms;
	zCom /= nAtoms;
}