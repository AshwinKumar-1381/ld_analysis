// irvingKirkwood.h

#ifndef IK_H
#define IK_H

#include "./library.h"
#include "analysis.h"
#include "../../LD/LD-cpp/src/interactions.h"

using namespace analysis;
using namespace program;

namespace analysis
{
	void computeInteractionPressure(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, pair_style *INTERACTION, Bin1D **Pin);
	void computeKineticPressure(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, Bin1D **Pkin);
	void computeSwimPressure(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, Bin1D *Pswim, float PeA, float PeB);
}

#endif /*IK_H*/