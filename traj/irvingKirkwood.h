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
	void computeInteractionPressure(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, pair_style *INTERACTION, particleBin **Pin);
	void computeSwimPressure(Trajectory *TRAJ, atom_style *ATOMS, System *BOX, particleBin *Pswim);
}

#endif /*IK_H*/