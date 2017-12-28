// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file src/protocols/hydrate/Hydrate.cc
/// @brief The Hydrate Protocol
/// @detailed
/// @author Joaquin Ambia

#ifndef INCLUDED_protocols_hydrate_hydrate_functions_HH
#define INCLUDED_protocols_hydrate_hydrate_functions_HH

// Protocols
#include <protocols/hydrate/Hydrate.hh>
#include <protocols/moves/Mover.hh>

// Core
#include <core/chemical/AA.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh> //yumeng 03/20/13
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>

// Basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

namespace protocols{
namespace hydrate{

using namespace core;

void
read_hyfile( 
	std::string const & hyfile, 
	utility::vector1< bool > & enforced_V, 
	utility::vector1< bool > & hydrate_V
);

void
hydrate_hyfile( 
	pose::Pose & pose, 
	utility::vector1< bool > const & hydrate_V 
);

void
set_water_info_and_add_de_novo_water( 
	pose::Pose & pose
);

bool
atom_is_hydratable( 
	pose::Pose const & pose, 
	Size const & residue, 
	std::string const & atom
);

bool
atom_is_hydratable( 
	pose::Pose const & pose, 
	Size const & residue, 
	Size const & atom
);

bool
is_inside(
	pose::Pose const & pose, 
	Vector const & xyz
);

// Moves a fraction of water molecules away from the protein, not to include during dew
//
void
set_dew_waters_not_to_be_included(
	pose::Pose & pose, 
	Real const partial_hydrate_dew
);

bool		// yumeng
residue_near_water(
	pose::Pose const & pose, 
	Size const & ii
);


void
set_task_and_movemap(
	pose::Pose const & pose,
	std::string const & protein_flexibility,
	pack::task::PackerTaskOP & task,
	kinematics::MoveMap & mm, 
	bool const & minimize_bb_where_packing
);

void
remove_high_energy_water_molecules( 
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn
);

void
calculate_water_overcoordinated_hb_correction( 
	pose::Pose const & pose, 
	utility::vector1< Real > & water_hb_correction
);

void
enforce_all_waters(
	pose::Pose & pose
);

void
get_ready_for_sew_packing(
	pose::Pose & pose, 
	pack::task::PackerTaskOP & task
);

void
remove_all_anchors(
	pose::Pose & pose
);

void // yumeng
set_bb_movemap(
	pose::Pose const & pose,
	std::string const & mini_backbone_file_name,
	kinematics::MoveMap & mm
);

void
remove_non_buried_wat(
	pose::Pose & pose
);

void
water_specific_scores( 
	pose::Pose & pose, 
	scoring::ScoreFunction const & scorefxn
);

void
show_water_hb_network(
	pose::Pose const & pose
);

/*
void
hydrate_cavities( 
	pose::Pose & pose, 
	utility::vector1< bool > & enforced_V
);

void
calculate_water_entropy_correction( 
	pose::Pose & pose, 
	utility::vector1< Real > & water_S_correction
);

void
remove_exposed_wat(
	pose::Pose & pose
);

void
slide_jumps_to_closest_residue(
	pose::Pose & pose
);

void
set_hyb_movemap(
	pose::Pose & pose,
	kinematics::MoveMap & movemap
);

void
set_two_regions_movemaps(
	pose::Pose & pose,
	kinematics::MoveMap & not_wat_movemap,
	kinematics::MoveMap & wat_movemap
);

void
set_task_with_de_novo_water_using_resfile(
  pose::Pose & pose,
  std::string resfile,
  pack::task::PackerTaskOP & task
);


void
make_hybrid_energy_map(
	pose::Pose & pose,
	scoring::EnergyMap & EM,
	scoring::ScoreFunction & scorefxn
);

void
account_for_entropy(
	pose::Pose & pose,
	scoring::EnergyMap & EM,
	scoring::ScoreFunction & scorefxn
);

*/
}	//	close hydrate
}	//	close hydrate

#endif
