// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file src/protocols/hydrate/Hydrate.cc
/// @brief A protocol to add explicit water molecules.
/// @detailed
/// @author Joaquin Ambia Garrido (jgarrido@bcm.edu)


// Protocols
#include <protocols/hydrate/Hydrate.hh>
#include <protocols/hydrate/hydrate_functions.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/AtomTreeDiffJobOutputter.hh>
#include <protocols/jd2/Job.hh>

// Core
//#include <core/chemical/AA.hh>
//#include <core/pack/rotamer_set/RotamerCouplings.hh>
//#include <core/conformation/Residue.hh>
//#include <core/conformation/ResidueMatcher.hh>
//#include <core/conformation/ResidueFactory.hh>
//#include <core/chemical/ResidueTypeSet.hh>
//#include <core/chemical/ResidueSelector.hh>
//#include <core/chemical/VariantType.hh>
//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/task/TaskFactory.hh>
//#include <core/pack/task/operation/TaskOperations.hh>
//#include <core/pack/task/ResfileReader.hh> //yumeng 03/15/13
//#include <core/pack/min_pack.hh>
#include <core/pack/pack_rotamers.hh>
//#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
//#include <core/pack/rotamer_set/WaterPackingInfo.hh>
//#include <core/pose/datacache/CacheableDataType.hh>
//#include <core/pose/Pose.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

//#include <core/kinematics/FoldTree.hh>
//#include <core/kinematics/MoveMap.hh>
//#include <core/import_pose/import_pose.hh>
//#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>

//#include <core/scoring/constraints/util.hh>
//#include <core/scoring/hbonds/HBondOptions.hh>
//#include <core/scoring/methods/EnergyMethodOptions.hh>
//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/Energies.hh>
//#include <core/scoring/EnergyGraph.hh>
//#include <core/types.hh>

// Basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/keys/constraints.OptionKeys.gen.hh>
//#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/hydrate.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// Numeric headers

// C++ headers
#include <iostream>
#include <string>

// Construct tracer.
static basic::Tracer TR("protocols.hydrate.Hydrate");


namespace protocols {
namespace hydrate {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
// Default constructor
/// @details  
Hydrate::Hydrate(): Mover(),
	score_fxn_(new core::scoring::ScoreFunction()),
	main_task_factory_(new core::pack::task::TaskFactory())
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::hydrate;
	using namespace pack::task;

	TR << "Initializing hydrate protocol" << std::endl;

	// Setting up task
	main_task_factory_->push_back( new operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory_->push_back( new operation::ReadResfile );
	} else {
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory_->push_back( rtrop );
	}

	// Setting up scorefxn and constraints
	score_fxn_ = scoring::getScoreFunction();
	if ( option[ in::file::centroid_input ].user() ) {
		scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	} else {
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	}

	// Setting up options
	if (option[ OptionKeys::packing::resfile ].user())
		protein_flexibility_ = "resfile";
	else
		protein_flexibility_ = option[ protein_flexibility ]();
	hydrate_all_ = option[ hydrate_all ]();
	partial_hydrate_dew_ = option[ partial_hydrate_dew ]();
	short_residence_time_mode_ = option[short_residence_time_mode ]();
	near_water_threshold_ = option[ near_water_threshold ](); 	//yumeng
	minimize_bb_where_packing_ = option[ minimize_bb_where_packing ]();		// yumeng

	type("Hydrate");
}

Hydrate::Hydrate( 
	scoring::ScoreFunctionOP scorefxn,
	std::string sf_type,
	std::string protein_flexibility
): Mover(),
	main_task_factory_(new core::pack::task::TaskFactory())
{	
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::hydrate;
	using namespace pack::task;

	TR << "Initializing hydrate protocol" << std::endl;

	// Setting up task
	main_task_factory_->push_back( new operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory_->push_back( new operation::ReadResfile );
	} else {
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory_->push_back( rtrop );
	}

	// Setting up scorefxn and constraints
	if ( sf_type == "input" ){
		score_fxn_ = scorefxn;
	} else if ( sf_type == "explicit"){
	} else if ( sf_type == "hybrid"){
	} else {
		TR << "WARNING!!!\n\n Invalid sf_type in Hydrate constructor. Valid options are: input, explicit and hybrid.\n";
		TR << "Using input sf as default." << std::endl;
		score_fxn_ = scorefxn;
	}

	if ( option[ in::file::centroid_input ].user() ) {
		scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	} else {
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	}

	// Setting up options
	protein_flexibility_ = protein_flexibility;
	hydrate_all_ = option[ hydrate_all ]();
	partial_hydrate_dew_ = option[ partial_hydrate_dew ]();
	short_residence_time_mode_ = option[short_residence_time_mode ]();
	near_water_threshold_ = option[ near_water_threshold ](); 	//yumeng
	minimize_bb_where_packing_ = option[ minimize_bb_where_packing ]();		// yumeng

	type("Hydrate");
}

// Copy constructor
Hydrate::Hydrate(Hydrate const & hyd): Mover(hyd)
{
	copy_data(*this, hyd);
}

// Assignment operator
Hydrate &
Hydrate::operator=(Hydrate const & hyd)
{
	// Abort self-assignment.
	if (this == &hyd) {
		return *this;
	}

	Mover::operator=(hyd);
	copy_data(*this, hyd);
	return *this;
}

// Destructor
Hydrate::~Hydrate() {}


// Mover methods
std::string
Hydrate::get_name() const
{
	return type();
}

protocols::moves::MoverOP
Hydrate::clone() const
{
	return new Hydrate(*this);
}

protocols::moves::MoverOP
Hydrate::fresh_instance() const
{
	return new Hydrate();
}


/// @details  
/// @param    <pose>: the structure to be moved
/// @remarks  
void
Hydrate::apply(Pose & pose)
{
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys::hydrate;

	TR << "Applying " << get_name() << std::endl;

	if ( option[ just_score ]() ){
		water_specific_scores( pose, *score_fxn_ );
		show_water_hb_network( pose );
		return;
	}

	set_water_info_and_add_de_novo_water( pose );

	///////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//		Packing (Double edge water rotamers)
	//
	//	Packing, during this first round, de novo waters generate rotamers with two optimized
	//	hydrogen bonds (double edge (dew)).
	// Some de novo water molecules will not be included in this dew packing round, they will
	// be included in the sew round.
	if ( short_residence_time_mode_  ) partial_hydrate_dew_ = 0;
	if ( partial_hydrate_dew_ < 1 ) set_dew_waters_not_to_be_included(pose, partial_hydrate_dew_);
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line(); // -ex1 -ex2  etc.
	kinematics::MoveMap mm;
	set_task_and_movemap( pose, protein_flexibility_, task, mm, minimize_bb_where_packing_ );
	if ( !short_residence_time_mode_ ){
		pack::pack_rotamers_loop( pose, *score_fxn_, task, 25 );	
		remove_high_energy_water_molecules( pose, *score_fxn_ );
	}

	// If running in short residence time mode, all waters are added as single edge (enforced), 
	// and then packed without anchor atom, like if they were input waters
	if ( short_residence_time_mode_ ) enforce_all_waters( pose );

	// Single Edge (anchor) Water (sew) rotamers packing
	//
	// Any de novo water molecule that could not find an energetically favorable position near the protein
	// in the previous packing, now attempts to find it by constructing rotamers with a single
	// optimized hydrogen bond (single edge (sew)) to its anchor. 
	// This does not imply that the rotamers will have only that one hb.
	pack::task::PackerTaskOP sew_task( pack::task::TaskFactory::create_packer_task( pose ));
	sew_task->initialize_from_command_line(); // -ex1 -ex2  etc.
	(*score_fxn_)(pose);
	get_ready_for_sew_packing( pose, sew_task);
	pack::pack_rotamers_loop( pose, *score_fxn_, sew_task, 25);
	if ( !short_residence_time_mode_ ) remove_high_energy_water_molecules( pose, *score_fxn_);
	// If running in short residence time mode, once water moleculer are added as sew, we pack them
	// concurrently with the rest of the protein considered flexible, using input water rotamers (anchorless)
	if ( short_residence_time_mode_ ){
		remove_all_anchors( pose );
		pack::task::PackerTaskOP task_strm( pack::task::TaskFactory::create_packer_task( pose ));
		task_strm->initialize_from_command_line(); // -ex1 -ex2  etc.
		kinematics::MoveMap mm_strm;
		set_task_and_movemap (pose, protein_flexibility_, task_strm, mm_strm, minimize_bb_where_packing_);
		pack::pack_rotamers_loop( pose, *score_fxn_, task_strm, 25);	
		remove_high_energy_water_molecules( pose, *score_fxn_);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//		Minimizing
	//
	if ( option[ min_backbone_file ].user() ) {
		set_bb_movemap( pose, option[ min_backbone_file ](), mm);
	}
	if ( option[ show_derivatives_check ]() ){
		TR << "Minimizing and showing derivatives" << std::endl;
		optimization::AtomTreeMinimizer().run( pose, mm, *score_fxn_, optimization::MinimizerOptions("dfpmin",0.001,true, 
			true, true));
	} else {
		TR << "Minimizing" << std::endl;
		optimization::AtomTreeMinimizer().run( pose, mm, *score_fxn_, optimization::MinimizerOptions("dfpmin",0.001,true));
	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	//
	//		Removing non buried water molecules
	//
	remove_non_buried_wat(pose);

	/////////////////////////////////////////////////////////////////////////
	//
	//		Adding over coordinated hydrogen bonds scores
	//		And if the hybrid score function is used, it is considered here
	//		[[This must be the last action of the apply function]]
	//
	water_specific_scores( pose, *score_fxn_ );
	show_water_hb_network( pose );
}


void
Hydrate::copy_data(
		Hydrate hyd_to,
		Hydrate hyd_from)
{
	hyd_to.score_fxn_ = hyd_from.score_fxn_;
	hyd_to.main_task_factory_ = hyd_from.main_task_factory_;
}

}  // namespace hydrate
}  // namespace protocols

