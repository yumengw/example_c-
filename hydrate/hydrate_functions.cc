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

// Protocols
#include <protocols/hydrate/Hydrate.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/AtomTreeDiffJobOutputter.hh>
#include <protocols/jd2/Job.hh>

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

#include <core/scoring/constraints/util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>
#include <core/kinematics/MoveMap.hh> //yuemng 03/20/13
#include <core/kinematics/FoldTree.hh> 
#include <core/pack/task/ResfileReader.hh>

// Basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/hydrate.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

// Numeric 
#include <numeric/NumericTraits.hh>

// C++
#include <fstream>
#include <iostream>
#include <string>
#include <iterator>

// Numeric
#include <numeric/random/random.hh>

static basic::Tracer TR("protocols.hydrate.hydrate_functions");

namespace protocols{
namespace hydrate{

static numeric::random::RandomGenerator RG(12321);

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys::hydrate;

// Read hyfile; it should be a list with two columns, 1 with the keyword "HYD" or "ENF" and 2 with the residue number
void
read_hyfile(
	std::string const & hyfile_name,
	utility::vector1< bool > & enforced_V,
	utility::vector1< bool > & hydrate_V
){

	TR << "Reading hyfile " << hyfile_name << std::endl;

	utility::io::izstream input( hyfile_name.c_str() );
	std::string line;
	if ( !input ) {
		TR << "Can't open file " << hyfile_name << "!" << std::endl;
		exit(0);
	}
	
	while ( getline( input, line ) ) {
		if ( line.substr(0,3).compare("ENF") == 0 ){	  // This comparison means it IS "ENF"
			TR << line << std::endl;
			std::stringstream ss(line);
			std::istream_iterator<std::string> begin(ss);
			std::istream_iterator<std::string> end;
			std::vector<std::string> vstrings(begin, end);
			Size rsd = atoi(vstrings[1].c_str());
			enforced_V[rsd] = true;
		}
		if ( line.substr(0,3).compare("HYD") == 0 ){	  // This comparison means it IS "HYD"
			TR << line << std::endl;
			std::stringstream ss(line);
			std::istream_iterator<std::string> begin(ss);
			std::istream_iterator<std::string> end;
			std::vector<std::string> vstrings(begin, end);
			Size rsd = atoi(vstrings[1].c_str());
			hydrate_V[rsd] = true;
		}
	}
}

// Determines if the location of the Vector is "inside" the pose or not
bool
is_inside(
	pose::Pose const & pose,
	Vector const & xyz
)
{
	const Real PI = numeric::NumericTraits<Real>::pi();

	Size sec_resolution =7;
	Real resolution = sec_resolution;
	Real sec_size_cos_phi = 2/resolution;
	Real sec_size_theta = PI/resolution;
	utility::vector1 < bool > neighbor_sections (2*sec_resolution*sec_resolution);

	for (Size ii = 1; ii<=pose.total_residue(); ++ii){
		for (Size jj = 1; jj<=pose.residue(ii).nheavyatoms(); ++jj){

			Vector neighbor_xyz (pose.residue(ii).xyz(jj));
			Real distance = neighbor_xyz.distance(xyz);

			if (distance <=11 && distance > 0.1){
				Vector direction = ( neighbor_xyz - xyz).normalized();
				Real cos_phi = direction[2];
				Real theta = atan2(direction[1],direction[0]);
				Size sec_cos_phi = floor((cos_phi + 1)/sec_size_cos_phi);
				Size sec_theta = floor((theta + PI)/sec_size_theta);
				Size section = sec_cos_phi*2*sec_resolution + sec_theta +1;
				if (section == 0) section = 1;
				neighbor_sections[ section ] = true;
    	}
		}
	}

	Real total_sec (0);
	for (Size ii = 1; ii<= neighbor_sections.size(); ++ii){
		total_sec += int( neighbor_sections[ii] );
	}
	Real ratio = float(total_sec/neighbor_sections.size());

	if (ratio >= 0.7)
		return true;
	else
		return false;
}

// This function defines an atom as hydratable when it has neighboring positions that are:
// large enough for a water molecule and they are inside.
bool
atom_is_hydratable(
    pose::Pose const & pose,
    Size const & residue,
    std::string const & atom
)
{
	Real resolution = 24;
	Real step_cos_phi = 2/resolution;
	Real step_theta = numeric::NumericTraits<Real>::pi()/resolution;
	Real r_cavity = 1.8;

	Vector ref_x (1,0,0);
	Vector ref_y (0,1,0);
	Vector ref_z (0,0,1);
	Vector atom_xyz (pose.residue(residue).xyz(atom));

	for ( Real ro = 2.8; ro <= 3.2; ro += 0.2)
		for ( Real cos_phi = -1; cos_phi < 0.99; cos_phi += step_cos_phi){
			Real phi = acos(cos_phi);
			for ( Real theta = 0; theta < 6.28; theta += step_theta){

				// Position of the hypothetical solvent atom (O)
				Vector hyp_O (atom_xyz  + ro*( sin(phi)*(ref_x*cos(theta) + ref_y*sin(theta)) + ref_z*cos_phi ));
				if ( atom[1] == 'H'){
					hyp_O = atom_xyz  + (ro - 1.0)*( sin(phi)*(ref_x*cos(theta) + ref_y*sin(theta)) + ref_z*cos_phi );
				}

				if ( !is_inside(pose, hyp_O) ) continue;

				bool has_room (true);
				for (Size ii = 1; ii<=pose.total_residue(); ++ii){
					for (Size jj = 1; jj<=pose.residue(ii).nheavyatoms(); ++jj){    // Ignore H because they might not be well located
						

						if ( ii == residue && jj == pose.residue(ii).atom_index( atom ) ) continue;

						r_cavity = 1.8;
						if (pose.residue(ii).atom_name(jj)[1] == 'C')     r_cavity = 1.7;
						if (pose.residue(ii).atom_name(jj)[1] == 'O')     r_cavity = 1.05;	// smaller than vdw
						if (pose.residue(ii).atom_name(jj)[1] == 'N')     r_cavity = 1.55;
						if (pose.residue(ii).atom_name(jj)[1] == 'F')     r_cavity = 1.47;
						r_cavity += 1.52 - 0.25;     // To account for the solvent radious minus some tolerance

						Vector neighbor_xyz (pose.residue(ii).xyz(jj));
						if ( neighbor_xyz.distance(hyp_O)  < r_cavity){
							has_room = false;
						}
					}
				}
				if (has_room)
					return true;

			}   // theta
		}   // phi
	return false;
}

bool
atom_is_hydratable(
	pose::Pose const & pose,
	Size const & residue,
	Size const & atom
)
{
	return atom_is_hydratable(pose, residue, pose.residue(residue).atom_name(atom));
}


// Add water molecules according to the hyfile
void
hydrate_hyfile(
	pose::Pose & pose,
	utility::vector1< bool > const & hydrate_V
){
	TR << "Hydrating according to hyfile." << std::endl;

	conformation::ResidueOP tp3( conformation::ResidueFactory::create_residue( 
		pose.residue(1).residue_type_set().name_map("TP3") ) );
	pack::rotamer_set::WaterPackingInfoOP water_info ( static_cast< pack::rotamer_set::WaterPackingInfo & >
		( pose.data().get( pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) );

	for ( Size ii = 1; ii <= hydrate_V.size(); ++ii){
		if ( hydrate_V[ii] == true ){

			TR << "Hydrating residue " << ii << std::endl;
			conformation::Residue const & rsd( pose.residue(ii) );

			for ( Size jj=1; jj<= Size(rsd.natoms()); ++jj ) {
				if ( (rsd.atom_type(jj).is_acceptor() || rsd.atom_type(jj).is_polar_hydrogen())){
					
					if ( !option[ OptionKeys::hydrate::attempt_all_polar ]() && !atom_is_hydratable(pose,ii,jj) ) continue;
					Size anchor_atom = jj;
					tp3->set_xyz( "O", rsd.xyz(jj) ); // for neighbor calculation
					pose.append_residue_by_jump( *tp3, 1 );
					Size const pos1( pose.total_residue() );
					TR << "Water: " << pos1 << " appended to residue: " << ii << " atom: ";
					TR << rsd.atom_name(anchor_atom) << std::endl;
					(*water_info)[ pos1 ].anchor_residue( ii );
					(*water_info)[ pos1 ].anchor_atom( rsd.atom_name(anchor_atom ));
					(*water_info)[ pos1 ].aa( rsd.aa() );
					(*water_info)[ pos1 ].nstep( 1 );
				}
			}

		}	// hydrate true
	}
	pose.data().set( pose::datacache::CacheableDataType::WATER_PACKING_INFO, water_info );
}

void
set_water_info_and_add_de_novo_water(
	pose::Pose & pose
){
	if ( pose.data().has( pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) return;

	TR << "Setting up water information." << std::endl;

	Size const & total_input_residues ( pose.total_residue() );
    utility::vector1< bool > enforced_V( total_input_residues, false );
    utility::vector1< bool > hydrate_V( total_input_residues, false );
	if ( option[ hyfile ].user() )	read_hyfile( option[ hyfile ](), enforced_V, hydrate_V);

	pack::rotamer_set::WaterPackingInfoOP water_info( new pack::rotamer_set::WaterPackingInfo() );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii){
		if ( pose.residue(ii).name() == "TP3" ) {
			(*water_info)[ ii ].anchor_atom( "NONE" );
			if ( enforced_V[ii] )	(*water_info)[ ii ].enforced( true );
		}
	}   
	pose.data().set( pose::datacache::CacheableDataType::WATER_PACKING_INFO, water_info );

	if (!option[ hyfile ].user()) return;
	TR << "Adding de novo water molecules." << std::endl;
	if (!option[ hyfile ].user()) hydrate_hyfile(pose, hydrate_V);
}

/*void
add_de_novo_water(
	pose::Pose & pose
){
	if (!option[ hyfile ].user()) return;

	TR << "Adding de novo water molecules." << std::endl;

	Size const & total_input_residues ( pose.total_residue() );
    utility::vector1< bool > enforced_V( total_input_residues, false );
    utility::vector1< bool > hydrate_V( total_input_residues, false );

	read_hyfile(option[ hyfile ](), enforced_V, hydrate_V);
	hydrate_hyfile(pose, hydrate_V);
}
*/
void
set_dew_waters_not_to_be_included(
    pose::Pose & pose, 
    Real const partial_hydrate_dew
){

	TR << "Considering " << partial_hydrate_dew*100 << "% of the de novo water molecules in dew packing" << std::endl;

	pack::rotamer_set::WaterPackingInfoOP water_info
			( static_cast< pack::rotamer_set::WaterPackingInfo & >
			( pose.data().get( pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) );
	conformation::ResidueOP tp3( core::conformation::ResidueFactory::create_residue(
			pose.residue(1).residue_type_set().name_map("TP3") ) );

	numeric::xyzVector < Real > awayO(0,0,0);
	numeric::xyzVector < Real > awayH1 (0.9572,0,0);
	numeric::xyzVector < Real > awayH2 (-0.2399872,0.92662721,0);
	numeric::xyzVector < Real > far_away(10000,10000,10000);

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii){
		if ( pose.residue(ii).name() != "TP3") continue;
		if ( (*water_info)[ii].anchor_atom() == "NONE") continue;

		if (RG.uniform() >= partial_hydrate_dew){
			TR << "Not including water " << ii << " in first round of hydration (dew)." << std::endl;
			tp3->set_xyz("O", ii*far_away + awayO);
			tp3->set_xyz("H1", ii*far_away + awayH1);
			tp3->set_xyz("H2", ii*far_away + awayH2);
			pose.replace_residue(ii,*tp3,false);
		}
	}
}

//yumeng 03/15/2013
//calculate whether the residue is near water
//considering all heavy atoms
bool
residue_near_water(
	pose::Pose const & pose, 
	Size const & ii
){
	Real near_water_threshold ( option[ OptionKeys::hydrate::near_water_threshold ]() );
	for (Size jj = 1; jj <= pose.total_residue(); ++jj){
		if (ii == jj) continue;
		if (pose.residue(jj).name() == "TP3") {
		Vector water_oxygen_xyz (pose.residue(jj).xyz(1));
			for (Size kk = 1; kk <= pose.residue(ii).nheavyatoms(); ++kk) {
				Vector atom_heavy_xyz (pose.residue(ii).xyz(kk));
				if (atom_heavy_xyz.distance(water_oxygen_xyz) <= near_water_threshold) 
					return true;
			}
		}
	}
	return false;
}

// yumeng
void
set_task_and_movemap( 
	pose::Pose const & pose, 
	std::string const & protein_flexibility, 
	pack::task::PackerTaskOP & task, 
	kinematics::MoveMap & mm, 
	bool const & minimize_bb_where_packing
){
	TR << "Setting task and movemap with protein flexibility: " << protein_flexibility << std::endl;

	kinematics::FoldTree f = pose.fold_tree();

	//	First deal with water, which is independent of protein_flexibility
	// All water near the protein will be packed at this stage
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii){
		mm.set_bb(ii, false);									// first set all bb moves as false (not water specific)
		if ( pose.residue(ii).name() != "TP3" ) continue;
		if ( pose.residue(ii).xyz(1).x() > 10000 ) {    // It's away and should not be packed at this stage
			task->nonconst_residue_task(ii).prevent_repacking();
		}
		else
			task->nonconst_residue_task(ii).restrict_to_repacking();
		mm.set_chi(ii, true);
	}
	// All jumps connecting water will always be flexible
	// For now, we don't expect any other jumps to be flexible
	for (Size jj = 1; jj <= f.num_jump() ; ++jj){
		if ( pose.residue(f.downstream_jump_residue(jj)).name() == "TP3" ){
			mm.set_jump(jj, true);
		}
	}

	if (protein_flexibility == "not") {
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii){
			if ( pose.residue(ii).name() == "TP3" ) continue;
			task->nonconst_residue_task(ii).prevent_repacking();
			mm.set_chi(ii,false);
		}
	}

	if (protein_flexibility == "all"){
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii){
			if ( pose.residue(ii).name() == "TP3" ) continue;
			task->nonconst_residue_task(ii).restrict_to_repacking();
			mm.set_chi(ii, true);
			if ( minimize_bb_where_packing ) mm.set_bb(ii, true);
		}
	}

	if ( protein_flexibility == "near_water" ){
		for (Size ii = 1; ii <= pose.total_residue(); ++ii){
			if ( pose.residue(ii).name() == "TP3" ) continue;
			if ( residue_near_water(pose, ii) ) {
				task->nonconst_residue_task(ii).restrict_to_repacking();
				mm.set_chi(ii, true);
				if ( minimize_bb_where_packing ) mm.set_bb(ii, true);
			}
			else {
				task->nonconst_residue_task(ii).prevent_repacking();
				mm.set_chi(ii, false);
			}
		}
	}

	if (protein_flexibility == "resfile"){
		pack::task::PackerTaskOP temp_task( pack::task::TaskFactory::create_packer_task( pose ));
		pack::task::parse_resfile(pose, *temp_task);				// Must be here to have water info in task
		for (Size ii = 1; ii <= pose.total_residue(); ++ii) {
			if ( temp_task->being_packed(ii) || pose.residue(ii).name() == "TP3" ) {
//				task->nonconst_residue_task(ii).restrict_to_repacking();	// To avoid problems with design, now possible
				mm.set_chi(ii, true); //check whether the resi is modified by packer
				if ( minimize_bb_where_packing ) mm.set_bb(ii, true);
			}
			else {
				mm.set_chi(ii, false);
				task->nonconst_residue_task(ii).prevent_repacking();
			}

		}
	}
}

void
calculate_water_overcoordinated_hb_correction(
	pose::Pose const & pose,
	utility::vector1< Real > & hb_correction
){
	using namespace core::scoring;

	Real const hbond_threshold( option[ OptionKeys::hydrate::hbond_threshold]() );

	Real bb_sc_weight = pose.energies().weights()[ hbond_bb_sc ];
	Real sc_weight = pose.energies().weights()[ hbond_sc ];
	if ( pose.energies().weights()[ hbond_wat ] != 0 ) {
		bb_sc_weight = pose.energies().weights()[ hbond_wat ];
		sc_weight = pose.energies().weights()[ hbond_wat ];
	}

	scoring::hbonds::HBondSet  hbond_set;
	hbond_set.clear();
	scoring::hbonds::fill_hbond_set( pose, false, hbond_set );

	for (Size ii= 1; ii <= pose.total_residue(); ++ii)  // We make sure we start from 0 for all residues
		hb_correction[ii] = 0;

	Size hb_count(0);
	utility::vector1< Real > hb_energy(10);
	utility::vector1< Size > hb_don_res(10); // we also keep track of the donor residue to calculate corrections
	utility::vector1< Size > hb_acc_res(10); // we also keep track of the acceptr residue to calculate corrections
	utility::vector1< Size > hb_involved(10);

	// For simplicity, everything is with respect to the acceptors
	for (Size ii= 1; ii <= pose.total_residue(); ++ii){ 			// Go over all residues
		for ( Size jj=1; jj<= pose.residue(ii).natoms(); ++jj){         // and all atoms

			if (!pose.residue(ii).atom_type(jj).is_acceptor()) continue;	// if not acceptor, move on

			hb_count = 0;

			// Go over all hydrogen bonds and find the ones where atom jj is involved in
			// We only correct for hydrogen bonds involving water
			for ( Size hb=1; hb<= Size(hbond_set.nhbonds()); ++hb ){
				if ( hbond_set.hbond(hb).acc_res() == ii && hbond_set.hbond(hb).acc_atm() == jj // atom of interest
						&& ( pose.residue(hbond_set.hbond(hb).acc_res()).name() == "TP3"
							|| pose.residue(hbond_set.hbond(hb).don_res()).name() == "TP3" )  // hb involves water
						&& hbond_set.hbond(hb).energy()*hbond_set.hbond(hb).weight() < hbond_threshold  // relevant hb
				){

					// 'weight' is the sf weight, hbond_set.hbond(hb).weight() is the hb_enviroment (burial dependent) weight
					Size weight = sc_weight;
					if ( hbond_set.hbond(hb).acc_atm_is_protein_backbone() ||hbond_set.hbond(hb).don_hatm_is_protein_backbone() ){
						weight = bb_sc_weight;
					}

					// For this atom we count the number of hbonds it is involved in (hb_count),
					// and keep track of their energy (hb_energy) and interacting residue (hb_don_res)
					++hb_count;
					hb_energy[hb_count] = hbond_set.hbond(hb).energy()*hbond_set.hbond(hb).weight()*weight;
					hb_don_res[hb_count] = hbond_set.hbond(hb).don_res();

					// We always keep the 2 strongest hbs in the first two places of hb_energy and hb_don_res
					// All others will be deleted
					if ( hb_energy[hb_count] < hb_energy[2] && hb_count >= 3 ){
						Real temp = hb_energy[hb_count];             // Here we make sure we have the first two hb energies
						hb_energy[hb_count] = hb_energy[2];		     // in the first two places 
						hb_energy[2] = temp;
						Real temp_don = hb_don_res[hb_count];                // also the doner residue
						hb_don_res[hb_count] = hb_don_res[2];
						hb_don_res[2] = temp_don;
					}
					if ( hb_energy[2] < hb_energy[1] && hb_count >= 2 ){
						Real temp = hb_energy[2];
						hb_energy[2] = hb_energy[1];
						hb_energy[1] = temp;
						Real temp_don = hb_don_res[2];
						hb_don_res[2] = hb_don_res[1];
						hb_don_res[1] = temp_don;
					}
				}
			}	// all hydrogen bonds; Finished collecting info on hb_count, hb_energy[] and hb_don_res[] for acceptor

			// Corrections are pertinent only if the acceptor makes two or more hb
			// They are divided in between the residues involved to be consistent with Rosetta standards for scoring
			if (hb_count >= 2){

				// The acceptor is N and can have only one hb
				if ( pose.residue(ii).atom_name(jj)[1] == 'N' ){
					for (Size kk = 2; kk<=hb_count; ++kk){
						hb_correction[ii] -= hb_energy[kk]/2.0;
						hb_correction[hb_don_res[kk]] -= hb_energy[kk]/2.0;
					}
				}
				else{				// The acceptor is O
					if ( hb_energy[2] < 0.2*hb_energy[1] ){ // Make sure the interference is comparable
						hb_correction[ii] -= 0.2*(hb_energy[1] + hb_energy[2])/2.0;
						hb_correction[hb_don_res[1]] -= 0.2*hb_energy[1]/2.0;
						hb_correction[hb_don_res[2]] -= 0.2*hb_energy[2]/2.0;
					}
					else{
						hb_correction[ii] -= hb_energy[2]/2.0;
						hb_correction[hb_don_res[2]] -= hb_energy[2]/2.0;
					}
					for (Size kk = 3; kk<=hb_count; ++kk){
						hb_correction[ii] -= hb_energy[kk]/2.0;
						hb_correction[hb_don_res[kk]] -= hb_energy[kk]/2.0;
					}
				}
			}

		}	// atoms
	}	// residues
}

void
remove_high_energy_water_molecules(
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn
){
	(scorefxn)(pose);	// needed to generate an energy object 
	TR << "Looking for high energy water molecules" << std::endl;

	Real water_energy_threshold ( option[OptionKeys::hydrate::water_energy_threshold ]() );

	pack::rotamer_set::WaterPackingInfoOP water_info
			( static_cast< pack::rotamer_set::WaterPackingInfo & >
			( pose.data().get( pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) );

	numeric::xyzVector < Real > awayO(0,0,0);
	numeric::xyzVector < Real > awayH1 (0.9572,0,0);
	numeric::xyzVector < Real > awayH2 (-0.2399872,0.92662721,0);
	numeric::xyzVector < Real > far_away(10000,10000,10000);

	utility::vector1 < Real > water_hb_correction (pose.total_residue(), 0 );
	calculate_water_overcoordinated_hb_correction(pose, water_hb_correction);

	for (Size ii = 1; ii<= pose.total_residue(); ++ii){
		if (pose.residue(ii).name() != "TP3") continue;		// Make sure it's water
		if ( (*water_info)[ii].enforced() ){ 		// and it's not enforced
			TR << "ENFORCED WATER MOLECULE " << ii << std::endl;
			continue;
		}

		Real water_energy ( pose.energies().residue_total_energy(ii) + water_hb_correction[ii] );
		if ( water_energy > water_energy_threshold){
			TR << "Removing high energy water: "<< ii << " with energy: " << water_energy << std::endl;
			conformation::ResidueOP temp_residue(pose.residue(ii).clone() );
			temp_residue->set_xyz("O", ii*far_away + awayO);
			temp_residue->set_xyz("H1", ii*far_away + awayH1);
			temp_residue->set_xyz("H2", ii*far_away + awayH2);
			pose.replace_residue(ii,*temp_residue,false);
		}
	}
}

void
enforce_all_waters(
	pose::Pose & pose
){
	pack::rotamer_set::WaterPackingInfoOP water_info
			( static_cast< pack::rotamer_set::WaterPackingInfo & >
			( pose.data().get( pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) );

	for (Size ii = 1; ii <= pose.total_residue(); ++ii){
		if (pose.residue(ii).name() == "TP3") {
			(*water_info)[ ii ].enforced( true );
		}
	}

	pose.data().set( pose::datacache::CacheableDataType::WATER_PACKING_INFO, water_info );
}

// We attempt to pack de novo water molecules that are away from the protein at this point by using just
// one optimized hydrogen bond to their anchor
void
get_ready_for_sew_packing(
	pose::Pose & pose,
	pack::task::PackerTaskOP & task
){
	TR << "Setting task for single edge water molecules" << std::endl;

	conformation::ResidueOP tp3( core::conformation::ResidueFactory::create_residue( 
		pose.residue(1).residue_type_set().name_map("TP3") ) );
	pack::rotamer_set::WaterPackingInfoOP water_info
		( static_cast< pack::rotamer_set::WaterPackingInfo & >
		( pose.data().get( pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) );

	for (Size ii=1; ii<=pose.total_residue(); ++ii){
		if (pose.residue(ii).name() == "TP3" && pose.energies().residue_total_energy(ii) == 0
				&& (*water_info)[ii].anchor_atom() != "NONE"){
			task->nonconst_residue_task( ii ).restrict_to_repacking();
			tp3->set_xyz( "O", pose.residue((*water_info)[ii].anchor_residue()).xyz( (*water_info)[ii].anchor_atom() ) );
			pose.replace_residue(ii,*tp3,false);
			(*water_info)[ii].rotamer_bonds( "SINGLE" );
			TR << "Will repack using single edge rotamers: " << ii << std::endl;
		}
		else {
			task->nonconst_residue_task( ii ).prevent_repacking();
		}
	}
	pose.data().set( pose::datacache::CacheableDataType::WATER_PACKING_INFO, water_info );
}

void
remove_all_anchors(
	pose::Pose & pose
){
	pack::rotamer_set::WaterPackingInfoOP water_info
			( static_cast< pack::rotamer_set::WaterPackingInfo & >
			( pose.data().get( pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) );

	for (Size ii = 1; ii <= pose.total_residue(); ++ii){
		if (pose.residue(ii).name() == "TP3") {
			(*water_info)[ ii ].anchor_atom( "NONE" );
		}
	}

	pose.data().set( pose::datacache::CacheableDataType::WATER_PACKING_INFO, water_info );
}

void 	// yumeng
read_header(
	Size const & total_res,
	std::string const & line,
	kinematics::MoveMap & mm
){
	if (line == "NONE") {
		//do nothing and use default setting false
		TR << "No backbone movement allowed" << std::endl;
	}
	if (line == "ALL") {
		for (Size ii = 1; ii <= total_res; ++ii) {
			mm.set_bb(ii, true);
		}
		TR << "Set backbone movement for all residues" << std::endl;
	}
}

void	// yumeng
read_body(
    Size const & total_res,
    std::string const line,
    kinematics::MoveMap & mm
){
	std::istringstream istr(line);
	std::string method; // action type
	std::string resid; // residue id
	istr >> method;
	istr >> resid;
	core::Size id = atoi(resid.c_str());

	if ( method.compare("MINIBB") == 0 ) {
		if (id > total_res || id == 0) {
			TR << "residue out of range" << id << std::endl;
			exit(0);
		}
		else mm.set_bb(id, true);
		//TR << "set backbone movement for residue: " << id << std::endl;
	}
	//add more options if needed
}

void		// yumeng
set_bb_movemap(
	pose::Pose const & pose,
	std::string const & mini_backbone_file_name,
	kinematics::MoveMap & mm
){
	
	TR << "Reading mini_backbone file " << mini_backbone_file_name << std::endl;

	utility::io::izstream input( mini_backbone_file_name.c_str() ); //get file name

	if ( !input ) {
		TR << "Can't open file " << mini_backbone_file_name << "!" << std::endl;
		exit(0);
	}

	bool not_start = true;
	std::string line;
	while ( getline( input, line ) ) {
		if (line == "START") {
			not_start = false;
			continue;
		}
		if (not_start) read_header( pose.total_residue(), line, mm);
		else read_body( pose.total_residue(), line, mm);
	}
}

void
remove_non_buried_wat(
	pose::Pose & pose
){
	if ( option[ keep_non_buried_waters ]() ) return;

    TR << "Removing non buried water molecules." << std::endl;

	numeric::xyzVector < Real > awayO(0,0,0);
	numeric::xyzVector < Real > awayH1 (0.9572,0,0);
	numeric::xyzVector < Real > awayH2 (-0.2399872,0.92662721,0);
	numeric::xyzVector < Real > far_away(100,100,100);

	utility::vector1 < bool > remove_wat ( pose.total_residue() );
	for ( Size ii = 1; ii<=pose.total_residue() ; ++ii ){
		remove_wat[ii] = false;
		if (pose.residue(ii).name() != "TP3") continue;
		if (is_inside( pose, pose.residue(ii).xyz("O") )) continue;
		remove_wat[ii] = true;
	}
	
	for ( Size ii = pose.total_residue(); ii>=1 ; --ii ){
		if ( remove_wat[ii] ){
			TR << "Removing non buried water " << ii << std::endl;
			conformation::ResidueOP temp_residue(pose.residue(ii).clone() );
			temp_residue->set_xyz("O", ii*far_away + awayO);
			temp_residue->set_xyz("H1", ii*far_away + awayH1);
			temp_residue->set_xyz("H2", ii*far_away + awayH2);
			pose.replace_residue(ii,*temp_residue,false);
		}
	}

}

void
water_specific_scores(
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn
){
	using namespace scoring;

	TR << "Adding water overcoordinated hydrogen bond scores" << std::endl;

	(scorefxn)(pose);
	utility::vector1 < Real > water_hb_correction (pose.total_residue(), 0 );
	calculate_water_overcoordinated_hb_correction(pose, water_hb_correction);

	Real water_overcoor_hb (0);
	TR << "Residue \t wat_overcoor_hb" << std::endl;
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii){
		water_overcoor_hb += water_hb_correction[ii];
		if ( water_hb_correction[ii] != 0)
			TR << ii << "_wat_corr\t" << water_hb_correction[ii] << std::endl;
	}

	std::map< std::string, Real > scores;
	import_pose::atom_tree_diffs::map_of_weighted_scores(pose, scorefxn, scores);

	scores["wat_overcoor_hb"] = water_overcoor_hb;
	scores["total_score"] = scores["total_score"] + water_overcoor_hb;

	// Update scores within the job distrubutor (will output)
	protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	for ( std::map< std::string, Real >::const_iterator curr( scores.begin() ), end( scores.end() );
		curr != end; ++curr) {
		curr_job->add_string_real_pair(curr->first, curr->second);
	}
}

void
show_water_hb_network(
	pose::Pose const & pose
){
	if ( !option[ display_water_hb_network ]() ) return;

	using namespace core::scoring;

	Real hb_sc_weight = pose.energies().weights()[hbond_sc];
	Real hb_bb_sc_weight = pose.energies().weights()[hbond_bb_sc];
	if ( pose.energies().weights()[ hbond_wat ] != 0 ){
		hb_sc_weight = hb_bb_sc_weight = pose.energies().weights()[ hbond_wat ];
	}

	hbonds::HBondSet  hbond_set;
	hbond_set.clear();
	hbonds::fill_hbond_set( pose, false, hbond_set );
	TR << "Water_Hydrogen_Bond_list \t acc acc_atom \t don \t don_atom \t Enegry \t HB_env_weight" << std::endl;
	for ( Size hb=1; hb<= Size(hbond_set.nhbonds()); ++hb ){
		if ( pose.residue( hbond_set.hbond(hb).acc_res() ).name() == "TP3" 
				|| pose.residue( hbond_set.hbond(hb).don_res() ).name() == "TP3"){

			Size weight = hb_sc_weight;
			if ( hbond_set.hbond(hb).acc_atm_is_protein_backbone() ||hbond_set.hbond(hb).don_hatm_is_protein_backbone() ){
				weight = hb_bb_sc_weight;
			}

			TR << "Hydrogen_Bond_Between \t ";
			if ( pose.residue(hbond_set.hbond(hb).acc_res()).name() == "TP3" )
				TR << "wat_" << hbond_set.hbond(hb).acc_res() << " \t ";
			else
				TR << "res_" << hbond_set.hbond(hb).acc_res() << " \t ";
			TR << pose.residue(hbond_set.hbond(hb).acc_res()).atom_name( hbond_set.hbond(hb).acc_atm() );

			if ( pose.residue(hbond_set.hbond(hb).don_res()).name() == "TP3" )
				TR << "\t wat_" << hbond_set.hbond(hb).don_res() << " \t ";
			else
				TR << "\t res_" << hbond_set.hbond(hb).don_res() << " \t ";
			TR << pose.residue(hbond_set.hbond(hb).don_res()).atom_name( hbond_set.hbond(hb).don_hatm() ) << " \t ";
			TR << hbond_set.hbond(hb).energy()*weight*hbond_set.hbond(hb).weight() << " \t ";
			TR << hbond_set.hbond(hb).weight() << std::endl;

		}
	}	// hb
}




/*


void
read_water_near_res_file(
	std::string const & wnnr_file_name,
	utility::vector1< bool > & wnnr_V
){

	TR << "Reading water_near_res file " << wnnr_file_name << std::endl;

	utility::io::izstream input( wnnr_file_name.c_str() );
	std::string line;
	if ( !input ) {
		TR << "Can't open file " << wnnr_file_name << "!" << std::endl;
		exit(0);
	}
	
	while ( getline( input, line ) ) {
		if ( line.substr(0,3).compare("RES") == 0 ){	  // This comparison means it IS "ENF"
			TR << line << std::endl;
			std::stringstream ss(line);
			std::istream_iterator<std::string> begin(ss);
			std::istream_iterator<std::string> end;
			std::vector<std::string> vstrings(begin, end);
			Size rsd = atoi(vstrings[1].c_str());
			wnnr_V[rsd] = true;
		}
	}

}




void
hydrate_cavities( 
	pose::Pose & pose, 
	utility::vector1< bool > & enforced_V
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys::hydrate;

	TR << "Hydrating cavities." << std::endl;

	scoring::hbonds::HBondSet  hbond_set;
	hbond_set.clear();
	Real const & hbond_threshold( option[ OptionKeys::hydrate::hbond_threshold ]() );
	scoring::hbonds::fill_hbond_set( pose, false, hbond_set );
	Size non_water_residues = pose.total_residue();
	for (Size ii = 1; ii <= pose.total_residue(); ++ii){
		if (pose.residue(ii).name() == "TP3"){
			non_water_residues = ii -1;
			break;
		}
	}
	Size old_pose_total_res = pose.total_residue();
	conformation::ResidueOP tp3( core::conformation::ResidueFactory::create_residue( 
			pose.residue(1).residue_type_set().name_map("TP3") ) );
	pack::rotamer_set::WaterPackingInfoOP water_info
			( static_cast< pack::rotamer_set::WaterPackingInfo & >
			( pose.data().get( pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) );

	for (Size ii = 1; ii <= non_water_residues; ++ii){
		conformation::Residue const & rsd( pose.residue(ii) );
		for ( Size jj=1; jj<= Size(rsd.natoms()); ++jj ) {  // Go over all polar atoms
			bool hydrate_atm = false;

			if (rsd.atom_type(jj).is_acceptor() ){
				hydrate_atm = true;         // Check if it has the all hbs satisfied
				Size bonds = 0;
				for ( Size hb=1; hb<= Size(hbond_set.nhbonds()); ++hb ){
					if ( option[ OptionKeys::hydrate::ignore_hb_depth_on_hydration ]() ){
						if (hbond_set.hbond(hb).energy() < hbond_threshold
								&& ii == hbond_set.hbond(hb).acc_res() && jj == hbond_set.hbond(hb).acc_atm() ){
							++bonds;
						}
					}
					else if (hbond_set.hbond(hb).energy()*hbond_set.hbond(hb).weight() < hbond_threshold
							&& ii == hbond_set.hbond(hb).acc_res() && jj == hbond_set.hbond(hb).acc_atm() ){
						++bonds;
					}
//std::cout << "weight = " << hbond_set.hbond(hb).weight() << "\t";
//std::cout << ".energy = " << hbond_set.hbond(hb).energy() << "\t";
//std::cout << "hbond_threshold = " << hbond_threshold << "\t";
//std::cout << "residue,atom = " << hbond_set.hbond(hb).acc_res() << "\t" << hbond_set.hbond(hb).acc_atm() << "\n";
				}
				if (bonds >= 2) hydrate_atm = false;
				if ( option[ OptionKeys::hydrate::dont_hydrate_hb_engaged_bb_O ]() && rsd.atom_name(jj) == " O  " ) {
					if ( bonds >=1 ) hydrate_atm = false;
				}
			}
			if ( rsd.atom_type(jj).is_polar_hydrogen()){
				hydrate_atm = true;
				for ( Size hb=1; hb<= Size(hbond_set.nhbonds()); ++hb ){
					if ( option[ OptionKeys::hydrate::ignore_hb_depth_on_hydration ]() ){
						if (hbond_set.hbond(hb).energy() < hbond_threshold
								&& ii == hbond_set.hbond(hb).don_res() && jj == hbond_set.hbond(hb).don_hatm() ){
							hydrate_atm = false;
						}
					}
					else if (hbond_set.hbond(hb).energy()*hbond_set.hbond(hb).weight() < hbond_threshold 
							&& ii == hbond_set.hbond(hb).don_res() && jj == hbond_set.hbond(hb).don_hatm() ){
						hydrate_atm = false;
					}
				}
			}
			// Check that the specific atom in not hydrated already
			for (Size kk = non_water_residues +1; kk <= old_pose_total_res; kk++){
				if ((*water_info)[kk].anchor_residue() == ii && (*water_info)[kk].anchor_atom() == rsd.atom_name(jj))
					hydrate_atm = false;
			}


			if (hydrate_atm == true){   // Now we know the atom is polar and has unsatisfied hb

				if ( !atom_is_hydratable(pose,ii,jj) ) continue;
				tp3->set_xyz( "O", rsd.xyz(jj) ); // for neighbor calculation
				if ( option[ OptionKeys::hydrate::water_jumps_to_root]() ){
					pose.append_residue_by_jump( *tp3, 1 );
				} else {
					pose.append_residue_by_jump( *tp3, ii );
				}
				Size const pos1( pose.total_residue() );

				Size anchor_atom = jj;
				if (rsd.name() == "HIS_D")     // jag Hack to be able to use HIS_D
					if (anchor_atom == 15) anchor_atom = 7;
				if (rsd.name() == "HIS")
					if (anchor_atom == 17) anchor_atom = 10;

				(*water_info)[ pos1 ].anchor_residue( ii );
				(*water_info)[ pos1 ].anchor_atom( rsd.atom_name(anchor_atom) );
				(*water_info)[ pos1 ].aa( rsd.aa() );
				(*water_info)[ pos1 ].nstep( 1 );
				enforced_V.push_back(false);
				TR << "Water: " << pos1 << " appended to residue: " << ii << " " << rsd.name();
				TR << " -> " << rsd.atom_name(anchor_atom) << std::endl;
				
			}
		}
	}
	pose.data().set( core::pose::datacache::CacheableDataType::WATER_PACKING_INFO, water_info );

}



Real
atom_ASA(
    pose::Pose & pose,
    Size & residue,
    std::string const & atom
)
{
	Real resolution = 24;
	Real step_cos_phi = 2/resolution;
	Real step_theta = 3.14159265/resolution;

	Vector ref_x (1,0,0);
	Vector ref_y (0,1,0);
	Vector ref_z (0,0,1);
	Vector atom_xyz (pose.residue(residue).xyz(atom));

	Real r_atom (1.8);
	if (atom[1] == 'H' )    r_atom = 1.2;
	if (atom[1] == 'C')     r_atom = 1.7;
	if (atom[1] == 'O')     r_atom = 1.52;
	if (atom[1] == 'N')     r_atom = 1.55;
	if (atom[1] == 'F')     r_atom = 1.47;
	r_atom += 1.52;     // To account for the solvent radious

	Size blocked_sections = 0;
	for ( Real cos_phi = -1; cos_phi < 0.99; cos_phi += step_cos_phi){
		Real phi = acos(cos_phi);
		for ( Real theta = 0; theta < 6.28; theta += step_theta){

			// Position of the hypothetical solvent (Oxygen)
			Vector hyp_O (atom_xyz  + r_atom*( sin(phi)*(ref_x*cos(theta) + ref_y*sin(theta)) + ref_z*cos_phi ));

			bool blocked (false);
			for (Size ii = 1; ii<=pose.total_residue(); ++ii){
				for (Size jj = 1; jj<=pose.residue(ii).natoms(); ++jj){

					if (ii == residue && jj == pose.residue(residue).atom_index(atom) ) continue;

					Real r_neighbor (1.8);
					if (pose.residue(ii).atom_name(jj)[1] == 'H')       r_neighbor = 1.2;
					if (pose.residue(ii).atom_name(jj)[1] == 'C')       r_neighbor = 1.7;
					if (pose.residue(ii).atom_name(jj)[1] == 'O')       r_neighbor = 1.52;
					if (pose.residue(ii).atom_name(jj)[1] == 'N')       r_neighbor = 1.55;
					if (pose.residue(ii).atom_name(jj)[1] == 'F')       r_neighbor = 1.47;
					r_neighbor += 1.52;     // To account for the solvent radious

					Vector neighbor_xyz (pose.residue(ii).xyz(jj));
					Real distance = neighbor_xyz.distance(hyp_O);

					if ( distance <= r_neighbor ){
						++blocked_sections;
						blocked = true;
						break;
					}
				}
				if (blocked) break;
			}
		}
	}

	// returns the Accessible Surface Area in A
	return (1 - blocked_sections/resolution/resolution/2)*4*3.1415926535*(r_atom - 1.52)*(r_atom - 1.52);
}


Real
atom_ASA(
	pose::Pose & pose,
	Size & residue,
	Size & atom
)
{
	return atom_ASA(pose, residue, pose.residue(residue).atom_name(atom));
}

void
remove_exposed_wat(
	pose::Pose & pose
){
	TR << "Removing exposed water molecules " << std::endl;

	numeric::xyzVector < Real > awayO(0,0,0);
	numeric::xyzVector < Real > awayH1 (0.9572,0,0);
	numeric::xyzVector < Real > awayH2 (-0.2399872,0.92662721,0);
	numeric::xyzVector < Real > far_away(100,100,100);

	utility::vector1< bool > exposed_wat( pose.total_residue() );
	for (Size ii = 1; ii <= pose.total_residue(); ++ii){
		if ( pose.residue(ii).name() == "TP3" ){
			exposed_wat[ii] = false;
			Real total_ASA (0);
			for ( Size jj=1; jj<= pose.residue(ii).natoms(); ++jj){
				total_ASA += atom_ASA(pose,ii,jj);
			}
			if ( total_ASA >= 0.1 ){
				exposed_wat[ii] = true;
			}
		}
	}
	for (Size ii = 1; ii <= pose.total_residue(); ++ii){
		if ( pose.residue(ii).name() == "TP3" && exposed_wat[ii] == true && pose.residue(ii).xyz(1).x() < 1000){
			TR << "Removing exposed water "<< ii << std::endl;
			conformation::ResidueOP temp_residue(pose.residue(ii).clone() );
			temp_residue->set_xyz("O", ii*far_away + awayO);
			temp_residue->set_xyz("H1", ii*far_away + awayH1);
			temp_residue->set_xyz("H2", ii*far_away + awayH2);
			pose.replace_residue(ii,*temp_residue,false);
		}
	}

}

void
slide_jumps_to_closest_residue(
	pose::Pose & pose
){
	TR << "Sliding jumps to closest residue instead of default root. " << std::endl;

	kinematics::FoldTree f = pose.fold_tree();
	Size protein_residues (0);
	for (Size ii=1; ii<=pose.total_residue(); ++ii){
		if ( pose.residue(ii).atom_name(1) != " N  " &&  pose.residue(ii).atom_name(2) != " CA "){	// Not protein

			Real jump_res(1);
			std::string jump_atom (" CA ");
			Real closest_res_d(100);
			Vector wat_xyz (pose.residue(ii).xyz(1));	// this can actually be water, ion or ligand

			for (Size jj=1; jj<= protein_residues; jj++){
				for (Size atom=1; atom<= pose.residue(jj).nheavyatoms(); ++atom){
					Vector atom_xyz(pose.residue(jj).xyz(atom));
					Real distance = atom_xyz.distance(wat_xyz);
					if ( distance < closest_res_d){
						closest_res_d = distance;
						jump_res = jj;
						jump_atom = pose.residue(jj).atom_name(atom);
					}
				}
			}
			Size jump = f.get_jump_that_builds_residue(ii);
			f.slide_jump(jump, jump_res, ii);
			f.set_jump_atoms(jump, jump_res, jump_atom, ii, pose.residue(ii).atom_name(1));
		}
		else {
			protein_residues = ii;
		}
	}
	pose.fold_tree(f);
}


void
set_hyb_movemap(
	pose::Pose & pose,
	kinematics::MoveMap & movemap
){
	// All movemap elements start being true, so we just need to restrict it.
	// In our case we just restrict water molecules away from protein
	for (Size ii = 1; ii <= pose.total_residue(); ++ii){
		if ( pose.residue(ii).name() == "TP3" && pose.energies().residue_total_energy(ii) == 0 )
			movemap.set_chi(ii, false);
	}
	
	// All jumps will be flexible in the near water region
	kinematics::FoldTree f = pose.fold_tree();
	for (Size jj = 1; jj <= f.num_jump() ; ++jj){
		movemap.set_jump(jj, true);
	}

}

void
set_two_regions_movemaps(
	pose::Pose & pose,
	core::kinematics::MoveMap & not_wat_movemap,
	core::kinematics::MoveMap & wat_movemap
){
	// All movemap elements start being true, so we just need to restrict it.

	for (Size ii = 1; ii <= pose.total_residue(); ++ii){

		// It's away water and should not be packed at this stage at all
		// Will be packed in second round with single edge water rotamers ( sew )
		if ( pose.residue(ii).name() == "TP3" && pose.energies().residue_total_energy(ii) == 0 ){
			not_wat_movemap.set_chi(ii, false);			
			wat_movemap.set_chi(ii, false);			
			continue;
		}

		// near water region
		if (residue_near_water(pose, ii) || pose.residue(ii).name() == "TP3") {
			not_wat_movemap.set_chi(ii, false);
			not_wat_movemap.set_bb(ii, false);
		}
		// far from water region
		else {
			wat_movemap.set_chi(ii, false);
			wat_movemap.set_bb(ii, false);
		}
	}

	// All jumps will be flexible in the near water region
	kinematics::FoldTree f = pose.fold_tree();
	for (Size jj = 1; jj <= f.num_jump() ; ++jj){
		wat_movemap.set_jump(jj, true);
	}
}


void
set_task_with_de_novo_water_using_resfile( 
	pose::Pose & pose, 
	std::string resfile, 
	pack::task::PackerTaskOP & task
){
	TR << "Setting task with de novo water, using resfile " << std::endl;

	pack::task::PackerTaskOP temp_task( pack::task::TaskFactory::create_packer_task( pose ));
	pack::task::parse_resfile(pose, *temp_task, resfile);				// Must be here to have water info in task
	for (Size ii = 1; ii <= pose.total_residue(); ++ii) {
		if ( pose.residue(ii).name() == "TP3" ) {
			task->nonconst_residue_task(ii).restrict_to_repacking();	// To avoid problems with design
		}
		else if ( !temp_task->being_packed(ii)) {
			task->nonconst_residue_task(ii).prevent_repacking();
		}
	}

}






void
calculate_water_entropy_correction(
	pose::Pose & pose,
	utility::vector1< Real > & S_correction
){

	Real const hbond_threshold( option[ OptionKeys::hydrate::hbond_threshold]() );

	scoring::hbonds::HBondSet  hbond_set;
	hbond_set.clear();
	scoring::hbonds::fill_hbond_set( pose, false, hbond_set );

	for (Size ii= 1; ii <= pose.total_residue(); ++ii)  // We make sure we start from 0 for all residues
		S_correction[ii] = 0;

	for (Size ww = 1; ww <= pose.total_residue(); ++ww){
		if ( pose.residue(ww).name() != "TP3") continue; 					// Go over all waters
		if ( pose.energies().residue_total_energy(ww) == 0 ) continue;		// That are near the protein

		Real hb_count (0);

			// Go over all hydrogen bonds and find the ones this water (ww) is involved in
			// Surely this can be done more efficiently...
		for ( Size hb=1; hb<= Size(hbond_set.nhbonds()); ++hb ){
			if ( ( hbond_set.hbond(hb).acc_res() == ww || hbond_set.hbond(hb).don_res() == ww ) // water or interest
					&& hbond_set.hbond(hb).energy()*hbond_set.hbond(hb).weight() < hbond_threshold  // relevant hb
			){
				if( pose.residue(hbond_set.hbond(hb).acc_res()).name() != "TP3"
						|| pose.residue(hbond_set.hbond(hb).don_res()).name() != "TP3" ){
				++hb_count;			// Only water-protein hbs contribute to the water entropy
				}
			}
		}
		// The water entropy is assumend proportional to the number of hbs it makes with the protein
		// See Steven W. Rick's paper (J. Phys. Chem. B, 2010, 114, 11552-11560)
		// The maximum entropic cost is 4 RU ~ 2 kcal/mol 
		// See Jack D. Dunitz paper (Science, vol 264, 670)
		if (hb_count > 4)	hb_count = 4;
		S_correction[ww] = hb_count*option[ OptionKeys::hydrate::water_entropy_weight]() ;

	}	// all waters

}


// All energies are not weighted
void
make_hybrid_energy_map(
	pose::Pose & pose,
	scoring::EnergyMap & EM,
  scoring::ScoreFunction & scorefxn
){

	// transfer energy on hyb terms

	scoring::hbonds::HBondSet  hbond_set;
	hbond_set.clear();
	scoring::hbonds::fill_hbond_set( pose, false, hbond_set );

	if (option[ OptionKeys::hydrate::separate_wat_wat_vs_wat_res_hb_optE]()){

		for ( Size hb=1; hb<= Size(hbond_set.nhbonds()); ++hb ){
			if ( pose.residue(hbond_set.hbond(hb).acc_res()).name() == "TP3" &&
						pose.residue(hbond_set.hbond(hb).don_res()).name() == "TP3" ){
				EM[ scoring::wat_wat_hb ] += hbond_set.hbond(hb).energy();
				EM[ scoring::hbond_sc ] -= hbond_set.hbond(hb).energy();
			}
			else if ( pose.residue(hbond_set.hbond(hb).acc_res()).name() == "TP3" ||
						pose.residue(hbond_set.hbond(hb).don_res()).name() == "TP3" ){
				EM[ scoring::wat_res_hb ] += hbond_set.hbond(hb).energy();
				if ( hbond_set.hbond(hb).acc_atm_is_protein_backbone() ||
							hbond_set.hbond(hb).don_hatm_is_protein_backbone() ){
					EM[ scoring::hbond_bb_sc ] -= hbond_set.hbond(hb).energy();
				} else { 
					EM[ scoring::hbond_sc ] -= hbond_set.hbond(hb).energy();
				}
			}
		}

	} else {

		for ( Size hb=1; hb<= Size(hbond_set.nhbonds()); ++hb ){
		  if ( pose.residue(hbond_set.hbond(hb).acc_res()).name() == "TP3" ||
						pose.residue(hbond_set.hbond(hb).don_res()).name() == "TP3" ){
				EM[ scoring::hbond_wat ] += hbond_set.hbond(hb).energy();
				if ( hbond_set.hbond(hb).acc_atm_is_protein_backbone() ||
							hbond_set.hbond(hb).don_hatm_is_protein_backbone() ){
					EM[ scoring::hbond_bb_sc ] -= hbond_set.hbond(hb).energy();
				} else { 
					EM[ scoring::hbond_sc ] -= hbond_set.hbond(hb).energy();
				}
			}
		}

	}
	
	//	account for wat_3b_hb
	utility::vector1 < Real > water_3b_hb_correction (pose.total_residue(), 0 );
	calculate_water_overcoordinated_hb_correction(pose, water_3b_hb_correction);
	for (Size ii=1; ii<=pose.total_residue(); ++ii){	// needs to be divided by weight to not be weighted
		EM [ scoring::wat_3b_hb ] += water_3b_hb_correction[ii]/pose.energies().weights()[ scoring::hbond_sc ];
	}

	//	account for wat_entropy		// jag, deprecated 
//	utility::vector1 < Real > water_S_correction (pose.total_residue(), 0 );
//	calculate_water_entropy_correction(pose, water_S_correction);
//	for (Size ii=1; ii<=pose.total_residue(); ++ii){
//		EM [ scoring::wat_entropy ] += water_S_correction[ii];
//	}
	
	// Separate fa_sol wat-res 
	// Generate a pose without any water and calculate wat_fa_sol
	pose::Pose no_wat_pose = pose;
	Real wat_fa_sol = 0;
	numeric::xyzVector < core::Real > awayO(0,0,0);
	numeric::xyzVector < core::Real > awayH1 (0.9572,0,0);
	numeric::xyzVector < core::Real > awayH2 (-0.2399872,0.92662721,0);
	numeric::xyzVector < core::Real > far_away(1000,1000,1000);
	for (Size ii = 1; ii <= no_wat_pose.total_residue(); ++ii){
		if ( no_wat_pose.residue(ii).name() == "TP3"){
			conformation::ResidueOP temp_residue(no_wat_pose.residue(ii).clone() );
			temp_residue->set_xyz("O", ii*far_away + awayO);
			temp_residue->set_xyz("H1", ii*far_away + awayH1);
			temp_residue->set_xyz("H2", ii*far_away + awayH2);
			no_wat_pose.replace_residue(ii,*temp_residue,false);

			wat_fa_sol += pose.energies().residue_total_energies(ii)[ scoring::fa_sol ];
		}
	}
	(scorefxn)(no_wat_pose);
	Real residues_wat_contribution = pose.energies().total_energies()[ scoring::fa_sol ]
																		- wat_fa_sol
																		- no_wat_pose.energies().total_energies()[ scoring::fa_sol ];

	// Since it's hybrid, fa sol should not have any water contribution
	EM [ scoring::fa_sol ] -= wat_fa_sol;

	// Move the residues' fa_sol due to water to a different score, so it can be optimized separately
	EM [ scoring::fa_sol ] -= residues_wat_contribution;
//	EM [ scoring::res_due_wat_fa_sol ] += residues_wat_contribution;

}	// All energies are unweighted

// All energies are not weighted
void
account_for_entropy(
	pose::Pose & pose,
	scoring::EnergyMap & EM,
  scoring::ScoreFunction & scorefxn
){

	//	account for wat_3b_hb
	utility::vector1 < Real > water_3b_hb_correction (pose.total_residue(), 0 );
	calculate_water_overcoordinated_hb_correction(pose, water_3b_hb_correction);
	for (Size ii=1; ii<=pose.total_residue(); ++ii){	// needs to be divided by weight to not be weighted
		EM [ scoring::wat_3b_hb ] += water_3b_hb_correction[ii]/pose.energies().weights()[ scoring::hbond_sc ];
	}

	//	account for wat_entropy
//	utility::vector1 < Real > water_S_correction (pose.total_residue(), 0 );
//	calculate_water_entropy_correction(pose, water_S_correction);
//	for (Size ii=1; ii<=pose.total_residue(); ++ii){
//		EM [ scoring::wat_entropy ] += water_S_correction[ii];
//	}

}


*/
}	//	close hydrate
}	//	close hydrate
