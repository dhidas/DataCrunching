#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <sstream>
#include <time.h>

#include "conf_gen_ag.h"
#include "fingerprint.h"
#include "master_score.h"
#include "simplex.h"
#include "trace.h"

class Bump_Filter;
class Master_Score;

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++
AG_Conformer_Search::AG_Conformer_Search()
{
    ie_vdwA = NULL;
    ie_vdwB = NULL;  
    write_fragment_libraries = false;
}


// +++++++++++++++++++++++++++++++++++++++++
AG_Conformer_Search::~AG_Conformer_Search()
{
    delete[]ie_vdwA;
    delete[]ie_vdwB; 
}


// +++++++++++++++++++++++++++++++++++++++++
void
AG_Conformer_Search::input_parameters(Parameter_Reader & parm)
{
    string         tmp;
    write_fragment_libraries = (parm.query_param("write_fragment_libraries", "no", "yes no") == "yes");
    if (write_fragment_libraries) {
        fragment_library_prefix = (parm.query_param("fragment_library_prefix", "fraglib"));
        fragment_library_freq_cutoff = atoi(parm.query_param("fragment_library_freq_cutoff", "1").c_str());
        fragment_library_sort_method = parm.query_param("fragment_library_sort_method", "freq", "freq | fingerprint");
        fragment_library_trans_origin = (parm.query_param("fragment_library_trans_origin", "no", "yes no") == "yes");
        user_specified_anchor = false;
        limit_max_anchors = false;
    }
    else {
    

    cout << "\nAnchor & Grow Parameters" << endl;
    cout <<
    "------------------------------------------------------------------------------------------" 
    << endl;

    user_specified_anchor = false;
    limit_max_anchors = false;
    anchor_size = 40;
    max_anchor_num = 10;

    user_specified_anchor = (parm.query_param("user_specified_anchor", "no", "yes no") == "yes");
    if (user_specified_anchor){
       //cout << "Specify atom name,number seperated by a ','" << endl;
       atom_in_anchor = parm.query_param("atom_in_anchor", "C1,1");
    } else {
        limit_max_anchors = (parm.query_param("limit_max_anchors", "no", "yes no") == "yes");
        if (limit_max_anchors){
        //    anchor_size = 3; // min_anchor_size is set to 3
            max_anchor_num = atoi(parm.query_param("max_anchor_num", "1").c_str());
        }

         anchor_size = atoi(parm.query_param("min_anchor_size", "5").c_str());
         if (anchor_size <= 0) {
             cout << "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                  << endl;
             exit(0);
         }
    } 


    cluster = (parm.query_param("pruning_use_clustering", "yes", "yes no") == "yes");
    if(cluster){
        num_anchor_poses =
            atoi(parm.query_param("pruning_max_orients", "1000").c_str());
        if (num_anchor_poses <= 0) {
            cout << "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        pruning_clustering_cutoff = atoi(parm.query_param(
                "pruning_clustering_cutoff", "100").c_str());
        if (pruning_clustering_cutoff <= 0) {
            cout << "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        
        anchor_score_cutoff = 1000.0;
        num_growth_poses = INT_MAX;
        // DTM-11-10-08
        //growth_cutoff = false;
        //growth_cutoff = true;
        growth_score_cutoff_begin = atof(parm.query_param("pruning_conformer_score_cutoff", "100.0").c_str());
        growth_score_cutoff = growth_score_cutoff_begin;
        growth_score_scaling_factor = atof(parm.query_param("pruning_conformer_score_scaling_factor", "1.0").c_str());

    } else { 
        num_anchor_poses =
            atoi(parm.query_param("pruning_max_orients", "1000").c_str());
        if (num_anchor_poses <= 0) {
            cout << "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        anchor_score_cutoff = atof(parm.query_param("pruning_orient_score_cutoff", "100.0").c_str());
        num_growth_poses = atoi(parm.query_param("pruning_max_conformers", "75").c_str());
        //growth_cutoff = true;
        growth_score_cutoff_begin = atof(parm.query_param("pruning_conformer_score_cutoff", "100.0").c_str());
        growth_score_cutoff = growth_score_cutoff_begin;
        growth_score_scaling_factor = atof(parm.query_param("pruning_conformer_score_scaling_factor", "1.0").c_str());
    }
    use_clash_penalty = (parm.query_param("use_clash_overlap", "no", "yes no") == "yes");

    if (use_clash_penalty) {
        clash_penalty = atof(parm.query_param("clash_overlap", "0.5").c_str());
        if (clash_penalty <= 0.0) {
            cout << "clash_overlap should be a positive number" << endl;
            exit(0);
        }
    }

    print_growth_tree = (parm.query_param("write_growth_tree", "no", "yes no") == "yes");
    if (print_growth_tree) {
        cout << "Warning: Writing the growth tree increases memory usage and can generate lots of large files." << endl
             << "Concatenating and compressing growth tree branches is recommended" << endl;
    }

/*
    // For de novo fragment library generation

    //write_fragment_libraries = (parm.query_param("write_fragment_libraries", "no", "yes no") == "yes");
    if (write_fragment_libraries) {
        fragment_library_prefix = (parm.query_param("fragment_library_prefix", "fraglib"));
        fragment_library_freq_cutoff = atoi(parm.query_param("fragment_library_freq_cutoff", "1").c_str());
        fragment_library_sort_method = parm.query_param("fragment_library_sort_method", "freq", "freq | fingerprint");
        fragment_library_trans_origin = (parm.query_param("fragment_library_trans_origin", "no", "yes no") == "yes");
    }
*/
    verbose = 0 != parm.verbosity_level();   // -v is for extra scoring info
}
}

// +++++++++++++++++++++++++++++++++++++++++
void
AG_Conformer_Search::initialize()
{
    // Separate input section from initializing section.
    cout << endl;
    // I.e., this is the first output for the initializing section.
    cout << "Initializing Conformer Generator Routines...\n";
}


// +++++++++++++++++++++++++++++++++++++++++
// Trent E Balius 2009-02-11
// this function gets the internal energy parms from Master_Conformer_Search.
void
AG_Conformer_Search::initialize_internal_energy_parms(bool uie, int rep_exp, int att_exp, float diel, float iec)
{
     use_internal_energy = uie;
     ie_rep_exp = rep_exp;
     ie_att_exp = att_exp;
     ie_diel    = diel;
     internal_energy_cutoff = iec;
}

// +++++++++++++++++++++++++++++++++++++++++
// Dwight 2016-06-16
// this function simply initializes the internal energy to False when a scoring 
//function is specified but internal energy is set to NO and it gets the internal
//energy value from Master_Conformer_Search
void
AG_Conformer_Search::initialize_internal_energy_null(bool uie)
{
 use_internal_energy = uie;

}

// +++++++++++++++++++++++++++++++++++++++++
void
AG_Conformer_Search::prepare_molecule(DOCKMol & mol)
{
        // DTM - 11-12-08
        int        i;

        copy_molecule(orig, mol);

    bond_list.clear();
    atom_seg_ids.clear();
    bond_seg_ids.clear();
    anchors.clear();
    layers.clear();
    pruned_confs.clear();
    orig_segments.clear();
    layer_segments.clear();

    current_anchor = 0;

    identify_rigid_segments(mol);
    id_anchor_segments();

    // DTM - 11-12-08 - copy segment assignments to new atom_segment_ids array in the DOCKMol object
    for(i=0;i<atom_seg_ids.size();i++) {
            orig.atom_segment_ids[i] = atom_seg_ids[i];
            mol.atom_segment_ids[i] = atom_seg_ids[i];
    }

    // CDS - 04-14 - copy segment assignments to new bond_segment_ids array in the DOCKMol object
    for(i=0;i<bond_seg_ids.size();i++) {
            orig.bond_segment_ids[i] = bond_seg_ids[i];
            mol.bond_segment_ids[i] = bond_seg_ids[i];
    }

    // Write fragments for de novo libraries
    if (write_fragment_libraries){
        count_fragments(mol);
    }

    anchor_positions.clear();
    anchor_positions.reserve(1000); // can store more than 1000 anchors dynamically
}


// +++++++++++++++++++++++++++++++++++++++++
void
AG_Conformer_Search::identify_rigid_segments(DOCKMol & mol)
{
    int             i,
                    j,
                    max_central,
                    central;
    int             a1,
                    a2,
                    a3,
                    a4;
    int             t1,
                    t2,
                    nbr;
    vector < int   >nbrs;

    BREADTH_SEARCH  bfs;

    // Identify rigid segments & which segment each atom belongs to
    for (i = 0; i < orig.num_atoms; i++)
        atom_seg_ids.push_back(-1);

    for (i = 0; i < orig.num_bonds; i++)
        bond_seg_ids.push_back(-1);

    extend_segments(0, 0, mol);

    for (i = 0; i < orig.num_atoms; i++) {
        if (atom_seg_ids[i] < 0) {
            cout << "ERROR: molecule seems to consist of multiple unconnected fragments.  Program will terminate." << endl;
            exit(1);
        }
    }

    for (i = 0; i < orig.num_bonds; i++) {

        t1 = orig.bonds_origin_atom[i];
        t2 = orig.bonds_target_atom[i];

        if (atom_seg_ids[t1] == atom_seg_ids[t2]) {
            bond_seg_ids[i] = atom_seg_ids[t1];
            orig_segments[bond_seg_ids[i]].bonds.push_back(i);
        }

    }

    // ID the inter-segment rot-bonds
    for (i = 0; i < bond_list.size(); i++) {

        a2 = orig.bonds_origin_atom[bond_list[i].bond_num];
        a3 = orig.bonds_target_atom[bond_list[i].bond_num];

        max_central = -1;
        nbrs = mol.get_atom_neighbors(a2);

        for (j = 0; j < nbrs.size(); j++) {
            nbr = nbrs[j];
            if (nbr != a3) {
                central = bfs.get_search_radius(mol, nbr, a2);
                if (central > max_central) {
                    a1 = nbr;
                    max_central = central;
                }
            }
        }

        max_central = -1;
        nbrs = mol.get_atom_neighbors(a3);

        for (j = 0; j < nbrs.size(); j++) {
            nbr = nbrs[j];
            if (nbr != a2) {
                central = bfs.get_search_radius(mol, nbr, a3);
                if (central > max_central) {
                    a4 = nbr;
                    max_central = central;
                }
            }
        }

        bond_list[i].atom1 = a1;
        bond_list[i].atom2 = a2;
        bond_list[i].atom3 = a3;
        bond_list[i].atom4 = a4;
        bond_list[i].seg1 = atom_seg_ids[bond_list[i].atom2];   // -1
        bond_list[i].seg2 = atom_seg_ids[bond_list[i].atom3];   // -1
        bond_list[i].initial_angle =
            orig.get_torsion(bond_list[i].atom1, bond_list[i].atom2,
                             bond_list[i].atom3, bond_list[i].atom4);
    }

}


// +++++++++++++++++++++++++++++++++++++++++
// is this function growing out a single segment ?
void
AG_Conformer_Search::extend_segments(int atom_num, int segment, DOCKMol & mol)
{
    int             nbr;
    int             bond;
    int             i;
    vector < int   >nbrs;

    SEGMENT         tmp_segment;
    ROT_BOND        tmp_bond;

    if (atom_num == 0)
        orig_segments.push_back(tmp_segment);

    atom_seg_ids[atom_num] = segment;
    orig_segments[segment].atoms.push_back(atom_num);

    if (orig.amber_at_heavy_flag[atom_num])
        orig_segments[segment].num_hvy_atoms++;

    nbrs = orig.get_atom_neighbors(atom_num);

    for (i = 0; i < nbrs.size(); i++) {

        nbr = nbrs[i];
        bond = orig.get_bond(atom_num, nbr);

        if (atom_seg_ids[nbr] == -1) {

            if (mol.bond_is_rotor(bond)) {

                bond_list.push_back(tmp_bond);
                bond_list[bond_list.size() - 1].bond_num = bond;

                orig_segments.push_back(tmp_segment);

                orig_segments[segment].neighbors.push_back(orig_segments.
                                                           size() - 1);
                orig_segments[segment].neighbor_bonds.push_back(bond_list.
                                                                size() - 1);
                orig_segments[segment].neighbor_atoms.push_back(nbr);

                orig_segments[orig_segments.size() -
                              1].neighbors.push_back(segment);
                orig_segments[orig_segments.size() -
                              1].neighbor_bonds.push_back(bond_list.size() - 1);
                orig_segments[orig_segments.size() -
                              1].neighbor_atoms.push_back(atom_num);

                extend_segments(nbr, orig_segments.size() - 1, mol);

            } else {

                extend_segments(nbr, segment, mol);

            }
        }
    }
}


// +++++++++++++++++++++++++++++++++++++++++
// TEB ADD 2010-01-23
// Idenify if atom in string is in the
// segment
// +++++++++++++++++++++++++++++++++++++++++
bool
AG_Conformer_Search::atom_in_anchor_segments(SEGMENT seg)
{
 bool flag = false;
 int i;
 stringstream s;
 string atomname, current_atomname,temp;
 int atomid, current_atomid ;
 s << atom_in_anchor;
 getline(s,atomname,',');
 s >> atomid;
// cout <<"lets see" <<atomname << " " << atomid << endl;

 cout << atom_in_anchor << " :: " << atomname << "," << atomid <<endl;

 for (i = 0; i < seg.atoms.size(); i++){
     // the atom num is not stored in dockmol.
     current_atomid      = seg.atoms[i];
//     cout << current_atomid << "," ;
     current_atomname = orig.atom_names[current_atomid];
//     cout << current_atomname << endl;

     if (atomname == current_atomname && (current_atomid+1) == atomid){
//         cout << "condition met!" << endl;
         flag = true;
         break;
     }
 }

 if (flag){
   for (i = 0; i < seg.atoms.size(); i++){
     int current_atomid      = seg.atoms[i]; // the atom num is not stored.
     string current_atomname = orig.atom_names[current_atomid];
//     cout << current_atomname << "," << current_atomid << endl;
   }
 }

return flag;
}


// +++++++++++++++++++++++++++++++++++++++++
void
AG_Conformer_Search::id_anchor_segments()
{
    int             i;
    INTPair         tmp;

//    cout << "I AM HERE in AG_Conformer_Search::id_anchor_segments()" << endl;

    for (i = 0; i < orig_segments.size(); i++) {
        tmp.first = orig_segments[i].num_hvy_atoms + orig_segments[i].neighbors.size(); 
        // problem with anchor defs this way
        //  sudipto: who put this comment here and why?
        tmp.second = i;
        anchors.push_back(tmp);
    }

    sort(anchors.begin(), anchors.end());
    reverse(anchors.begin(), anchors.end());

    // TEB ADD 2010-01-23
    // pick only the user spesified anchor for placement in the array.
    if (user_specified_anchor){
//        cout << "user_specified_anchor is true" << endl;
        vector <INTPair> temp_anchors;
        temp_anchors = anchors;
        anchors.clear();
        for (i = 0; i < temp_anchors.size(); i++) {
            if (atom_in_anchor_segments(orig_segments[temp_anchors[i].second])){
//               cout << " user anchor defined." <<endl;
               //anchors.push_back(temp_anchors[i]);
               tmp.first = temp_anchors[i].first;
               tmp.second = temp_anchors[i].second;
               anchors.push_back(tmp);
            }
        }
        if (anchors.size() == 0){
            cout << "atom_in_anchor (" << atom_in_anchor <<") not found no segment" <<endl;
            exit(0);
        }
        if (anchors.size() > 1){
            cout << "atom_in_anchor (" << atom_in_anchor <<") is shared amoug multiple anchors" <<endl;
            cout << "pick new atom and try again" <<endl;
        }
    }

    // TEB ADD 2010-01-25
    // pick only the top N (spesified by the user) anchors
//  if (limit_max_anchors)
    if (limit_max_anchors){
        vector <INTPair> temp_anchors;
        temp_anchors = anchors;
        anchors.clear();
        int num_anchors_val;

        cout << max_anchor_num  << " " << temp_anchors.size() << endl;

        if (max_anchor_num > temp_anchors.size() ){
           cout << "max_anchor_nunum of segments." << endl;
           cout << "Only "<< temp_anchors.size() << " anchors are used." << endl;
           num_anchors_val = temp_anchors.size();
        } else {
           cout << max_anchor_num << " anchors are used." << endl;
           num_anchors_val = max_anchor_num;
        }
        for (i = 0; i < num_anchors_val ; i++) {
           cout << i << endl;
           tmp.first = temp_anchors[i].first;
           tmp.second = temp_anchors[i].second;
           anchors.push_back(tmp);
        }

    }

}


// +++++++++++++++++++++++++++++++++++++++++
bool
AG_Conformer_Search::next_anchor(DOCKMol & mol)
{
    DOCKMol         tmp_mol,
                    blank_mol;
    CONFORMER       tmp_conf;
    vector < CONFORMER > conf_vec;
    int             i;

    // trent & sudipto Jan 15, 09
    // this makes ure we no longer have the orients of the previous anchor
    anchor_positions.clear();

    // If there are no structs in anchor_confs, generate a new anchor and
    // expand it
    if (anchor_confs.size() == 0) {

        // if there are no more anchors, or you have reached the end of the
        // list
        if ((anchors.size() == 0) || (current_anchor == anchors.size())
            || ((anchors[current_anchor].first < anchor_size)
                && (current_anchor > 0))){
            if (verbose) cout << "No more anchor fragments to be docked." << endl;
            return false;
        }

        copy_molecule(tmp_mol, orig);

        // clear bond directionality during minization
        bond_tors_vectors.clear();
        bond_tors_vectors.resize(orig.num_bonds, -1);
        //bond_tors_vectors.resize(mol.num_bonds, -1);

        // reset the assigned atoms list
        assigned_atoms.clear();
        assigned_atoms.resize(orig.num_atoms, false);
        //assigned_atoms.resize(mol.num_atoms, false);

        // clear and resize the layer_segment list to be the same size as the
        // orig_segments list
        layer_segments.clear();
        layer_segments.resize(orig_segments.size());

        layers.clear();
        extend_layers(anchors[current_anchor].second,
                      anchors[current_anchor].second, 0);

        // cout << "@@@\t" << orig_segments.size() << "\t" << layers.size() <<
        // "\t" << layer_segments.size() << endl;

        for (i = 0; i < layers[0].segments.size(); i++) {
            activate_layer_segment(tmp_mol, 0, i);
        }

        current_anchor++;

        // copy the anchor to the anchor_conf list
        anchor_confs.clear();
        anchor_confs.resize(anchor_confs.size() + 1);
        copy_molecule(anchor_confs[anchor_confs.size() - 1], tmp_mol);

    }

/**
        int i,j,k;
        for(i=0;i<layers.size();i++) {
                for(j=0;j<layers[i].segments.size();j++) {
                        for(k=0;k<layer_segments[layers[i].segments[j]].atoms.size();k++) {
                                cout << i << "\t" << j << "\t" << layer_segments[layers[i].segments[j]].atoms[k]+1 << endl;
                        }
                }
        }
**/

    // copy last mol from anchor_conf to mol
    copy_molecule(mol, anchor_confs[anchor_confs.size() - 1]);
    anchor_confs.pop_back();

    if (anchor_confs.size() > 0)
        last_conformer = false;
    else
        last_conformer = true;

    return true;


}


// +++++++++++++++++++++++++++++++++++++++++
void
AG_Conformer_Search::extend_layers(int previous_segment, int current_segment,
                                   int current_layer)
{
    int             i;
    INTPair         tmp_pair;
    LAYER           tmp_layer;
    LAYER_SEGMENT   tmp_segment;

    if (current_layer >= layers.size())
        layers.push_back(tmp_layer);

    layers[current_layer].segments.push_back(current_segment);
    layers[current_layer].num_segments = layers[current_layer].segments.size();

    // loop over the atoms in each orig_segment and add them to the
    // layer_segment
    // if the atom has not already been assigned to a previous segment
    for (i = 0; i < orig_segments[current_segment].atoms.size(); i++) {
        if (!assigned_atoms[orig_segments[current_segment].atoms[i]]) {
            layer_segments[current_segment].atoms.
                push_back(orig_segments[current_segment].atoms[i]);
            assigned_atoms[orig_segments[current_segment].atoms[i]] = true;
        }
    }

    // copy the bonds from orig_segment to layer_segment
    layer_segments[current_segment].bonds =
        orig_segments[current_segment].bonds;

    for (i = 0; i < orig_segments[current_segment].neighbors.size(); i++) {

        if (orig_segments[current_segment].neighbors[i] != previous_segment) {
            // if this bond leads to the next layer, then add the leading atom
            // and bond to the current segment

            layer_segments[current_segment].bonds.
                push_back(bond_list
                          [orig_segments[current_segment].neighbor_bonds[i]].
                          bond_num);

            if (!assigned_atoms
                [orig.
                 bonds_origin_atom[bond_list
                                   [orig_segments[current_segment].
                                    neighbor_bonds[i]].bond_num]]) {
                layer_segments[current_segment].atoms.push_back(orig.
                                                                bonds_origin_atom
                                                                [bond_list
                                                                 [orig_segments
                                                                  [current_segment].
                                                                  neighbor_bonds
                                                                  [i]].
                                                                 bond_num]);
                assigned_atoms[orig.
                               bonds_origin_atom[bond_list
                                                 [orig_segments
                                                  [current_segment].
                                                  neighbor_bonds[i]].
                                                 bond_num]] = true;

                // set bond direction vector
                bond_tors_vectors[bond_list
                                  [orig_segments[current_segment].
                                   neighbor_bonds[i]].bond_num] =
                    orig.
                    bonds_target_atom[bond_list
                                      [orig_segments[current_segment].
                                       neighbor_bonds[i]].bond_num];
            }

            if (!assigned_atoms
                [orig.
                 bonds_target_atom[bond_list
                                   [orig_segments[current_segment].
                                    neighbor_bonds[i]].bond_num]]) {
                layer_segments[current_segment].atoms.push_back(orig.
                                                                bonds_target_atom
                                                                [bond_list
                                                                 [orig_segments
                                                                  [current_segment].
                                                                  neighbor_bonds
                                                                  [i]].
                                                                 bond_num]);
                assigned_atoms[orig.
                               bonds_target_atom[bond_list
                                                 [orig_segments
                                                  [current_segment].
                                                  neighbor_bonds[i]].
                                                 bond_num]] = true;

                // set bond direction vector
                bond_tors_vectors[bond_list
                                  [orig_segments[current_segment].
                                   neighbor_bonds[i]].bond_num] =
                    orig.
                    bonds_origin_atom[bond_list
                                      [orig_segments[current_segment].
                                       neighbor_bonds[i]].bond_num];
            }

            extend_layers(current_segment,
                          orig_segments[current_segment].neighbors[i],
                          current_layer + 1);

        } else {
            // if this bond comes from the previous layer, then add it as the
            // rot_bond in the segment

            layer_segments[current_segment].rot_bond =
                orig_segments[current_segment].neighbor_bonds[i];
            layer_segments[current_segment].origin_segment = previous_segment;
        }
    }

}


// +++++++++++++++++++++++++++++++++++++++++
bool
AG_Conformer_Search::submit_anchor_orientation(DOCKMol & mol, bool more_orients)
{
    Trace trace("AG_Conformer_Search::submit_anchor_orientation");
    SCOREMol        tmp_mol;
    float           mol_score;
    int             insert_point;
    int             i;

    mol_score = mol.current_score;

    // Loop over existing list of anchors
    insert_point = -1;
    for (i = 0; i < anchor_positions.size(); i++) {

        // finds the first instance of a lower score than the mol, and selects
        // that as the insert point
        if (insert_point == -1) {
            if (anchor_positions[i].first > mol_score)
                insert_point = i;
        }

    }

    //insert mol into list, and pop last member off
    if (insert_point == -1) {
        // if insert_point = -1 (mol to be appended to end of list)
        // only append mol if list is less than number of confs/layer
        if (anchor_positions.size() < num_anchor_poses) {
                anchor_positions.push_back(tmp_mol);
                copy_molecule(anchor_positions[anchor_positions.size() - 1].
                              second, mol);
                anchor_positions[anchor_positions.size() - 1].first = mol_score;
        }

    } else {
        // if molecule is to be inserted within the list

        // if list is smaller than limit, append a blank mol on the end
        if (anchor_positions.size() < num_anchor_poses) {
            anchor_positions.push_back(tmp_mol);
        }
        // start at end of list, and shift all mols one spot to the end of
        // the list, until the insert point is reached
        for (i = anchor_positions.size() - 1; i > insert_point; i--) {
            copy_molecule(anchor_positions[i].second,
                              anchor_positions[i - 1].second);
            anchor_positions[i].first = anchor_positions[i - 1].first;
        }

        // copy mol into insert point
        copy_molecule(anchor_positions[insert_point].second, mol);
        anchor_positions[insert_point].first = mol_score;
   
    }

    // enter into loop if no more orients...
    if (! more_orients){
    // trent e balius 2011-05-30
        sort(anchor_positions.begin(), anchor_positions.end(), less_than_pair);
        if (verbose) cout << "# of anchor positions: " << anchor_positions.size() << endl;
        return true;
    } else {
        //if (verbose) cout << "there are no more orients and not last_conformer" << endl;
        return false;
    }

}


// +++++++++++++++++++++++++++++++++++++++++
float
AG_Conformer_Search::calc_layer_rmsd(CONFORMER & a, CONFORMER & b)
{
    int             i,
                    j,
                    k;
    int             atom;
    float           rmsd;
    int             heavy_total;

    rmsd = 0.0;
    heavy_total = 0;

    for (i = 0; i < layers.size(); i++) {
        for (j = 0; j < layers[i].segments.size(); j++) {
            for (k = 0; k < layer_segments[layers[i].segments[j]].atoms.size();
                 k++) {

                atom = layer_segments[layers[i].segments[j]].atoms[k];

                if (a.structure.atom_active_flags[atom]
                    && a.structure.amber_at_heavy_flag[atom]) {

                    rmsd +=
                        (pow((a.structure.x[atom] - b.structure.x[atom]), 2) +
                         pow((a.structure.y[atom] - b.structure.y[atom]),
                             2) + pow((a.structure.z[atom] -
                                       b.structure.z[atom]), 2)) * (i + 1);

                    heavy_total += i + 1;
                }
            }
        }
    }

    if (heavy_total > 0)
        rmsd = sqrt(rmsd / (float) heavy_total);
    else
        rmsd = 0.0;

    return rmsd;
}


// +++++++++++++++++++++++++++++++++++++++++
float
AG_Conformer_Search::calc_active_rmsd(CONFORMER & ref, CONFORMER & conf)
{
// this function calculates rmsd between the active atoms of the ref structure and the same atoms of conf.
// the ref structure is usually the anchor
// only heavy atom rmsd is reported
// the rmsd of the active atoms in the reference is reported
    int    i;
    float  rmsd = 0.0;
    int    atom_num_total = 0;

    for (i = 0; i < ref.structure.num_atoms; i++) {
          if (ref.structure.atom_active_flags[i] && ref.structure.amber_at_heavy_flag[i]){
                    rmsd +=
                        (pow((ref.structure.x[i] - conf.structure.x[i]), 2) +
                         pow((ref.structure.y[i] - conf.structure.y[i]), 2) + 
                         pow((ref.structure.z[i] - conf.structure.z[i]), 2));

                    atom_num_total += 1;
          }
    }

    if (atom_num_total > 0)
        rmsd = sqrt(rmsd / (float) atom_num_total);
    else
        rmsd = 0.0;

    return rmsd;
}

// +++++++++++++++++++++++++++++++++++++++++
void
AG_Conformer_Search::grow_periphery(Master_Score & score,
                                    Simplex_Minimizer & simplex, Bump_Filter & bump)
{ Trace trace("AG_Conformer_Search::grow_periphery()");
    //cout << "AG_Conformer_Search::grow_periphery" << endl;
    int             i,
                    j,
                    k,
                    l;
    CONFORMER       tmp_conf;
    SCOREMol        tmp_scoremol;

    vector < SCOREMol > anchor_positions_b4min;
    vector < CONFORMER > seeds, exp_seeds;
    vector < CONFORMER > b4min_seeds;
    vector < CONFORMER > all_gen_seeds; // store tree with each conformer pointing to its parent
    vector < CONFORMER > all_gen_b4min_seeds; // all unminimized conformers
    float           rmsd;

    // initialize
    bool valid_orient = false;
    bool ie_prune = false; // prune conformers using internal_energy_cutoff
    int num_layers = layers.size();
    clock_t start = clock();

    //reserve space in vectors
    anchor_positions_b4min.reserve(anchor_positions.size());
    all_gen_seeds.reserve(2000);
    all_gen_b4min_seeds.reserve(2000);

    // check if a scoring function is being used.
    if (score.use_primary_score) {

        //initialize internal energy
        trace.note("AG_conformer::grow_periphery:reinitialize energy ");
        score.primary_score->use_internal_energy = use_internal_energy;
        if (use_internal_energy) {
           score.primary_score->ie_att_exp = ie_att_exp;
           score.primary_score->ie_rep_exp = ie_rep_exp;
           score.primary_score->ie_diel = ie_diel;
           ie_prune=true;

           // need to use DOCKMol with radii and segments assigned
           // it does not matter which atoms are labeled active
           score.primary_score->initialize_internal_energy(anchor_positions[0].second);
        }
    }// else do nothing to internal energy (because it is not used when scoring function is not
     // used

    // /////////////////////////// ANCHORS //////////////////////////////

    if (verbose) cout << "-----------------------------------" << endl 
         << "VERBOSE GROWTH STATS : ANCHOR #" <<  current_anchor << endl << endl;

    // add anchors to seed list

    // moved from main dock loop for flexible docking only
    // sudipto & trent Jan-16-08
    for (i = 0; i < anchor_positions.size(); i++){
       anchor_positions_b4min.push_back(anchor_positions[i]);  // save unmin anchor in anchor_positions_b4min
       simplex.minimize_rigid_anchor(anchor_positions[i].second, score); //minimize anchor
    }

    for (i = 0; i < anchor_positions.size(); i++) {
        // compute the score for the molecule
        if (score.use_primary_score) {
            valid_orient = score.compute_primary_score(anchor_positions[i].second);

            // code to check if there was an error in the scoring function: sudipto 
            // the scoring function stores an error description in dockmol.current_data 
            if (!valid_orient && verbose) cout << "Error scoring Anchor #" << current_anchor 
                 << " Orient #" << i << ": " << anchor_positions[i].second.current_data << endl;

        } else valid_orient = true;

        if (valid_orient) {
            exp_seeds.push_back(tmp_conf);
            exp_seeds[exp_seeds.size() - 1].layer_num = 1;
            exp_seeds[exp_seeds.size() - 1].score = anchor_positions[i].second.current_score;
            exp_seeds[exp_seeds.size() - 1].used = false;
            copy_molecule(exp_seeds[exp_seeds.size() - 1].structure, anchor_positions[i].second);

            b4min_seeds.push_back(tmp_conf);
            b4min_seeds[b4min_seeds.size() - 1].layer_num = 1;

            //b4_min_seeds[b4_min_seeds.size() - 1].score = anchor_positions_b4min[i].second.current_score;
            // we want the b4min score to match the after min score so sorting and pruning the anchors remain in step

            b4min_seeds[b4min_seeds.size() - 1].score = anchor_positions[i].second.current_score;
            b4min_seeds[b4min_seeds.size() - 1].used = false;
            copy_molecule(b4min_seeds[b4min_seeds.size() - 1].structure, anchor_positions_b4min[i].second);
        }
    }

    //assert(exp_seeds.size() == b4min_seeds.size());

    // Pruning section
    seeds.clear(); seeds.reserve(500);
    sort(exp_seeds.begin(), exp_seeds.end(), conformer_less_than);
    sort(b4min_seeds.begin(), b4min_seeds.end(), conformer_less_than);
  
    // make sure we are not missing any unminimized structures
    // remove once we are sure growth tree is bug-free: sudipto
    //assert(exp_seeds.size() == b4min_seeds.size());

    // filter confs by scores
    for (j = 0; j < exp_seeds.size(); j++) {
        if (exp_seeds[j].score > anchor_score_cutoff)  // +1000 for clustering, +25 otherwise
            exp_seeds[j].used = true;
    }

    // remove confs that fail the rank/rmsd test
    if(cluster){
        float RMSD_CUTOFF = pruning_clustering_cutoff;
        for (j = 0; j < exp_seeds.size(); j++) {
               if (!exp_seeds[j].used) {
                    for (k = j + 1; k < exp_seeds.size(); k++) {
                       if (!exp_seeds[k].used) {

                            rmsd = calc_layer_rmsd(exp_seeds[j], exp_seeds[k]);
                            rmsd = MAX(rmsd, 0.001);

                            if ((float) k / rmsd > RMSD_CUTOFF) {
                                exp_seeds[k].used = true;
                            }
                       }        
                    }
              }
        }
    }

    // End pruning section
    int anchor_count = 0; // for assigning anchor number 
    conf_anchors.clear(); // remove old anchor orients from previous calls
    conf_anchors.reserve(anchor_positions.size()); //anchor_positions.size() is an upper limit for # of clustered anchors

    // copy pruned confs to seed list
    for (j = 0; j < exp_seeds.size(); j++) {

            // This section is required for the growth tree.
            exp_seeds[j].conformer_num = count_conf_num;
            exp_seeds[j].parent_num   = -1; // when parent is -1, we are at the top of the tree
            b4min_seeds[j].conformer_num = count_conf_num;
            b4min_seeds[j].parent_num   = -1; 
            
            // print anchors to file
            // print_conformer(exp_seeds[j]); 

            if (!exp_seeds[j].used) {
           
                exp_seeds[j].anchor_num   = anchor_count; //anchor position in array conf_anchors 
                b4min_seeds[j].anchor_num = anchor_count; 
                conf_anchors.push_back(b4min_seeds[j]); //unminimized anchors
           
                // Calculate header string for anchor. 
                // The conf_before_min parameter in this function is not used for anchors
                // This section is required for the growth tree.
                // When the growth tree is turned off, computing rmsd for all the conformers
                // is not needed.
                if (print_growth_tree) conf_header(exp_seeds[j], "Minimized Anchor", score);
                if (print_growth_tree) conf_header(b4min_seeds[j], "Unminimized Anchor", score);
                if (print_growth_tree) b4min_seeds[j].structure.simplex_text = "";
           
                seeds.push_back(exp_seeds[j]);
                if (print_growth_tree) all_gen_seeds.push_back(exp_seeds[j]); 
                if (print_growth_tree) all_gen_b4min_seeds.push_back(b4min_seeds[j]); 
                anchor_count++;
            } 
            count_conf_num++;   // this counter also includes pruned anchors
    }
    // /////////////////////////// END ANCHORS //////////////////////////

    // print anchor pruning staticstics
    // why are no anchors getting pruned??
    if (verbose) cout << seeds.size() << "/" <<  anchor_positions.size() 
               << " anchor orients retained (max " << num_anchor_poses << ")" 
               << " t=" << ((double) (clock() - start)) / CLOCKS_PER_SEC << "s" << endl;

    // /////////////////////////// GROWTH ///////////////////////////////

    // counters for printing pruning stats
    int confs_pruned_bad_score = 0;
    int confs_pruned_clash_overlap = 0;
    int confs_pruned_outside_grid = 0;
    int confs_pruned_bump_filter = 0;
    int confs_pruned_clustered = 0;

    if (verbose) {
      for (int layer_count = 1; layer_count < num_layers; layer_count++)
        cout << "Lyr " << layer_count << "-" << layers[layer_count].segments.size() << " Segs|";
      cout << endl; 
    } 

    // loop over layers
    for (i = 1; i < num_layers; i++) {  // 

        // loop over rot bonds in segment
        for (l = 0; l < layers[i].segments.size(); l++) {

            exp_seeds.clear(); exp_seeds.reserve(1000);
            b4min_seeds.clear(); b4min_seeds.reserve(1000);

            // print details of rotatable bond sampled
            if (verbose && seeds.size() > 0) {
              int bond = layer_segments[layers[seeds[0].layer_num].segments[l]].rot_bond;  
              cout << "Lyr:" << i << " Seg:" << l << " Bond:" << bond_list[bond].bond_num 
                   << " : Sampling " <<  seeds[0].structure.amber_bt_torsion_total[bond_list[bond].bond_num]
                   << " dihedrals " 
                   << seeds[0].structure.atom_names[bond_list[bond].atom1] << "("
                   << seeds[0].structure.atom_types[bond_list[bond].atom1] << ")  "
                   << seeds[0].structure.atom_names[bond_list[bond].atom2] << "("
                   << seeds[0].structure.atom_types[bond_list[bond].atom2] << ")  "
                   << seeds[0].structure.atom_names[bond_list[bond].atom3] << "("
                   << seeds[0].structure.atom_types[bond_list[bond].atom3] << ")  "
                   << seeds[0].structure.atom_names[bond_list[bond].atom4] << "("
                   << seeds[0].structure.atom_types[bond_list[bond].atom4] << ")"
                   //<< " t=" << ((double) (clock() - start)) / CLOCKS_PER_SEC << "s"
                   << endl;
            } 

            // loop over seed conformers
            for (j = 0; j < seeds.size(); j++) {
                
                // Required for supporting multiple grid docking
                // Duplicates of each exp_seed is made for each grid with the corresponding grid_num
                if (score.ir_ensemble) { 
                        segment_torsion_drive(seeds[j], l, exp_seeds, score.c_mg_nrg.numgrids);
                        //cout << "grow_periphery: " << score.c_mg_nrg.numgrids << " grids" << endl;
                }
                else 
                        segment_torsion_drive(seeds[j], l, exp_seeds, 1);
                // num_grids = 1 when multiple grid docking is turned off
     
                // score/minimize the new confs
                for (k = 0; ((k < exp_seeds.size())&&(j==seeds.size()-1)); k++) {

                    // sudipto & trent - 11-14-08 - save current conf before minimizing
                    CONFORMER conf_before_min = exp_seeds[k];
                    b4min_seeds.push_back(conf_before_min);
                             
                    if (segment_clash_check(exp_seeds[k].structure, i, l)) {    // i+1
                    // Note: clash-checking not compatible with ir_resecore
                        
                        if (bump.check_growth_bumps(exp_seeds[k].structure)) {
                        // Note: bump-checking not compatible with ir_resecore

                             //print location in growth loop and minimize conformer 
                           //POSSIBLE PLACE TO  REINITIALIZE THE ANCHORS  
                             simplex.minimize_flexible_growth(exp_seeds[k].structure, score, bond_tors_vectors);
                             
                            // compute the score for the molecule
                            // now internal energy is computed within compute_primary_score
                            if (score.use_primary_score){
                                valid_orient = score.compute_primary_score(exp_seeds[k].structure);
                                valid_orient = score.compute_primary_score(b4min_seeds[k].structure);
                            } else
                                valid_orient = true;
 
                            if (valid_orient) {

                                 // compute internal energy for pruning
                                 //score.primary_score->compute_ligand_internal_energy(exp_seeds[k].structure);
                                 //score.primary_score->compute_ligand_internal_energy(b4min_seeds[k].structure);
                                         
                                 // internal_energy here is just the repulsive component of the vdw energy
                                 ostringstream text;
                                 text <<  i << ":" << l << ":" << j << ":" << k; 
                                 
                                 // This section is required for the growth tree.
                                 // When the growth tree is turned off, computing rmsd for all the conformers
                                 // is not needed.
                                 if (print_growth_tree) conf_header(exp_seeds[k],   text.str(), score);
                                 if (print_growth_tree) conf_header(b4min_seeds[k], text.str(), score);
                                 if (print_growth_tree) b4min_seeds[k].structure.simplex_text = "";

                                 // exp_seeds.score is used for sorting by energy. Only the top #max_conformers will
                                 // be reported, so internal energy should not be use for sorting. 
                                 // Sorting should be done only on the grid score
                                 exp_seeds[k].score = exp_seeds[k].structure.current_score + 
                                                      exp_seeds[k].structure.internal_energy;

                                 // set the score for the b4min structure to be the same as the minimized structure
                                 // this is a hack to make sure b4min_seeds gets sorted in the same order as exp_seeds
                                 b4min_seeds[k].score = exp_seeds[k].score;
 
                                 // this prints every conformer ever generated in the mol.conf_no.mol2 format
                                 // superseded by branch_* and pruned_* mol2 files
                                 /*if (verbose) cout << "Writing conformer #" << exp_seeds[k].conformer_num 
                                      << " k=" << k << " count_conf_num=" << count_conf_num << endl;
                                 print_conformer(exp_seeds[k]); 
                                 */
                                         growth_score_cutoff = growth_score_cutoff_begin  / ( pow (growth_score_scaling_factor, i) );

                                         if (exp_seeds[k].structure.current_score <= growth_score_cutoff) {
                                             if (ie_prune) {
                                                 if (exp_seeds[k].structure.internal_energy > internal_energy_cutoff) {
                                                     exp_seeds[k].used = true;
                                                     confs_pruned_bad_score++;
                                                 }
                                             }
                                         } else {
                                               exp_seeds[k].used = true;
                                               confs_pruned_bad_score++;
                                         }

                            } else {
                                exp_seeds[k].used = true;  // scoring failed
                                confs_pruned_outside_grid++;
                                //if (verbose) cout << "Error scoring Conformer #" << exp_seeds[k].conformer_num << endl;
                            } 
                        } else { 
                             exp_seeds[k].used = true;  // failed bump filter 
                             confs_pruned_bump_filter++;
                             //if (verbose) cout << " Conformer #" << exp_seeds[k].conformer_num << " failed bump filter" << endl;
                        }
                    } else { 
                        exp_seeds[k].used = true;  // clash_overlap check 
                        confs_pruned_clash_overlap++;
                        //if (verbose) cout << " Conformer #" << exp_seeds[k].conformer_num << " failed clash check" << endl;
                    }

                    // count_conf_num should be incremented inside the torsion function only
                    // no new conformers are produced outside of segment_torsion_drive()
                } // Loop over new conformers (k)

            }  // Loop over seed conformers (j)
 
            // Pruning section
            seeds.clear(); seeds.reserve(500);
            sort(exp_seeds.begin(), exp_seeds.end(), conformer_less_than);
            sort(b4min_seeds.begin(), b4min_seeds.end(), conformer_less_than);

            // DTM 11-11-08 - code to prune bad energies has been moved into the inner loop where confs are generated
          
            // remove confs that fail the rank/rmsd test
            if(cluster) {
                float RMSD_CUTOFF = pruning_clustering_cutoff;
                for (j = 0; j < exp_seeds.size(); j++) {
                    if (!exp_seeds[j].used) {
                        for (k = j + 1; k < exp_seeds.size(); k++) {
                            if (!exp_seeds[k].used) {

                                rmsd = calc_layer_rmsd(exp_seeds[j], exp_seeds[k]);
                                rmsd = MAX(rmsd, 0.001);

                                if ((float) k / rmsd > RMSD_CUTOFF) {
                                    exp_seeds[k].used = true;
                                    confs_pruned_clustered++;
                                }

                            }
                        }
                    }
                }
             }
    
            // End pruning section

            // DTM 11-11-08 - modified this to cut off the list of new seeds at the limit num_growth_poses
            // copy pruned confs to seed list
            for (j = 0; ((j < exp_seeds.size())&&(seeds.size()<num_growth_poses)); j++) {

                if (!exp_seeds[j].used) {
                    seeds.push_back(exp_seeds[j]);
                    if (print_growth_tree) all_gen_seeds.push_back(exp_seeds[j]);
                    if (print_growth_tree) all_gen_b4min_seeds.push_back(b4min_seeds[j]);
                } else {
                    // print pruned branch
                    // note that the current pruned conformer is not present in the all_gen_seeds list,
                    // but that does not matter because the function prints the grown ligand (i.e. current pruned conf)
                    // first and then looks for its parents conformers which are present in this list
                }
            }  // End loop copy pruned confs to seed list (j)

            // print the sizes of vectors to show how many were pruned
            if (verbose) {  
                cout << "Lyr:" << i << " Seg:" << l << " "; 
                cout << seeds.size() << "/" << exp_seeds.size() << " retained, Pruning: "; 
                if (confs_pruned_outside_grid>0) cout << confs_pruned_outside_grid << "-outside grid ";
                if (confs_pruned_bad_score>0) cout << confs_pruned_bad_score << "-score ";
                if (cluster && confs_pruned_clustered>0)           cout << confs_pruned_clustered << "-clustered ";
                if (use_clash_penalty && confs_pruned_clash_overlap>0) cout << confs_pruned_clash_overlap << "-clash overlap ";
                if (bump.bump_filter && confs_pruned_bump_filter>0)  cout << confs_pruned_bump_filter << "-bump filter "; 
                if (print_growth_tree) cout << " (" << all_gen_seeds.size() << " in growth tree)";
                cout << " t=" << ((double) (clock() - start)) / CLOCKS_PER_SEC << "s";
                cout << endl;
            }

            //reset pruning stats counters
            confs_pruned_bad_score = 0;
            confs_pruned_clash_overlap = 0;
            confs_pruned_outside_grid = 0;
            confs_pruned_bump_filter = 0;
            confs_pruned_clustered = 0;

        }  // End loop over segments (l)

    }   // End loop over layers (i)

    if (print_growth_tree)
      for (i = 0; i < seeds.size(); i++)
        print_branch(all_gen_seeds, all_gen_b4min_seeds, seeds[i],score);

    // ///////////////////////// END GROWTH /////////////////////////////

    // copy mols to pruned confs
    pruned_confs.clear(); pruned_confs.reserve(seeds.size());
    for (i = 0; i < seeds.size(); i++) {
        pruned_confs.push_back(tmp_scoremol);
        pruned_confs[i].first = seeds[i].score;
        copy_molecule(pruned_confs[i].second, seeds[i].structure);
    }
    count_conf_num++;// to make sure that the conf_num does not repeat.

    //clear all the vectors: do we really need to do this? (sudipto)
    anchor_positions_b4min.clear();
    seeds.clear();
    exp_seeds.clear();
    b4min_seeds.clear();
    all_gen_seeds.clear(); 
    all_gen_b4min_seeds.clear(); 
}

// +++++++++++++++++++++++++++++++++++++++++
void
AG_Conformer_Search::conf_header(CONFORMER & conf, string text, Master_Score & score)
{
    // sudipto & trent - 11-14-08 a print header for score,etc inside growth loop
    // text = "Anchor" for the anchor and lyr:seg:cnf:tor for growth conformers
    ostringstream conformer_text;
    // For consistency use the formatting for score outputting.
    const string DELIMITER = Base_Score::DELIMITER;

    // Print conformer and parent number afmin 
    conformer_text << DELIMITER << "Conformer #:\t\t" << conf.conformer_num << endl 
                   << DELIMITER << "Parent #:\t\t" << conf.parent_num << endl 
                   << DELIMITER << "Anchor #:\t" << current_anchor << endl; 

    // calculate grid score for the current conformer
    // also, calculate internal score for the current conformer
    score.compute_primary_score(conf.structure);

    // String dockmol.current_data contains the output score summary
    conformer_text << conf.structure.current_data; 

    if (conf.parent_num==-1)  // Anchors 
      conformer_text << DELIMITER << text  << endl;
    else                      // Growth Conformers
      conformer_text << DELIMITER << "lyr:seg:cnf:tor:\t" 
                     << fixed << setprecision(3) << text << endl;

    conformer_text << DELIMITER << "RMSD2Anc:\t" // calculates RMSD to unminimized anchor
                     << calc_active_rmsd(conf_anchors[conf.anchor_num], conf) << endl;
    
    // create a new conformer object for orig as calc_active_rmsd needs a CONFORMER
    // DOCKMol orig = xtal input ligand
    CONFORMER    orig_conf;
    copy_molecule(orig_conf.structure, orig);

    // this RMSD considers only the active atoms in the partially gorwn conformer
    // which is why conf must precede orig when calling calc_active_rmsd
    conformer_text << DELIMITER << "RMSD2Orig:\t" // calculates RMSD to xtal input ligand
                     << calc_active_rmsd(conf, orig_conf) << endl;

    // for optimization we should really stop creating a new conformer for orig
    // every time this function is called. We can either store orig_conf in the class
    // or overload calc_active_rmsd() to work with a DOCKMol as well: sudipto
    
    // save the header text in the current conformer
    conf.header = conformer_text.str();
}


// +++++++++++++++++++++++++++++++++++++++++
void
AG_Conformer_Search::segment_torsion_drive(CONFORMER & conf, int current_bond, vector < CONFORMER > &return_list, int num_rec)
{
    int             num_torsions,
                    bond;
    CONFORMER       new_conf,
                    tmp_conf;
    float           new_angle;
    int             seg_id;

    // bond =
    // layer_segments[layers[conf.layer_num-1].segments[current_bond]].rot_bond;
    bond = layer_segments[layers[conf.layer_num].segments[current_bond]].rot_bond;      // -1

    num_torsions = conf.structure.amber_bt_torsion_total[bond_list[bond].bond_num];

    // sudipto & trent Dec 09, 2008
    // print num of torsions sampled
    //if (verbose) cout << num_torsions << " torsions sampled for bond_num=" << bond_list[bond].bond_num 
    //    << " Conf_#=" << conf.conformer_num << " Score=" << conf.score << endl;

    for (int i = 0; i < num_torsions; i++) {

        copy_molecule(new_conf.structure, conf.structure);
        new_conf.layer_num = conf.layer_num;
        new_conf.score = conf.score;
        new_conf.used = conf.used;
        new_conf.anchor_num = conf.anchor_num;// trent balius 2008-11-17
        new_conf.parent_num = conf.conformer_num; //trent balius 2008-12-03 
        new_conf.conformer_num = count_conf_num; //trent balius 2008-12-03
//        cout << new_conf.conformer_num << endl;

        new_angle =
            new_conf.structure.amber_bt_torsions[bond_list[bond].bond_num][i];
        new_angle = (PI / 180) * new_angle;

        if (bond_list[bond].seg1 == layer_segments[layers[new_conf.layer_num].segments[current_bond]].origin_segment) { // -1
            new_conf.structure.set_torsion(bond_list[bond].atom1,
                                           bond_list[bond].atom2,
                                           bond_list[bond].atom3,
                                           bond_list[bond].atom4, new_angle);
            seg_id = bond_list[bond].seg2;
        }

        if (bond_list[bond].seg2 == layer_segments[layers[new_conf.layer_num].segments[current_bond]].origin_segment) { // -1
            new_conf.structure.set_torsion(bond_list[bond].atom4,
                                           bond_list[bond].atom3,
                                           bond_list[bond].atom2,
                                           bond_list[bond].atom1, new_angle);
            seg_id = bond_list[bond].seg1;
        }

        if (layer_segments[layers[new_conf.layer_num].segments[current_bond]].origin_segment != bond_list[bond].seg1)   // -1
            if (layer_segments[layers[new_conf.layer_num].segments[current_bond]].origin_segment != bond_list[bond].seg2)       // -1
                cout << "Layer growth error!" << endl;

        activate_layer_segment(new_conf.structure, new_conf.layer_num, current_bond);   // -1

        // Write_Mol2(new_conf.structure,
        // cout);////////////////////////////////////////////////////////////////

        if (current_bond == layers[new_conf.layer_num].segments.size() - 1) {   // -1
            new_conf.layer_num++;
        }

      // Generate a copy for each grid
      // Generate copies of each segment_torsion_drive seed for each individual receptor grid_num
      for (int mg_num=0; mg_num < num_rec; mg_num++) {
              return_list.push_back(tmp_conf);        // ///////////////////////////////////////////////////////////////
              return_list[return_list.size() - 1].layer_num = new_conf.layer_num;
              return_list[return_list.size() - 1].score = new_conf.score;
              return_list[return_list.size() - 1].used = new_conf.used;
              return_list[return_list.size() - 1].anchor_num = new_conf.anchor_num;// trent balius 2008-11-17
              return_list[return_list.size() - 1].conformer_num = new_conf.conformer_num;// trent balius 2008-12-03
              return_list[return_list.size() - 1].parent_num = new_conf.parent_num;// trent balius 2008-12-03
              copy_molecule(return_list[return_list.size() - 1].structure, new_conf.structure);
              return_list[return_list.size() - 1].structure.grid_num = mg_num;
              count_conf_num++; //trent balius 2008-12-03
      }
    }

}


// +++++++++++++++++++++++++++++++++++++++++
// trent balius 2008-12-03
void
AG_Conformer_Search::print_conformer(CONFORMER & conf)
{
       ostringstream file;
       file << "mol." << conf.conformer_num << "_anchor_" << current_anchor << ".mol2";
       fstream fout;
       fout.open (file.str().c_str(), fstream::out|fstream::app);
                            
       fout << conf.header;
       Write_Mol2(conf.structure, fout);
       return;
}


// +++++++++++++++++++++++++++++++++++++++++
// trent balius 2008-12-03
void
AG_Conformer_Search::print_branch(vector <CONFORMER> &list, vector <CONFORMER> &b4min, const CONFORMER & grown_lig, Master_Score & score)
{
  // this function takes vector of conformers with all levels of growth and prints a branch of the tree
  //cout << "start AG_Conformer_Search::print_branch" << endl;
  int parent = grown_lig.parent_num;
  vector < CONFORMER > new_list;

  // set name for branch output file
  fstream foutbranch;
  char fname[80];
  if (!grown_lig.used) sprintf(fname,"%s_anchor%d_branch%d.mol2",
             grown_lig.structure.title.c_str(),current_anchor,grown_lig.conformer_num);
    else sprintf(fname,"pruned_%d_anchor_%d.mol2",grown_lig.conformer_num,current_anchor);
  foutbranch.open (fname, fstream::out);

  // grown_lig is the start of the growth tree  
  new_list.push_back(grown_lig);

  // add the b4min grown lig conformer
  for (int i=b4min.size()-1; i>=0; i--) 
    if (b4min[i].conformer_num == grown_lig.conformer_num) {
      new_list.push_back(b4min[i]);
      break;
    }

  // store the tree in a vector from grown ligand to anchor
  int n = list.size();  // list of seeds from all generations (parent conformers)

  // count_conf_num backwards in parent conformers list searching for parents till anchor is reached
  // when parent is found append to new_list and set new parent
  for (int i = n - 1; i >= 0; i--){
       if (parent == list[i].conformer_num)
       {
           parent = list[i].parent_num;
           new_list.push_back(list[i]);
           new_list.push_back(b4min[i]);
       }
  }

  // Print from anchor to grown ligand
  n = new_list.size();
  for (int i = n - 1; i >=0; i--)
  {
       // calculate anchor rmsd
       float rmsd_anchor_intial=0, rmsd_prev=0;

       // new_list[n-1].structure is the the anchor for this branch
       // calc_active_rmsd() uses the first dockmol object to fix the number of active atoms
       // so the anchor must be the first one
       rmsd_anchor_intial = calc_active_rmsd(new_list[n-1],new_list[i]);

       // calculates the rmsd between the previous growth phase and the current one
       // when i=0 i.e. anchor, previous rmsd is meaningless and so set to zero
       if (i == n-1) rmsd_prev = 0;
       else rmsd_prev = calc_active_rmsd(new_list[i+1],new_list[i]);

       //foutbranch << "########## loop num:\t\t" << i << endl;
       //foutbranch << "########## RMSD Anchor:\t\t" << rmsd_anchor_intial << endl;
       //           << "########## RMSD Previous:\t" << rmsd_prev << endl;
       ostringstream conformer_text;
       conformer_text << new_list[i].header << "########## RMSD2Prev:\t" << rmsd_prev << endl;
       conformer_text << new_list[i].structure.simplex_text << endl;

       foutbranch << conformer_text.str();

       Write_Mol2(new_list[i].structure, foutbranch);
  }
   foutbranch.close();
}


// +++++++++++++++++++++++++++++++++++++++++
void
AG_Conformer_Search::activate_layer_segment(DOCKMol & mol, int layer_num,
                                            int curr_seg)
{
    int             i,
                    j,
                    k;
    int             seg_num;

    reset_active_lists(mol);

    // activate all layers but last
    for (i = 0; i < layer_num; i++) {
        for (j = 0; j < layers[i].segments.size(); j++) {
            for (k = 0; k < layer_segments[layers[i].segments[j]].atoms.size();
                 k++) {
                mol.atom_active_flags[layer_segments[layers[i].segments[j]].
                                      atoms[k]] = true;
            }

            for (k = 0; k < layer_segments[layers[i].segments[j]].bonds.size();
                 k++) {
                mol.bond_active_flags[layer_segments[layers[i].segments[j]].
                                      bonds[k]] = true;
            }
        }
    }

    // activate last layer up to current segment
    for (i = 0; i <= curr_seg; i++) {

        seg_num = layers[layer_num].segments[i];

        for (j = 0; j < layer_segments[seg_num].atoms.size(); j++)
            mol.atom_active_flags[layer_segments[seg_num].atoms[j]] = true;

        for (j = 0; j < layer_segments[seg_num].bonds.size(); j++)
            mol.bond_active_flags[layer_segments[seg_num].bonds[j]] = true;

    }

    mol.num_active_atoms = 0;
    mol.num_active_bonds = 0;

    for (i = 0; i < mol.num_atoms; i++)
        if (mol.atom_active_flags[i])
            mol.num_active_atoms++;

    for (i = 0; i < mol.num_bonds; i++)
        if (mol.bond_active_flags[i])
            mol.num_active_bonds++;

}


// +++++++++++++++++++++++++++++++++++++++++
void
AG_Conformer_Search::reset_active_lists(DOCKMol & mol)
{
    int             i;

    for (i = 0; i < mol.num_atoms; i++)
        mol.atom_active_flags[i] = 0;

    for (i = 0; i < mol.num_bonds; i++)
        mol.bond_active_flags[i] = 0;

    mol.num_active_atoms = 0;
    mol.num_active_bonds = 0;

}


// +++++++++++++++++++++++++++++++++++++++++
// Return true if the conf is not rejected by a vDW clash check.
bool
AG_Conformer_Search::segment_clash_check(DOCKMol & conf, int layer_num,
                                         int current_bond)
{
    int             a1,
                    a2;
    float           dist_squared,
                    ref;

    int             i,
                    j,
                    k,
                    l;
    int             seg_num,
                    new_seg;

    bool            skip_flag = false;

    if (use_clash_penalty) {

        // ID current segment
        seg_num = layers[layer_num].segments[current_bond];

        // loop over atoms in newest segment
        for (i = 0; i < layer_segments[seg_num].atoms.size(); i++) {
            a1 = layer_segments[seg_num].atoms[i];

            // loop from last layer to all inner layers
            for (j = layer_num; j >= 0; j--) {

                // loop over segments in current layer
                for (k = ((j == layer_num) ? (current_bond - 1) :
                         (layers[j].segments.size() - 1)); k >= 0; k--) {

                    new_seg = layers[j].segments[k];

                    // loop over atoms in this segment
                    for (l = 0; l < layer_segments[new_seg].atoms.size(); l++) {
                        a2 = layer_segments[new_seg].atoms[l];

                        if ((conf.get_bond(a1, a2) == -1)
                            && (!conf.atoms_are_one_three(a1, a2))) {

                            dist_squared = pow((conf.x[a1] - conf.x[a2]), 2) +
                                           pow((conf.y[a1] - conf.y[a2]), 2) +
                                           pow((conf.z[a1] - conf.z[a2]), 2);
                            ref = clash_penalty * (conf.amber_at_radius[a1] +
                                                   conf.amber_at_radius[a2]);
                            skip_flag = dist_squared < ref * ref;

                            if ((!conf.atom_active_flags[a1])
                                || (!conf.atom_active_flags[a2]))
                                cout << "Layer ERROR!!!!" << endl;
                        }
                    }
                }
            }
        }
    }

    return !skip_flag;
}


// +++++++++++++++++++++++++++++++++++++++++
bool
AG_Conformer_Search::next_conformer(DOCKMol & mol)
{



    if (pruned_confs.size() > 0) {
        copy_molecule(mol, pruned_confs[pruned_confs.size() - 1].second);
        pruned_confs.pop_back();
        return true;
    } else {
        return false;
    }



}



// +++++++++++++++++++++++++++++++++++++++++
// Store fragment fingerprints, frequencies, and in certain cases the DOCKMol
// object in a map as molecules are dissected at rotatable bonds
void
AG_Conformer_Search::count_fragments(DOCKMol & mol)
{

//TODO
// Output some statistics about the input library, e.g. max number of layers, max number of
// scaffolds per layer, etc. Stuff that will help us constrict growth space during dn design

    // Make a temporary DOCKMol and Fingerprint object
    DOCKMol tmp_mol;
    Fingerprint finger;

    // Send DOCKMol to fingerprint function and return a list of allowable
    // torsion environments for this molecule only
    vector <string> tmp_torsions;
    finger.write_torsion_environments(mol, tmp_torsions);

    // Add the allowable torsions to a map
    for (int i=0; i<tmp_torsions.size(); i++){
        if (!torsions_map[tmp_torsions[i]]){
            torsions_map[tmp_torsions[i]] = 1;
        } else {
            torsions_map[tmp_torsions[i]]++;
        }
    }

    tmp_torsions.clear();


    // Loop through each of the segments in tmp_mol and activate those atoms/
    // bond plus neighbor_atoms and neighbor_bonds for the given segment
    for (int i=0; i<orig_segments.size(); i++) {

        // Copy the instance of the molecule from the mol2 (orig) onto tmp_mol
        copy_molecule(tmp_mol, orig);

        // And activate the segment + neighbor atoms / bonds
        activate_fragment(tmp_mol, i);

        // Compute a fingerprint for just the active segment of the molecule
        string tmp_string;
        tmp_string = finger.compute_unique_string_active(tmp_mol);


        // Check to see if the tmp_string has already been seen as part of the
        // following hash (if it has not yet been seen):
        if (!segment_fingerprints[tmp_string].second){

            // Then add it to the hash with a frequency of '1'
            segment_fingerprints[tmp_string].second = 1;

        } else { // Else it has already been added to the hash,

            // So increment the frequency counter
            segment_fingerprints[tmp_string].second++;
        }


        // If the frequency counter is exactly the number of the 
        // fragment_library_freq_cutoff, then remember the mol object
        // Note: 
        if (segment_fingerprints[tmp_string].second == fragment_library_freq_cutoff){

            // Translate one of the dummy atoms to the origin
            if (fragment_library_trans_origin){

                // Declare a dockvector for the translation
                DOCKVector trans_vec;
    
                // Look through the atoms to identify the first that is active 
                // and a dummy
                for (int j=0; j<tmp_mol.num_atoms; j++){
    
                    if (tmp_mol.atom_active_flags[j] && tmp_mol.atom_types[j].compare("Du")==0){
                        trans_vec.x = -tmp_mol.x[j];
                        trans_vec.y = -tmp_mol.y[j];
                        trans_vec.z = -tmp_mol.z[j];
                        break;    
                    }
                }
    
                // Translate the whole dockmol based on the coordinates of that 
                // dummy atom
                tmp_mol.translate_mol(trans_vec);
            }

            // Add the dockmol object to the hash
            copy_molecule(segment_fingerprints[tmp_string].first, tmp_mol);
        }
    }


    return;

} // end AG_Conformer_Search::count_fragments()



// +++++++++++++++++++++++++++++++++++++++++
// Set active atom & active bond flags true for each atom and bond in a segment,
// as well as each neighbor_atom and neighbor_bond of that segment
// (Neighbors are remembered as Du - dummy)
void
AG_Conformer_Search::activate_fragment(DOCKMol & mol, int curr_seg)
{
    // Reset all active atoms and bonds
    reset_active_lists(mol);


    // Activate all of the atoms in the current segment
    for (int i=0; i < orig_segments[curr_seg].atoms.size(); i++) {
          mol.atom_active_flags[orig_segments[curr_seg].atoms[i]] = true;
    }

    // Activate all of the atoms that are neighbors to the current segment
    // change the atom type to 'Du' and atom name to 'DuX' where X = 1...N
    int dummy_counter = 1;
    for (int i=0; i<orig_segments[curr_seg].neighbor_atoms.size(); i++){
        mol.atom_active_flags[orig_segments[curr_seg].neighbor_atoms[i]] = true;
        mol.atom_types[orig_segments[curr_seg].neighbor_atoms[i]] = "Du";

        // Also change the partial charge of dummy atoms to 0.0
        mol.charges[orig_segments[curr_seg].neighbor_atoms[i]] = 0.0;

        stringstream ss;
        ss << "Du" << dummy_counter;
        mol.atom_names[orig_segments[curr_seg].neighbor_atoms[i]] = ss.str();
        ss.clear();

        dummy_counter++;
    }


    // Loop over all bonds and check whether the origin and target atoms of
    // those bonds are both active (tried looping over just SEGMENT.neighbor_
    // bonds, but it was not activating the correct bonds - bug?)
    for (int i=0; i<mol.num_bonds; i++){

        // For a given bond, if both target and origin atoms are active...
        if (mol.atom_active_flags[mol.bonds_origin_atom[i]] &&
            mol.atom_active_flags[mol.bonds_target_atom[i]]   ){

            // ...then make the bond active as well.
            mol.bond_active_flags[i] = true;
        }
    }


    // Write_Mol2 needs these to be defined
    mol.num_active_atoms = 0;
    mol.num_active_bonds = 0;

    for (int i=0; i<mol.num_atoms; i++)
        if (mol.atom_active_flags[i])
            mol.num_active_atoms++;

    for (int i=0; i<mol.num_bonds; i++)
        if (mol.bond_active_flags[i])
            mol.num_active_bonds++;

    return;

} // end AG_Conformer_Search::activate_fragment()



// +++++++++++++++++++++++++++++++++++++++++
// Post-processes fragment libraries to make them unique
void
AG_Conformer_Search::write_unique_fragments()
{
   // These are the same as those in Base_Score.
   const string DELIMITER    = Base_Score::DELIMITER;
   const int    FLOAT_WIDTH  = Base_Score::FLOAT_WIDTH;
   const int    STRING_WIDTH = Base_Score::STRING_WIDTH;

    // Make a vector for sorting the torsion environments
    vector <string> tmp_torsions;
    map <string, int>::iterator iter;
    for ( iter = torsions_map.begin(); iter != torsions_map.end(); iter++ ){
        ostringstream tmp_string;
        tmp_string <<iter->first <<"-" <<iter->second;
        tmp_torsions.push_back(tmp_string.str());
    }
    sort(tmp_torsions.begin(), tmp_torsions.end());

    // Write the torsion environments to file
    ostringstream file_torenv;
    file_torenv <<fragment_library_prefix <<"_torenv.dat";
    fstream fout_torenv;
    fout_torenv.open (file_torenv.str().c_str(), fstream::out|fstream::app);
    for (int i=0; i<tmp_torsions.size(); i++){
        fout_torenv <<tmp_torsions[i] <<"\n";
    }

    // Clear some memory
    fout_torenv.close();
    tmp_torsions.clear();
    torsions_map.clear();



    // Make a vector for sorting the fragments
    vector <pair <string, int> > tmp_vector;

    // given the map, remove elements that don't have the dockmol object attached
    map < string, pair < DOCKMol , int > >::iterator it;
    for( it = segment_fingerprints.begin(); it != segment_fingerprints.end(); it++){
        if (it->second.second >= fragment_library_freq_cutoff){
            pair <string, int> tmp_pair;
            tmp_pair.first = it->first;
            tmp_pair.second = it->second.second;
            tmp_vector.push_back(tmp_pair);
        }
    }

    // Sort the vector by frequency or by fingerprint
    if (fragment_library_sort_method.compare("freq") == 0){
        sort( tmp_vector.begin(), tmp_vector.end(), frequency_sort );
    }
    else if (fragment_library_sort_method.compare("fingerprint") == 0){
        sort( tmp_vector.begin(), tmp_vector.end(), fingerprint_sort );
    }
    else {
        cout <<"Note :: Unrecognized fragment library sorting method" <<endl;
    }


    // Create the output filenames
    ostringstream file_scaffold, file_linker, file_sidechain, file_rigid;

    file_scaffold << fragment_library_prefix << "_scaffold.mol2";
    file_linker << fragment_library_prefix << "_linker.mol2";
    file_sidechain << fragment_library_prefix << "_sidechain.mol2";
    file_rigid << fragment_library_prefix << "_rigid.mol2";

    // Use those filenames to open output filestreams
    fstream fout_scaffold, fout_linker, fout_sidechain, fout_rigid;

    fout_scaffold.open (file_scaffold.str().c_str(), fstream::out|fstream::trunc);
    fout_linker.open (file_linker.str().c_str(), fstream::out|fstream::trunc);
    fout_sidechain.open (file_sidechain.str().c_str(), fstream::out|fstream::trunc);
    fout_rigid.open (file_rigid.str().c_str(), fstream::out|fstream::trunc);
    
    global_frag_index=0;//initialize global_counter for fragments

    // Iterate over the strings in the vector, to print the dockmols from the map
    for (int i=0; i<tmp_vector.size(); i++) {

        int counter = 0;
        for(int j=0; j<segment_fingerprints[tmp_vector[i].first].first.num_atoms; j++) {

            if(segment_fingerprints[tmp_vector[i].first].first.atom_types[j].compare("Du") == 0){
                counter++;
            }
        }

        // And print to the correct filehandle
        if (counter > 2) {

            ostringstream fragment_name;
            fragment_name <<"scf." <<global_frag_index;
            segment_fingerprints[tmp_vector[i].first].first.energy = fragment_name.str();
            global_frag_index++;

            fout_scaffold <<DELIMITER<< setw(STRING_WIDTH) <<"TYPE:" << setw(FLOAT_WIDTH) << "Scaffold" <<endl;
            fout_scaffold <<DELIMITER<<setw(STRING_WIDTH)<<"FR_NAME:"<<setw(FLOAT_WIDTH)<<fragment_name.str() <<endl;
            fout_scaffold <<DELIMITER<<setw(STRING_WIDTH)<<"FREQ:"<<setw(FLOAT_WIDTH)<<fixed<<tmp_vector[i].second <<endl;
            fout_scaffold <<DELIMITER<<setw(STRING_WIDTH)<<"CONN_PTS:"<<setw(FLOAT_WIDTH)<<fixed << counter <<endl;
            fout_scaffold <<DELIMITER<<endl;
            Write_Mol2(segment_fingerprints[tmp_vector[i].first].first, fout_scaffold);

        } else if (counter == 2){

            ostringstream fragment_name;
            fragment_name <<"lnk." <<global_frag_index;
            segment_fingerprints[tmp_vector[i].first].first.energy = fragment_name.str();
            global_frag_index++;

            fout_linker<<DELIMITER<<setw(STRING_WIDTH)<<"TYPE:"<<setw(FLOAT_WIDTH) <<"Linker"<<endl;
            fout_linker<<DELIMITER<<setw(STRING_WIDTH)<<"FR_NAME:"<<setw(FLOAT_WIDTH) << fragment_name.str()<<endl;
            fout_linker<<DELIMITER<<setw(STRING_WIDTH)<<"FREQ:"<<setw(FLOAT_WIDTH)<<fixed<<tmp_vector[i].second<<endl;
            fout_linker<<DELIMITER<<setw(STRING_WIDTH)<<"CONN_PTS:"<<setw(FLOAT_WIDTH)<<fixed<<counter <<endl;
            fout_linker<<DELIMITER<<endl;
            

            Write_Mol2(segment_fingerprints[tmp_vector[i].first].first, fout_linker);

        } else if (counter == 1){

            ostringstream fragment_name;
            fragment_name <<"sid." <<global_frag_index;
            segment_fingerprints[tmp_vector[i].first].first.energy = fragment_name.str();
            global_frag_index++;

            fout_sidechain <<DELIMITER<<setw(STRING_WIDTH)<<"TYPE:"<<setw(FLOAT_WIDTH) <<"Sidechain"<<endl;
            fout_sidechain <<DELIMITER<<setw(STRING_WIDTH)<<"FR_NAME:"<<setw(FLOAT_WIDTH) << fragment_name.str()<<endl;
            fout_sidechain <<DELIMITER<<setw(STRING_WIDTH)<<"FREQ:"<<setw(FLOAT_WIDTH)<<fixed<<tmp_vector[i].second <<endl;
            fout_sidechain <<DELIMITER<<setw(STRING_WIDTH)<<"CONN_PTS:"<<setw(FLOAT_WIDTH)<<fixed<<counter <<endl;
            fout_sidechain <<DELIMITER<<endl;
            Write_Mol2(segment_fingerprints[tmp_vector[i].first].first, fout_sidechain);

        } else if (counter == 0){

            ostringstream fragment_name;
            fragment_name <<"rig."  <<global_frag_index;
            segment_fingerprints[tmp_vector[i].first].first.energy = fragment_name.str();
            global_frag_index++;

            fout_rigid <<DELIMITER<<setw(STRING_WIDTH)<<"TYPE:"<<setw(FLOAT_WIDTH)<<"Rigid"<<endl;
            fout_rigid <<DELIMITER<<setw(STRING_WIDTH)<<"FR_NAME:"<<setw(FLOAT_WIDTH)<< fragment_name.str()<<endl;
            fout_rigid <<DELIMITER<<setw(STRING_WIDTH)<<"FREQ:"<<setw(FLOAT_WIDTH)<<fixed<<tmp_vector[i].second <<endl;
            fout_rigid <<DELIMITER<<setw(STRING_WIDTH)<<"CONN_PTS:"<<setw(FLOAT_WIDTH)<<fixed<<counter <<endl;
            fout_rigid <<DELIMITER<<endl;
            Write_Mol2(segment_fingerprints[tmp_vector[i].first].first, fout_rigid);

        } else {
            cout <<"Warning :: This type of molecule should not exist " << endl;
        }
    }


    // Clear some memory
    fout_scaffold.close();
    fout_linker.close();
    fout_sidechain.close();
    fout_rigid.close();
    tmp_vector.clear();
    segment_fingerprints.clear();

    return;

} // end AG_Conformer_Search::write_unique_fragments()



// +++++++++++++++++++++++++++++++++++++++++
// Sort functions called in write_unique_fragments()
int
frequency_sort(pair <string, int> a, pair <string, int> b)
{
    return (a.second > b.second);
}

int
fingerprint_sort(pair <string, int> a, pair <string, int> b)
{
    return (a.first > b.first);
}

