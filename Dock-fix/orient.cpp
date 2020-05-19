#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include "orient.h"

using namespace std;


/*** matfit stuff ***/
#define SMALL  1.0e-20
#define SMALSN 1.0e-10
#define ABS(x)   (((x)<0)   ? (-(x)) : (x))

/************************************************************************/
static void
minimized_fit(double umat[3][3], double rm[3][3])
{
    // Sudipto: could someone add comments to explain what each varible
    // here describes, and more details about the algorithm used?

    double          rot[3][3],
                    turmat[3][3],
                    c[3][3],
                    coup[3],
                    dir[3],
                    step[3],
                    v[3],
                    rtsum,
                    rtsump,
                    rsum,
                    stp,
                    stcoup,
                    ud,
                    tr,
                    ta,
                    cs,
                    sn,
                    ac,
                    delta,
                    deltap,
                    gfac,
                    cle,
                    clep;
    int             i,j,k,l,m,  //loop counters
                    jmax,
                    ncyc,
                    nsteep,
                    nrem;

    /*
     * Rotate repeatedly to reduce couple about initial direction to zero.
     * Clear the rotation matrix
     */
    for (l = 0; l < 3; l++) {
        for (m = 0; m < 3; m++)
            rot[l][m] = 0.0;
        rot[l][l] = 1.0;
    }

    /*
     * Copy vmat[][] (sp) into umat[][] (dp) 
     */
    jmax = 30;
    rtsum = umat[0][0] + umat[1][1] + umat[2][2];
    delta = 0.0;

    for (ncyc = 0; ncyc < jmax; ncyc++) {
        /*
         * Modified CG. For first and every NSTEEP cycles, set previous step as 
         * zero and do an SD step 
         */
        nsteep = 3;
        nrem = ncyc - nsteep * (int) (ncyc / nsteep);

        if (!nrem) {
            for (i = 0; i < 3; i++)
                step[i] = 0.0;
            clep = 1.0;
        }

        /*
         * Couple 
         */
        coup[0] = umat[1][2] - umat[2][1];
        coup[1] = umat[2][0] - umat[0][2];
        coup[2] = umat[0][1] - umat[1][0];
        cle = sqrt(coup[0] * coup[0] + coup[1] * coup[1] + coup[2] * coup[2]);

        /*
         * Gradient vector is now -coup 
         */
        gfac = (cle / clep) * (cle / clep);

        /*
         * Value of rtsum from previous step 
         */
        rtsump = rtsum;
        deltap = delta;
        clep = cle;
        if (cle < SMALL){
            //cout <<ncyc<< ": cle < SMALL:" << cle << "<" << SMALL << endl;
            break;
        }

        /*
         * Step vector conjugate to previous 
         */
        stp = 0.0;
        for (i = 0; i < 3; i++) {
            step[i] = coup[i] + step[i] * gfac;
            stp += (step[i] * step[i]);
        }
        stp = 1.0 / sqrt(stp);

        /*
         * Normalised step 
         */
        for (i = 0; i < 3; i++)
            dir[i] = stp * step[i];

        /*
         * Couple resolved along step direction 
         */
        stcoup = coup[0] * dir[0] + coup[1] * dir[1] + coup[2] * dir[2];

        /*
         * Component of UMAT along direction 
         */
        ud = 0.0;
        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++)
                ud += umat[l][m] * dir[l] * dir[m];


        tr = umat[0][0] + umat[1][1] + umat[2][2] - ud;
        ta = sqrt(tr * tr + stcoup * stcoup);
        cs = tr / ta;
        sn = stcoup / ta;

        /*
         * If cs<0 then posiiton is unstable, so don't stop 
         */
        if ((cs > 0.0) && (ABS(sn) < SMALSN)){
            //cout <<ncyc <<": (cs > 0.0) && (ABS(sn) < SMALSN):"<< cs << ">"<< 0.0 << " && " << ABS(sn) << "<" << SMALSN << endl;
            break;
        }

        /*
         * Turn matrix for correcting rotation:
         * 
         * Symmetric part 
         */
        ac = 1.0 - cs;
        for (l = 0; l < 3; l++) {
            v[l] = ac * dir[l];
            for (m = 0; m < 3; m++)
                turmat[l][m] = v[l] * dir[m];
            turmat[l][l] += cs;
            v[l] = dir[l] * sn;
        }

        /*
         * Asymmetric part 
         */
        turmat[0][1] -= v[2];
        turmat[1][2] -= v[0];
        turmat[2][0] -= v[1];
        turmat[1][0] += v[2];
        turmat[2][1] += v[0];
        turmat[0][2] += v[1];

        /*
         * Update total rotation matrix 
         */
        for (l = 0; l < 3; l++) {
            for (m = 0; m < 3; m++) {
                c[l][m] = 0.0;
                for (k = 0; k < 3; k++)
                    c[l][m] += turmat[l][k] * rot[k][m];
            }
        }

        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++)
                rot[l][m] = c[l][m];

        /*
         * Update umat tensor 
         */
        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++) {
                c[l][m] = 0.0;
                for (k = 0; k < 3; k++)
                    c[l][m] += turmat[l][k] * umat[k][m];
            }

        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++)
                umat[l][m] = c[l][m];

        rtsum = umat[0][0] + umat[1][1] + umat[2][2];
        delta = rtsum - rtsump;

        /*
         * If no improvement in this cycle then stop 
         */
        if (ABS(delta) < SMALL){
            //cout <<ncyc<< ": ABS(delta) < SMALL: " << ABS(delta) << "<" << SMALL << endl;
            break;
        }

        /*
         * Next cycle 
         */
    }

    rsum = rtsum;

    /*
     * Copy rotation matrix for output 
     */

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            rm[i][j] = rot[i][j];       // can be transposed
}

/*************************************************************************/
int
compute_rot_matrix(XYZVec & x1, XYZVec & x2, double rm[3][3], int n)
{
    int             i,
                    j;
    double          umat[3][3];


    if (n < 2) {
        return (0);
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            umat[i][j] = 0.0;
    }

    for (j = 0; j < n; j++) {
        umat[0][0] += x1[j].x * x2[j].x;
        umat[1][0] += x1[j].y * x2[j].x;
        umat[2][0] += x1[j].z * x2[j].x;

        umat[0][1] += x1[j].x * x2[j].y;
        umat[1][1] += x1[j].y * x2[j].y;
        umat[2][1] += x1[j].z * x2[j].y;

        umat[0][2] += x1[j].x * x2[j].z;
        umat[1][2] += x1[j].y * x2[j].z;
        umat[2][2] += x1[j].z * x2[j].z;
    }

    minimized_fit(umat, rm);

    return (1);
}


/************************************************/
/************************************************/
/************************************************/
Orient::Orient()
{

    sph_dist_mat = NULL;
    lig_dist_mat = NULL;
    residual_mat = NULL;

}
/************************************************/
Orient::~Orient()
{

    delete[]sph_dist_mat;
    delete[]lig_dist_mat;
    delete[]residual_mat;

}
/************************************************/
void
Orient::input_parameters(Parameter_Reader & parm)
{
    cout << "\nOrient Ligand Parameters\n";
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    use_chemical_matching = false;
    orient_ligand = parm.query_param("orient_ligand", "yes", "yes no") == "yes";

    if (orient_ligand) {
        automated_matching = 
            parm.query_param("automated_matching", "yes", "yes no") == "yes";

        if (automated_matching) {
            tolerance = 0.25;
            dist_min = 2.0;
            min_nodes = 3;
            max_nodes = 10;
        } else {
            tolerance =
                atof(parm.query_param("distance_tolerance", "0.25").c_str());
            if (tolerance <= 0.0) {
                cout <<
                    "ERROR: Parameter must be a float greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
            dist_min =
                atof(parm.query_param("distance_minimum", "2.0").c_str());
            if (dist_min <= 0.0) {
                cout <<
                    "ERROR: Parameter must be a float greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
            min_nodes = atoi(parm.query_param("nodes_minimum", "3").c_str());
            if (min_nodes <= 0) {
                cout <<
                    "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
            max_nodes = atoi(parm.query_param("nodes_maximum", "10").c_str());
            if (max_nodes <= 0) {
                cout <<
                    "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }

        }

        Active_Site_Spheres :: set_sphere_file_name( parm );

        max_orients = atoi(parm.query_param("max_orientations", "1000").c_str());
        if (max_orients <= 0) {
            cout <<
                "ERROR: Parameter must be a integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }

        critical_points =
            parm.query_param("critical_points", "no", "yes no") == "yes";

        use_chemical_matching =
            parm.query_param("chemical_matching", "no", "yes no") == "yes";
        if (use_chemical_matching) {
            chem_match_tbl_fname =
                parm.query_param("chem_match_tbl", "chem_match.tbl");
        }

        use_ligand_spheres =
            parm.query_param("use_ligand_spheres", "no", "yes no") == "yes";
        if (use_ligand_spheres) {
            lig_sphere_filename =
                parm.query_param("ligand_sphere_file", "ligand.sph");
        }

        verbose = 0 != parm.verbosity_level();   // -v is for verbose flag
        //TODO: Print more info about orienting progress on verbose
    }
}

/************************************************/
void
Orient::initialize(int argc, char **argv)
{

    if (orient_ligand) {
        cout << "Initializing Orienting Routines...\n";

        orig_tolerance = tolerance;
        prepare_receptor();

        // check validity of manual matching
        if (max_nodes < min_nodes) {
            cout << endl <<
                "ERROR:  Invalid range of nodes.  Program will terminate." <<
                endl;
            exit(0);
        }
    }
}

/************************************************/
void
Orient::prepare_receptor()
{
    get_spheres();
    calculate_sphere_distance_matrix();
    if (use_chemical_matching)
        read_chem_match_tbl();
}

/************************************************/
void
Orient::clean_up()
{
    current_clique = 0; //TEB 2012-02-12
    delete[]lig_dist_mat;
    delete[]residual_mat;
    lig_dist_mat = NULL;
    residual_mat = NULL;
}

/************************************************/
void
Orient::get_spheres()
{
    bool            found_cluster;
    CRITICAL_CLUSTER tmp_cluster;

    spheres = Active_Site_Spheres :: get_instance();
    num_spheres = spheres.size();

    if (verbose) cout << "Read in " << num_spheres 
                      << " spheres for orienting." << endl;

    // process critical points
    if (critical_points) {
        receptor_critical_clusters.clear();

        // loop over spheres
        for (int i = 0; i < spheres.size(); i++) {

            // if sphere belongs to a critical cluster
            if (spheres[i].critical_cluster > 0) {

                found_cluster = false;

                for (int j = 0; j < receptor_critical_clusters.size(); j++) {
                    // if the cluster has been identified previously
                    if (spheres[i].critical_cluster ==
                        receptor_critical_clusters[j].index) {
                        // add the sphere to the existing cluster
                        receptor_critical_clusters[j].spheres.push_back(i);
                        found_cluster = true;
                        break;
                    }
                }

                // else create new cluster entry
                if (!found_cluster) {
                    tmp_cluster.index = spheres[i].critical_cluster;
                    tmp_cluster.spheres.clear();
                    tmp_cluster.spheres.push_back(i);
                    receptor_critical_clusters.push_back(tmp_cluster);
                }
            }
        }
    }
}

/************************************************/
void
Orient::calculate_sphere_distance_matrix()
{

    sph_dist_mat = new double[num_spheres * num_spheres];

    for (int i = 0; i < num_spheres; i++)
        for (int j = 0; j < num_spheres; j++) {
            sph_dist_mat[num_spheres * i + j] = spheres[i].distance(spheres[j]);
        }
}

/************************************************/
// Function will use active heavy atoms or dummy atoms as spheres
void
Orient::get_centers(DOCKMol & mol)
{
    centers.clear();
    num_centers = 0;

    if (!use_ligand_spheres) {

        Sphere tmp;
        for (int atom = 0; atom < mol.num_atoms; atom++) {
            // CDS-09/26/16: added dummy atom
            if (mol.amber_at_heavy_flag[atom] || mol.atom_types[atom] == "Du") {
                if (mol.atom_active_flags[atom]) {
                    tmp.crds.x = mol.x[atom];
                    tmp.crds.y = mol.y[atom];
                    tmp.crds.z = mol.z[atom];
                    tmp.radius = 0.0;
                    tmp.surface_point_i = 0;
                    tmp.surface_point_j = 0;
                    tmp.critical_cluster = 0;
                    if (use_chemical_matching)
                        tmp.color = mol.chem_types[atom];
                    centers.push_back(tmp);
                    num_centers++;
                }
            }
        }

    } else {
        //read in ligand spheres 
        num_centers = read_spheres( lig_sphere_filename, centers );
    }
    // print the number of anchor heavy atoms
    if (verbose) cout << "Orienting " << num_centers 
                      << " anchor heavy atom centers" << endl;
}

/************************************************/
// Function for ligand to ligand matching, which will prepare the reference ligand spheres
// similar to get_spheres().  Added by CDS - 09/24/16
// If Critical Points is true, it will make the dummy spheres critical only! 
void
Orient::get_lig_reference_spheres(DOCKMol & mol)
{
   orig_tolerance = tolerance;
   // Use the ligand reference function to identify the heavy atoms/spheres
   // Clear spheres
   spheres.clear();
   num_spheres = 0;

   if (!use_ligand_spheres) {

      Sphere tmp;
      int critical_cluster = 1; // place holder to update critical cluster
      for (int atom = 0; atom < mol.num_atoms; atom++) {
          // CDS-09/26/16: added dummy atom
          if (mol.amber_at_heavy_flag[atom] || mol.atom_types[atom] == "Du") {
              if (mol.atom_active_flags[atom]) {
                 tmp.crds.x = mol.x[atom];
                 tmp.crds.y = mol.y[atom];
                 tmp.crds.z = mol.z[atom];
                 tmp.radius = 0.0;
                 tmp.surface_point_i = 0;
                 tmp.surface_point_j = 0;
                 // If there is a dummy atom, make it a critical cluster
                 if ( mol.atom_types[atom] == "Du" ){
                    tmp.critical_cluster = critical_cluster;
                    critical_cluster++;
                 }
                 else{
                    tmp.critical_cluster = 0;
                 }
                 if (use_chemical_matching)
                    tmp.color = mol.chem_types[atom];
                 spheres.push_back(tmp);
                 num_spheres++;
              }
          }
      }

    } else {
        //read in ligand spheres 
        num_spheres = read_spheres( lig_sphere_filename, centers );
    }
    // print the number of anchor heavy atoms
    if (verbose) cout << "Orienting " << num_spheres 
                      << " anchor heavy atom centers" << endl;



   // The code below is from get_spheres()
   bool            found_cluster;
   CRITICAL_CLUSTER tmp_cluster;

   if (verbose) cout << "Read in " << num_spheres 
                     << " spheres for orienting." << endl;

   // process critical points
   if (critical_points) {
       receptor_critical_clusters.clear();

       // loop over spheres
       for (int i = 0; i < spheres.size(); i++) {

           // if sphere belongs to a critical cluster
           if (spheres[i].critical_cluster > 0) {

               found_cluster = false;

               for (int j = 0; j < receptor_critical_clusters.size(); j++) {
                   // if the cluster has been identified previously
                   if (spheres[i].critical_cluster ==
                       receptor_critical_clusters[j].index) {
                       // add the sphere to the existing cluster
                       receptor_critical_clusters[j].spheres.push_back(i);
                       found_cluster = true;
                       break;
                   }
               }

               // else create new cluster entry
               if (!found_cluster) {
                   tmp_cluster.index = spheres[i].critical_cluster;
                   tmp_cluster.spheres.clear();
                   tmp_cluster.spheres.push_back(i);
                   receptor_critical_clusters.push_back(tmp_cluster);
               }
           }
       }
   }
   calculate_sphere_distance_matrix();
   if (use_chemical_matching)
      read_chem_match_tbl();
}

/************************************************/
void
Orient::calculate_ligand_distance_matrix()
{

    lig_dist_mat = new double[num_centers * num_centers];

    for (int i = 0; i < num_centers; i++)
        for (int j = 0; j < num_centers; j++) {
            lig_dist_mat[num_centers * i + j] = centers[i].distance(centers[j]);
        }
}

/************************************************/
// This fuction is called in the main dock loop.

void
Orient::match_ligand(DOCKMol & mol)
{
    int             s1,
                    s2,
                    c1,
                    c2;
    int             i,
                    j,
                    k;
    INTVec          cmt_sphere_idx,
                    cmt_center_idx;

    if (orient_ligand) {

        if (verbose) cout << "-----------------------------------" << endl 
             << "VERBOSE ORIENTING STATS :" << endl << endl;

        //cached_orient.clear_molecule();
        original.clear_molecule();
        centers.clear();
        clique_spheres.clear();
        clique_centers.clear();

        num_orients = 0;
        last_orient_flag = false;

        copy_molecule(original, mol);

        get_centers(original);

        // calc chem matching table if chem matching is used
        if (use_chemical_matching) {

            // idx = center*num_spheres + sphere
            chem_match_align_tbl.clear();
            chem_match_align_tbl.resize((spheres.size() * centers.size()), 0);

            // assign chemical match type indices to the spheres
            cmt_sphere_idx.clear();
            cmt_sphere_idx.resize(spheres.size(), 0);

            // loop over spheres
            for (i = 0; i < spheres.size(); i++) {
                // loop over match types
                for (j = 0; j < chem_match_tbl_labels.size(); j++) {
                    // if types match, make the assignment
                    if (spheres[i].color == chem_match_tbl_labels[j]) {
                        cmt_sphere_idx[i] = j;
                    }
                }
            }

            // assign chemical match type indices to the centers
            cmt_center_idx.clear();
            // bug fix from revision 1.32 and 1.33; srb
            cmt_center_idx.resize( centers.size(), 0 );

            // loop over spheres
            for (i = 0; i < centers.size(); i++) {
                // loop over match types
                for (j = 0; j < chem_match_tbl_labels.size(); j++) {
                    // if types match, make the assignment
                    if (centers[i].color == chem_match_tbl_labels[j]) {
                        cmt_center_idx[i] = j;
                    }
                }
            }

            // generate table of legal sphere/center pairings
            for (i = 0; i < spheres.size(); i++) {
                for (j = 0; j < centers.size(); j++) {
                    chem_match_align_tbl[j * spheres.size() + i] =
                        chem_match_tbl_matrix[cmt_sphere_idx[i] *
                                              chem_match_tbl_labels.size() +
                                              cmt_center_idx[j]];
                }
            }

        }
        // end chemical matching code

        calculate_ligand_distance_matrix();

	// num_centers is the number of anchor heavy/dummy atoms
        num_nodes = num_spheres * num_centers;

        // initialize matrix full of node-node residuals
        if (num_nodes > sqrt(INT_MAX)) {
            // The new below will be hit by an integer overflow
            // and behave in an undefined way.
            //
            // TEST shows a case that produces a postive but wrong value when using
            //      (int) causing the new to possibly succeed with the wrong amount
            //      of memory. The overflow can also produce negative values, but
            //      new would catch that case.
            //cout << "TEST: " << (int)((int)65600)*((int)65600) << endl;
            //cout << "TEST: " << (long)((long)65600)*((long)65600) << endl;
            cout << "Too many spheres or atoms: num_nodes = num_spheres * num_centers = " << num_nodes << endl
                 << "num_spheres = " << num_spheres << endl
                 << "num_centers = " << num_centers << endl;
            cerr << "Too many spheres or atoms: num_nodes = num_spheres * num_centers = " << num_nodes << endl
                 << "num_spheres = " << num_spheres << endl
                 << "num_centers = " << num_centers << endl;
            throw bad_alloc();
        }
        residual_mat = new double[num_nodes * num_nodes];

        k = 0;

	// populate residual matrix
        for (i = 0; i < num_nodes; i++) {
            s1 = i % num_spheres;
            c1 = i / num_spheres;

            for (j = 0; j < num_nodes; j++) {
                s2 = j % num_spheres;
                c2 = j / num_spheres;

                residual_mat[k] =
                    fabs(sph_dist_mat[s1 * num_spheres + s2] -
                         lig_dist_mat[c1 * num_centers + c2]);

                k++;
            }
        }

        level = 0;
        am_iteration_num = 1;

        // perform clique detection
        id_all_cliques();

    } else
        last_orient_flag = false;
}

/************************************************/
void
Orient::calculate_translations()
{

    double sph_com_x = 0.0;
    double sph_com_y = 0.0;
    double sph_com_z = 0.0;
    double cen_com_x = 0.0;
    double cen_com_y = 0.0;
    double cen_com_z = 0.0;

    for (int i = 0; i < clique_size; i++) {
        cen_com_x += clique_centers[i].x;
        cen_com_y += clique_centers[i].y;
        cen_com_z += clique_centers[i].z;
        sph_com_x += clique_spheres[i].x;
        sph_com_y += clique_spheres[i].y;
        sph_com_z += clique_spheres[i].z;
    }

    sph_com_x = sph_com_x / clique_size;
    sph_com_y = sph_com_y / clique_size;
    sph_com_z = sph_com_z / clique_size;
    cen_com_x = cen_com_x / clique_size;
    cen_com_y = cen_com_y / clique_size;
    cen_com_z = cen_com_z / clique_size;

    spheres_com.x = sph_com_x;
    spheres_com.y = sph_com_y;
    spheres_com.z = sph_com_z;
    centers_com.x = cen_com_x;
    centers_com.y = cen_com_y;
    centers_com.z = cen_com_z;
}

/************************************************/
void
Orient::translate_clique_to_origin()
{

    for (int i = 0; i < clique_size; i++) {

        clique_centers[i].x = clique_centers[i].x - centers_com.x;
        clique_centers[i].y = clique_centers[i].y - centers_com.y;
        clique_centers[i].z = clique_centers[i].z - centers_com.z;

        clique_spheres[i].x = clique_spheres[i].x - spheres_com.x;
        clique_spheres[i].y = clique_spheres[i].y - spheres_com.y;
        clique_spheres[i].z = clique_spheres[i].z - spheres_com.z;
    }
}

/************************************************/
void
Orient::calculate_rotation()
{
    int n = clique_size;
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            rotation_matrix[i][j] = 0.0;

    compute_rot_matrix(clique_centers, clique_spheres, rotation_matrix, n);
}

/************************************************/
bool
Orient::more_orientations()
{
    if (orient_ligand) {
        return ! last_orient_flag;
    } else
        return false;
}

/************************************************/
// called in match_ligand()

void
Orient::id_all_cliques()
{
    int             i,
                    j,
                    k,
                    size;
    int             next_cand,
                    index;
    // int notcount;
    CLIQUE          tmp_clique;
    double          *new_level_residuals;
    double           tmp_resid;

    int limit_cliques = 0;
    am_iteration_num = 1;
    cliques.clear();

    // init arrays
    size = num_centers * num_spheres;
    bool *new_candset = new bool[size * max_nodes];
    bool *new_notset = new bool[size * max_nodes];
    int *new_state = new int[max_nodes];
    new_level_residuals = new double[max_nodes];


    if (verbose){
           cout <<"Sphere Center Matching Parameters:" << endl
                <<"   tolerance: "<<tolerance 
                <<"; dist_min: "<< dist_min
                <<"; min_nodes: "<<min_nodes
                <<"; max_nodes: "<<max_nodes<< endl;
    }
    // Loop over automated matching loop &&(am_iteration_num < 2)
    while ((cliques.size() <= max_orients)
           && ((automated_matching) || (am_iteration_num == 1))
           && (am_iteration_num < 10)) {

        tolerance = am_iteration_num * orig_tolerance;

        for (i = 0; i < size * max_nodes; i++) {
            if (i < size)
                new_candset[i] = true;
            else
                new_candset[i] = false;

            new_notset[i] = false;

            if (i < max_nodes) {
                new_state[i] = -1;
                new_level_residuals[i] = 0.0;
            }
        }

        int new_level = 0;

        // main loop over levels
        while (new_level > -1) {
            // find next true in candset
            next_cand = -1;
            for (i = new_state[new_level] + 1; (i < size) && (next_cand == -1);
                 i++) {
                index = new_level * size + i;

                if ((new_candset[index] == true)
                    && (new_notset[index] == false))
                    next_cand = i;

                // compute residuals
                if (next_cand != -1) {

                    // compute residuals
                    if (new_level > 0)
                        new_level_residuals[new_level] =
                            new_level_residuals[new_level - 1];
                    else if (new_level == 0)
                        new_level_residuals[new_level] = 0.0;
                    else 
                        cout << "ERROR new_level is negative." << endl;

                    for (j = new_level - 1; j > -1; j--) {
                        k = next_cand * size + new_state[j];
                        // sum total resid method
                        new_level_residuals[new_level] += residual_mat[k];

                        // single max resid method
                        // if(residual_mat[k] > new_level_residuals[new_level])
                        // new_level_residuals[new_level] = residual_mat[k];
                    }

                    if (new_level_residuals[new_level] > tolerance) {
                        next_cand = -1;
                        new_candset[index] = false;
                        new_notset[index] = true;
                    }

                }
                // end compute residuals
            }

            // if a candidate node is found at the current level
            if (next_cand != -1) {

                new_state[new_level] = next_cand;
                new_candset[new_level * size + new_state[new_level]] = false;

                new_level++;
                if (new_level < max_nodes)
                    new_state[new_level] = new_state[new_level - 1];
                if (new_level >= max_nodes) { 
                    // add state to cliques if tolerance is proper
                    if ((new_level_residuals[new_level - 1] >
                         orig_tolerance * (am_iteration_num - 1))
                        && (new_level_residuals[new_level - 1] <= tolerance)) {
                        limit_cliques++;

                        tmp_clique.nodes.clear();
                        tmp_clique.residual =
                            new_level_residuals[new_level - 1];
                        for (i = 0; i < new_level; i++)
                            tmp_clique.nodes.push_back(new_state[i]);

                        if (check_clique_critical_points(tmp_clique))
                            if (check_clique_chemical_match(tmp_clique))
                                cliques.push_back(tmp_clique);

                    }
                    // end clique add code

                    // if you've reached the end of the tree and still have
                    // nodes to add
                    new_level--;

                    if (new_level >= 0)
                        new_notset[new_level * size + new_state[new_level]] =
                            true;

                } else {

                    // recompute candset and notset
                    for (i = 0; i < size; i++) {
                        index = new_state[new_level - 1] * size + i;
                        new_candset[new_level * size + i] =
                            ((new_candset[(new_level - 1) * size + i]));
                        new_notset[new_level * size + i] =
                            ((new_notset[(new_level - 1) * size + i]));

                        // compute residuals
                        if ((new_candset[new_level * size + i])
                            || (new_notset[new_level * size + i])) {

                            tmp_resid = new_level_residuals[new_level - 1];

                            for (j = new_level - 1; j > -1; j--) {
                                k = i * size + new_state[j];

                                // sum total resid method
                                tmp_resid += residual_mat[k];

                                // single max resid method
                                // if(residual_mat[k] > tmp_resid)
                                // tmp_resid = residual_mat[k];

                            }

                            if (tmp_resid > tolerance) {
                                new_candset[new_level * size + i] = false;
                                new_notset[new_level * size + i] = false;
                            }
                        }
                        // End residual comp
                    }

                }

            } else {            // if no candidates are found
                //cout << "a";
                if (new_level >= min_nodes) {

                    // add state to cliques
                    if ((new_level_residuals[new_level - 1] >
                         orig_tolerance * (am_iteration_num - 1))
                        && (new_level_residuals[new_level - 1] <= tolerance)) {

                        tmp_clique.nodes.clear();
                        tmp_clique.residual =
                            new_level_residuals[new_level - 1];
                        for (i = 0; i < new_level; i++)
                            tmp_clique.nodes.push_back(new_state[i]);

                        // perform critical point checking
                        if (check_clique_critical_points(tmp_clique))
                            if (check_clique_chemical_match(tmp_clique))
                                cliques.push_back(tmp_clique);

                    }
                    // end clique add code
                }

                new_state[new_level] = -1;
                new_level--;
            }

            // end if 100X the # of orients are found.  This is mostly for the
            // dense tree cases.
            if (limit_cliques >= 100 * max_orients) {
                cout << "Warning:  Match Search Truncated due to too many " <<
                    max_nodes << " cliques." << endl;
                break;
            }

        }

        am_iteration_num++;

    } // End automated matching loop

    // sort by residuals
    sort(cliques.begin(), cliques.end());

    if (verbose){
        cout.precision(4);
        cout << fixed;
        cout << "Num of cliques generated: " << cliques.size() << endl;
        cout << " Residual Info:"  << endl;
        cout << "   min residual:    " << cliques[0].residual << endl;
        cout << "   median residual: " << cliques[(int)(cliques.size() / 2)].residual << endl;
        cout << "   max residual:    " << cliques[cliques.size() - 1].residual << endl;
        double temp_resid_sum=0;
        double temp_resid2_sum=0;
        int max_node_size = 0;
        int min_node_size = 999999;
        double node_size_sum = 0;
        for (int i = 0; i < cliques.size(); i++){
            temp_resid_sum += cliques[i].residual;
            temp_resid2_sum += pow(cliques[i].residual,2);
            if (min_node_size > cliques[i].nodes.size())
                 min_node_size = cliques[i].nodes.size();
            if (max_node_size < cliques[i].nodes.size())
                 max_node_size = cliques[i].nodes.size();
            node_size_sum += cliques[i].nodes.size();
        }
        double mean = temp_resid_sum / cliques.size();
        double std  = sqrt(temp_resid2_sum/ cliques.size() - pow(mean,2));
        cout << "   mean residual:   " << mean << endl;
        cout << "   std residual:    "  << std  << endl;
        cout << " Node Sizes:"  << endl;
        cout << "   min nodes:    "  << min_node_size << endl;
        cout << "   max nodes:    "  << max_node_size << endl;
        cout << "   mean nodes:   " << node_size_sum / cliques.size() << endl;
        cout.unsetf ( ios_base::fixed  ); 
             
    }

    current_clique = 0;

    // clean up arrays
    delete[]new_candset;
    delete[]new_notset;
    delete[]new_state;
    delete[]new_level_residuals;
    new_candset = NULL;
    new_notset = NULL;
    new_state = NULL;
    new_level_residuals = NULL;

}

/************************************************/
void
Orient::new_extract_coords_from_clique(CLIQUE & clique)
{
    // XYZCRD tmp;
    Sphere          tmp;
    int             index,
                    center,
                    sphere;
    int             i;

    clique_spheres.clear();
    clique_centers.clear();
    clique_size = clique.nodes.size();

    for (i = 0; i < clique_size; i++) {
        index = clique.nodes[i];
/*
        sphere = index / num_centers;
        center = index % num_centers;
*/
        // replace it with the proper sphere/center indexing
        sphere = index % num_spheres;
        center = index / num_spheres;
        // DTM - End removal of faulty code - 1/30/07

        tmp = spheres[sphere];
        clique_spheres.push_back(tmp.crds);

        tmp = centers[center];
        clique_centers.push_back(tmp.crds);
    }

}

/************************************************/
bool
Orient::check_clique_chemical_match(CLIQUE & clique)
{
    int             i;
    int             idx,
                    sphere_idx,
                    center_idx;

    if (use_chemical_matching) {

        // loop over clique nodes
        for (i = 0; i < clique.nodes.size(); i++) {

            // extract sphere & center indices
            idx = clique.nodes[i];
/*
            sphere_idx = idx / num_centers;
            center_idx = idx % num_centers;
*/
            // correct orinting: sudipto & DTM
            // replace it with the proper sphere/center indexing
            sphere_idx = idx % num_spheres;
            center_idx = idx / num_spheres;
            // DTM - End removal of faulty code - 1/30/07

            // check the chem_match_align_tbl - return false if match is
            // illegal
            idx = center_idx * num_spheres + sphere_idx;

            if (chem_match_align_tbl[idx] == 0)
                return false;

        }

        // else return true if match is legal
        return true;

    } else
        return true;

}

/************************************************/
bool
Orient::check_clique_critical_points(CLIQUE & clique)
{
    Sphere          tmp;
    int             index,
                    sphere;
    int             i,
                    j,
                    k;
    INTVec          tmp_clique_spheres;
    INTVec          hits;
    int             hit_sum;

    if (critical_points) {

        tmp_clique_spheres.clear();

        for (i = 0; i < clique.nodes.size(); i++) {
            index = clique.nodes[i];
            sphere = index % num_spheres;    // correct orienting: sudipto & DTM
            tmp_clique_spheres.push_back(sphere);
        }

        hits.clear();
        hits.resize(receptor_critical_clusters.size(), 0);

        // loop over the clusters
        for (i = 0; i < hits.size(); i++) {

            // loop over the cluster spheres
            for (j = 0; j < receptor_critical_clusters[i].spheres.size(); j++) {

                // loop over the clique spheres
                for (k = 0; k < tmp_clique_spheres.size(); k++) {

                    if (receptor_critical_clusters[i].spheres[j] ==
                        tmp_clique_spheres[k]) {
                        hits[i] = 1;
                        break;
                    }
                }

                if (hits[i] == 1)
                    break;
            }
        }

        hit_sum = 1;

        for (i = 0; i < hits.size(); i++) {
            hit_sum *= hits[i];
        }

        if (hit_sum == 0)
            return false;
        else
            return true;

    } else {
        return true;
    }
}

/************************************************/
// Called in main loop in dock.cpp.  Is a condition in while loop.

bool
Orient::new_next_orientation(DOCKMol & mol)
{
    //cout << "new_next_orientation" << endl;
    if (orient_ligand) {

        // in case no cliques could be found
        if (cliques.size() == 0) {
	    if (verbose) cout << "No orients found for current anchor" << endl;
            return false;
	}

        if (last_orient_flag) {
	    //if (verbose) cout << "Current clique:" << current_clique << endl;
	    //if (verbose) cout << "(last_orient_flag==true)Current clique:" << current_clique << endl;
            clean_up();
	    //if (verbose) cout << "(last_orient_flag==true)Current clique:" << current_clique << endl;
            return false;
        }

        new_extract_coords_from_clique(cliques[current_clique]);
        calculate_translations();
        translate_clique_to_origin();
        calculate_rotation();

        copy_molecule(mol, original);
        mol.translate_mol(-centers_com);
        mol.rotate_mol(rotation_matrix);
        mol.translate_mol(spheres_com);

        current_clique++;

        // DTM - change this line to ensure all cliques are examined as orients (stop skipping the last one) - 1/30/07
        if ((current_clique == max_orients) || (current_clique == cliques.size())){
	    //if (verbose) cout << "Current clique:" << current_clique << endl;
            last_orient_flag = true;
        } else{
            last_orient_flag = false;
            
        }

        return true;

    } else {
        // code to return the input mol once, and then return false after that
        // (allow one pass through the loop)
        if (!last_orient_flag) {
            last_orient_flag = true;
            return true;
        } else {
            return false;
        }
    }

}

/************************************************/
void
Orient::read_chem_match_tbl()
{
    FILE           *ifp;
    char            line[100],
                    tmp[100];
    string          tmp_string;
    int             i,
                    j,
                    z;

    if (use_chemical_matching) {

        // open match table file and loop over lines
        ifp = fopen(chem_match_tbl_fname.c_str(), "r");

        if (ifp == NULL) {
            cout << "\n\nCould not open " << chem_match_tbl_fname <<
                " for reading.  Program will terminate." << endl << endl;
            exit(0);
        }

        while (fgets(line, 100, ifp) != NULL) {

            // read in match table chemical labels
            if (!strncmp(line, "label", 5)) {
                sscanf(line, "%*s %s", tmp);
                tmp_string = tmp;
                chem_match_tbl_labels.push_back(tmp_string);
            }

            chem_match_tbl_matrix.
                resize((chem_match_tbl_labels.size() *
                        chem_match_tbl_labels.size()), 0);

            // read in match table matrix
            if (!strncmp(line, "table", 5)) {
                for (i = 0; i < chem_match_tbl_labels.size(); i++) {
                    for (j = 0; j <= i; j++) {
                        fscanf(ifp, "%d", &z);
                        chem_match_tbl_matrix[i * chem_match_tbl_labels.size() +
                                              j] = z;
                        chem_match_tbl_matrix[j * chem_match_tbl_labels.size() +
                                              i] = z;
                    }
                }
            }

        }

        fclose(ifp);
    }
    /*
     * for(i=0;i<chem_match_tbl_labels.size();i++) {
     * for(j=0;j<chem_match_tbl_labels.size();j++) { cout <<
     * chem_match_tbl_matrix[i*chem_match_tbl_labels.size() + j] << " "; } cout 
     * << endl; } 
     */
}
