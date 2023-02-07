#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <chrono>  // timer
#include <sstream>
#include <fstream>  // for the file

// include helper 
#include "experimentsHelper.h"

// include types
#include "NViewsTypes.h"
// include utils
#include "NViewsUtils.h"
// triangulation solver
#include "NViewsClass.h"
// include certifier
#include "NViewsCertifier.h"


// data generation
#include "../utils/generatePointCloud.h"

// for eigendecomposition
#include <eigen3/Eigen/Dense>
#include <Eigen/Eigenvalues> 


#include <chrono>  // timer


using namespace std;
using namespace NPlanarTrian;
using namespace std::chrono;

#define SAVEPROB true

enum class methodGenPose{
        ORBITAL = 0, 
        LATERAL,
        GENERAL 
};  // end of enum class



Vector3 returnThisTranslation( double max_parallax, const Vector3 & dir_parallax)
{
        return (dir_parallax);

}; 


double distEuc(const Vector3 & X, const Vector3 & Y)
{
        double dist = 110.0; 
        dist = (X-Y).norm();
        return dist;
};

double evalConstr(const Matrix3 & H, const Vector2 & p1, const Vector2 & p2)
{
        Vector3 o1, o2; 
        o1 << p1, 1; 
        o2 << p2, 1; 

        Vector3 q = H * o1; 
        q /= q(2); 
        
        double res = (o2 - q).squaredNorm(); 
        return res;
};

int main(int argc, char** argv)
{
    std::cout << "Generic test:\nN views triangulation comparison\n"; 



    /* Read params from input */ 

    string name_in_params = "basic_params.txt"; 

    /* Read the name of the file */
    if (argc > 1)
        name_in_params = argv[1]; 


    SceneOptions options; 

    std::cout << "Generic test file !\n";
    std::cout << "Input for the test: " << name_in_params << std::endl; 


    // read params from file
    bool valid_options = readOptionsFromFile(name_in_params, options);   


    std::srand(std::time(nullptr));

    double width = 3; 
    double height = 2; 
    double n_rows = 3; 
    double n_cols = 2; 
                                   
    int size_img_i = 512;
          
    for (size_t noise_id=0; noise_id < options.n_noise; noise_id++)
    {
        double noise_i = options.noise[noise_id];

        for (size_t cam_id = 0; cam_id < options.n_arr_cams; cam_id++)
        {
                int M_cam_i = options.arr_cams[cam_id];
        
        
                for (size_t par_id = 0; par_id < options.n_parallax; par_id++)
                {
                        double par_i = options.max_parallax[par_id];    
        
        
                        for (size_t focal_id = 0; focal_id < options.n_focal; focal_id++)
                        {
                                double focal_i = options.focal_length[focal_id];
                                
                               for (size_t d_c_id = 0; d_c_id < options.n_dist_centers; d_c_id++)
                               {
                                        double dist_center_i = options.dist_centers[d_c_id]; 
                                         
                                         // This file saves all our resutls
                                         auto name_f_sol = "res/sol_noise_" + std::to_string(noise_i) 
                                                                     + "_cams_" + std::to_string((int)M_cam_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fsol(name_f_sol);
                                         
                                         auto name_f_lin_3d = "res/lin_3D_noise_" + std::to_string(noise_i) 
                                                                     + "_cams_" + std::to_string((int)M_cam_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream flin3d(name_f_lin_3d);
                                         
                                         auto name_f_3d = "res/err_3D_noise_" + std::to_string(noise_i) 
                                                                     + "_cams_" + std::to_string((int)M_cam_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream f3d(name_f_3d);
                                         
                                         auto name_f_l2 = "res/err_l2_noise_" + std::to_string(noise_i) 
                                                                     + "_cams_" + std::to_string((int)M_cam_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fl2(name_f_l2);
         
                                         
                                         
                                         auto name_f_l1 = "res/err_l1_noise_" + std::to_string(noise_i) 
                                                                     + "_cams_" + std::to_string((int)M_cam_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fl1(name_f_l1);
                                         
                                         auto name_f_linfty = "res/err_linfty_noise_" + std::to_string(noise_i) 
                                                                     + "_cams_" + std::to_string((int)M_cam_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream flinfty(name_f_linfty);
                                         
                                         auto name_f_times = "res/times_noise_" + std::to_string(noise_i) 
                                                                     + "_cams_" + std::to_string((int)M_cam_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream ftime(name_f_times);
                                         
                                         auto name_f_exp = "res/times_exp_noise_" + std::to_string(noise_i) 
                                                                     + "_cams_" + std::to_string((int)M_cam_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream ft_exp(name_f_exp);
                                         
                                         auto name_diff = "res/diff_noise_" + std::to_string(noise_i) 
                                                                     + "_cams_" + std::to_string((int)M_cam_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fdiff(name_diff);
                                         
                                         
                                         auto name_cert = "res/cert_noise_" + std::to_string(noise_i) 
                                                                     + "_cams_" + std::to_string((int)M_cam_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fcert(name_cert);
                                         
                                         
                                         #if SAVEPROB
                                         
                                                 auto name_prob = "res/prob_" + std::to_string(noise_i) 
                                                                             + "_cams_" + std::to_string((int)M_cam_i) 
                                                                             + "_par_" + std::to_string(par_i) 
                                                                             + "_focal_" + std::to_string((int)focal_i) 
                                                                             + "_center_" + std::to_string(dist_center_i) 
                                                                             + ".txt";
                                                 std::ofstream fprob(name_prob);
                                         
                                         #endif
                                         auto name_constr = "res/constr_" + std::to_string(noise_i) 
                                                                     + "_cams_" + std::to_string((int)M_cam_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fconstr(name_constr);
                                         
                                                                                
                                         
                                       for (size_t n_iter = 0; n_iter < options.max_iter; n_iter++)
                                       {
                                         
                                               // define struct with params
                                               PCParams str_in = PCParams(); 
                                               
                                               str_in.focal_length = focal_i; 
                                               str_in.noise = noise_i; 
                                               str_in.size_img = size_img_i;                                                  
                                               str_in.width = width; 
                                               str_in.height = height; 
                                               str_in.n_rows = n_rows; 
                                               str_in.n_cols = n_cols; 
                                               str_in.d_plane = (double) dist_center_i; 
                                               str_in.max_parallax = par_i;
                                               str_in.max_angle = options.max_rotation; 
                                               str_in.M_cameras = M_cam_i;

                                               // Select pose generation 
                                               // param: options.method_trans in {1, 2, 3}
                                               methodGenPose m_gen_pose = static_cast<methodGenPose>(options.method_trans);

                                               
                                               // generate problem
                                               PCRes str_out = PCRes(); 
                                               std::cout << "Selection method for pose generation\n"; 
                                               
                                               switch (m_gen_pose)
                                               {
                                                case methodGenPose::ORBITAL:
                                                {
                                                        std::cout << "[ORBITAL CAMERA]\n";
                                                        str_in.max_parallax = 5;  // radius circle
                                                        
                                                        Vector3 dist_vector; 
                                                        dist_vector << dist_center_i, 0, 0; 
                                                        str_in.dir_parallax = dist_vector;
                                                        
                                                        str_in.max_angle = dist_center_i; 
                                                        // generateOrbitalRotation takes the distance to center and the translation vector
                                                        
                                                        // In this configuration, rotation and translation depends on dist_center_i 
                                                        // and radius = max_parallax = 5 units
                                                        str_in.noise_rot = options.noise_rot; 
                                                        str_in.noise_trans = options.noise_trans;
                                                        
                                                        str_out = generatePointCloud(str_in, generateOrbitalTranslation, generateOrbitalRotation); 
                                                }
                                                break;
                                                
                                                case methodGenPose::LATERAL: 
                                                        std::cout << "[LATERAL CAMERA]\n";
                                                        // In this configuration, 
                                                        // translation has form X, 0 0
                                                        // and rotation is identity
                                                        
                                                        str_in.max_angle = 0.0;                 
                                                        str_in.max_parallax = par_i;                                       
                                                        str_in.noise_trans = options.noise_trans; 
                                                        str_in.noise_rot = options.noise_rot; 
                                                        str_out = generatePointCloud(str_in, generateTranslationStereo);                                                 
                                                break;
                                                
                                                /*
                                                case methodGenPose::GENERAL: 
                                                        std::cout << "[GENERAL CAMERA]\n";
                                                        str_out = generatePointCloud(str_in); 
                                                break;
                                                */
        
                                                default:
                                                        std::cout << "[GENERAL CAMERA]\n";
                                                        str_in.max_angle = options.max_rotation;                 
                                                        str_in.max_parallax = par_i;                                       
                                                        str_in.noise_trans = options.noise_trans; 
                                                        str_in.noise_rot = options.noise_rot; 
                                                        str_out = generatePointCloud(str_in);                                      
                                                break;
                                               
                                               }
                                               
                                               Matrix3 K = Matrix3::Identity(); 
                                               K(0,0) = focal_i; 
                                               K(1,1) = focal_i; 
                                               K(0,2) = (double) size_img_i; 
                                               K(1,2) = (double) size_img_i;
                                               Matrix3 iK = K.inverse();

                                                
                                                // extract data
                                                // generate full graph (M 2) combinations
                
                                                Eigen::MatrixXd idx_matrix;
                                                int n_comb = -1;
                                                if (options.use_linear)
                                                        n_comb = generateLinGraph(M_cam_i, idx_matrix);
                                                else
                                                        n_comb = generateM2Comb(M_cam_i, idx_matrix);
                                                                                          
                                                // generate correspondences
                                                std::vector<PairObj> set_corr; 
                                                set_corr.empty();
                                                for (int jj=0; jj<n_comb; jj++)
                                                {
                                                        
                                                        int id1 = idx_matrix(0, jj); 
                                                        int id2 = idx_matrix(1, jj); 
                                                        Matrix3 R1 = str_out.set_rot[id1]; 
                                                        Matrix3 R2 = str_out.set_rot[id2]; 
                                                        Vector3 t1 = str_out.set_trans[id1]; 
                                                        Vector3 t2 = str_out.set_trans[id2]; 
                                                        Matrix4 P1 = Matrix4::Identity(); 
                                                        Matrix4 P2 = Matrix4::Identity(); 
                                                        P1.block<3,3>(0,0) = R1; 
                                                        P1.block<3,1>(0,3) = t1;
                                                        P2.block<3,3>(0,0) = R2; 
                                                        P2.block<3,1>(0,3) = t2;
                                                        
                                                        Matrix4 Prel = P2 * P1.inverse(); 
                                                        Matrix3 Rrel = Prel.block<3,3>(0,0); 
                                                        Vector3 trel = Prel.block<3,1>(0,3);                 
                                                        
                                                                                                          
                                                        
                                                        Vector3 n; 
                                                        n << 0, 0, 1; 
                                                        Vector3 m = R1 * n; 
                                                        double d1 = dist_center_i + m.dot(t1); 
                                                        Matrix3 Hi = Rrel + trel * m.transpose() / d1; 

                                                        PairObj corr_i; 
                                                        corr_i.id1 = id1; 
                                                        corr_i.id2 = id2; 
                                                        corr_i.H = Hi;                 
                                                        corr_i.p1 = iK * str_out.obs[0].col(id1); 
                                                        corr_i.p2 = iK * str_out.obs[0].col(id2);     
                                                        set_corr.push_back(corr_i);       
                                                                                             
                                                }
                                                
                                                
                                                
                                                // Save pose 
                                                #if SAVEPROB
                                                {
                                                        for (int jj=0; jj < M_cam_i; jj++)
                                                        {
                                                                Matrix3 R = str_out.set_rot[jj]; 
                                                                Vector3 t = str_out.set_trans[jj]; 
                                                                fprob << t(0) << "," << t(1) << "," << t(2) << ",";
                                                                
                                                                for (int id_r = 0; id_r < 3; id_r++)
                                                                {
                                                                        for (int id_c=0; id_c < 3; id_c++)
                                                                                fprob << R(id_r, id_c) << ",";
                                                                } 
                                                                fprob << std::endl;
                                                        }
                                                }
                                                #endif
                                              
                                               // for each point
                                               int idx = 0; 
                                               int n_points = str_out.points_3D.cols(); 
                                               for (idx = 0; idx < n_points; idx++)
                                               {
                                                
                                                        // 1. Loop through cameras
                                                        for (int jj=0; jj<n_comb; jj++)
                                                        {
                                                                int id1 = idx_matrix(0, jj); 
                                                                int id2 = idx_matrix(1, jj); 
                                                                set_corr[jj].p1 = iK * str_out.obs[idx].col(id1); 
                                                                set_corr[jj].p2 = iK * str_out.obs[idx].col(id2);   
                                                             
                                                        }
                                                        
                                                        Vector3 X = str_out.points_3D.col(idx);
                                                        
                                                        // run correction method
                                                        NViewsClass corr_N_view; 
                                                        // 1. Create constraint matrices
                                                        corr_N_view.createProblemMatrices(set_corr, M_cam_i); 
                                                       
                                                        // 2. Run correction
                                                        NViewsOptions options_corr;                                                         
                                                        options_corr.save_val_constr = false;
                                                        options_corr.debug = false; 
                                                        
                                                        auto start_t_ours = high_resolution_clock::now();
                                                        NViewsResult res_corr = corr_N_view.correctObservations(options_corr);
                                                        auto time_ours = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_ours);
                                                        
                                                        
                                                        // check optimality of the solution
                                                        std::vector<Constr2View> constr_red_1, constr_red_2;
                                                        
                                                        corr_N_view.getConstrRed(constr_red_1, constr_red_2);
                                                        bool debug_cert = false; 
                                                        NViewCertClass cert_obj(M_cam_i, constr_red_1, constr_red_2, debug_cert); 
                                                                                                        
                                                        // call function
                                                        auto start_t_cert = high_resolution_clock::now();
                                                        NCertRes res_cert = cert_obj.checkOptimality(res_corr.sol_final); 
                                                        auto time_cert = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_cert);
                                                        
                                                                      
                                        
                                                        // c. Reconstruct 3D point by linear method
                                                    
                                                        // Matrix projections
                                                        std::vector<Matrix4> proj_s; 
                                                        proj_s.empty();
                                                        std::vector<Vector3> obs_s, obs_init, obs_ref; 
                                                        obs_s.empty(); 
                                                        obs_init.empty(); 
                                                        obs_ref.empty(); 
                                                        
                                                        for (int jc=0; jc<M_cam_i;jc++)
                                                        {
                                                                // matrix projection for this camera
                                                                Matrix3 R = str_out.set_rot[jc]; 
                                                                Vector3 t = str_out.set_trans[jc]; 
                                                                Matrix4 P1 = Matrix4::Identity(); 
                                                                P1.block<3,3>(0,0) = R; 
                                                                P1.block<3,1>(0,3) = t;
                                                                proj_s.push_back(P1); 
                                                                
                                                                // observations
                                                                
                                                                Vector3 pt = iK * str_out.obs[idx].col(jc); 
                                                                obs_s.push_back(pt);
                                                                
                                                                // update observation init
                                                                Vector3 delta_init; 
                                                                delta_init << res_corr.sol_init( jc*2), res_corr.sol_init(jc*2 + 1), 0;
                                                                obs_init.push_back(pt + delta_init);
                                                                
                                                                // update observation refinenement  
                                                                Vector3 delta_ref; 
                                                                delta_ref << res_corr.sol_final( jc*2), res_corr.sol_final(jc*2 + 1), 0; 
                                                                obs_ref.push_back(pt + delta_ref);
                                                        }
                                                        
                                                                                                          
                                                        // triangulate point
                                                        Vector3 P_lin; 
                                                        Eigen::VectorXd depths_lin; 
                                                        double error_lin = triangulateNPoint(proj_s, obs_s, P_lin, depths_lin);
                                                        
                                                        Vector3 P_init; 
                                                        Eigen::VectorXd depths_init; 
                                                        double error_init = triangulateNPoint(proj_s, obs_init, P_init, depths_init);
                                                        
                                                        Vector3 P_ref; 
                                                        Eigen::VectorXd depths_ref; 
                                                        double error_ref = triangulateNPoint(proj_s, obs_ref, P_ref, depths_ref);
                                                        
                                                                                                                                                                
                                                
                                                        
                                                        /* Save results */ 
                                                        // solution of problem 
                                                        fsol << P_lin(0) << "," << P_lin(1) << "," << P_lin(2) << "," ; 
                                                        fsol << P_init(0) << "," << P_init(1) << "," << P_init(2) << "," ; 
                                                        fsol << P_ref(0) << "," << P_ref(1) << "," << P_ref(2) << "," ; 
                                                        fsol << std::endl; 
                                                        
                                                        // Euclidean distance 3D point
                                                        f3d << (X-P_lin).norm() << ",";                 
                                                        f3d << (X-P_init).norm() << ",";       
                                                        f3d << (X-P_ref).norm(); 
                                                        f3d << std::endl;
                                                        
                                                        // Error for the linear method
                                                        flin3d << error_lin << ","; 
                                                        flin3d << error_init << ",";
                                                        flin3d << error_ref << ","; 
                                                        flin3d << std::endl;                                                         
                                                        
                                                        // l2 norm for observations 
                                                        fl2 << res_corr.sol_init.squaredNorm()   << ","; 
                                                        fl2 << res_corr.sol_final.squaredNorm()   << ","; 
                                                        fl2 << std::endl;
                                                        
                                                        // L1 norm for observations 
                                                        fl1 << res_corr.sol_init.lpNorm<1>() << ","; 
                                                        fl1 << res_corr.sol_final.lpNorm<1>()<< ","; 
                                                        fl1 << std::endl;
                                                        
                                                        // Linfty norm for observations 
                                                        flinfty << res_corr.sol_init.lpNorm<Eigen::Infinity>() << ","; 
                                                        flinfty << res_corr.sol_final.lpNorm<Eigen::Infinity>() << ",";
                                                        flinfty  << std::endl;
                                                        
                                                        
                                                        ftime << res_corr.time_init << ",";
                                                        ftime << (double) time_ours.count() << ",";  
                                                        ftime << M_cam_i << ","; 
                                                        ftime << (double) time_cert.count() << ",";
                                                        ftime << std::endl;        
                                                        
                                                        ft_exp <<  res_corr.time_init << ","; 
                                                        ft_exp <<  res_corr.time_ref << ","; 
                                                        ft_exp <<  res_corr.n_iters << ","; 
                                                        ft_exp <<  res_cert.time_mult << ","; 
                                                        ft_exp <<  res_cert.time_hess << ",";
                                                        ft_exp << std::endl;        
                                                            
                                                                                                                                                                       
                                                                                                                
                                                        /* Certifier */
                                                        fcert << res_cert.min_eig << std::endl;
                                                        
                                                        /* epipolar constraint */   
                                                        // compute and save epipolar error 
                                                        fconstr << res_corr.max_constr_init << ","; 
                                                        fconstr << res_corr.max_constr << ","; 
                                                        fconstr << res_corr.tot_constr_init << ","; 
                                                        fconstr << res_corr.tot_constr << ",";                                                         
                                                        fconstr << res_corr.sq_constr_init << ",";
                                                        fconstr << res_corr.sq_constr << ",";
                                                        fconstr << res_corr.error_lin << ",";
                                                        fconstr << std::endl;  
                                                        
                                                        
                                                        
                                                }  // end of each point             
                                                } // end of each iteration 
                                                fsol.close();    // solution
                                                flin3d.close();  // error linear method
                                                f3d.close();     // error 3d point
                                                fl2.close();     // error l2 obs
                                                fl1.close();     // error l1 obs
                                                flinfty.close(); // error linfty obs
                                                ftime.close();   // time 
                                                ft_exp.close();   // time dual diagonal
                                                
                                                #if SAVEPROB
                                                        fprob.close();  // close the problem file if opened
                                                #endif
                                                fcert.close();   // min eigenvalue hessian
                                                fconstr.close(); // epipolar constraint
                                                
                                }  // end for dist center
      
                         }  // enf for focal
       
                 }  // end for parallax
      
        }  // end for M_cameras
      
      }  // end for noise              
       
  return 0;

}  // end of main fcn
