#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>


// include types
#include "NViewsTypes.h"
// include utils
#include "NViewsUtils.h"
// triangulation solver
#include "NViewsClass.h"
// certifier 
#include "NViewsCertifier.h"

// data generation
#include "../utils/generatePointCloud.h"



// for eigendecomposition
#include <eigen3/Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include <chrono>  // timer

using namespace std::chrono;


using namespace std;
using namespace NPlanarTrian;


size_t generateRandomIndex(size_t N)
{
        size_t n = N; 
        n = (((double)std::rand() - 0.5) / (double) RAND_MAX) * 2.0 * N;   
        return n; 
}

double distEuc(const Vector3 & X, const Vector3 & Y)
{
        double dist = 110.0; 
        dist = (X-Y).norm();
        return dist;
} 

double evalConstr(const Matrix3 & H, const Vector2 & p1, const Vector2 & p2)
{
        Vector3 o1, o2; 
        o1 << p1, 1; 
        o2 << p2, 1; 

        Vector3 Hp1 = H * o2;
        Hp1 /= Hp1(2); 
        
        double res = (Hp1 - o1).squaredNorm(); 
        return res;
}

int main(int argc, char** argv)
{

        std::cout << "Example N view triangulation\n"; 
    
        
        // parameters for estimation
        double noise = 1.0;
        double max_parallax = 2.0;  // in meters
        double focal_length = 512; 
        size_t size_img = 512;
       
        double width = 3; 
        double height = 2; 
        double n_rows = 1; 
        double n_cols = 1; 
        double d_plane = 5.0;
        
        double max_rot = 0.5;
        int M_cameras = 4;
          
        std::srand(std::time(nullptr));

                                       
        // generate problem
        PCRes str_out(n_cols * n_rows);
        PCParams str_in; 
        str_in.noise = noise; 
        str_in.focal_length = focal_length; 
        str_in.size_img = size_img; 
        str_in.M_cameras = M_cameras;
        str_in.width = width; 
        str_in.height = height; 
        str_in.d_plane = d_plane; 
        str_in.n_rows = n_rows; 
        str_in.n_cols = n_cols; 
        
        str_out = generatePointCloud(str_in); //, UtilsTwoView::generateTranslationStereo);    
       
        Matrix3 K = Matrix3::Identity(); 
        K(0,0) = focal_length; 
        K(1,1) = focal_length; 
        K(0,2) = (double) size_img; 
        K(1,2) = (double) size_img;
        Matrix3 iK = K.inverse();

        // generate full graph (M 2) combinations
        std::cout << "Number cameras: " << M_cameras << std::endl;
        std::cout << "Number points: " << str_out.points_3D.cols() << std::endl;
        Eigen::MatrixXd idx_matrix;
        int n_comb = generateLinGraph(M_cameras, idx_matrix);
        // int n_comb = generateM2Comb(M_cameras, idx_matrix);
        // std::cout << "Graph matrix:\n" << idx_matrix << std::endl;
        
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
                double d1 = d_plane + m.dot(t1); 
                
                
                Matrix3 Hi = Rrel + trel * m.transpose() / d1; 
                // check constraint
                Vector3 q2 = Hi * iK * str_out.obs[0].col(id2); 
                q2 /= q2(2); 
                
                
                q2 = Hi * iK * str_out.obs[0].col(id1); 
                q2 /= q2(2); 
                
                /* std::cout << "Index i: " << id1 << " and index j: " << id2 << std::endl; 
                std::cout << "p2:\n" << iK * str_out.obs[0].col(id2) << std::endl; 
                std::cout << "H * p1:\n" << q2 << std::endl; 
                */
                
                
                
                PairObj corr_i; 
                corr_i.id1 = id1; 
                corr_i.id2 = id2; 
                corr_i.H = Hi;                 
                corr_i.p1 = iK * str_out.obs[0].col(id1); 
                corr_i.p2 = iK * str_out.obs[0].col(id2);     
                set_corr.push_back(corr_i);        
                
                        
        }
        
        // run correction method
        NViewsClass corr_N_view; 
        // std::cout << "Creating matrices\n";
        // 1. Create constraint matrices
        corr_N_view.createProblemMatrices(set_corr, M_cameras); 
       
       
        // 2. Run correction
        // std::cout << "Setting options\n";
        NViewsOptions options_corr; 
        options_corr.save_val_constr = false;
        options_corr.debug = false; 
        // std::cout << "Correcting points\n";
        NViewsResult res_corr = corr_N_view.correctObservations(options_corr);
        
        // Show results
        corr_N_view.printResult(res_corr); 
        
        
        // check optimality of the solution
        // std::vector<Eigen::MatrixXd> constr; 
        std::vector<Constr2View> constr_red_1, constr_red_2;
        // corr_N_view.getConstr(constr); 
        corr_N_view.getConstrRed(constr_red_1, constr_red_2);
        bool debug_cert = false; 
        NViewCertClass cert_obj(M_cameras, constr_red_1, constr_red_2, debug_cert); 
        
        // call function
        // std::cout << "Calling certifier\n\n"; 
        NCertRes res_cert = cert_obj.checkOptimality(res_corr.sol_final); 
        // std::cout << "Minimum eigenvalue Hessian: " << res_cert.min_eig << std::endl; 
                
        cert_obj.printResult(res_cert);
        // std::cout << "Hess:\n" << res_cert.Hess << std::endl;

    
    
        // Matrix projections
        std::vector<Matrix4> proj_s; 
        proj_s.empty();
        std::vector<Vector3> obs_s, obs_init, obs_ref; 
        obs_s.empty(); 
        obs_init.empty(); 
        obs_ref.empty(); 
        
        for (int jc=0; jc<M_cameras;jc++)
        {
                // matrix projection for this camera
                Matrix3 R = str_out.set_rot[jc]; 
                Vector3 t = str_out.set_trans[jc]; 
                Matrix4 P1 = Matrix4::Identity(); 
                P1.block<3,3>(0,0) = R; 
                P1.block<3,1>(0,3) = t;
                proj_s.push_back(P1); 
                
                // observations
                
                Vector3 pt = iK * str_out.obs[0].col(jc); 
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
                                      
                                   
        
        std::cout << "\n\n\n\nChecking solutions\n";                              
                                     
                                                
        std::cout << "Number constraints: " << n_comb << std::endl;
        std::cout << "Error linear: " << error_lin << std::endl; 
        std::cout << "error init: " << error_init << std::endl; 
        std::cout << "Error final: " << error_ref << std::endl;

        
        std::cout << "P3d:\n" << str_out.points_3D.col(0) << std::endl;
        std::cout << "P linear:\n" << P_lin << std::endl;
        std::cout << "P init:\n" << P_init << std::endl;
        std::cout << "P ref:\n" << P_ref << std::endl;
        
       
                                
        return 0;
    

}  // end of main fcn
