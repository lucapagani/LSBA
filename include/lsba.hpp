//  Copyright (C) 2014 Luca Pagani
// 
// This file is part of LSBA.
// 
// LSBA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// LSBA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with LSBA.  If not, see <http://www.gnu.org/licenses/>.

#ifndef LSBA_H__
#define LSBA_H__

#include <MatSparse.h>
#include <SmoothMatrix.h>
#include <GenMatrix.h>
#include <UCBtypedef.h>
#include <UCBsplines.h>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <cubature/cubature.h>

#define lsba_real float

using namespace Eigen;

class LSBA {

public:

  //! Default constructor
  LSBA () = default;

  /** Constructor with boost shared pointers to scattered data
   * @param u: u-values
   * @param v: v-values
   * @param z: z-values
   * @param phi: spline coefficients
   * @param smoothingFac: Smoothing factor; a weight that determines the smoothness of the surface
   */
  LSBA ( boost::shared_ptr<std::vector<double> > u,
         boost::shared_ptr<std::vector<double> > v,
         boost::shared_ptr<std::vector<double> > z,
         boost::shared_ptr<GenMatrix<lsba_real> > phi, // solution vector
         double smoothingFac = 0.0
       );

  /** Constructor with shared pointers to scattered data
   * @param u: u-values
   * @param v: v-values
   * @param z: z-values
   * @param w: weights for weighted least square
   * @param phi: spline coefficients
   * @param smoothingFac: Smoothing factor; a weight that determines the smoothness of the surface
   */
  LSBA ( boost::shared_ptr<std::vector<double> > u,
         boost::shared_ptr<std::vector<double> > v,
         boost::shared_ptr<std::vector<double> > z,
         boost::shared_ptr<std::vector<double> > w,
         boost::shared_ptr<GenMatrix<lsba_real> > phi, // solution vector
         double smoothingFac = 0.0
       );

  //! Destructor
  ~LSBA() {}

  /** Set the weights for weithed least square
   * @param weights: vector of weights
   */
  void
  SetWeights ( boost::shared_ptr<std::vector<double> > weights );

  /** Compute the spline coefficients with the Conjugate Gradient method
   * @param no_iterations: maximum number of iterations
   * @param compute_smoothing_matrix: if the smoothing matrix must be computed; if true the algorithom can bacome slow
   */
  void
  Compute ( int no_iterations,
            bool compute_smoothing_matrix = false
          );

  /** Compute the spline coefficients with the direct sparse cholesky factorization
   *
   * @param compute_smoothing_matrix: if the smoothing matrix must be computed; if true the algorithom can bacome slow
   */
  void
  ComputeDirect ( bool compute_smoothin_matrix = false );

  /** Set the domain over which the surface is to be defined.
    * The default is the xy-range of the scattered data.
    * If used, this must be done before creating the surface.
    *
    * @param umin: minimum value of the u-coordinate
    * @param vmin: minimum value of the v-coordinate
    * @param umax: maximum value of the u-coordinate
    * @param vmax: maximum value of the v-coordinate
    *
    * \note This function can only be used to expand the domain beyond the uv-range
    *       of the scattered data. It is the users responsibility to check that
    *       no scattered data falls outside the domain.
    *       (use std::min_element and std::max_element to find range of data
    *       for std::vector)
    */
  void
  SetDomain ( double umin,
              double vmin,
              double umax,
              double vmax
            );

  /** Get the domain over which the surface is defined
    *
    * @param umin: minimum value of the u-coordinate
    * @param vmin: minimum value of the v-coordinate
    * @param umax: maximum value of the u-coordinate
    * @param vmax: maximum value of the v-coordinate
    */
  void
  GetDomain ( double& umin,
              double& vmin,
              double& umax,
              double& vmax
            ) const;

  /** Number of unknowns in the equation system, i.e. B-spline coefficients
   *
   * @return number of unknowns spline coefficients
   */
  int
  NumberOfUnknowns() const {
    return x_->noX() * x_->noY(); // Which is the same as n1_*n2_
  }

  /** Adjust the thin plate spline energy.
   *  By default this value is 1.0. To get a surface that is smoother (but less accurate)
   *  than the default, a value grater than 1.0 should be given and vice versa.
   *  This will typically be done after relaxation has been run and one want a
   *  a surface that is smoother, or less smooth and more accurate, than obtained in
   *  from previous relaxation. Relaxation must then be run afterwards
   *
   * @param smoothing_fac: new value of the smoothing factor
   */
  void
  SetSmoothingFactor ( double smoothing_fac );

  //! Compute F'F and F'z
  void
  ComputeParameters ();

  //! Build sparse structure
  void
  BuildSparseStructure ();

  //! Compute sqrt(W) * F, sqrt(W) * z, F' * W * F and F'z
  void
  ComputeParametersW ();

  /** Prediction in the location (u,v)
   *
   * @param u: value of the u coordinate of the point
   * @param v: value of the v coordinate of the point
   *
   * @return value of the predicted point
   */
  double
  Predict ( double u,
            double v
          ) const;

  /** Predict the mean and the variance in the location (u,v)
   *
   * @param u: value of the u coordinate of the point
   * @param v: value of the v coordinate of the point
   * @param mean: value of the predicted mean
   * @param variance: value of the variance of the predicted point in the location (u,v)
   */
  void
  Predict ( double u,
            double v,
            double& mean,
            double& variance,
            size_t max_iterations = 0
          ) const;

  /** Predict the mean and the variance in the location (u,v)
   *
   * @param u: value of the u coordinate of the point
   * @param v: value of the v coordinate of the point
   * @param mean: value of the predicted mean
   * @param variance: value of the variance of the predicted point in the location (u,v)
   * @param max_iteration: maximum number of iteration for inverse the matrix
   */
  void
  Predict ( boost::shared_ptr< std::vector<double> > u,
            boost::shared_ptr< std::vector<double> > v,
            std::vector<double>& mean,
            std::vector<double>& variance,
            int max_iteration
          ) const;

  /** Predict the mean and the variance in the location (u,v)
   *
   * @param u: value of the u coordinate of the point
   * @param v: value of the v coordinate of the point
   * @param mean: value of the predicted mean
   * @param variance: value of the variance of the predicted point in the location (u,v)
   * @param max_iteration: maximum number of iteration for inverse the matrix
   * @param rhs_low: right-hand side term in order to compute the prediction variance of the linkage model
   */
  void
  Predict ( double u,
            double v,
            double& mean,
            double& variance,
            Matrix<lsba_real, Dynamic, 1> f0,
            Matrix<lsba_real, Dynamic, 1>& rhs_low,
            size_t max_iterations = 0
          ) const;

  /** Compute estimation of variance. It can be slow if the smoothing term is not 0.0
   */
  void
  ComputeVariance ();

  /** Compute estimation of variance as the population variance of the residuals
   */
  void
  ComputeVarianceFast ();

  /** Compute the mean variance in the uv domain [u_min, u_max] x [v_min, v_max]
   *
   */
  double
  ComputeMeanVariance ( double u_min,
                        double v_min,
                        double u_max,
                        double v_max,
                        size_t max_eval = 0,
                        double req_abs_error = 1e-4,
                        double req_rel_error = 1e-4,
                        size_t max_iterations = 0
                      ) const;

  /** Compute the mean standard deviation in the uv domain [u_min, u_max] x [v_min, v_max]
   *
   */
  double
  ComputeMeanSd ( double u_min,
                  double v_min,
                  double u_max,
                  double v_max,
                  size_t max_eval = 0,
                  double req_abs_error = 1e-4,
                  double req_rel_error = 1e-4,
                  size_t max_iterations = 0
                ) const;

  /** Set the variance of the model
   *
   * @param value of the estimated variance
   */
  void
  set_variance ( double sigma2 ) {
    sigma2_ = sigma2;
  };                
                
  /** Return the variance of the model
   *
   * @return value of the estimated variance
   */
  double
  get_variance () const {
    return sigma2_;
  };

  /** Return the spline size in u and v direction
   *
   * @param n1 spline size in direction u
   * @param n2 spline size in direction v
   */
  void
  get_spline_size ( int& n1, int& n2 ) const {
    n1 = n1_;
    n2 = n2_;
  }

  /** Return the smoothing factor
   *
   * @return value of the smoothing factor lambda
   */
  double
  get_smoothing_factor () const {
    return lambda_;
  }

//   boost::shared_ptr<SparseMatrix<lsba_real, ColMajor> >
//   get_FF_matrix () {
//     boost::shared_ptr<SparseMatrix<lsba_real, ColMajor> > FF_ptr = boost::make_shared<SparseMatrix<lsba_real, ColMajor> > ( FF_ );
//
//     return FF_ptr;
//   }


  /** Calculate the l2 norm of the current solution (Frobenius) */
//   double CurrentSolutionNorm() const;

  /** The discrete l_2 norm of the residual ||b-Ax||_2 possibly scaled by area of grid cell */
//   double Residual_l2 ( bool scaled=false ) const;
  /** The max norm of the residual (unscaled) : ||b-Ax||_inf */
//   double Residual_linf() const;

protected:

  //! Vector of the u coordinates
  boost::shared_ptr<std::vector<double> > u_;
  //! Vector of the v coordinates
  boost::shared_ptr<std::vector<double> > v_;
  //! Vector of the z coordinates
  boost::shared_ptr<std::vector<double> > z_;
  //! Vector of the weights
  boost::shared_ptr<std::vector<double> > w_;

  //! Model matrix
  SparseMatrix<lsba_real, RowMajor> F_;
  //! sqrt(W) * F_
  SparseMatrix<lsba_real, RowMajor> Fw_;
  //! The unknown B-spline coefficient matrix (corresponds to PHI_ in MBA)
  boost::shared_ptr<GenMatrix<lsba_real> > x_;
  //! Vector of unkown B-spline coefficient
  Matrix<lsba_real, Dynamic, 1> beta_;
  //! sqrt(W) * zz_
  Matrix<lsba_real, Dynamic, 1> zz_;
  //! Vector of observed points
  Matrix<lsba_real, Dynamic, 1> zzw_;
  //! Model variance
  double sigma2_;

  // Useful parameters
  //! F'F
  SparseMatrix<lsba_real, ColMajor> FF_;
  //! F'z
  Matrix<lsba_real, Dynamic, 1> Fz_;
  //! FF_ + lambda_ * E
  SparseMatrix<lsba_real, ColMajor> S_;
  //! Smoothing matrix
  SparseMatrix<lsba_real, ColMajor> Sl_;

  //! Solver for prediction variance
//   ConjugateGradient<SparseMatrix<lsba_real> > solver_;

  //! Size of spline coefficient matrix: x_.noX(), x_.noY() in A_ x_= b_
  int n1_, n2_;
  //! Order of spline surface
  int K_,   L_;

  //! Smoothing factor calculated from the l_2 norms of the matrices:
  double lambda_;
  //! Factor multiplied with lambda (default = 0, least square solution)
  double lambda_fac_;

  //! Domain of the surface
  std::vector<double> domain_;

  //! Compute the domain of the surface
  void
  CalculateDomain();

  /** Return the minimum u value
   *
   * @return minimum u value
   */
  double
  umin() const {
    return domain_[0];
  }

  /** Return the minimum v value
   *
   * @return minimum v value
   */
  double
  vmin() const {
    return domain_[1];
  }

  /** Return the maximum u value
   *
   * @return maximum u value
   */
  double
  umax() const {
    return domain_[2];
  }

  /** Return the maximum v value
   *
   * @return maximum v value
   */
  double
  vmax() const {
    return domain_[3];
  }

  //! Fill the vector beta with the initial values in row-major order (phi)
  void FillBeta ();

  //! Compute the matrix of the smoothing term (E)
  void
  ComputeSmoothingTerm ();

  /** Add a point to the model matrix (F)
   *
   * @param k: number of the row of the model matix
   * @param u: value of the u coordinate of the point
   * @param v: value of the v coordinate of the point
   */
  void
  AddPointF ( int k,
              double u,
              double v
            );

//   void
//   SetLimits();

  //! Compute lambda as ||F'F||_F / ||E||_F
  void
  ComputeLambda(); // for smoothing term

  static int
  PredictVariance ( unsigned n_dim,
                    const double* uv,
                    void* util,
                    unsigned f_dim,
                    double* f_val
                  );

  static int
  PredictSd ( unsigned n_dim,
              const double* uv,
              void* util,
              unsigned f_dim,
              double* f_val
            );  
  
  /** Compute the value of the spline function in location (u,v)
   *
   * @param u: value of the u coordinate
   * @param v: value of the v coordinate
   * @param f: vector of spline function in location (u,v)
   */
  void
  PointF ( double u,
           double v,
           SparseVector<lsba_real>& f
         ) const;

  /** Compute the value of the spline function in location (u,v)
   *
   * @param u: value of the u coordinate
   * @param v: value of the v coordinate
   * 
   * @return f: vector of spline function in location (u,v)
   */
  Matrix<lsba_real, Dynamic, 1>
  PointF ( double u,
           double v
         ) const;

  void
  AddPointFw ( int k,
               double u,
               double v
             );

};

#endif
