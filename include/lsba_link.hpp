//  Copyright (C) 2014-2015 Luca Pagani
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

#ifndef LSBA_LINK__
#define LSBA_LINK__

#include "lsba.hpp"

using std::vector;
using namespace Eigen;

class LSBA_EXPORT LSBALink: public LSBA {

public:

  //! Default constructor
  LSBALink () = default;

  /*! Constructor with dimension 1
   *
   * @param u: vector of u values
   * @param v: vector of v values
   * @param z: vector of z values
   * @param phi: spline coefficients
   * @param smoothing_fac: Smoothing factor; a weight that determines the smoothness of the surface
   * @param lsba_low: lsba object of the low fidelity model
   * @param compute_weights: true if the weights must be computed as the inverse of the variance of prediction
   * @param build_structure: true if the matrix F and the system terms must be computed
   */
  LSBALink ( boost::shared_ptr< vector<double> > u,
             boost::shared_ptr< vector<double> > v,
             boost::shared_ptr< vector<double> > z,
             boost::shared_ptr<GenMatrix<lsba_real> > phi, // solution vector
             boost::shared_ptr<LSBA> lsba_low,
             double smoothing_fac = 0.0,
             bool compute_weights = true,
             bool build_structure = true
           );

  //! Desctructor
  ~LSBALink () {};

  /** Build the matrix F and the system terms
   */
  void
  BuildStructure () {
   LSBA::BuildSparseStructure ();
  };

  /** Set the domain over which the surface is to be defined.
    * The default is the xy-range of the scattered data.
    * If used, this must be done before creating the surface.
    *
    * @param umin: minimum value of the u-coordinate
    * @param vmin: minimum value of the v-coordinate
    * @param umax: maximum value of the u-coordinate
    * @param vmax: maximum value of the v-coordinate
    * @param build_structure: true if the matrix F and the system terms must be computed
    *
    * \note This function can only be used to expand the domain beyond the uv-range
    *       of the scattered data. It is the users responsibility to check that
    *       no scattered data falls outside the domain.
    *       (use std::min_element and std::max_element to find range of data
    *       for std::vector)
    */
  void
  set_domain ( double umin,
               double vmin,
               double umax,
               double vmax,
               bool build_structure = true
             )
  {
    LSBA::set_domain ( umin, vmin, umax, vmax, build_structure );
  };

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
            int max_iterations = 0
          ) const;

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

  /** Set the variance of the low fidelity points
   *
   * @param var_low variace of the low fidelity points
   */
  void
  set_variance_low ( boost::shared_ptr< vector<double> > var_low ) {
    if ( var_low->size() != z_->size() ) {
      std::cerr << "The length of the vector must be equal to the length of the points." << std::endl;
      exit ( EXIT_FAILURE );
    }
    var_low_ = var_low;
  };

  /** Set the variance of the model
   *
   * @param value of the estimated variance
   */
  void
  set_variance ( double sigma2 ) {
    sigma2_ = sigma2;
  };

private:

  //! Compute residual between Hi-Fi points and Lo-Fi prediction
  void
  ComputeResidual ();

  //! Vector of z points
  boost::shared_ptr< vector<double> > z_;
  //! Vector of z residuals points
  boost::shared_ptr< vector<double> > z_res_;
  //! Vector of weights
  boost::shared_ptr< vector<double> > weights_;

  //! Low fidelity lsba object
  boost::shared_ptr<LSBA> lsba_low_;

  //! true if the weights must be computed as the inverse of the variance of prediction
  bool compute_weights_;

  //! Variance of the low fidelity points
  boost::shared_ptr< vector<double> > var_low_;

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

  /** Add a point to the model matrix (F1)
   *
   * @param F: model matrix
   * @param k: number of the row of the model matix
   * @param u: value of the u coordinate of the point
   * @param v: value of the v coordinate of the point
   */
  void
  AddPointF0 ( SparseMatrix<lsba_real, RowMajor>& F,
               int k,
               double u,
               double v
             ) const;

};
#endif
