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

#include "lsba_link.hpp"
#include <boost/make_shared.hpp>

/*! Constructor with dimension 1
 *
 * @param u: vector of u values
 * @param v: vector of v values
 * @param z: vector of z values
 * @param phi: spline coefficients
 * @param smoothingFac: Smoothing factor; a weight that determines the smoothness of the surface
 * @param lsba_low: lsba object of the low fidelity model
 */
LSBALink::LSBALink ( boost::shared_ptr< vector<double> > u,
                     boost::shared_ptr< vector<double> > v,
                     boost::shared_ptr< vector<double> > z,
                     boost::shared_ptr<GenMatrix<lsba_real> > phi, // solution vector
                     boost::shared_ptr<LSBA> lsba_low,
                     double smoothingFac,
                     bool compute_weights
                   ) : z_ ( z ), lsba_low_ ( lsba_low ), compute_weights_ ( compute_weights )
{
  u_ = u;
  v_ = v;

  z_res_.reset ( new vector<double> );

  ComputeResidual ();

  LSBA::z_ = z_res_;

  x_ = phi;
  lambda_fac_ = smoothingFac;

  n1_ = x_->noX();
  n2_ = x_->noY(); // This is also done when building equation system
  K_ = L_ = 4;

  CalculateDomain();

  // Resize model matrix
  F_.resize ( u_->size (), n1_ * n2_ );
  F_.reserve ( u_->size () * 16 );
  // Fill model matrix
  for ( size_t i = 0; i < u_->size (); ++i )
    AddPointF ( i, ( *u_ ) [i], ( *v_ ) [i] );
  F_.finalize();
  // Resize vector of unknown
  beta_.resize ( n1_ * n2_ );
  // Resize zz_
  zz_.resize ( z_->size () );
  // Fill zz_
  for ( size_t i = 0; i < z_->size(); ++i )
    zz_ ( i ) = static_cast<lsba_real> ( ( *z_res_ ) [i] );

  // Compute parameters
  ComputeParameters ();

  if ( lambda_fac_ > 0. ) {
    ComputeLambda ();
    lambda_ *= lambda_fac_;
    BuildSparseStructure();
    ComputeSmoothingTerm ();
  } else {
    lambda_ = 0.;
  }

  if ( compute_weights_ == true ) {
    LSBA::set_weights ( weights_ );
  }

}

void
LSBALink::set_domain ( double umin,
                      double vmin,
                      double umax,
                      double vmax
                    )
{
  if ( domain_.size() != 4 )
    domain_.resize ( 4 );

  domain_[0] = umin;
  domain_[1] = vmin;
  domain_[2] = umax;
  domain_[3] = vmax;

  // Resize model matrix
  F_.resize ( u_->size (), n1_ * n2_ );
  F_.reserve ( u_->size () * 16 );
  // Fill model matrix
  for ( size_t i = 0; i < u_->size (); ++i )
    AddPointF ( i, ( *u_ ) [i], ( *v_ ) [i] );
  F_.finalize();

  // Compute parameters
  ComputeParameters ();

  if ( lambda_fac_ > 0. ) {
    ComputeLambda ();
    lambda_ *= lambda_fac_;
    BuildSparseStructure();
    ComputeSmoothingTerm ();
  } else {
    lambda_ = 0.;
  }

  if ( compute_weights_ == true ) {
    LSBA::set_weights ( weights_ );
  }
}

//! Compute residual between Hi-Fi points and Lo-Fi prediction
void
LSBALink::ComputeResidual ()
{

  double mean ( 0. ), var ( 0. );
  z_res_->resize ( u_->size () );

  lsba_low_->ComputeVarianceFast();

  if ( compute_weights_ == true ) {

    weights_.reset ( new vector<double> );
    weights_->resize ( u_->size() );
    #pragma omp parallel for private ( mean, var )
    for ( size_t i = 0; i < u_->size (); ++i ) {
      lsba_low_->Predict ( ( *u_ ) [i], ( *v_ ) [i], mean, var, 1e2 );
      ( *z_res_ ) [i] = ( *z_ ) [i] - mean;
      ( *weights_ ) [i] = 1. / var;
    }

  } else {

    #pragma omp parallel for private ( mean, var )
    for ( size_t i = 0; i < u_->size (); ++i ) {
      mean = lsba_low_->Predict ( ( *u_ ) [i], ( *v_ ) [i] );
      ( *z_res_ ) [i] = ( *z_ ) [i] - mean;
//       ( *weights_ ) [i] = 1.;
    }

  }

}


/** Prediction in the location (u,v)
 *
 * @param u: value of the u coordinate of the point
 * @param v: value of the v coordinate of the point
 *
 * @return value of the predicted point
 */
double
LSBALink::Predict ( double u,
                    double v
                  ) const
{
  double value = LSBA::Predict ( u, v ) + lsba_low_->Predict ( u, v );
  return value;
}

/** Predict the mean and the variance in the location (u,v)
 *
 * @param u: value of the u coordinate of the point
 * @param v: value of the v coordinate of the point
 * @param mean: value of the predicted mean
 * @param variance: value of the variance of the predicted point in the location (u,v)
 */
void
LSBALink::Predict ( double u,
                    double v,
                    double& mean,
                    double& variance,
                    int max_iterations
                  ) const
{
  mean = variance = 0.;

  // Mean and variance of low-fidelity model
  double mean_low, variance_low;
  SparseVector<lsba_real> f;
  Matrix<lsba_real, Dynamic, 1> fh;

  PointF ( u, v, f );

  fh.resize ( f.rows() );

  for ( size_t i = 0; i < fh.rows(); ++i )
    fh ( i ) = f.coeff ( i );

//   for ( SparseVector<lsba_real>::InnerIterator it ( f ); it; ++it )
//     fh ( it.row () ) = it.value();

  int n1, n2;

  lsba_low_->get_spline_size ( n1, n2 );

  // Compute F0 matrix
  SparseMatrix<lsba_real, RowMajor> F0;
  F0.resize ( u_->size (), n1 * n2 );
  F0.reserve ( u_->size () * 16 );
  // Fill model matrix
  for ( size_t i = 0; i < u_->size (); ++i )
    AddPointF0 ( F0, i, ( *u_ ) [i], ( *v_ ) [i] );
  F0.makeCompressed();
  F0.finalize();

  double sigma2_low = lsba_low_->get_variance();

  Matrix<lsba_real, Dynamic, 1> ff;
  ConjugateGradient<SparseMatrix<lsba_real> > solver;

  if ( max_iterations > 0 )
    solver.setMaxIterations ( max_iterations );

  if ( lambda_ == 0. ) {
    solver.compute ( FF_ );
  } else {
    solver.compute ( S_ );
  }

  ff = solver.solve ( fh );

  if ( w_ == NULL && var_low_ == NULL ) {
    lsba_low_->Predict ( u, v, mean_low, variance_low );
    LSBA::Predict ( u, v, mean, variance );

    mean += mean_low;
    variance += variance_low;

  } else {
    Matrix<lsba_real, Dynamic, 1> rhs_cov;

    ff = F_ * ff ;

    if ( w_ == NULL && var_low_ != NULL ) {
      for ( size_t i = 0; i < ff.rows(); ++i )
        variance += ff ( i ) * ff ( i ) * ( ( *var_low_ ) [i] + sigma2_ );
    } else {
      for ( size_t i = 0; i < ff.rows(); ++i )
        variance += ff ( i ) * ff ( i ) * ( 1. + sigma2_ ) * ( *w_ ) [i];
    }

//     std::cout << "aaa" << std::endl;
//     Matrix<lsba_real, Dynamic, Dynamic> F0f = Matrix<lsba_real, Dynamic, Dynamic> ( F0 );

//     std::cout << "aaa" << std::endl;    
    Matrix<lsba_real, Dynamic, 1> f0 = F0.transpose() * ff;

    lsba_low_->Predict ( u, v, mean_low, variance_low, f0, rhs_cov, max_iterations );

    variance += variance_low + 2 * sigma2_low * f0.dot ( rhs_cov );

    mean = fh.dot ( beta_ ) + mean_low;
  }
//   else {
//     Matrix<lsba_real, Dynamic, 1> rhs_cov;
// //         SparseMatrix<lsba_real, RowMajor> Fww;
// //         Fww.resize ( Fw_.rows (), Fw_.cols () );
// //         Fww.reserve ( Fw_.nonZeros() );
// //
// //         for ( int k = 0; k < Fw_.outerSize(); ++k ) {
// //             Fww.startVec ( k );
// //             for ( SparseMatrix<lsba_real, RowMajor>::InnerIterator it ( F_, k ); it; ++it ) {
// //                 Fww.insertBack ( it.row (), it.col () ) = it.value () * ( 1. + sigma2_ ) * ( *w_ ) [it.row ()];
// //             }
// //         }
// //         Fww.finalize();
// 
//     ff = F_ * ff ;
// 
//     for ( size_t i = 0; i < ff.rows(); ++i )
//       variance += ff ( i ) * ff ( i ) * ( 1. + sigma2_ ) * ( *w_ ) [i];
// 
//     Matrix<lsba_real, Dynamic, Dynamic> F0f = Matrix<lsba_real, Dynamic, Dynamic> ( F0 );
// 
// //         Matrix<lsba_real, Dynamic, 1> f0 = F0f.transpose() * ( F_ * ff );
//     Matrix<lsba_real, Dynamic, 1> f0 = F0f.transpose() * ff;
// 
//     lsba_low_->Predict ( u, v, mean_low, variance_low, f0, rhs_cov, max_iterations );
// 
// //         Matrix<lsba_real, Dynamic, 1> temp = F_.transpose() * ( Fww * ff );
// //         variance = ff.dot ( temp ) + variance_low;
// 
//     variance += variance_low + 2 * sigma2_low * f0.dot ( rhs_cov );
// 
//     mean = fh.dot ( beta_ ) + mean_low;
// 
// //     std::cout << mean << " " << variance << std::endl;
//   }
}

/** Compute the mean variance in the uv domain [u_min, u_max] x [v_min, v_max]
 *
 */
double
LSBALink::ComputeMeanVariance ( double u_min,
                                double v_min,
                                double u_max,
                                double v_max,
                                size_t max_eval,
                                double req_abs_error,
                                double req_rel_error,
                                size_t max_iterations
                              ) const
{

  double range_min[]= {u_min, v_min};
  double range_max[]= {u_max, v_max};

  double value, err;

  std::pair<const LSBALink*, size_t>* util ( new std::pair<const LSBALink*, size_t> );
  *util = std::make_pair ( this, max_iterations );

  hcubature ( 1, PredictVariance, util, 2, range_min, range_max, max_eval, req_abs_error, req_rel_error, ERROR_INDIVIDUAL, &value, &err );

//   delete util->second;
//   delete util;

  return value / ( u_max - u_min ) / ( v_max - v_min );
}

/** Compute the mean variance in the uv domain [u_min, u_max] x [v_min, v_max]
 *
 */
double
LSBALink::ComputeMeanSd ( double u_min,
                          double v_min,
                          double u_max,
                          double v_max,
                          size_t max_eval,
                          double req_abs_error,
                          double req_rel_error,
                          size_t max_iterations
                        ) const
{

  double range_min[]= {u_min, v_min};
  double range_max[]= {u_max, v_max};

  double value, err;

  std::pair<const LSBALink*, size_t>* util ( new std::pair<const LSBALink*, size_t> );
  *util = std::make_pair ( this, max_iterations );

  hcubature ( 1, PredictSd, util, 2, range_min, range_max, max_eval, req_abs_error, req_rel_error, ERROR_INDIVIDUAL, &value, &err );

//   delete util->second;
//   delete util;

  return value / ( u_max - u_min ) / ( v_max - v_min );
}

int
LSBALink::PredictVariance ( unsigned n_dim,
                            const double* uv,
                            void* util,
                            unsigned f_dim,
                            double* f_val
                          )
{
  std::pair<const LSBALink*, size_t>* util_ = static_cast<std::pair<const LSBALink*, size_t>* > ( util );
  const LSBALink* lsba_ = util_->first;

  // Mean and variance of low-fidelity model
  double mean_low, variance_low;
  Matrix<lsba_real, Dynamic, 1> rhs_cov;
//     Matrix<lsba_real, Dynamic, 1> fh = lsba_->PointF ( uv[0], uv[1] );

  SparseVector<lsba_real> f;
  Matrix<lsba_real, Dynamic, 1> fh;

  lsba_->PointF ( uv[0], uv[1], f );

  fh.resize ( f.rows() );

  for ( size_t i = 0; i < fh.rows(); ++i )
    fh ( i ) = f.coeff ( i );

  int n1, n2;

  lsba_->lsba_low_->get_spline_size ( n1, n2 );

  // Compute F0 matrix
  SparseMatrix<lsba_real, RowMajor> F0;
  F0.resize ( lsba_->u_->size (), n1 * n2 );
  F0.reserve ( lsba_->u_->size () * 16 );
  // Fill model matrix
  for ( size_t i = 0; i < lsba_->u_->size (); ++i )
    lsba_->AddPointF0 ( F0, i, ( *lsba_->u_ ) [i], ( *lsba_->v_ ) [i] );
  F0.makeCompressed();
  F0.finalize();

  double sigma2_low = lsba_->lsba_low_->get_variance();

  Matrix<lsba_real, Dynamic, 1> ff;
  ConjugateGradient<SparseMatrix<lsba_real> > solver;

  if ( util_->second > 0 )
    solver.setMaxIterations ( util_->second );

  if ( lsba_->lambda_ == 0. ) {
    solver.compute ( lsba_->FF_ );
  } else {
    solver.compute ( lsba_->S_ );
  }

  ff = solver.solve ( fh );

  double mean ( 0. ), variance ( 0. );
  f_val[0] = 0.;
   
  if ( lsba_->w_ == NULL && lsba_->var_low_ == NULL ) {
    lsba_->lsba_low_->Predict ( uv[0], uv[1], mean_low, variance_low );
    lsba_->LSBA::Predict ( uv[0], uv[1], mean, variance );

    mean += mean_low;
    f_val[0] = variance + variance_low;

  } else {
    Matrix<lsba_real, Dynamic, 1> rhs_cov;

    ff = lsba_->F_ * ff ;

    if ( lsba_->w_ == NULL && lsba_->var_low_ != NULL ) {
      for ( size_t i = 0; i < ff.rows(); ++i )
        f_val[0] += ff ( i ) * ff ( i ) * ( ( *(lsba_->var_low_) ) [i] + lsba_->sigma2_ );
    } else {
      for ( size_t i = 0; i < ff.rows(); ++i )
        f_val[0] += ff ( i ) * ff ( i ) * ( 1. + lsba_->sigma2_ ) * ( *(lsba_->w_) ) [i];
    }

    Matrix<lsba_real, Dynamic, 1> f0 = F0.transpose() * ff;

    lsba_->lsba_low_->Predict ( uv[0], uv[1], mean_low, variance_low, f0, rhs_cov, util_->second );

    f_val[0] += variance_low + 2 * sigma2_low * f0.dot ( rhs_cov );

  }  
  
//   double sigma2_low = lsba_->lsba_low_->get_variance();
// 
//   Matrix<lsba_real, Dynamic, 1> ff;
//   ConjugateGradient<SparseMatrix<lsba_real> > solver;
// 
//   if ( util_->second > 0 )
//     solver.setMaxIterations ( util_->second );
// 
//   if ( lsba_->lambda_ == 0. ) {
//     solver.compute ( lsba_->FF_ );
//   } else {
//     solver.compute ( lsba_->S_ );
//   }
// 
//   ff = solver.solve ( fh );
// 
//   double mean ( 0. ), variance ( 0. );
// 
//   if ( lsba_->w_ == NULL && lsba_->var_low_ == NULL ) {
//     lsba_->lsba_low_->Predict ( uv[0], uv[1], mean_low, variance_low );
//     lsba_->LSBA::Predict ( uv[0], uv[0], mean, variance );
// 
//     f_val[0] = variance + variance_low;
// 
//   } else if ( lsba_->w_ == NULL && lsba_->var_low_ != NULL ) {
//     Matrix<lsba_real, Dynamic, 1> rhs_cov;
// 
//     ff = lsba_->F_ * ff ;
// 
//     for ( size_t i = 0; i < ff.rows(); ++i )
//       variance += ff ( i ) * ff ( i ) * ( ( * ( lsba_->var_low_ ) ) [i] + lsba_->sigma2_ );
// 
//     Matrix<lsba_real, Dynamic, Dynamic> F0f = Matrix<lsba_real, Dynamic, Dynamic> ( F0 );
// 
//     Matrix<lsba_real, Dynamic, 1> f0 = F0f.transpose() * ff;
// 
//     lsba_->lsba_low_->Predict ( uv[0], uv[1], mean_low, variance_low, f0, rhs_cov, util_->second );
// 
//     f_val[0] = variance + variance_low + 2 * sigma2_low * f0.dot ( rhs_cov );
//   } else {
//     Matrix<lsba_real, Dynamic, 1> rhs_cov;
//     SparseMatrix<lsba_real, RowMajor> Fww;
//     Fww.resize ( lsba_->Fw_.rows (), lsba_->Fw_.cols () );
//     Fww.reserve ( lsba_->Fw_.nonZeros() );
// 
//     for ( int k = 0; k < lsba_->Fw_.outerSize(); ++k ) {
//       Fww.startVec ( k );
//       for ( SparseMatrix<lsba_real, RowMajor>::InnerIterator it ( lsba_->F_, k ); it; ++it ) {
//         Fww.insertBack ( it.row (), it.col () ) = it.value () * ( *lsba_->w_ ) [it.row ()] * ( 1. + lsba_->sigma2_ );
//       }
//     }
//     Fww.finalize();
// 
//     Matrix<lsba_real, Dynamic, Dynamic> F0f = Matrix<lsba_real, Dynamic, Dynamic> ( F0 );
// 
//     Matrix<lsba_real, Dynamic, 1> f0 = F0f.transpose() * ( lsba_->F_ * ff );
// 
//     lsba_->lsba_low_->Predict ( uv[0], uv[1], mean_low, variance_low, f0, rhs_cov, util_->second );
// 
//     Matrix<lsba_real, Dynamic, 1> temp = lsba_->F_.transpose() * ( Fww * ff );
//     f_val[0] = ff.dot ( temp ) + variance_low;
// 
//     f_val[0] += 2 * sigma2_low * f0.dot ( rhs_cov );
//   }

  if ( f_val[0] < 0. ) {
    f_val[0] = 0.;
  }

  return 0;
}

int
LSBALink::PredictSd ( unsigned n_dim,
                      const double* uv,
                      void* util,
                      unsigned f_dim,
                      double* f_val
                    )
{
  std::pair<const LSBALink*, size_t>* util_ = static_cast<std::pair<const LSBALink*, size_t>* > ( util );
  const LSBALink* lsba_ = util_->first;

  // Mean and variance of low-fidelity model
  double mean_low, variance_low;
  Matrix<lsba_real, Dynamic, 1> rhs_cov;
//     Matrix<lsba_real, Dynamic, 1> fh = lsba_->PointF ( uv[0], uv[1] );

  SparseVector<lsba_real> f;
  Matrix<lsba_real, Dynamic, 1> fh;

  lsba_->PointF ( uv[0], uv[1], f );

  fh.resize ( f.rows() );

  for ( size_t i = 0; i < fh.rows(); ++i )
    fh ( i ) = f.coeff ( i );

  int n1, n2;

  lsba_->lsba_low_->get_spline_size ( n1, n2 );

  // Compute F0 matrix
  SparseMatrix<lsba_real, RowMajor> F0;
  F0.resize ( lsba_->u_->size (), n1 * n2 );
  F0.reserve ( lsba_->u_->size () * 16 );
  // Fill model matrix
  for ( size_t i = 0; i < lsba_->u_->size (); ++i )
    lsba_->AddPointF0 ( F0, i, ( *lsba_->u_ ) [i], ( *lsba_->v_ ) [i] );
  F0.makeCompressed();
  F0.finalize();

  double sigma2_low = lsba_->lsba_low_->get_variance();

  Matrix<lsba_real, Dynamic, 1> ff;
  ConjugateGradient<SparseMatrix<lsba_real> > solver;

  if ( util_->second > 0 )
    solver.setMaxIterations ( util_->second );

  if ( lsba_->lambda_ == 0. ) {
    solver.compute ( lsba_->FF_ );
  } else {
    solver.compute ( lsba_->S_ );
  }

  ff = solver.solve ( fh );

  double mean ( 0. ), variance ( 0. );
  f_val[0] = 0.;
   
  if ( lsba_->w_ == NULL && lsba_->var_low_ == NULL ) {
    lsba_->lsba_low_->Predict ( uv[0], uv[1], mean_low, variance_low );
    lsba_->LSBA::Predict ( uv[0], uv[1], mean, variance );

    mean += mean_low;
    f_val[0] = variance + variance_low;

  } else {
    Matrix<lsba_real, Dynamic, 1> rhs_cov;

    ff = lsba_->F_ * ff ;

    if ( lsba_->w_ == NULL && lsba_->var_low_ != NULL ) {
      for ( size_t i = 0; i < ff.rows(); ++i )
        f_val[0] += ff ( i ) * ff ( i ) * ( ( *(lsba_->var_low_) ) [i] + lsba_->sigma2_ );
    } else {
      for ( size_t i = 0; i < ff.rows(); ++i )
        f_val[0] += ff ( i ) * ff ( i ) * ( 1. + lsba_->sigma2_ ) * ( *(lsba_->w_) ) [i];
    }

    Matrix<lsba_real, Dynamic, 1> f0 = F0.transpose() * ff;

    lsba_->lsba_low_->Predict ( uv[0], uv[1], mean_low, variance_low, f0, rhs_cov, util_->second );

    f_val[0] += variance_low + 2 * sigma2_low * f0.dot ( rhs_cov );

  }

  if ( f_val[0] >= 0. ) {
    f_val[0] = sqrt ( f_val[0] );
  } else {
    f_val[0] = 0.;
  }

  return 0;
}

void
LSBALink::AddPointF0 ( SparseMatrix<lsba_real, RowMajor>& F,
                       int k,
                       double u,
                       double v
                     ) const
{
  // Add point to model matrix
  int n1, n2, m, n;

  lsba_low_->get_spline_size ( n1, n2 );
  m = n1 - 3;
  n = n2 - 3;

  double domain[4];

  lsba_low_->get_domain ( domain[0], domain[1], domain[2], domain[3] );

  // Map to the half open domain Omega = [0,m) x [0,n)
  // The mapped uc and vc must be (strictly) less than m and n respectively
  double uc = ( u - domain[0] ) / ( domain[2] - domain[0] ) * ( double ) m;
  double vc = ( v - domain[1] ) / ( domain[3] - domain[1] ) * ( double ) n;

  int i, j;
  double s, t;
  UCBspl::ijst ( m, n, uc, vc, i, j, s, t ); // i and j from -1

  double w_kl[4][4];

  // All the 16 tensors in the 4x4 neighbourhood of a point
  UCBspl::WKL ( s, t, w_kl ); // substituted by Odd Andersen, 16 dec. 2003

  // Fill model matrix
  if ( w_ == NULL ) {
    // Fill model matrix
    F.startVec ( k );
    for ( size_t ii = 0; ii < 4; ++ii )
      for ( size_t jj = 0; jj < 4; ++jj ) {
        F.insertBack ( k, ( i + 1 ) * n2 + j + 1 + ii * n2 + jj ) = w_kl[ii][jj];
      }
  } else {
    F.startVec ( k );
    for ( size_t ii = 0; ii < 4; ++ii )
      for ( size_t jj = 0; jj < 4; ++jj ) {
        F.insertBack ( k, ( i + 1 ) * n2 + j + 1 + ii * n2 + jj ) = w_kl[ii][jj] * ( *w_ ) [k];
      }
  }
//   return true;
}
