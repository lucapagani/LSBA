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

#include "lsba.hpp"
#include <stdexcept>

#include <mba/UCBsplines.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <chrono>

using namespace std;

LSBA::LSBA ( boost::shared_ptr<std::vector<double> > u,
             boost::shared_ptr<std::vector<double> > v,
             boost::shared_ptr<std::vector<double> > z,
             boost::shared_ptr<GenMatrix<lsba_real> > x,
             double lambda_fac,
             bool build_structure
           )
    : u_ ( u ), v_ ( v ), z_ ( z ), x_ ( x ), lambda_fac_ ( lambda_fac )
{
    n1_ = x_->noX();
    n2_ = x_->noY(); // This is also done when building equation system
    K_ = L_ = 4;

    CalculateDomain();

    // Resize vector of unknown
    beta_.resize ( n1_ * n2_ );
    // Resize zz_
    zz_.resize ( z_->size () );
    // Fill zz_
    for ( size_t i = 0; i < z_->size(); ++i )
        zz_ ( i ) = static_cast<lsba_real> ( ( *z_ ) [i] );

    if ( build_structure == true )
        BuildStructure ();

}

LSBA::LSBA ( boost::shared_ptr<std::vector<double> > u,
             boost::shared_ptr<std::vector<double> > v,
             boost::shared_ptr<std::vector<double> > z,
             boost::shared_ptr<std::vector<double> > w,
             boost::shared_ptr<GenMatrix<lsba_real> > x, // solution vector
             double lambda_fac,
             bool build_structure
           )
    : u_ ( u ), v_ ( v ), z_ ( z ), w_ ( w ), x_ ( x ), lambda_fac_ ( lambda_fac )
{
    n1_ = x_->noX() ;
    n2_ = x_->noY(); // This is also done when building equation system
    K_ = L_ = 4;

    CalculateDomain();

    // Resize vector of unknown
    beta_.resize ( n1_ * n2_ );
    // Resize zz_
    zz_.resize ( z_->size () );
    // Fill zz_
    for ( size_t i = 0; i < z_->size(); ++i )
        zz_ ( i ) = static_cast<lsba_real> ( ( *z_ ) [i] );

    if ( build_structure == true )
        BuildStructure ();

}

void LSBA::BuildStructure()
{
    // Resize model matrix
    F_.resize ( u_->size (), n1_ * n2_ );
    F_.reserve ( u_->size () * 16 );
    // Fill model matrix
    for ( size_t i = 0; i < u_->size (); ++i )
        AddPointF ( i, ( *u_ ) [i], ( *v_ ) [i] );
    F_.finalize();

    // Compute parameters
    if ( w_ == NULL )
        ComputeParameters ();
    else
        ComputeParametersW();

    if ( lambda_fac_ > 0. ) {
        ComputeLambda ();
        lambda_ *= lambda_fac_;
        BuildSparseStructure();
        ComputeSmoothingTerm ();
    } else {
        lambda_ = 0.;
    }
}

void LSBA::CalculateDomain()
{

    if ( domain_.size() != 4 )
        domain_.resize ( 4 );

    domain_[0] = *std::min_element ( u_->begin(), u_->end() );
    domain_[1] = *std::min_element ( v_->begin(), v_->end() );
    domain_[2] = *std::max_element ( u_->begin(), u_->end() );
    domain_[3] = *std::max_element ( v_->begin(), v_->end() );
}

void LSBA::set_domain ( double umin,
                        double vmin,
                        double umax,
                        double vmax,
                        bool build_structure
                      )
{
    if ( domain_.size() != 4 )
        domain_.resize ( 4 );

    domain_[0] = umin;
    domain_[1] = vmin;
    domain_[2] = umax;
    domain_[3] = vmax;

    if ( build_structure == true )
      BuildStructure();
}

void LSBA::get_domain ( double& umin,
                        double& vmin,
                        double& umax,
                        double& vmax
                      ) const
{

    umin = domain_[0];
    vmin = domain_[1];
    umax = domain_[2];
    vmax = domain_[3];
}

void
LSBA::AddPointF ( int k,
                  double u,
                  double v
                )
{
    // Add point to model matrix
    int m_ = n1_ - 3;
    int n_ = n2_ - 3;

    // Map to the half open domain Omega = [0,m) x [0,n)
    // The mapped uc and vc must be (strictly) less than m and n respectively
    double uc = ( u - umin() ) / ( umax()-umin() ) * ( double ) m_;
    double vc = ( v - vmin() ) / ( vmax()-vmin() ) * ( double ) n_;

//   if ( uc < 0.0 || vc < 0.0 || uc > m_ || vc > n_ ) {
//     throw runtime_error ( "ERROR in LSBA::AddPoint: A point was mapped "
//                           "to outside domain" );
//   }

    int i, j;
    double s, t;
    UCBspl::ijst ( m_, n_, uc, vc, i, j, s, t ); // i and j from -1

    double w_kl[4][4];

    // All the 16 tensors in the 4x4 neighbourhood of a point
    UCBspl::WKL ( s, t, w_kl ); // substituted by Odd Andersen, 16 dec. 2003

    // Fill model matrix
    F_.startVec ( k );
    for ( size_t ii = 0; ii < 4; ++ii )
        for ( size_t jj = 0; jj < 4; ++jj ) {
            F_.insertBack ( k, ( i + 1 ) * n2_ + j + 1 + ii * n2_ + jj ) = w_kl[ii][jj];
        }

//   return true;
}

void
LSBA::FillBeta ()
{
    for ( int i = 0; i < n1_; ++i )
        for ( int j = 0; j < n2_; ++j )
            beta_ ( i * n2_ + j ) = ( *x_ ) ( i - 1, j - 1 );
}

void
LSBA::Compute ( int MaxIteration,
                bool compute_smoothing_matrix
              )
{
    // Fill beta_
    FillBeta ();

    ConjugateGradient<SparseMatrix<lsba_real> > solver;
    solver.setMaxIterations ( MaxIteration );

    if ( lambda_ == 0. ) {
        solver.compute ( FF_ );
    } else {
        solver.compute ( S_ );
    }

    if ( compute_smoothing_matrix == false ) {
//     cout << beta_.rows() << " " << Fz_.rows() << " " << Fz_.cols() << endl;
//     cout << beta_ << endl;
//     cout << Fz_ << endl;
        beta_ = solver.solveWithGuess ( Fz_, beta_ );
    } else {

        SparseMatrix<lsba_real, ColMajor> F_t;

        if ( w_ == NULL ) {
            F_t  = F_.transpose();
        } else {
            F_t = Fw_.transpose();
        }

        SparseMatrix<lsba_real, ColMajor> temp = solver.solve ( F_t );

        if ( w_ != NULL ) {
            for ( int k = 0; k < temp.outerSize(); ++k ) {
                for ( SparseMatrix<lsba_real, RowMajor>::InnerIterator it ( F_, k ); it; ++it ) {
                    temp.coeffRef ( it.col (), it.row () ) *= sqrt ( ( *w_ ) [it.row ()] );
                }
            }
        }

        Sl_ = F_ * temp;
        beta_ = temp * zz_;
    }

    size_t index = 0;
    for ( int i = 0; i < n1_; ++i )
        for ( int j = 0; j < n2_; ++j, ++index )
            ( *x_ ) ( i - 1, j - 1 ) = beta_ ( index );

}

void
LSBA::ComputeDirect ( bool compute_smoothing_matrix )
{
    SimplicialLDLT <SparseMatrix<lsba_real> > solver;

    if ( lambda_ == 0. ) {
        solver.compute ( FF_ );
//     solver_.compute ( FF_ );
    } else {
        solver.compute ( S_ );
//     solver_.compute ( S_ );
    }

    if ( compute_smoothing_matrix == false ) {
//   cout << FF_ << endl << endl;
        beta_ = solver.solve ( Fz_ );
    } else {

        SparseMatrix<lsba_real, ColMajor> F_t;

        if ( w_ == NULL ) {
            F_t  = F_.transpose();
        } else {
            F_t = Fw_.transpose();
        }

        SparseMatrix<lsba_real, ColMajor> temp = solver.solve ( F_t );

        if ( w_ != NULL ) {
            for ( int k = 0; k < temp.outerSize(); ++k ) {
                for ( SparseMatrix<lsba_real, RowMajor>::InnerIterator it ( F_, k ); it; ++it ) {
                    temp.coeffRef ( it.col (), it.row () ) *= sqrt ( ( *w_ ) [it.row ()] );
                }
            }
        }

        Sl_ = F_ * temp;
        beta_ = temp * zz_;
    }

    if ( solver.info() != Eigen::Success ) {
        // decomposition failed
        cerr << "Decomposition failed" << endl;
    } else {
        size_t index = 0;
        for ( int i = 0; i < n1_; ++i )
            for ( int j = 0; j < n2_; ++j, ++index )
                ( *x_ ) ( i - 1, j - 1 ) = beta_ ( index );
    }
}

//! Compute F'F and F'z
void
LSBA::ComputeParameters ()
{
    FF_ = F_.transpose () * F_;
    FF_.makeCompressed();

    Fz_ = F_.transpose () * zz_;
}

//! Compute sqrt(W) * F, sqrt(W) * z, F' * W * F and F'z
void
LSBA::ComputeParametersW ()
{

    Fw_.resize ( F_.rows (), F_.cols () );
    Fw_.reserve ( F_.nonZeros() );

    for ( int k = 0; k < F_.outerSize(); ++k ) {
        Fw_.startVec ( k );
        for ( SparseMatrix<lsba_real, RowMajor>::InnerIterator it ( F_, k ); it; ++it ) {
            Fw_.insertBack ( it.row (), it.col () ) = it.value () * sqrt ( ( *w_ ) [it.row ()] );
        }
    }
    Fw_.finalize();

    zzw_ = zz_;
    for ( size_t i = 0; i < zzw_.rows (); ++i )
        zzw_[i] *= sqrt ( ( *w_ ) [i] );

    FF_ = Fw_.transpose () * Fw_;
    FF_.makeCompressed();

    Fz_ = Fw_.transpose () * zzw_;
}

void LSBA::ComputeLambda()
{

    // Calculate lambda from the norm of the
    // matrices: lambda = Norm_A/norm_E where E is the smoothing matrix.
    // This value will be multiplied with lambdaFac_ (if lamdaFac_ is different
    // from 1.0, which is the default).
    // Must  be done after A_ is calculated (from the scattered data)
    // (Note that after addPoint, the smoothing term is not correct.)

    // Calculate lambda for smoothing term
    double ff_norm = sqrt ( FF_.cwiseAbs2().sum() );

    int m_ = n1_-3;
    int n_ = n2_-3;

    double du = ( umax() - umin() ) / ( double ) m_;
    double dv = ( vmax() - vmin() ) / ( double ) n_;

    SmoothMatrix E ( n1_, n2_, du, dv );

    double e_norm = E.norm_l2();
//   lambda_ = 1.0;
    if ( e_norm != 0.0 ) {
        lambda_ = ff_norm / e_norm;
    } else {
        throw runtime_error ( "LSBA::ComputeLambda(): e_norm is zero when "
                              "calculating lambda, exit now in Debug mode" );
    }
}

void LSBA::ComputeSmoothingTerm ()
{
    // where lambdaVal is either the original lambda_*lambdaFac_,
    // which is given when the system is built, or it is adjustment
    // if we want to smooth/recover the surface later.

    // Note that the smoothing term does not depend on the scattered data,
    // Only the system matrix is affected.

    int m_ = n1_ - 3;
    int n_ = n2_ - 3;
    double du = ( umax() - umin() ) / ( double ) m_;
    double dv = ( vmax() - vmin() ) / ( double ) n_;

    SmoothMatrix E ( m_ + 3, n_ + 3, du, dv );

    int n_row, n_col;
    for ( int k = 0; k < FF_.outerSize(); ++k ) {
        for ( SparseMatrix<lsba_real, ColMajor>::InnerIterator it ( FF_, k ); it; ++it ) {
            S_.coeffRef ( it.row(), it.col() ) = it.value();
        }
    }

    for ( int k = 0; k < S_.outerSize(); ++k ) {
        int i = k % n1_;
        int j = k / n1_;
        for ( SparseMatrix<lsba_real, ColMajor>::InnerIterator it ( S_, k ); it; ++it ) {
            n_row = it.row();
            n_col = it.col();
            int r = n_row % n1_;
            int s = n_row / n1_;
            S_.coeffRef ( n_row, n_col ) += lambda_ * E ( i, j, r, s );
        }
    }

    S_.makeCompressed();
//   S_.finalize();

}

//! Compute the number of non zeros entries of the matrix E
static int
NoNonZeros ( int n1,
             int n2,
             int K,
             int L
           )
{
    // Count non-zeros
    int numNonZeros = 0;
    int imin,imax,jmin,jmax;
    for ( int j = 0; j < n2; j++ ) {
        jmin = max ( 0,j-L+1 );
        jmax = min ( n2-1,j+L-1 );
        for ( int i = 0; i < n1; i++ ) {
            imin = max ( 0,i-K+1 );
            imax = min ( n1-1,i+K-1 );
            numNonZeros += ( imax - imin + 1 ) * ( jmax - jmin +1 );
        }
    }
    return numNonZeros;
}

void
LSBA::BuildSparseStructure()
{
    // Must be run after n1_ and n2_ have been set

    if ( n1_ < 1 || n2_ < 1 ) {
        throw runtime_error ( "ERROR in buildSparseStructure, inproper initialization" );
    }

    // Build the data structure for FEl_ and init with zeros

    // Count non-zeros
    int numNonZeros = NoNonZeros ( n1_, n2_, K_, L_ );

    S_.resize ( n1_*n2_, n1_*n2_ );
    S_.reserve ( numNonZeros );

    for ( int j = 0; j < n1_; j++ ) {
        int jmin = max ( 0,j-K_+1 );
        int jmax = min ( n1_-1,j+K_-1 );
        for ( int i = 0; i < n2_; i++ ) {
            int ii = j * n2_ + i; // row number
            int imin = max ( 0,i-L_+1 );
            int imax = min ( n2_-1,i+L_-1 );
            S_.startVec ( ii );
            for ( int s = jmin; s <= jmax; s++ ) {
                for ( int r = imin; r <= imax; r++ ) {
                    S_.insertBack ( s * n1_ + r, ii ) = 0.0;
                }
            }
        }
    }

}

void
LSBA::set_weights ( boost::shared_ptr<std::vector<double> > weights,
  bool build_structure
)
{
    if ( weights->size() != u_->size() ) {
        std::cerr << "Weights size must be equal to points size" << std::endl;
        exit ( EXIT_FAILURE );
    }
    w_ = weights;

    if ( build_structure == true )
      BuildStructure ();
}

void
LSBA::set_smoothing_factor ( double smoothingFac )
{
    if ( smoothingFac < 0. ) {
        cerr << "Smoothing factor must be positive or null" << endl;
        exit ( EXIT_FAILURE );
    }
    if ( lambda_fac_ == 0. && smoothingFac > 0. ) {
        lambda_fac_ = smoothingFac;
        ComputeLambda ();
        lambda_ *= lambda_fac_;
        BuildSparseStructure();
        ComputeSmoothingTerm ();
    } else if ( smoothingFac > 0. ) {
        lambda_ *= smoothingFac / lambda_fac_;
        lambda_fac_ = smoothingFac;
        ComputeSmoothingTerm ();
    } else {
        lambda_ = 0.;
    }
}

//! Prediction in the location (u,v)
double
LSBA::Predict ( double u,
                double v
              ) const
{
    SparseVector<lsba_real> f;
//     Matrix<lsba_real, Dynamic, 1> f0;

    PointF ( u, v, f );

//     f0.resize ( f.rows(), f.cols() );
//     for ( size_t i = 0; i < f0.rows(); ++i )
//         f0 ( i ) = f.coeff( i );

    double val = f.dot ( beta_ );

//   f.resize ( 0 );

    return val;
}

//! Predict the mean and the variance in the location (u,v)
void
LSBA::Predict ( double u,
                double v,
                double& mean,
                double& variance,
                size_t max_iterations
              ) const
{
    mean = variance = 0.;
    SparseVector<lsba_real> f0;
    Matrix<lsba_real, Dynamic, 1> f;

    PointF ( u, v, f0 );

    f.resize ( f0.rows() );

    for ( size_t i = 0; i < f0.rows(); ++i )
        f ( i ) = f0.coeff( i );
//     f.setZero();

//     for ( SparseVector<lsba_real>::InnerIterator it ( f0 ); it; ++it )
//         f ( it.row () ) = it.value();

    mean = f.dot ( beta_ );

//   SparseVector<lsba_real> ff;
    Matrix<lsba_real, Dynamic, 1> ff;
    ConjugateGradient<SparseMatrix<lsba_real> > solver;
    if ( max_iterations > 0 )
        solver.setMaxIterations ( max_iterations );

    if ( lambda_ == 0. ) {
        solver.compute ( FF_ );
        ff = solver.solve ( f );
        variance = sigma2_ *  f.dot ( ff );
    } else {
        solver.compute ( S_ );
        ff = solver.solve ( f );
        Matrix<lsba_real, Dynamic, 1> temp = FF_ * ff;
        variance = sigma2_ * ff.dot ( temp );
    }

}

/** Predict the mean and the variance in the location (u,v)
 *
 * @param u: value of the u coordinate of the point
 * @param v: value of the v coordinate of the point
 * @param mean: value of the predicted mean
 * @param variance: value of the variance of the predicted point in the location (u,v)
 * @param max_iteration: maximum number of iteration for inverse the matrix
 */
void
LSBA::Predict ( boost::shared_ptr< std::vector<double> > u,
                boost::shared_ptr< std::vector<double> > v,
                std::vector<double>& mean,
                std::vector<double>& variance,
                int max_iterations
              ) const
{
    size_t size = u->size();

    if ( size != v->size() ) {
        cerr << "Vector u and v must have the same size." << endl;
        exit ( EXIT_FAILURE );
    }

//   mean.clear();
//   variance.clear();
    mean.resize ( size );
    variance.resize ( size );

    Matrix<lsba_real, Dynamic, 1> f, ff;
    ConjugateGradient<SparseMatrix<lsba_real> > solver;
    SparseMatrix<lsba_real, ColMajor> m_inv ( FF_.cols(), FF_.cols() ), identity ( FF_.cols(), FF_.cols() );

    identity.reserve ( identity.rows() );
    for ( size_t i = 0; i < identity.rows(); ++i ) {
        identity.startVec ( i );
        identity.insertBack ( i, i ) = 1.;
    }
    identity.finalize();

    if ( lambda_ == 0. ) {
        m_inv.reserve ( FF_.nonZeros() );

        solver.compute ( FF_ );
        solver.setMaxIterations ( max_iterations );

        Matrix<lsba_real, Dynamic, 1> base ( size ), invCol ( size );

        m_inv = solver.solve ( identity );

        for ( size_t i = 0; i < size; ++i ) {
            ff = m_inv * f;
            f = PointF ( ( *u ) [i], ( *v ) [i] );
            mean [i] = f.dot ( beta_ );
            variance [i] = sigma2_ *  f.dot ( ff );
        }
    } else {
        m_inv.reserve ( S_.nonZeros() );

        solver.compute ( S_ );
        solver.setMaxIterations ( max_iterations );

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        m_inv = solver.solve ( identity );

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        cout << "Inversion time: " << elapsed_seconds.count() << endl;


        Matrix<lsba_real, Dynamic, 1> temp;

        for ( size_t i = 0; i < size; ++i ) {
            f = PointF ( ( *u ) [i], ( *v ) [i] );

            ff = m_inv * f;
            temp = FF_ * ff;

            mean [i] = f.dot ( beta_ );
            variance [i] = sigma2_ *  ff.dot ( temp );
        }
    }
}

/** Predict the mean and the variance in the location (u,v)
 *
 * @param u: value of the u coordinate of the point
 * @param v: value of the v coordinate of the point
 * @param mean: value of the predicted mean
 * @param variance: value of the variance of the predicted point in the location (u,v)
 * @param rhs_low: right-hand side term in order to compute the prediction variance of the linkage model
 * @param max_iteration: maximum number of iteration for inverse the matrix
 */
void
LSBA::Predict ( double u,
                double v,
                double& mean,
                double& variance,
                Matrix<lsba_real, Dynamic, 1> f0,
                Matrix<lsba_real, Dynamic, 1>& rhs_low,
                size_t max_iterations
              ) const
{
    SparseVector<lsba_real> f1;
    Matrix<lsba_real, Dynamic, 1> f;

    PointF ( u, v, f1 );

    f.resize ( f1.rows() );

    for ( size_t i = 0; i < f.rows(); ++i )
        f ( i ) = f1.coeff( i );

    mean = f.dot ( beta_ );

    f -= f0;

//   SparseVector<lsba_real> ff;
//   Matrix<lsba_real, Dynamic, 1> ff;
    ConjugateGradient<SparseMatrix<lsba_real> > solver;
    if ( max_iterations > 0 )
        solver.setMaxIterations ( max_iterations );


    if ( lambda_ == 0. ) {
        solver.compute ( FF_ );
        rhs_low = solver.solve ( f );
        variance = sigma2_ *  f.dot ( rhs_low );
    } else {
        solver.compute ( S_ );
        rhs_low = solver.solve ( f );
        Matrix<lsba_real, Dynamic, 1> temp = FF_ * rhs_low;
        variance = sigma2_ * rhs_low.dot ( temp );
        rhs_low = solver.solve ( temp );

    }
}

int
LSBA::PredictVariance ( unsigned n_dim,
                        const double* uv,
                        void* util,
                        unsigned f_dim,
                        double* f_val
                      )
{

//   size_t max_iterations = * ( static_cast<size_t*> ( p_max_iterations ) );
    std::pair<const LSBA*, size_t>* util_ = static_cast<std::pair<const LSBA*, size_t>* > ( util );
//   size_t max_iterations = 0;
//   LSBA* lsba_ = static_cast<LSBA*> ( lsba );
    const LSBA* lsba_ = util_->first;

    SparseVector<lsba_real> f0;
    Matrix<lsba_real, Dynamic, 1> f;

    lsba_->PointF ( uv[0], uv[1], f0 );

    f.resize ( f0.rows() );

    for ( size_t i = 0; i < f.rows(); ++i )
        f ( i ) = f0.coeff( i );

    Matrix<lsba_real, Dynamic, 1> ff;
    ConjugateGradient<SparseMatrix<lsba_real> > solver;
    if ( util_->second > 0 )
        solver.setMaxIterations ( util_->second );


    if ( lsba_->lambda_ == 0. ) {
        solver.compute ( lsba_->FF_ );
        ff = solver.solve ( f );
        f_val[0] = lsba_->sigma2_ *  f.dot ( ff );
    } else {
        solver.compute ( lsba_->S_ );
        ff = solver.solve ( f );
        Matrix<lsba_real, Dynamic, 1> temp = lsba_->FF_ * ff;
        f_val[0] = lsba_->sigma2_ * ff.dot ( temp );
    }

    return 0;
}


int
LSBA::PredictSd ( unsigned n_dim,
                  const double* uv,
                  void* util,
                  unsigned f_dim,
                  double* f_val
                )
{

    std::pair<const LSBA*, size_t>* util_ = static_cast<std::pair<const LSBA*, size_t>* > ( util );
    const LSBA* lsba_ = util_->first;

    SparseVector<lsba_real> f0;
    Matrix<lsba_real, Dynamic, 1> f;

    lsba_->PointF ( uv[0], uv[1], f0 );

    f.resize ( f0.rows() );

    for ( size_t i = 0; i < f.rows(); ++i )
        f ( i ) = f0.coeff( i );

    Matrix<lsba_real, Dynamic, 1> ff;
    ConjugateGradient<SparseMatrix<lsba_real> > solver;
    if ( util_->second > 0 )
        solver.setMaxIterations ( util_->second );


    if ( lsba_->lambda_ == 0. ) {
        solver.compute ( lsba_->FF_ );
        ff = solver.solve ( f );
        f_val[0] = lsba_->sigma2_ *  f.dot ( ff );
    } else {
        solver.compute ( lsba_->S_ );
        ff = solver.solve ( f );
        Matrix<lsba_real, Dynamic, 1> temp = lsba_->FF_ * ff;
        f_val[0] = lsba_->sigma2_ * ff.dot ( temp );
    }

    f_val[0] = sqrt ( f_val[0] );

    return 0;
}

/** Compute the mean variance in the uv domain [u_min, u_max] x [v_min, v_max]
 *
 */
double
LSBA::ComputeMeanVariance ( double u_min,
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

    std::pair<const LSBA*, size_t>* util ( new std::pair<const LSBA*, size_t> );
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
LSBA::ComputeMeanSd ( double u_min,
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

    std::pair<const LSBA*, size_t>* util ( new std::pair<const LSBA*, size_t> );
    *util = std::make_pair ( this, max_iterations );

    hcubature ( 1, PredictSd, util, 2, range_min, range_max, max_eval, req_abs_error, req_rel_error, ERROR_INDIVIDUAL, &value, &err );

//   delete util->second;
//   delete util;

    return value / ( u_max - u_min ) / ( v_max - v_min );
}

void
LSBA::PointF ( double u,
               double v,
               SparseVector<lsba_real>& f
             ) const
{
    // Add point to model matrix
    int m_ = n1_ - 3;
    int n_ = n2_ - 3;

    // Map to the half open domain Omega = [0,m) x [0,n)
    // The mapped uc and vc must be (strictly) less than m and n respectively
    double uc = ( u - umin() ) / ( umax()-umin() ) * ( double ) m_;
    double vc = ( v - vmin() ) / ( vmax()-vmin() ) * ( double ) n_;

//   if ( uc < 0.0 || vc < 0.0 || uc > m_ || vc > n_ ) {
//     throw runtime_error ( "ERROR in LSBA::AddPoint: A point was mapped "
//                           "to outside domain" );
//   }

    int i, j;
    double s, t;
    UCBspl::ijst ( m_, n_, uc, vc, i, j, s, t ); // i and j from -1

    double w_kl[4][4];

    // All the 16 tensors in the 4x4 neighbourhood of a point
    UCBspl::WKL ( s, t, w_kl ); // substituted by Odd Andersen, 16 dec. 2003

    // Fill model matrix
    f.resize ( F_.cols() );
    for ( size_t ii = 0; ii < 4; ++ii )
        for ( size_t jj = 0; jj < 4; ++jj ) {
            f.insertBack ( ( i + 1 ) * n2_ + j + 1 + ii * n2_ + jj ) = w_kl[ii][jj];
        }

}

Matrix<lsba_real, Dynamic, 1>
LSBA::PointF ( double u,
               double v
             ) const
{
    // Add point to model matrix
    int m_ = n1_ - 3;
    int n_ = n2_ - 3;

    // Map to the half open domain Omega = [0,m) x [0,n)
    // The mapped uc and vc must be (strictly) less than m and n respectively
    double uc = ( u - umin() ) / ( umax()-umin() ) * ( double ) m_;
    double vc = ( v - vmin() ) / ( vmax()-vmin() ) * ( double ) n_;

//   if ( uc < 0.0 || vc < 0.0 || uc > m_ || vc > n_ ) {
//     throw runtime_error ( "ERROR in LSBA::AddPoint: A point was mapped "
//                           "to outside domain" );
//   }

    int i, j;
    double s, t;
    UCBspl::ijst ( m_, n_, uc, vc, i, j, s, t ); // i and j from -1

    double w_kl[4][4];

    // All the 16 tensors in the 4x4 neighbourhood of a point
    UCBspl::WKL ( s, t, w_kl ); // substituted by Odd Andersen, 16 dec. 2003

    // Fill model matrix
    Matrix<lsba_real, Dynamic, 1> f;
    f.resize ( F_.cols() );
    f.setZero ();
    for ( size_t ii = 0; ii < 4; ++ii )
        for ( size_t jj = 0; jj < 4; ++jj ) {
            f ( ( i + 1 ) * n2_ + j + 1 + ii * n2_ + jj ) = w_kl[ii][jj];
        }

    return f;
}

//! Compute estimation of variance
void
LSBA::ComputeVariance ()
{
    sigma2_ = 0.;
    size_t n = z_->size ();

    if ( w_ == NULL ) {
        for ( size_t i = 0; i < n; ++i ) {
            sigma2_ += pow ( ( *z_ ) [i] - Predict ( ( *u_ ) [i], ( *v_ ) [i] ), 2 );
        }
    } else {
        for ( size_t i = 0; i < n; ++i ) {
            sigma2_ += pow ( ( *z_ ) [i] - Predict ( ( *u_ ) [i], ( *v_ ) [i] ), 2 ) * ( *w_ ) [i];
        }
    }

    if ( lambda_ == 0. ) {
        sigma2_ /= ( n - F_.cols () );
    } else {
        // Check if smoothing matrix exesists
        if ( Sl_.rows() == 0 ) {
            ConjugateGradient<SparseMatrix<lsba_real> > solver;
            solver.compute ( S_ );
            solver.setMaxIterations ( 1e2 );

            SparseMatrix<lsba_real, ColMajor> F_t;

            if ( w_ == NULL ) {
                F_t  = F_.transpose();
            } else {
                F_t = Fw_.transpose();
            }

            SparseMatrix<lsba_real, ColMajor> temp = solver.solve ( F_t );

            if ( w_ != NULL ) {
                for ( int k = 0; k < temp.outerSize(); ++k ) {
                    for ( SparseMatrix<lsba_real, RowMajor>::InnerIterator it ( F_, k ); it; ++it ) {
                        temp.coeffRef ( it.col (), it.row () ) *= sqrt ( ( *w_ ) [it.row ()] );
                    }
                }
            }

            Sl_ = F_ * temp;
        }

        double df = Sl_.diagonal().sum();

        if ( df < n ) {
            sigma2_ /= ( n - df );
        } else {
            cerr << "\033[01;33mDegree of freedom to high -> call ComputeVarianceFast()\033[0m" << endl;
            ComputeVarianceFast();
        }
    }
}

//! Compute estimation of variance
void
LSBA::ComputeVarianceFast ()
{
    sigma2_ = 0.;
    size_t n = z_->size ();

    if ( w_ == NULL ) {
        for ( size_t i = 0; i < n; ++i ) {
            sigma2_ += pow ( ( *z_ ) [i] - Predict ( ( *u_ ) [i], ( *v_ ) [i] ), 2 );
        }
    } else {
        for ( size_t i = 0; i < n; ++i ) {
            sigma2_ += pow ( ( *z_ ) [i] - Predict ( ( *u_ ) [i], ( *v_ ) [i] ), 2 ) * ( *w_ ) [i];
        }
    }

    sigma2_ /= n;
}

void
LSBA::AddPointFw ( int k,
                   double u,
                   double v
                 )
{
    // Add point to model matrix
    int m_ = n1_ - 3;
    int n_ = n2_ - 3;

    // Map to the half open domain Omega = [0,m) x [0,n)
    // The mapped uc and vc must be (strictly) less than m and n respectively
    double uc = ( u - umin() ) / ( umax()-umin() ) * ( double ) m_;
    double vc = ( v - vmin() ) / ( vmax()-vmin() ) * ( double ) n_;

    if ( uc < 0.0 || vc < 0.0 || uc > m_ || vc > n_ ) {
        throw runtime_error ( "ERROR in LSBA::AddPoint: A point was mapped "
                              "to outside domain" );
    }

    int i, j;
    double s, t;
    UCBspl::ijst ( m_, n_, uc, vc, i, j, s, t ); // i and j from -1

    double w_kl[4][4];

    // All the 16 tensors in the 4x4 neighbourhood of a point
    UCBspl::WKL ( s, t, w_kl ); // substituted by Odd Andersen, 16 dec. 2003

    // Fill model matrix
    Fw_.startVec ( k );
    for ( size_t ii = 0; ii < 4; ++ii )
        for ( size_t jj = 0; jj < 4; ++jj ) {
            Fw_.insertBack ( k, ( i + 1 ) * n2_ + j + 1 + ii * n2_ + jj ) = w_kl[ii][jj] * sqrt ( ( *w_ ) [k] );
        }

}
