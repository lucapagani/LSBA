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

#include <lsba/lsba_link.hpp>
#include <iostream>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <mba/MBA.h>
#include <ctime>

using namespace std;

typedef vector<double> VecD;
typedef boost::shared_ptr<vector<double> > VecDPtr;

int RANK = 1;

// Nominal function
double
peaks ( double u,
        double v
      )
{
  return ( 3. * pow ( 1. - u, 2 ) * exp ( - ( pow ( u, 2 ) ) - pow ( v + 1., 2 ) ) - 10. * ( u / 5. - pow ( u, 3 ) - pow ( v, 5 ) ) * exp ( - pow ( u, 2 ) - pow ( v, 2 ) ) - 1. / 3. * exp ( - pow ( u + 1., 2 ) - pow ( v, 2 ) ) );
}


// Estimate parameters LSBA
boost::shared_ptr<LSBA>
EstimateParameters ( VecDPtr u,
                     VecDPtr v,
                     VecDPtr z,
                     int h,
                     double lambda = 1.,
                     VecDPtr w = NULL
                   )
{

  MBA mba ( u, v, z );
  mba.setDomain ( -3., -3., 3., 3. );
  mba.MBAalg ( 1, 1, h );

  UCBspl::SplineSurface surf = mba.getSplineSurface();
  boost::shared_ptr<GenMatrixType> phi = surf.getCoefficients();

  boost::shared_ptr<LSBA> lsba ( new LSBA ( u, v, z, phi, lambda, false ) );
  lsba->set_domain ( -3., -3., 3., 3., false );

  if ( w != NULL )
    lsba->set_weights ( w );

  lsba->BuildStructure();

  lsba->Compute ( 1e2 );

  return lsba;

}

// Estimate paremeters LSBA linking model
boost::shared_ptr<LSBALink>
EstimateParameters ( VecDPtr u_h,
                     VecDPtr v_h,
                     VecDPtr z_h,
                     int h,
                     boost::shared_ptr<LSBA> lsba_low,
                     bool compute_var_fast = false,
                     double lambda = 1.,
                     bool compute_w = true,
                     VecDPtr var_low = 0
                   )
{

  int n = u_h->size();

  VecDPtr z_res ( new VecD );
  z_res->reserve ( n );

  for ( size_t i = 0; i < n; ++i ) {
    z_res->push_back ( ( *z_h ) [i] - lsba_low->Predict ( ( *u_h ) [i], ( *v_h ) [i] ) );
  }

  MBA mba ( u_h, v_h, z_res );
  mba.setDomain ( -3., -3., 3., 3. );
  mba.MBAalg ( 1, 1, h );

  UCBspl::SplineSurface surf = mba.getSplineSurface();
  boost::shared_ptr<GenMatrixType> phi = surf.getCoefficients();

  boost::shared_ptr<LSBALink> lsba ( new LSBALink ( u_h, v_h, z_h, phi, lsba_low, lambda, compute_w ) );
  lsba->set_domain ( -3., -3., 3., 3. );
  lsba->Compute ( 1e2 );

  if ( compute_w == false ) {
    if ( var_low != 0 ) {
      var_low->resize ( n );
      #pragma omp parallel for
      for ( size_t i = 0; i < n; ++i ) {
        double mean ( 0. ), variance ( 0. );
        lsba_low->Predict ( ( *u_h ) [i], ( *v_h ) [i], mean, variance );
        ( *var_low ) [i] = variance;
      }
      lsba->set_variance_low ( var_low );
    } else {
      double variance = lsba_low->ComputeMeanVariance ( -3., -3., 3., 3., 1e2, 1e-4, 1e-4, 1e2 );
      VecDPtr var ( new VecD ( n,  variance ) );
      lsba->set_variance_low ( var );
    }
  }

  compute_var_fast == false ? lsba->ComputeVariance() : lsba->ComputeVarianceFast();

  return lsba;
}

// Predict z points
template<typename model>
void
Predict ( boost::shared_ptr<model> lsba,
          VecDPtr u,
          VecDPtr v,
          VecDPtr z_pred
        )
{
  size_t n = u->size();

  z_pred->clear();
  z_pred->reserve ( n );

  for ( size_t i = 0; i < n; ++i ) {
    z_pred->push_back ( lsba->Predict ( ( *u ) [i], ( *v ) [i] ) );
  }
}

// Simulate Lo-Fi data
void
SimulateLo ( VecDPtr u,
             VecDPtr v,
             VecDPtr z,
             std::default_random_engine &generator
           )
{
  z->resize ( u->size () );

  std::normal_distribution<double> rnorm ( 0., 0.4 );

  for ( auto it_u = u->begin (), it_v = v->begin (), it_z = z->begin(); it_u != u->end (); ++it_u, ++it_v, ++it_z ) {
    *it_z = peaks ( *it_u, *it_v ) + rnorm ( generator ) + pow ( *it_u, 2 ) / 10. + pow ( *it_v, 2 ) / 10.;
  }
}

// Simulate Hi-Fi data
void
Simulate ( VecDPtr u,
           VecDPtr v,
           VecDPtr z,
           std::default_random_engine &generator
         )
{
  z->resize ( u->size () );

  std::normal_distribution<double> rnorm ( 0., 0.1 );

  for ( auto it_u = u->begin (), it_v = v->begin (), it_z = z->begin (); it_u != u->end (); ++it_u, ++it_v, ++it_z ) {
    *it_z = peaks ( *it_u, *it_v ) + rnorm ( generator );
  }
}

int main ( int argc, char **argv )
{
  // Initialize vectors
  VecDPtr u_low ( new VecD ), v_low ( new VecD ), z_low ( new VecD );
  VecDPtr u_high ( new VecD ), v_high ( new VecD ), z_high ( new VecD );
  VecDPtr u_test ( new VecD ), v_test ( new VecD ), z_test ( new VecD );

  // Number of Lo-Fi points
  int n_low = .5e3;
  u_low->reserve ( n_low );
  v_low->reserve ( n_low );

  double u0, v0;

  double d = 6. / ( n_low - 1. );
  for ( size_t i = 0; i < n_low; ++i )
    for ( size_t j = 0; j < n_low; ++j ) {
      u0 = -3 + i * d;
      v0 = -3 + j * d;

      u_low->push_back ( u0 );
      v_low->push_back ( v0 );
    }

  // Number of Hi-Fi points
  int n_high = 20;
  u_high->reserve ( n_high );
  v_high->reserve ( n_high );

  d = 6. / ( n_high - 1. );
  for ( size_t i = 0; i < n_high; ++i )
    for ( size_t j = 0; j < n_high; ++j ) {
      u_high->push_back ( -3. + i * d );
      v_high->push_back ( -3. + j * d );
    }

  // Number of test points for prediction
  int n_test = 100;
  u_test->reserve ( n_test );
  v_test->reserve ( n_test );

  d = 6. / ( n_test - 1. );
  for ( size_t i = 0; i < n_test; ++i )
    for ( size_t j = 0; j < n_test; ++j ) {
      u_test->push_back ( -3. + i * d );
      v_test->push_back ( -3. + j * d );
    }

  // Choose number of levels of MBA algorithm
  int h_low = 5, h_high = 5, h_link = 4;

  // lambda of for the LSBA algorithm
  double lambda_lo = 10.;
  double lambda = 1.;
  
  // Random number generator
  std::default_random_engine generator;
  
  std::clock_t start;
  double duration;

  // Simulate Lo-Fi data
  SimulateLo ( u_low, v_low, z_low, generator );
 
  // Estimate Lo-Fi parameters
  start = std::clock();  
  boost::shared_ptr<LSBA> lsba_low = EstimateParameters ( u_low, v_low, z_low, h_low, lambda_lo );
  duration = ( std::clock() - start ) / ( double ) CLOCKS_PER_SEC;
  cout << "Lo-Fi estimation time: " << duration << " s" << endl;

  lsba_low->ComputeVarianceFast();

  // Simulate Hi-Fi data
  Simulate ( u_high, v_high, z_high, generator );
  
  // Estimate Hi-Fi parameters
  start = std::clock();
  boost::shared_ptr<LSBA> lsba_hi = EstimateParameters ( u_high, v_high, z_high, h_high, lambda );
  duration = ( std::clock() - start ) / ( double ) CLOCKS_PER_SEC;
  cout << "Hi-Fi estimation time: " << duration << " s" << endl;

  // Estimate link parameters
  boost::shared_ptr<LSBALink> lsba_link = EstimateParameters ( u_high, v_high, z_high, h_link, lsba_low, true, lambda, false );
  duration = ( std::clock() - start ) / ( double ) CLOCKS_PER_SEC;
  cout << "Fusion estimation time: " << duration << " s" << endl;

  // Predict in the test locations
  VecDPtr pred_low ( new vector<double> ), pred_hi ( new vector<double> ), pred_link ( new vector<double> );
  Predict<LSBA> ( lsba_low, u_test, v_test, pred_low );
  Predict<LSBA> ( lsba_hi, u_test, v_test, pred_hi );
  Predict<LSBALink> ( lsba_link, u_test, v_test, pred_link );

  // Write output files
  ofstream fout;

  fout.precision ( 6 );
  fout << fixed;

  fout.open ( "lo_fi_data.asc" );
  for ( size_t i = 0; i < u_low->size(); ++i )
    fout << ( *u_low ) [i] << " " << ( *v_low ) [i] << " " << ( *z_low ) [i] << "\n";
  fout.close();

  fout.open ( "hi_fi_data.asc" );
  for ( size_t i = 0; i < u_high->size(); ++i )
    fout << ( *u_high ) [i] << " " << ( *v_high ) [i] << " " << ( *z_high ) [i] << "\n";
  fout.close();

  fout.open ( "lo_fi_prediction.asc" );
  for ( size_t i = 0; i < u_test->size(); ++i )
    fout << ( *u_test ) [i] << " " << ( *v_test ) [i] << " " << ( *pred_low ) [i] << "\n";
  fout.close();

  fout.open ( "hi_fi_prediction.asc" );
  for ( size_t i = 0; i < u_test->size(); ++i )
    fout << ( *u_test ) [i] << " " << ( *v_test ) [i] << " " << ( *pred_hi ) [i] << "\n";
  fout.close();

  fout.open ( "fusion_prediction.asc" );
  for ( size_t i = 0; i < u_test->size(); ++i )
    fout << ( *u_test ) [i] << " " << ( *v_test ) [i] << " " << ( *pred_link ) [i] << "\n";
  fout.close();

  return 0;
}