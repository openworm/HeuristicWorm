/*
 *  WJordanSolver.cpp
 *  Created on: Oct 18, 2011
 *  Author: Alexander Dibert
 */
#include <cmath>
#include <iostream>
#include <fstream>

#include "WJordanSolver.h"

/*
IKs 0
IKf 1
ICa 2
Il 3
NN 0
PP 1
QQ 2
EE 3
FF 4
HH 5
*/
//-----------------------------------------------------------------------
WJordanSolver::~WJordanSolver() {}

//-----------------------------------------------------------------------
WJordanSolver::WJordanSolver() {}

//-----------------------------------------------------------------------
void WJordanSolver::Init( const WGenome< float > & iGenome )
{
	m_totalcost = 0;
	m_Istim[ 0 ] = -100e-12;
	m_Istim[ 1 ] = -400e-12;
	m_Istim[ 2 ] = -700e-12;
	m_Vhold[ 0 ] = 0;
	m_Vhold[ 1 ] = 20e-3;
	m_Vhold[ 2 ] = 40e-3;

	m_deltat = 0.1e-3;
	m_numpoints = floor( 0.05 / m_deltat );

	m_Cm = iGenome.GetGene( 0 );
	for( int i = 0; i < 4; ++i )
		m_gxs[ i ] = m_Cm * iGenome.GetGene( 1 + i );
	for( int i = 0; i < 4; ++i )
		m_vrevs[ i ] = iGenome.GetGene( 5 + i );
	for( int i = 0; i < 6; ++i )
		m_vhalf_x[ i ] = iGenome.GetGene( 9 + i );
	for( int i = 0; i < 6; ++i )
		m_k_x[ i ] = iGenome.GetGene( 15 + i );
	for( int i = 0; i < 5; ++i )
		m_T_x[ i ] = iGenome.GetGene( 21 + i );
	m_TCa = iGenome.GetGene( 26 );
	m_alphaCa = iGenome.GetGene( 27 );
	m_thiCa = 6.1e-6 / ( m_gxs[ 2 ] * m_TCa );
	m_vhalf_x[ 5 ] *= 1e-9;
	m_k_x[ 5 ] *= 1e-9;
}

//-----------------------------------------------------------------------
void WJordanSolver::Train( WGenome< float > const & iGenome, std::vector< float > & oResult )
{
	this->Init( iGenome );
	oResult.clear();

	// I100
	for( unsigned int i = 0; i < m_numpoints; ++i )
		m_I[ i ] = ( m_gxs[ 3 ] * ( 75e-3 + m_vrevs [3 ] ) );
	for( int unsigned i = m_numpoints / 5; i <= ( m_numpoints / 5 ) * 3; ++i )
		m_I[ i ] += m_Istim[ 0 ]; // '0' means I100

	this->IClamp();
	for( unsigned int i = ( int )floor( m_numpoints / 5); i < m_numpoints; i += 20 )
		oResult.push_back( m_V[ i ] );

	// I700
	for( unsigned int i = 0; i < m_numpoints; ++i )
		m_I[ i ] = ( m_gxs[ 3 ] * ( 75e-3 + m_vrevs[ 3 ] ) );
	for( int unsigned i = m_numpoints / 5; i <= ( m_numpoints / 5 ) * 3; ++i )
		m_I[ i ] += m_Istim[ 2 ]; // '2' means I700

	this->IClamp();
	for( unsigned int i = (int)floor( m_numpoints / 5); i < m_numpoints; i += 20 )
		oResult.push_back( m_V[i] );

	// V0
	for( unsigned int i = 0; i < m_numpoints; ++i )
		m_V[i] = -70e-3;
	for(int i = m_numpoints / 5; i <= 3 * (int)floor( m_numpoints / 5); ++i )
		m_V[i] = m_Vhold[0]; // '0' means V0

	this->VClamp();
	for( int i = (int)floor( m_numpoints / 5 ) + 10; i <= 3 * (int)floor( m_numpoints / 5 ); i += 10) //we add 2 to count so we get 20 pts
		oResult.push_back( m_I[i] );

	// V40
	for( unsigned int i = 0; i < m_numpoints; ++i )
		m_V[i] = -70e-3;
	for(int i = m_numpoints / 5; i <= 3 * (int)floor( m_numpoints / 5); ++i )
		m_V[i] = m_Vhold[2]; // '2' means V40

	this->VClamp();
	for( int i = (int)floor( m_numpoints / 5 ) + 10; i <= 3 * (int)floor( m_numpoints / 5 ); i += 10) //we add 2 to count so we get 20 pts
		oResult.push_back( m_I[i] );
}

//-----------------------------------------------------------------------
void WJordanSolver::Test( WGenome< float > const & iGenome, std::vector< float > & oResult )
{
	this->Init( iGenome );
	oResult.clear();

  // I400
	for( unsigned int i = 0; i < m_numpoints; ++i )
		m_I[ i ] = ( m_gxs[ 3 ] * ( 75e-3 + m_vrevs [3 ] ) );
	for( int unsigned i = m_numpoints / 5; i <= ( m_numpoints / 5 ) * 3; ++i )
		m_I[ i ] += m_Istim[ 1 ]; // '1' means I400

	this->IClamp();
	for( unsigned int i = ( int )floor( m_numpoints / 5); i < m_numpoints; i += 20 )
		oResult.push_back( m_V[ i ] );

	// V20
	for( unsigned int i = 0; i < m_numpoints; ++i )
		m_V[i] = -70e-3;
	for(int i = m_numpoints / 5; i <= 3 * (int)floor( m_numpoints / 5); ++i )
		m_V[i] = m_Vhold[1]; // '1' means V20

	this->VClamp();
	for( int i = (int)floor( m_numpoints / 5 ) + 10; i <= 3 * (int)floor( m_numpoints / 5 ); i += 10) //we add 2 to count so we get 20 pts
		oResult.push_back( m_I[i] );
}

//----------------------------------------------------------
/*///////////////////////Counters/////////////////////////*/

//-----------------------------------------------------------------------
void WJordanSolver::IClamp()
{
		m_V[0] = -0.075;
		m_Ca[0] = 0;
		for(int i = 0; i <= 4; ++i)
			m_gatevars[i] = 1 / ( 1 + exp( ( m_vhalf_x[i] - m_V[0] ) / m_k_x[i] ) );

		m_gatevars[5] = 1 / ( 1 + exp( ( m_vhalf_x[5] - m_Ca[0] ) / m_k_x[5] ) );

		// Now do main loop
		for( unsigned  int i = 1; i < m_numpoints; ++i )
		{
			float k1[7], k2[7], k3[7], k4[7];
			float h_2 = m_deltat / 2;
			float tempgate[6];
			float tempCa, tempV;

			//step 1 - get k1 values
			for(int j = 0; j <= 4; ++j){
				k1[j] = dx(m_gatevars[j], m_V[i-1], m_T_x[j], m_vhalf_x[j], m_k_x[j]);
			}
			k1[5] = dCa( m_gxs, m_gatevars, m_vrevs, m_alphaCa, m_V[i-1], m_TCa, m_thiCa, m_Ca[i-1] );
			k1[6] = dV( m_gxs, m_gatevars, m_vrevs, m_alphaCa, m_I[i], m_Cm, m_V[i-1]);

			//step 2  - get k2 values
			tempV = m_V[i-1] + k1[6] * h_2;
			tempCa = m_Ca[i-1] + k1[5] * h_2;
			for(int j = 0; j <= 4; ++j){
				tempgate[j] =  (m_gatevars[j] + k1[j] * h_2);
				k2[j] = dx(tempgate[j], tempV, m_T_x[j], m_vhalf_x[j], m_k_x[j]);
			}
			tempgate[5] = 1 / ( 1 + exp( ( m_vhalf_x[5] - tempCa) / m_k_x[5] ) );
			k2[5] = dCa( m_gxs, tempgate, m_vrevs, m_alphaCa, tempV, m_TCa, m_thiCa, tempCa );
			k2[6] = dV (m_gxs, tempgate, m_vrevs, m_alphaCa, m_I[i], m_Cm, tempV );

			//step 3  - get k3 values
			tempV = m_V[i-1] + k2[6] * h_2;
			tempCa = m_Ca[i-1] + k2[5] * h_2;
			for(int j = 0; j <= 4; ++j){
				tempgate[j] =  (m_gatevars[j] + k2[j] * h_2);
				k3[j] = dx(tempgate[j], tempV, m_T_x[j], m_vhalf_x[j], m_k_x[j]);
			}
			tempgate[5] = 1 / ( 1 + exp( ( m_vhalf_x[5] - tempCa ) / m_k_x[5] ) );
			k3[5] = dCa( m_gxs, tempgate, m_vrevs, m_alphaCa, tempV, m_TCa, m_thiCa, tempCa );
			k3[6] = dV( m_gxs, tempgate, m_vrevs, m_alphaCa, m_I[i], m_Cm, tempV );

			//step 4  - get k4 values
			tempV = m_V[i-1] + k3[6] * h_2 * 2;
			tempCa = m_Ca[i-1] + k3[5] * h_2 * 2;
			for(int j = 0; j <= 4; ++j){
				tempgate[j] =  ( m_gatevars[j] + k3[j] * h_2 * 2 );
				k4[j] = dx(tempgate[j], tempV, m_T_x[j], m_vhalf_x[j], m_k_x[j]);
			}
			tempgate[5] = 1 / ( 1 + exp( ( m_vhalf_x[5] - tempCa ) / m_k_x[5] ) );
			k4[5] = dCa( m_gxs, tempgate, m_vrevs, m_alphaCa, tempV, m_TCa, m_thiCa, tempCa );
			k4[6] = dV( m_gxs, tempgate, m_vrevs, m_alphaCa, m_I[i], m_Cm, tempV );

			// now we must combine them to get final value
			for(int j = 0; j <= 4; ++j){
				m_gatevars[j] += ( m_deltat / 6 ) * ( k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j] );
			}
			m_Ca[i] = m_Ca[i-1] + ( m_deltat / 6 ) * ( k1[5] + 2 * k2[5] + 2 * k3[5] + k4 [5] );
			m_gatevars[5] = 1 / ( 1 + exp( ( m_vhalf_x[5] - tempCa ) / m_k_x[5] ) );
			m_V[i] = m_V[i-1] + ( m_deltat / 6 ) * ( k1[6] + 2 * k2[6] + 2 * k3[6] + k4[6] );
		}
		return;
}

//-----------------------------------------------------------------------
void WJordanSolver::VClamp()
{
	// First inita3aise gating variables to steady state values
	m_Ca[0] = 0;
	for( int i = 0; i <= 4; ++i ) {
		m_gatevars[i] = 1 / ( 1 + exp( ( m_vhalf_x[i] - m_V[0] ) / m_k_x[i] ) );
	}
	m_gatevars[5] = 1 / ( 1 + exp( ( m_vhalf_x[5] - m_Ca[0] ) / m_k_x[5] ) );
	m_I[0] = m_gxs[0] * m_gatevars[0] * ( m_V[0] - m_vrevs[0] ) + m_gxs[1] * pow( m_gatevars[1], 4 ) * m_gatevars[2] *
			( m_V[0] - m_vrevs[1] ) + m_gxs[2] * pow( m_gatevars[3], 2 ) * m_gatevars[4] * ( 1 + ( m_gatevars[5] - 1 ) * m_alphaCa ) *
			( m_V[0]-m_vrevs[2] ) + m_gxs[3] * ( m_V[0] - m_vrevs[3] );

	// Now do main loop
	for( unsigned int i = 1; i < m_numpoints; ++i ) {
		float k1[6], k2[6], k3[6], k4[6];
		float h_2 = m_deltat / 2;
		float tempgate[6];
		float tempCa;

		//step 1 - get k1 values
		for( int j = 0; j <= 4; ++j ){
			k1[j] = dx( m_gatevars[j], m_V[i], m_T_x[j], m_vhalf_x[j], m_k_x[j] );
		}
		k1[5] = dCa( m_gxs, m_gatevars, m_vrevs, m_alphaCa, m_V[i], m_TCa, m_thiCa, m_Ca[i-1] );

		//step 2  - get k2 values
		tempCa = m_Ca[i-1] + k1[5] * h_2;
		for( int j = 0; j <= 4; ++j ){
			tempgate[j] =  ( m_gatevars[j] + k1[j] * h_2 );
			k2[j] = dx( tempgate[j], m_V[i], m_T_x[j], m_vhalf_x[j], m_k_x[j] );
		}
		tempgate[5] = 1 / ( 1 + exp( ( m_vhalf_x[5] - tempCa ) / m_k_x[5] ) );
		k2[5] = dCa( m_gxs, tempgate, m_vrevs, m_alphaCa, m_V[i], m_TCa, m_thiCa, tempCa );

		//step 3  - get k3 values
		tempCa = m_Ca[i-1] + k2[5] * h_2;
		for( int j = 0; j <= 4; ++j ){
			tempgate[j] =  ( m_gatevars[j] + k2[j] * h_2 );
			k3[j] = dx( tempgate[j], m_V[i], m_T_x[j], m_vhalf_x[j], m_k_x[j] );
		}
		tempgate[5] = 1 / ( 1 + exp( ( m_vhalf_x[5] - tempCa ) / m_k_x[5] ) );
		k3[5] = dCa( m_gxs, tempgate, m_vrevs, m_alphaCa, m_V[i], m_TCa, m_thiCa, tempCa );

		//step 4  - get k4 values
		tempCa = m_Ca[i-1] + k3[5] * h_2 * 2;
		for( int j = 0; j <= 4; ++j ){
			tempgate[j] =  ( m_gatevars[j] + k3[j] * h_2 * 2 );
			k4[j] = dx( tempgate[j], m_V[i], m_T_x[j], m_vhalf_x[j], m_k_x[j] );
		}
		tempgate[5] = 1 / ( 1 + exp( ( m_vhalf_x[5] - tempCa ) / m_k_x[5] ) );
		k4[5] = dCa( m_gxs, tempgate, m_vrevs, m_alphaCa, m_V[i], m_TCa, m_thiCa, tempCa );

		// now we must combine them to get final value
		for( int j = 0; j <= 4; ++j ){
			m_gatevars[j] += ( m_deltat / 6 ) * ( k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j] );
		}
		m_Ca[i] = m_Ca[i-1] + ( m_deltat / 6 ) * ( k1[5] + 2 * k2[5] + 2 * k3[5] + k4[5] );

		m_gatevars[5] = 1 / ( 1 + exp( ( m_vhalf_x[5] - m_Ca[i] ) / m_k_x[5] ) );
		m_I[i] = m_gxs[0] * m_gatevars[0] * ( m_V[i] - m_vrevs[0] ) + m_gxs[1] * pow( m_gatevars[1], 4 ) * m_gatevars[2] *
				( m_V[i] - m_vrevs[1] ) + m_gxs[2] * pow( m_gatevars[3], 2 ) * m_gatevars[4] * ( 1 + ( m_gatevars[5] - 1 ) * m_alphaCa ) *
				( m_V[i]-m_vrevs[2] );// + m_gxs[3]*(m_V[i] - m_vrevs[3]);  NOTE: 3 removed for consistency with Jospin

		//if(isnan(m_I[i])) cerr << "Vclamp has returned a NAN";
	}
}

//-----------------------------------------------------------------------
float WJordanSolver::IVCost( WTestData<float> & iData )
{
	float cost = 0;
	float n_inf, p_inf, q_inf, e_inf, f_inf, h_inf;
	float K_ref, Ca_ref;
	float K_points[10], Ca_points[10];
	float V;

	for( int i = 0; i < 10; ++i ) {
		V = -50e-3 + i * 10e-3;
		n_inf = 1 / ( 1 + exp( ( m_vhalf_x[0] - V ) / m_k_x[0] ) );
		p_inf = 1 / ( 1 + exp( ( m_vhalf_x[1] - V ) / m_k_x[1] ) );
		q_inf = 1 / ( 1 + exp( ( m_vhalf_x[2] - V ) / m_k_x[2] ) );
		K_points[i] = m_gxs[1] * pow( p_inf, 4 ) * q_inf * ( V - m_vrevs[1] ) + m_gxs[0] * n_inf * ( V - m_vrevs[0] );
		K_points[i] /= m_Cm;

		//NOTE: in order to avoid problems with Ca2+ inactivation, I have decided to use the peak IV relationship for ICa
		V = -30e-3 + i * 10e-3;
		e_inf = 1 / ( 1 + exp( ( m_vhalf_x[3] - V ) / m_k_x[3] ) );
		f_inf = 1;
		h_inf = 1;

		Ca_points[i] = m_gxs[2] * pow( e_inf, 2 ) * f_inf * h_inf * ( V - m_vrevs[2] );
		Ca_points[i] /= m_Cm;
	}

	for( int i = 0; i < 10 ; ++i ) {
		K_ref = iData.Get( 6, i );
		cost +=	pow( ( K_ref - K_points[i] ) *4,2) * 50;
	}

	for( int i = 0; i < 10 ; ++i ) {
		Ca_ref = iData.Get( 7, i );
		cost +=	pow( ( Ca_ref - Ca_points[i] ) * 4,2 ) *100;
	}
	if(isnan( cost) ) {
		//cerr << "IVcost.hh has returned NAN!";
		cost = 1e60;
	}
	//cerr << "\n Good cost! = " << cost;
	return cost;
}

//----------------------------------------------------------
/*///////////////////////Derivatives//////////////////////*/


//-----------------------------------------------------------------------
float WJordanSolver::dx( float iX, float iV, float iT, float iVhalf, float iK ) const {
	float x_inf = 1 / (1 + exp( ( iVhalf - iV ) / iK ) );
	return ( x_inf - iX ) / iT;
}

//-----------------------------------------------------------------------
float WJordanSolver::dV(float iGxs[4], float iGatevars[6], float iVrevs[4], float iAlphaCa, float iIin, float iCmem, float iV) const {
	float Imem = iGxs[0] * iGatevars[0] * ( iV - iVrevs[0]) + iGxs[1] * pow( iGatevars[1], 4 ) * iGatevars[2] *
			( iV - iVrevs[1] ) + iGxs[2] * pow( iGatevars[3], 2) * iGatevars[4] * ( 1 + ( iGatevars[5] - 1 ) * iAlphaCa ) *
			( iV - iVrevs[2] ) + iGxs[3] * ( iV - iVrevs[3] );
	return -( iIin + Imem ) / iCmem;
}

//-----------------------------------------------------------------------
float WJordanSolver::dCa (float iGxs[4], float iGatevars[6], float iVrevs[4], float iAlphaCa,  float iV, float iTCa, float iThiCa, float iCa) const {
	float Ca_flow = iGxs[2] * pow( iGatevars[3], 2 ) * iGatevars[4] * ( 1 + ( iGatevars[5] - 1 ) * iAlphaCa ) * ( iV - iVrevs[2] );
	return -( iCa / iTCa + Ca_flow * iThiCa );
}
