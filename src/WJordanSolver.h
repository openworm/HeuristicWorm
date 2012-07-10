/*
 * WJordanSolver.h
 *  Created on: Feb 07, 2012
 *  Author: Alexander Dibert
 *  Uses code by Jordan Boyle
 */

#ifndef WJordanSolver_H_
#define WJordanSolver_H_

#include "WSolver.h"
#include "WTestData.h"

/*
 * Solver based on code by Jordan Boyle.
 * Solves simple HH equations.
 */

class WJordanSolver : public WSolver
{
public:
	WJordanSolver();
	virtual ~WJordanSolver();
	void Init( const WGenome< float > & iGenome );
	void Train( WGenome< float > const & iGenome, std::vector< float > & oResult );
	void Test( WGenome< float > const & iGenome, std::vector< float > & oResult );

protected:

	//Counters
	void IClamp();
	void VClamp();
	float IVCost( WTestData<float> & iData );
	float CountScore(int iTestNum, WTestData<float> & iData, bool iLearn );

	//Derivatives
	float dx( float iX, float iV, float iT, float iVhalf, float iK ) const;
	float dV( float iGxs[4], float iGatevars[6], float iVrevs[4], float iAlphaCa, float iIin, float iCmem, float iV ) const;
	float dCa (float iGxs[4], float iGatevars[6], float iVrevs[4], float iAlphaCa,  float iV, float TiCa, float iThiCa, float iCa ) const;

protected:
	double m_totalcost;
	float  m_Istim[3];
	float  m_Vhold[3];
	float  m_deltat;
	unsigned int    m_numpoints;
	float  m_Cm;
	float  m_gxs[4];
	float  m_vrevs[4];
	float  m_vhalf_x[6];
	float  m_k_x[6];
	float  m_T_x[5];
	float  m_TCa;
	float  m_alphaCa;
	float  m_thiCa;
	float  m_I[500];
	float  m_V[500];
	float  m_Ca[500];
	float  m_gatevars[6];
};

#endif
