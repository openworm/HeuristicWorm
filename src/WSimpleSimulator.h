/*
 * WSimpleSimulator.h
 *  Created on: Feb 07, 2012
 *  Author: Alexander Dibert
 */

#ifndef WSimpleSimulator_h_
#define WSimpleSimulator_h_

#include "WSimulator.h"

/*
 * @class Simulator
 * Computes fitness function as a square of weighted distance between target vector and the one computed by solver.
 */

class WSimpleSimulator : public WSimulator
{
public:
	WSimpleSimulator( WSolver * iSolver, std::istream & iTrainDataStream, std::istream & iTestDataStream );
	virtual ~WSimpleSimulator();

private:
	WSimpleSimulator();

	void Init( WSolver * iSolver, std::istream & iTrainDataStream, std::istream & iTestDataStream ); // TODO make it work with data of arbitrary length
	float Evaluate( WTestData< float > const & iTarget );
	void Parse( WGenome< float > const & iGenome );

private:
	float  m_ICoefficients[ 20 ]; // TODO make it not hard-coded.
	float  m_UCoefficients[ 20 ];
};

#endif