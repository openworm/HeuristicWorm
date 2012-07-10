/*
 *  WSimulator.h
 *  Created on: Oct 12, 2011
 *  Author: Alexander Dibert
 */

#ifndef WSimulator_h_
#define WSimulator_h_

#include "WPopulation.hpp"
#include "WTestData.h"
#include "WSolver.h"

/*
 * Simulator Interface
 */
class WSimulator
{
public:
	virtual ~WSimulator(){};

	virtual void Init( WSolver * iSolver, std::istream & iTrainDataStream, std::istream & iTestDataStream );
	virtual float Train( WGenome< float > const & iGenome );
	virtual float Test( WGenome< float > const & iGenome );
	virtual std::vector< float > & GetSimulationResult( WGenome< float > const & iGenome );
	virtual std::vector< float > & GetTestSimulationResult( WGenome< float > const & iGenome );

private:
	virtual float Evaluate( WTestData< float > const & iTarget )  = 0; //TODO make it const
	virtual void Parse( WGenome< float > const & iGenome ) = 0;

protected:
	WSolver * m_solver;

	WTestData< float > m_trainData;
	WTestData< float > m_testData;
	std::vector< float > m_result;

};

#endif /* WSimulator_h_ */
