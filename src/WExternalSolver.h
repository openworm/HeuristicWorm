/*
 *  WExternalSolver.h
 *  Created on: Jan 18, 2011
 *  Author: Alexander Dibert
 */

#ifndef WExternalSolver_h
#define WExternalSolver_h

#include "WSolver.h"

/*
 * Handles a connection with external simulator. // TODO Think about it!
 */
class WExternalSolver : public WSolver
{
public:
	WExternalSolver( std::istream & iTrainStream, std::istream & iTestStream );
	void Init( WGenome< float > const & iGenome );

	float Train( WGenome< float > const & iGenome );
	float Test( WGenome< float > const & iGenome );

protected:
	void Parse( WGenome< float > const & iGenome );
};

#endif
