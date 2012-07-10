/*
 * WEvolutionStrategyModified.h
 *
 *  Created on: Apr 8, 2012
 *  Author: Alexander Dibert
 */

#ifndef WEvolutionStrategyModified_h
#define WEvolutionStrategyModified_h

#include "WEvolutionStrategy.h"

/*
 * Modified uncorrelated Evolution strategy class: not muated previous population is included into selection pool.
 */

class WEvolutionStrategyModified: public WEvolutionStrategy
{
public:
	WEvolutionStrategyModified( WSimulator * iSimulator, bool iIsTestRun = false );
	WEvolutionStrategyModified( std::istream & iInput, WSimulator * iSimulator, bool iIsTestRun = false );
	virtual ~WEvolutionStrategyModified();

protected:
	virtual void Crossover();
	virtual void Mutation();

private:
	WEvolutionStrategyModified();
};

#endif /* WEvolutionStrategyModified_h */