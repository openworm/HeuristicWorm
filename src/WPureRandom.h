/*
 *  WPureRandom.h
 *  Created on: Apr 8, 2012
 *  Author: Alexander Dibert
 */

#ifndef WPureRandom_h
#define WPureRandom_h

#include "WHeuristicAlg.hpp"

/*
 * Picks up a new population by random guess.
 */

class WPureRandom : public WHeuristicAlg< float, float > {
public:
	WPureRandom( WSimulator * iSimulator, bool iIsTestRun = false );
	WPureRandom( std::istream & iInput, WSimulator * iSimulator, bool iIsTestRun = false );
	virtual ~WPureRandom();

	virtual void Start();

protected:
	virtual void ReadData( std::istream & );

	virtual void Mutation();
	virtual void Crossover() {};
	virtual void Selection();


private:
	WPureRandom();
};

#endif
