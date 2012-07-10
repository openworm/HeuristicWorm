/*
 * WDiffEvolution.h
 *
 *  Created on: Oct 8, 2011
 *  Author: Alexander Dibert
 */

#ifndef WDIFFEVOLUTION_H_
#define WDIFFEVOLUTION_H_

#include "WHeuristicAlg.hpp"

/*
 * Simple Differential Evolution class
 */

class WDiffEvolution : public WHeuristicAlg< float, float > {
public:
	WDiffEvolution( WSimulator * iSimulator, bool iIsTestRun = false );
	WDiffEvolution( std::istream & iInput, WSimulator * iSimulator, bool iIsTestRun = false );
	virtual ~WDiffEvolution();

	virtual void Start();

protected:
	virtual void ReadData( std::istream & );

	virtual void Mutation();
	virtual void Crossover();
	virtual void Selection();

	virtual void SetCR( int i ){} // used only by descendant
	virtual void SetF( int i ){} // used only by descendant

private:
	WDiffEvolution();

protected:
	// Differential Evolution parameters.
	float m_F;
	float m_CR;
};

#endif /* WDIFFEVOLUTION_H_ */
