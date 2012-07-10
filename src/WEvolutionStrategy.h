/*
 * WEvolutionStrategy.h
 *
 *  Created on: Apr 8, 2012
 *  Author: Alexander Dibert
 */

#ifndef WEvolutionStrategy_h
#define WEvolutionStrategy_h

const float nextPopulationSizeCoefficient = 7;

#include "WHeuristicAlg.hpp"

/*
 * Simple uncorrelated Evolution strategy class
 */

class WEvolutionStrategy : public WHeuristicAlg< float, float > {
public:
	WEvolutionStrategy( WSimulator * iSimulator, bool iIsTestRun = false );
	WEvolutionStrategy( std::istream & iInput, WSimulator * iSimulator, bool iIsTestRun = false );
	virtual ~WEvolutionStrategy();

	virtual void GeneratePopulation( WPopulation< float > * iPopulation = NULL );
	virtual void Start();

protected:
	virtual void ReadData( std::istream & );
	virtual void Init( std::istream & iInput );
	virtual void InitParams();

	virtual void Mutation();
	virtual void Crossover();
	virtual void Selection();

private:
	WEvolutionStrategy();

protected:
	// Evolution Strategy parameters.
	float m_learningRate;
	float m_mainLearningRate;
	std::vector< std::pair< double, int > > m_scoreIndex;

};

#endif /* WEvolutionStrategy_h */
