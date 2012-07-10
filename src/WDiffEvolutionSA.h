/*
 * WDiffEvolutionSA.h
 *
 *  Created on: Oct 18, 2011
 *  Author: Alexander Dibert
 */

#ifndef WDIFFEVOLUTIONSA_H_
#define WDIFFEVOLUTIONSA_H_

#include "WDiffEvolution.h"

/*
 * Self-adaptive Differential Evolution class
 */

	class WDiffEvolutionSA : public WDiffEvolution {
	public:
		WDiffEvolutionSA( WSimulator * iSimulator, bool iTestRun = false );
		WDiffEvolutionSA( std::istream & iInput, WSimulator * iSimulator, bool iTestRun = false );
		void Start();
		void Init( std::istream & iInput);

	private:
		void Iteration();
		void Selection();
		void InitEvolutionParams();
		void EvolveEvolutionParams();
		void SetCR( int i );
		void SetF( int i );

	private:
		std::vector< std::pair< float, float >  > m_evolutionParams;
		std::vector< std::pair< float, float >  > m_newEvolutionParams;

		std::ofstream m_paramTrace;   // Just for tracing
	};

#endif
