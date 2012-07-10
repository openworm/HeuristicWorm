/*
 *  WEvolutionStrategyModified.cpp
 *
 *  Created on: Apr 14, 2012
 *  Author: Alexander Dibert
 */

#include "WEvolutionStrategyModified.h"
#include "WRandom.h"

//-----------------------------------------------------------------------
WEvolutionStrategyModified::WEvolutionStrategyModified( WSimulator * iSimulator, bool iIsTestRun ): WEvolutionStrategy( iSimulator, iIsTestRun )
{
    m_prefix = "ESMod";
}

//-----------------------------------------------------------------------
WEvolutionStrategyModified::WEvolutionStrategyModified( std::istream & iInput, WSimulator * iSimulator, bool iIsTestRun ): WEvolutionStrategy( iInput, iSimulator, iIsTestRun )
{
	m_prefix = "ESMod";
}

//-----------------------------------------------------------------------
WEvolutionStrategyModified::~WEvolutionStrategyModified() {}

//-----------------------------------------------------------------------
void WEvolutionStrategyModified::Crossover()
{
	int f, s = 0;
	m_nextPop.Clear();
	::srand( ::time( NULL ) );
	for( int iter = 0; iter < ( nextPopulationSizeCoefficient - 1 ) * m_populationSize; ++iter )
	{
		WGenome< float > genome;
		genome.Resize( m_genomeLength * 2 );
		for( int g = 0; g < m_genomeLength; ++g )
		{
			f = ::rand() % m_populationSize;
			do {
				s = ::rand() % m_populationSize;

			} while( s == f );

			if( rand() % 2 ) {
				genome[ g ] = m_currPop[ f ].GetGene( g );
				genome[ g + m_genomeLength ] = m_currPop[ f ].GetGene( g + m_genomeLength );
			}
			else {
				genome[ g ] = m_currPop[ s ].GetGene( g );
				genome[ g + m_genomeLength ] = m_currPop[ s ].GetGene( g + m_genomeLength );
			}
		}
 		m_nextPop.AddGenome( genome );
	}
}

//-----------------------------------------------------------------------
void WEvolutionStrategyModified::Mutation()
{
	::WEvolutionStrategy::Mutation();
	for(int  i = 0; i < m_populationSize; ++i ) {
		m_nextPop.AddGenome( m_currPop.GetGenome( i ) );
	}
}