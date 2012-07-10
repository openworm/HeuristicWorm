/*
 * WDiffEvolutionSA.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: shamdor
 */

#include <fstream>
#include <cstdlib>

#include "WDiffEvolutionSA.h"
#include "Helpers.h"

namespace DiffEvolutionSA {
	const float Fl = 0.1;
	const float Fu = 0.9;
	const float tao1 = 0.1;
	const float tao2 = 0.1;
}

//-----------------------------------------------------------------------
WDiffEvolutionSA::WDiffEvolutionSA( WSimulator * iSimulator, bool iIsTestRun ) : WDiffEvolution( iSimulator, iIsTestRun )
{
	m_prefix = "DESA";
}

//-----------------------------------------------------------------------
WDiffEvolutionSA::WDiffEvolutionSA( std::istream & iInput, WSimulator * iSimulator, bool iIsTestRun ) : WDiffEvolution( iInput, iSimulator, iIsTestRun )
{
	m_prefix = "DESA";
}

//-----------------------------------------------------------------------
void WDiffEvolutionSA::InitEvolutionParams()
{
	m_evolutionParams.clear();
	m_newEvolutionParams.clear();
	for( int i = 0 ; i < m_populationSize; ++i ) {
		m_evolutionParams.push_back( std::pair< float, float >( DiffEvolutionSA::Fl, m_CR )  );
		m_newEvolutionParams.push_back( std::pair< float, float >( 0, 0 )  );
	}
}

//-----------------------------------------------------------------------
void WDiffEvolutionSA::Start() {
	if( !m_isTestRun ) {
		for( short i = 0; i < m_runNum; ++i )
		{
			m_testTrace.open( ( "./Results/" + m_startTimeStr + "_" + m_prefix + "_Test_Fitness_trace_" + Helpers::ITOA( m_iterationNum ) + "_" + Helpers::ITOA( m_populationSize ) + "_" + Helpers::ITOA( i ) +  ".txt" ).c_str() );   // Just for tracing
			m_fitnessTrace.open( ( "./Results/" + m_startTimeStr + "_" + m_prefix + "_Fitness_trace_" + Helpers::ITOA( m_iterationNum ) + "_" + Helpers::ITOA( m_populationSize ) + "_" + Helpers::ITOA( i ) +  ".txt" ).c_str() );     // Just for tracing
			m_paramTrace.open( ( "./Results/" + m_startTimeStr + "_" + m_prefix + "_Param_trace_" + Helpers::ITOA( m_iterationNum ) + "_" + Helpers::ITOA( m_populationSize ) + "_" + Helpers::ITOA( i ) +  ".txt" ).c_str() );         // Just for tracing

			this->Fitness( m_currPop );
			this->Evolve();
			this->Test();
			this->WriteData( i );

			m_currPop.Clear();
			m_nextPop.Clear();
			this->GeneratePopulation();

			m_testTrace.close();
			m_fitnessTrace.close();
			m_paramTrace.close();
		}
	}
	else {
		this->LoadPopulation();
		this->Test();
		std::cout << "Genome test fitness is " << m_currPop.GetGenome( 0 ).GetTestScore() << std::endl;
	}
}

//-----------------------------------------------------------------------
void WDiffEvolutionSA::Init( std::istream & iInput ) {
	::WDiffEvolution::Init( iInput );
	this->InitEvolutionParams();
}

//-----------------------------------------------------------------------
void WDiffEvolutionSA::Iteration() {
	this->EvolveEvolutionParams();
	WDiffEvolution::Iteration();
}

//-----------------------------------------------------------------------
void WDiffEvolutionSA::EvolveEvolutionParams() {
	::srand( ::time( NULL ) );

	for( int i = 0; i < m_populationSize; ++i ) {
		if( ( float )rand() / RAND_MAX < DiffEvolutionSA::tao1 )
			m_newEvolutionParams[ i ].first = DiffEvolutionSA::Fl + DiffEvolutionSA::Fu * ( float )rand() / RAND_MAX;
		else
			m_newEvolutionParams[ i ].first = m_evolutionParams[ i ].first;
		if( ( float )rand() / RAND_MAX < DiffEvolutionSA::tao2 )
			m_newEvolutionParams[ i ].second = ( float )rand() / RAND_MAX;
		else
			m_newEvolutionParams[ i ].second = m_evolutionParams[ i ].second;
	}
}

//-----------------------------------------------------------------------
void WDiffEvolutionSA::SetCR( int i ) {
	m_CR = m_newEvolutionParams[ i ].second;
}

//-----------------------------------------------------------------------
void WDiffEvolutionSA::SetF( int i ) {
	m_F = m_newEvolutionParams[ i ].first;
}

//-----------------------------------------------------------------------
void WDiffEvolutionSA::Selection() {
	this->Fitness( m_nextPop );

// 	std::cout << "Selection: start\n";
	size_t bestNum = m_currPop.GetBestNum();
	double bestScore = m_currPop.GetBestScore();

	for( int i = 0; i < m_populationSize; ++i )
	{
		float const nextScore = m_nextPop[ i ].GetScore();
		if( nextScore < m_currPop[ i ].GetScore() ) {
			m_currPop[ i ] = m_nextPop[ i ];
			m_evolutionParams[ i ] = m_newEvolutionParams[ i ];
			if( nextScore < bestScore ) {
				bestNum = i;
				bestScore = nextScore;
			}
		}
	}
	m_currPop.SetBestNum( bestNum );
	m_currPop.SetBestScore( m_currPop[ bestNum ].GetScore() );

	std::cout << "Selection: done; Best F = " << m_evolutionParams[m_currPop.GetBestNum()].first << "; Best CR = " << m_evolutionParams[m_currPop.GetBestNum()].second << "\n";
	this->Test(); // Just for tracing
	m_paramTrace << m_evolutionParams[m_currPop.GetBestNum()].first << "\n" ;
	m_testTrace << m_currPop.GetBestGene().GetTestScore() << "\n";
	m_fitnessTrace << m_currPop.GetBestScore() << "\n" ;
	m_paramTrace.flush();
	m_testTrace.flush();
	m_fitnessTrace.flush();
    m_nextPop.Clear();
}