/*
 * WDiffEvolution.cpp
 *
 *  Created on: Oct 9, 2011
 *      Author: Alexander Dibert
 */

#include <fstream>
#include <string>
#include <cmath>

#include "WDiffEvolution.h"
#include "WRandom.h"
#include "Helpers.h"

//-----------------------------------------------------------------------
WDiffEvolution::WDiffEvolution( WSimulator * iSimulator, bool iIsTestRun ): WHeuristicAlg< float, float >( iSimulator, iIsTestRun )
{
    m_prefix = "DE";
}

//-----------------------------------------------------------------------
WDiffEvolution::WDiffEvolution( std::istream & iInput, WSimulator * iSimulator, bool iIsTestRun ): WHeuristicAlg< float, float >( iInput, iSimulator, iIsTestRun )
{
	m_prefix = "DE";
}

//-----------------------------------------------------------------------
WDiffEvolution::~WDiffEvolution() {}

//-----------------------------------------------------------------------
void WDiffEvolution::Start() {
	if( !m_isTestRun ) {
		for( short i = 0; i < m_runNum; ++i )
		{
			m_testTrace.open( ( "./Results/" + m_startTimeStr + "_" + m_prefix + "_Test_Fitness_trace_" + Helpers::ITOA( m_iterationNum ) + "_" + Helpers::ITOA( m_populationSize ) + "_" + Helpers::ITOA( i ) +  ".txt" ).c_str() );   // Just for tracing
			m_fitnessTrace.open( ( "./Results/" + m_startTimeStr + "_" + m_prefix + "_Fitness_trace_" + Helpers::ITOA( m_iterationNum ) + "_" + Helpers::ITOA( m_populationSize ) + "_" + Helpers::ITOA( i ) +  ".txt" ).c_str() );     // Just for tracing

			this->Fitness( m_currPop );
			this->Evolve();
			this->Test();
			this->WriteData( i );

			m_currPop.Clear();
			m_nextPop.Clear();
			this->GeneratePopulation();

			m_testTrace.close();
			m_fitnessTrace.close();
		}
	}
	else {
		this->LoadPopulation();
		this->Test();
		std::cout << "Genome test fitness is " << m_currPop.GetGenome( 0 ).GetTestScore() << std::endl;
	}
}

//-----------------------------------------------------------------------
void WDiffEvolution::ReadData( std::istream & iInput ) {
	std::string fileName;
	std::string buffer;

	while( !iInput.eof() ) {
		std::getline( iInput, buffer );
		std::getline( iInput, fileName );

		if( !buffer.compare( "Genome Structure" ) )
		{
			std::ifstream input( fileName.c_str() );
			if( !input.is_open() ) {
				std::cerr << "ERROR: can't open " << fileName.c_str() << std::endl;
				exit( 0 );
			}

			std::string name;
			float max, min;
			short genomeLength = 0;

			std::getline( input, name, ' ' );
			input >> max >> min;
			while( !input.eof() ) {
				WGenome<float>::m_genomeStruct.AddParam( min, max, name );
				++genomeLength;
				std::getline( input, name );
				std::getline( input, name, ' ' );
				input >> max >> min;
			}
			m_genomeLength = genomeLength;
            input.close();
            continue;
		}

		if( !buffer.compare( "Genome Template" ) )
		{
			std::ifstream input( fileName.c_str() );
			if( !input.is_open() ) {
				std::cerr << "ERROR: can't open " << fileName.c_str() << std::endl;
				exit( 0 );
			}

			std::string symb, str;

			std::getline( input, str, '%' );
			std::getline( input, symb, '%' );
			while( !input.eof() ) {
				WGenome<float>::m_genomeStruct.AddTemplateLine( str, symb );
				std::getline( input, str, '%' );
				std::getline( input, symb, '%' );
			}
			if( !str.empty() )
				WGenome<float>::m_genomeStruct.m_templateStrings.push_back( str );

            input.close();
            continue;
		}

		if( !buffer.compare( "Main Constants" ) )
		{
			std::string buff;
			std::ifstream input( fileName.c_str() );
			if( !input.is_open()) {
				std::cerr << "ERROR: can't open " << fileName.c_str() << std::endl;
				exit( 0 );
			}

			std::getline( input, buff, ' ' );
			while( !input.eof() ) {
				if( !buff.compare( "Population_size" ) ) {
					input >> m_populationSize;
					std::getline( input, buff );
					std::getline( input, buff, ' ' );
					continue;
				}
				if( !buff.compare( "F" ) ) {
					input >> m_F; // belongs to [0,2]
					if(m_F > 2 || m_F < 0)
						m_F = 0.1;
					std::getline( input, buff );
					std::getline( input, buff, ' ' );
					continue;
				}
				if( !buff.compare( "CR" ) ) {
					input >> m_CR; // belongs to [0,1]
					if( m_CR > 1 || m_CR < 0)
						m_CR = 0.2;
					std::getline( input, buff );
					std::getline( input, buff, ' ' );
					continue;
				}
				if( !buff.compare( "Runs" ) ) {
					input >> m_runNum;
					std::getline( input, buff );
					std::getline( input, buff, ' ' );
					continue;
				}
                if( !buff.compare( "Output" ) ) {
					std::getline( input, m_output );
					std::getline( input, buff, ' ' );
                	continue;
                }
                if( !buff.compare( "Iteration_Limit" ) ) {
                	input >> m_iterationNum;
					std::getline( input, buff );
					std::getline( input, buff, ' ' );
                    continue;
                }
                if( !buff.compare( "Tolerated_Score" ) ) {
                   	input >> m_triggerError;
					std::getline( input, buff );
					std::getline( input, buff, ' ' );
                    continue;
                }
                std::getline( input, buff );
				std::getline( input, buff, ' ' );
			}
            input.close();
		}
	}
}

//-----------------------------------------------------------------------
void WDiffEvolution::Crossover() {
	srand( ( unsigned )time( NULL ) );
	int a, b, c;
// 	std::cout << "Mutation: start\n";
	for(int i = 0; i < m_populationSize; ++i)
	{
		this->SetF( i );
		//std::cout << "Mutation[" << i <<"]: start\n";
		do {
			a = rand() % m_populationSize;
		}
		while( a == i );
		do {
			b = rand() % m_populationSize;
		}
		while( b == i || b == a );
		do {
			c = rand() % m_populationSize;
		}
		while( c == i || c == a || c == b );

		//std::cout << "Mutation[" << i <<"]: indexes are generated\n";

		m_nextPop.AddGenome( m_currPop[c] + ( m_currPop[a] - m_currPop[b] ) * m_F );
		m_nextPop[i].Normalize();

		//std::cout << "Crossover[" << i <<"]: done\n";
	}
// 	std::cout << "Mutation: done\n";
}

//-----------------------------------------------------------------------
void WDiffEvolution::Mutation() {
// 	std::cout << "Crossover: start\n";
	srand( ( unsigned )time( NULL ) );
	for( int i = 0; i < m_populationSize; ++i ) {
		int forcedReplacementIndex = rand() % m_genomeLength;
		for( int j = 0; j < m_genomeLength; ++j ) {
			this->SetCR( i );
			if( ( float )rand() / RAND_MAX > m_CR && j != forcedReplacementIndex ) {
				m_nextPop[ i ].SetGene( j, ( m_currPop[ i ] ).GetGene( j ) );
			}
		}
	}
// 	std::cout << "Crossover: done\n";
}

//-----------------------------------------------------------------------
void WDiffEvolution::Selection()
{
	this->Fitness( m_nextPop );

// 	std::cout << "Selection: start\n";
	size_t bestNum = m_currPop.GetBestNum();
	double bestScore = m_currPop.GetBestScore();

	for( int i = 0; i < m_populationSize; ++i )
	{
		float const nextScore = m_nextPop[ i ].GetScore();
		if( nextScore < m_currPop[ i ].GetScore() ) {
			m_currPop[ i ] = m_nextPop[ i ];
			if( nextScore < bestScore ) {
				bestNum = i;
				bestScore = nextScore;
			}
		}
	}
	m_currPop.SetBestNum( bestNum );
	m_currPop.SetBestScore( m_currPop[ bestNum ].GetScore() );

// 	std::cout << "Selection: done\n";
	this->Test(); // Just for tracing
	m_testTrace << m_currPop.GetBestGene().GetTestScore() << "\n";
	m_fitnessTrace << m_currPop.GetBestScore() << "\n" ;
	m_testTrace.flush();
	m_fitnessTrace.flush();
    m_nextPop.Clear();
}
