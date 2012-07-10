/*
 *  WEvolutionStrategy.cpp
 *
 *  Created on: Apr 9, 2012
 *  Author: Alexander Dibert
 */

#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>

#include "WEvolutionStrategy.h"
#include "WRandom.h"
#include "Helpers.h"

//-----------------------------------------------------------------------
WEvolutionStrategy::WEvolutionStrategy( WSimulator * iSimulator, bool iIsTestRun ): WHeuristicAlg< float, float >( iSimulator, iIsTestRun )
{
    m_prefix = "ES";
}

//-----------------------------------------------------------------------
WEvolutionStrategy::WEvolutionStrategy( std::istream & iInput, WSimulator * iSimulator, bool iIsTestRun ): WHeuristicAlg< float, float >( iInput, iSimulator, iIsTestRun )
{
	m_prefix = "ES";
}

//-----------------------------------------------------------------------
WEvolutionStrategy::~WEvolutionStrategy() {}

//-----------------------------------------------------------------------
void WEvolutionStrategy::Start() {
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
void WEvolutionStrategy::GeneratePopulation( WPopulation< float > * iPopulation )
{
	if( iPopulation == NULL )
		iPopulation = &m_currPop;

	WHeuristicAlg<float, float>::GeneratePopulation( iPopulation );
	for( int i = 0; i < m_populationSize; ++i )
	{
		WGenome< float > & gene = ( *iPopulation )[ i ];
		for( int j = 0; j < m_genomeLength; ++j )
			gene.AddGene( WGenome< float >::m_genomeStruct.m_minMax[ 2 ][ j ] / 64 );
	}
}

//-----------------------------------------------------------------------
void WEvolutionStrategy::ReadData( std::istream & iInput ) { //TODO Bad style!
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
void WEvolutionStrategy::Init( std::istream & iInput )
{
	WHeuristicAlg< float, float >::Init( iInput );
	this->InitParams();
	m_scoreIndex.resize( nextPopulationSizeCoefficient * m_populationSize );
}

//-----------------------------------------------------------------------
void WEvolutionStrategy::InitParams()
{
	m_learningRate = 1 / sqrt( 2 * sqrt( m_genomeLength ) );
	m_mainLearningRate = 1 / ::sqrt( 2 * m_genomeLength );
}

//-----------------------------------------------------------------------
void WEvolutionStrategy::Mutation()
{
	for( int g = 0; g < m_nextPop.GetSize(); ++g )
	{
		for( int i = 0; i < m_genomeLength; ++i )
		{
			float & value = this->m_nextPop[ g ][ m_genomeLength + i ];
			value *= exp( m_mainLearningRate * WNormRand( 1 ) * m_learningRate * WNormRand( 1 ) ); // TODO Optimize
			if( value <  WGenome< float >::m_genomeStruct.m_minMax[ 2 ][ i ] / 256 )
				value =  WGenome< float >::m_genomeStruct.m_minMax[ 2 ][ i ] / 256;
			m_nextPop[ g ][ i ] +=  value * WNormRand( 1 );
		}
		m_nextPop[ g ].Normalize();
	}
}

//-----------------------------------------------------------------------
void WEvolutionStrategy::Crossover()
{
	int f, s = 0;
	m_nextPop.Clear();
	::srand( ::time( NULL ) );
	for( int iter = 0; iter < ( nextPopulationSizeCoefficient ) * m_populationSize; ++iter )
	{
		f = ::rand() % m_populationSize;
		do {
			s = ::rand() % m_populationSize;
		} while( s == f );

 		WGenome< float > genome;
		genome.Resize( 2 * m_genomeLength );
 		for( int i = 0; i < m_genomeLength; ++i )
 		{
			bool rnd = ::rand() % 2;
			if( rnd )
			{
				genome[ i ] = m_currPop[ f ].GetGene( i );
				genome[ i + m_genomeLength ] = m_currPop[ f ].GetGene( i + m_genomeLength );
			}
			else
			{
				genome[ i ] = m_currPop[ s ].GetGene( i );
				genome[ i + m_genomeLength ] = m_currPop[ s ].GetGene( i + m_genomeLength );
			}
 		}
 		m_nextPop.AddGenome( genome );
	}
}

//-----------------------------------------------------------------------
void WEvolutionStrategy::Selection()
{
	this->Fitness( m_nextPop );

//  	std::cout << "Selection: start\n";

	for( int i = 0; i < nextPopulationSizeCoefficient * m_populationSize; ++i )
		m_scoreIndex[ i ] = std::make_pair< double, int >( m_nextPop[ i ].GetScore(), i );

	std::sort( m_scoreIndex.begin(), m_scoreIndex.end() );

	for( int i = 0; i < m_genomeLength; ++i )
	{
		m_currPop[ i ] = m_nextPop[ m_scoreIndex[ i ].second ];
	}

	m_currPop.SetBestNum( 0 );
	m_currPop.SetBestScore( m_scoreIndex[ 0 ].first );

// 	std::cout << "Selection: done\n";
	this->Test(); // Just for tracing
	m_testTrace << m_currPop.GetBestGene().GetTestScore() << "\n";
	m_fitnessTrace << m_currPop.GetBestScore() << "\n" ;
	m_testTrace.flush();
	m_fitnessTrace.flush();
    m_nextPop.Clear();
}