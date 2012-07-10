/*
 * WIHeuristicAlg.h
 *
 *  Created on: Oct 8, 2011
 *  Author: Alexander Dibert
 */

#ifndef WIHEURISTICALG_H_
#define WIHEURISTICALG_H_

#include <iostream>
#include <string>

#include "Helpers.h"
#include "WRandom.h"
#include "WPopulation.hpp"
#include "WTestData.h"
#include "WSimulator.h"

/*
 * Heuristic algorithm interface
 */

template< typename Type, typename TestType > class WHeuristicAlg {
public:
	//-----------------------------------------------------------------------
	WHeuristicAlg< Type, TestType >( WSimulator * iSimulator, bool iIsTestRun = false )
	{
		m_simulator = iSimulator;
		m_isTestRun = iIsTestRun;
	}

	//-----------------------------------------------------------------------
	WHeuristicAlg< Type, TestType >( std::istream & iInput, WSimulator * iSimulator, bool iIsTestRun )
	{
		m_simulator = iSimulator;
		this->Init( iInput );
		m_isTestRun = iIsTestRun;
	}

	//-----------------------------------------------------------------------
	virtual void Start() = 0;

	//-----------------------------------------------------------------------
	virtual void SerializePopulation( bool iuseTemplate = false ) {
		m_currPop.Serialize( "./CurrPopulation/Genome", iuseTemplate );
	}

	//-----------------------------------------------------------------------
	virtual void DeserializePopulation() {
		m_currPop.Resize( m_populationSize );
		m_currPop.Deserialize( "./CurrPopulation/Genome" );
	}

	//-----------------------------------------------------------------------
	virtual void Init( std::istream & iInput )
	{
		this->ReadData( iInput );

		m_startTime = ::time( NULL );
		m_startTimeStr= Helpers::ITOA( m_startTime );


		m_currPop.Reserve( m_populationSize );
		m_currPop.SetGenomeLength( m_genomeLength );
		m_nextPop.Reserve( m_populationSize );
		m_nextPop.SetGenomeLength( m_genomeLength );
	}

	//-----------------------------------------------------------------------
	virtual void GeneratePopulation( WPopulation< Type > * iPopulation = NULL )
	{
		if( iPopulation == NULL )
			iPopulation = &m_currPop;

		srand( ( unsigned )time( NULL ) );
		for( int i = 0; i < m_populationSize; ++i ) {
			WGenome< float > genome( m_genomeLength );
			for( int j = 0; j < m_genomeLength; ++j ) {
				genome.AddGene( WRand( WGenome< float >::m_genomeStruct.m_minMax[ 1 ][ j ], WGenome< float >::m_genomeStruct.m_minMax[ 0 ][ j ] ) );
			}
			iPopulation->AddGenome( genome );
		}
	}

	//-----------------------------------------------------------------------
	void Test() {
		this->Emulate( m_currPop[ m_currPop.GetBestNum() ], false );
	}

	//-----------------------------------------------------------------------
	void Fitness()
	{
		this->Fitness( m_currPop );
	}

protected:
	//-----------------------------------------------------------------------
	virtual void ReadData( std::istream & ) = 0;

	//-----------------------------------------------------------------------
	virtual void LoadPopulation()
	{
		m_populationSize = 1;
		float tmp;
		std::ifstream in("testGenome.txt");
		if( !in ) {
			std::cerr << "ERROR: Can't open test genome data\n";
			exit( 0 );
		}
		WGenome< float > genome( m_genomeLength );
		for(int j = 0; j < m_genomeLength; ++j){
			in >> tmp;
			genome.AddGene( tmp );
		}
		m_currPop.AddGenome( genome );
		in.close();
	}

	//-----------------------------------------------------------------------
	virtual void Evolve()
	{
	    int iterCount = 0;
		while ( ( ++iterCount <= m_iterationNum ) && ( m_currPop.GetBestScore() > m_triggerError ) ) {
			this->Iteration();
			std::cout << "Iteration number " << iterCount << " ; Fitness = " << m_currPop.GetBestScore() << "\n";
		}
	}

	//-----------------------------------------------------------------------
	virtual void Iteration()
	{
		this->Crossover();
		this->Mutation();
		this->Selection();
	}

	//-----------------------------------------------------------------------
	void Emulate( WGenome <float> & iGenome, bool iLearn )
	{
		if( iLearn )
			iGenome.SetScore( m_simulator->Train( iGenome ) );
		else
			iGenome.SetTestScore( m_simulator->Test( iGenome ) );
	}

	//-----------------------------------------------------------------------
	void Fitness( WPopulation< float > & iPopulation )
	{
// 		std::cout << "Fitness: start\n";
		for( unsigned int i = 0; i < iPopulation.GetSize(); ++i  )
			this->Emulate( iPopulation[ i ], true );
// 		std::cout << "Fitness: done\n";
	}

	//-----------------------------------------------------------------------
	virtual void Mutation() = 0;
	virtual void Crossover() = 0;
	virtual void Selection() = 0;


	//-----------------------------------------------------------------------
	void WriteData( short i ) {

		time_t timePassed = ::time( NULL ) - m_startTime;

		std::ofstream output( ( "./Results/" + m_startTimeStr + "_" + m_prefix + "_" + m_output + "_" + Helpers::ITOA( m_iterationNum ) + "_" + Helpers::ITOA( m_populationSize ) + "_" + Helpers::ITOA( i ) +  ".txt" ).c_str() );
		std::ofstream simpleOutput( ( "./Results/" + m_startTimeStr + "_" + m_prefix + "_" + m_output + "_Simple_" + Helpers::ITOA( m_iterationNum ) + "_" + Helpers::ITOA( m_populationSize ) + "_" + Helpers::ITOA( i ) + ".txt" ).c_str() );
		std::ofstream traceOutput( ( "./Results/" + m_startTimeStr + "_" + m_prefix + "_" + m_output + "_Trace_" + Helpers::ITOA( m_iterationNum ) + "_" + Helpers::ITOA( m_populationSize ) + "_" + Helpers::ITOA( i ) + ".txt" ).c_str() );
		std::ofstream testTraceOutput( ( "./Results/" + m_startTimeStr + "_" + m_prefix + "_" + m_output + "_TestTrace_" + Helpers::ITOA( m_iterationNum ) + "_" + Helpers::ITOA( m_populationSize ) + "_" + Helpers::ITOA( i ) + ".txt" ).c_str() );

		output << "Running time: \t" << timePassed << std::endl;
		output << "Best score: \t " << m_currPop.GetBestScore() << std::endl;
		output << "Test score: \t " << m_currPop.GetBestGene().GetTestScore() << std::endl << std::endl;
		for( int i = 0; i < m_genomeLength; ++i ) {
			output << WGenome<float>::m_genomeStruct.m_name[i] << " \t " << m_currPop.GetGenome( m_currPop.GetBestNum() ).GetGene( i ) << std::endl;
			simpleOutput << m_currPop.GetGenome( m_currPop.GetBestNum() ).GetGene( i ) << std::endl;
		}

		std::vector< float > & trace = m_simulator->GetSimulationResult( m_currPop.GetBestGene() );
		for( std::vector< float >::const_iterator iter = trace.begin(); iter != trace.end(); ++iter )
			traceOutput << *iter << std::endl;

		trace = m_simulator->GetTestSimulationResult( m_currPop.GetBestGene() );
		for( std::vector< float >::const_iterator iter = trace.begin(); iter != trace.end(); ++iter )
			testTraceOutput << *iter << std::endl;

		output.close();
		simpleOutput.close();
		traceOutput.close();
	}

protected:
	std::string m_startTimeStr;
	time_t m_startTime;

	bool m_isTestRun;

	int m_populationSize;
	int m_genomeLength;

	short m_runNum;
	int m_iterationNum;
	float m_triggerError;

	WPopulation< Type > m_currPop;
	WPopulation< Type > m_nextPop;
	std::string m_output;
	std::string m_prefix;        // Output files prefix
	std::ofstream m_testTrace;   // Just for tracing
	std::ofstream m_fitnessTrace;

	WSimulator * m_simulator;
};

#endif /* WIHEURISTICALG_H_ */
