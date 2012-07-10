/*
 *  WLeMassonSimulator.cpp
 *  Created on: Mar 29, 2012
 *  Author: Alexander Dibert
 */

#include "WLeMassonSimulator.h"

//----------------------------------------------------------------------
WLeMassonSimulator::WLeMassonSimulator( WSolver * iSolver, std::istream & iTrainDataStream, std::istream & iTestDataStream )
	: m_I100Score( "I" )
	, m_I400Score( "I" )
	, m_I700Score( "I" )
	, m_V0Score( "V" )
	, m_V20Score( "V" )
	, m_V40Score( "V" )
{
	this->Init( iSolver, iTrainDataStream, iTestDataStream );
}

//----------------------------------------------------------------------
void WLeMassonSimulator::Init( WSolver * iSolver, std::istream & iTrainDataStream, std::istream & iTestDataStream )
{
	m_deltaT = 10e-3; // TODO Do we need it?
	WSimulator::Init( iSolver, iTrainDataStream, iTestDataStream );
	m_I100Score.Init( WSimulator::m_trainData.Get( 0 ), m_deltaT );
	m_I400Score.Init( WSimulator::m_testData.Get( 0 ), m_deltaT );
	m_I700Score.Init( WSimulator::m_trainData.Get( 1 ), m_deltaT );
	m_V0Score.Init( WSimulator::m_trainData.Get( 2 ), m_deltaT, 1000, 0.000001 );
	m_V20Score.Init( WSimulator::m_testData.Get( 1 ), m_deltaT, 1000, 0.000001 );
	m_V40Score.Init( WSimulator::m_trainData.Get( 3 ), m_deltaT, 1000, 0.000001 );
}

//----------------------------------------------------------------------
float WLeMassonSimulator::Evaluate( WTestData< float > const & iTarget )  //TODO Make it MORE flexible.
{
	double result = 0;
	float prevValue, currValue;
	unsigned int j;
	int passedValuesCounter = 0;
	WScoreContainer * container;
	std::string name;
	for( unsigned int i = 0; i < iTarget.GetSize(); ++i  )
	{
		currValue = m_result[ passedValuesCounter ];
		WTest< float > const & test = iTarget.Get( i );

		name = test.GetName();
		if( !name.compare( "I100" ) ) // TODO bad style
			container = &m_I100Score;
		else if( !name.compare( "I400" ) )
			container = &m_I400Score;
		else if( !name.compare( "I700" ) )
			container = &m_I700Score;
		else if( !name.compare( "V0" ) )
			container = &m_V0Score;
		else if( !name.compare( "V20" ) )
			container = &m_V20Score;
		else if( !name.compare( "V40" ) )
			container = &m_V40Score;
		else {
			std::cerr << "ERROR: Wrong test type. Aborting.";
			::exit( 0 );
		}

		container->ResetScore();
		for( j = 1; j < test.GetSize(); ++j ) {
				prevValue = currValue;
				currValue = m_result[ passedValuesCounter + j ];
				container->UpdateScore( prevValue, currValue, m_deltaT * 25 );
			}
		passedValuesCounter += test.GetSize();
		result += container->ComputeScore();
	}
	return result;
}
//-----------------------------------------------------------------------
void WLeMassonSimulator::Parse( WGenome< float > const & iGenome ) {}

//-----------------------------------------------------------------------
WLeMassonSimulator::~WLeMassonSimulator() {}
