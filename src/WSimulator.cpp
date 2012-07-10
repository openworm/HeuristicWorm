/*
 *  WSimulator.cpp
 *  Created on: Mar 29, 2012
 *  Author: Alexander Dibert
 */
#include "WSimulator.h"

//-----------------------------------------------------------------------
void WSimulator::Init( WSolver * iSolver, std::istream & iTrainDataStream, std::istream & iTestDataStream )
{
	this->m_solver = iSolver;

	std::string buff;
	std::getline( iTrainDataStream, buff );
	while( !iTrainDataStream.eof() )
	{
		WTest< float > learningData = WTest< float >();
		learningData.SetType( buff.substr( 0, 1 ) );
		learningData.SetName( buff );
		std::getline( iTrainDataStream, buff );
		while( buff != "---" ) {
			learningData.Add( atof( buff.c_str() ) );
			std::getline( iTrainDataStream, buff );
		}
		this->m_trainData.Add( learningData );
		std::getline( iTrainDataStream, buff );
	}

	std::getline( iTestDataStream, buff );
	while( !iTestDataStream.eof() )
	{
		WTest< float > test = WTest< float >();
		test.SetType( buff.substr( 0, 1 ) );
		test.SetName( buff );
		std::getline( iTestDataStream, buff );
		while( buff != "---" ){
			test.Add( atof( buff.c_str() ) );
			std::getline( iTestDataStream, buff );
		}
		this->m_testData.Add( test );
		std::getline( iTestDataStream, buff );
	}
}

//-----------------------------------------------------------------------
float WSimulator::Train( WGenome< float > const & iGenome )
{
	m_solver->Train( iGenome, m_result );
	return this->Evaluate( m_trainData );
}

//-----------------------------------------------------------------------
float WSimulator::Test( WGenome< float > const & iGenome )
{
	m_solver->Test( iGenome, m_result );
	return this->Evaluate( m_testData );
}

//-----------------------------------------------------------------------
std::vector< float > & WSimulator::GetSimulationResult( WGenome< float > const & iGenome )
{
	m_solver->Train( iGenome, m_result );
	return m_result;
}

//-----------------------------------------------------------------------
std::vector< float > & WSimulator::GetTestSimulationResult( WGenome< float > const & iGenome )
{
	m_solver->Test( iGenome, m_result );
	return m_result;
}