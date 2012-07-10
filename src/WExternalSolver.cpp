/*
 *  WExternalSimulator.cpp
 *  Created on: Jan 18, 2011
 *  Author: Alexander Dibert
 */

#include"WExternalSolver.h"

//-----------------------------------------------------------------------
WExternalSolver::WExternalSolver( std::istream & iTrainStream, std::istream & iTestStream )
{
/*	std::string buf;
	while( !iTrainStream.eof() ){
		WTest< float > learningData = WTest< float >();
		iTrainStream >> buf;
		learningData.SetName(buf);
		iTrainStream >> buf;
		while( buf != "---" ){
			learningData.Add( atof( buf.c_str() ) );
			iTrainStream >> buf;
		}
		m_trainData.Add(learningData);
	}

	while( !iTestStream.eof() ){
		WTest< float > test = WTest< float >();
		iTestStream >> buf;
		test.SetName(buf);
		iTestStream >> buf;
		while( buf != "---" ){
			test.Add( atof( buf.c_str() ) );
			iTestStream >> buf;
		}
		m_testData.Add(test);
	} */
}

//-----------------------------------------------------------------------
void WExternalSolver::Init( const WGenome< float > & iGenome )
{
	//this->Parse( iGenome );
}

//-----------------------------------------------------------------------
void WExternalSolver::Parse( const WGenome< float > & iGenome ) {}

//-----------------------------------------------------------------------
/*float WExternalSimulator::Simulate( const WGenome< float > & iGenome, bool iLearn ) //TODO
{
	this->Init( iGenome );

	float totalcost = 0;
	WTestData<float> & data = m_learningData;

	if( iLearn ) {
		unsigned int init = 0, bound = 3;
		for( unsigned int j = init; j < bound; j += 2 )
			totalcost += this->CountScore( j, data, iLearn );
		for( unsigned int j = init; j < bound; j += 2 )
			totalcost += this->CountScore( 3 + j, data, iLearn );
	}
	return 0;
}*/

//-----------------------------------------------------------------------
float WExternalSolver::Train( WGenome< float > const & iGenome ) {}

//-----------------------------------------------------------------------
float WExternalSolver::Test( WGenome< float > const & iGenome ) {}
