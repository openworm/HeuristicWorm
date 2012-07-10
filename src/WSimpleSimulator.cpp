/*
 *  WSimpleSimulator.hpp
 *  Created on: Mar 29, 2012
 *  Author: Alexander Dibert
 */

#include "WSimpleSimulator.h"

//-----------------------------------------------------------------------
WSimpleSimulator::WSimpleSimulator( WSolver * iSolver, std::istream & iTrainDataStream, std::istream & iTestDataStream )
{
	this->Init( iSolver, iTrainDataStream, iTestDataStream );
}

//-----------------------------------------------------------------------
void WSimpleSimulator::Init( WSolver * iSolver, std::istream & iTrainDataStream, std::istream & iTestDataStream ) // TODO make it work with data of arbitrary length
{
	WSimulator::Init( iSolver, iTrainDataStream, iTestDataStream );

	std::string ICoefficients = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ";
	std::stringstream IStream( ICoefficients );
	for( int i = 0; i < 20; ++i ) // TODO make it not hard-codeed.
		IStream >> m_ICoefficients[ i ];

	std::string UCoefficients = "0 2000000000000 2000000000000 2000000000000 2000000000000 2000000000000 1000000000000 1000000000000 1000000000000 1000000000000 1000000000000 1000000000000 1000000000000 1000000000000 1000000000000 1000000000000 1000000000000 1000000000000 1000000000000 1000000000000";
	std::stringstream testStream( UCoefficients );
	for( int i = 0; i < 20; ++i )
		testStream >> m_UCoefficients[ i ];
}

//-----------------------------------------------------------------------

float WSimpleSimulator::Evaluate( WTestData< float > const & iTarget )
{
	double result = 0;
	int indx = 0;
	float * coeefs ;

	for( size_t i = 0; i < iTarget.GetSize(); ++i )
	{
		WTest< float > const & test = iTarget.Get( i );
		char testType = ( test.GetType() )[0];
		switch( testType ) {
			case 'I':
				coeefs = m_ICoefficients;
				break;
			case 'V':
				coeefs = m_UCoefficients;
				break;
			default:
				std::cerr << "Error:Wrong test type! Aborting";
				::exit( 0 );
				break;
		}
		for( size_t j = 0; j < test.GetSize() ; ++j ) {
			result += pow( ( double )( test.Get( j ) - this->m_result[ indx] ),  2) * coeefs[ j ] * 1000;
			++indx;
		}
	}
	return result;
}

//-----------------------------------------------------------------------
void WSimpleSimulator::Parse( WGenome< float > const & iGenome ) {}

//-----------------------------------------------------------------------
WSimpleSimulator::~WSimpleSimulator() {}