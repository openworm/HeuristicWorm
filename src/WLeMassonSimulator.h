/*
 *  WLeMassonSimulator.h
 *  Created on: Mar 27, 2012
 *  Author: Alexander Dibert
 */

#ifndef WLeMassonSimulator_h_
#define WLeMassonSimulator_h_

#include <map>
#include <limits>
#include <cmath>

#include "WSimulator.h"

static const unsigned short SCOREARRAYSIZE = 100; //50

class WLeMassonSimulator : public WSimulator
{
public:
	WLeMassonSimulator( WSolver * iSolver, std::istream & iTrainDataStream, std::istream & iTestDataStream );
	virtual ~WLeMassonSimulator();

private:
	void Init( WSolver * iSolver, std::istream & iTrainDataStream, std::istream & iTestDataStream );
	float Evaluate( WTestData< float > const & iTarget );
	void Parse( WGenome< float > const & iGenome );

private:
	struct WScoreContainer
	{
		WScoreContainer( std::string const & iString )
			: m_type( iString )
			, m_initialized( false )
		{}

		void Init( WTest< float > const & iTestData, float iDeltaT, float xMultiplier = 1, float dxMultiplier = 1 )
		{
			if( !m_initialized )
			{
				if( !m_type.compare( "I" ) )
				{
					m_minX = -80e-3;
					m_maxX = 120e-3;
					m_minDx = -.15; // * dxMultiplier; //-150
					m_maxDx = .25; // * dxMultiplier; //250
				}

				if( !m_type.compare( "V" ) )
				{
					m_minX = -80e-9;
					m_maxX = 120e-9;
					m_minDx = -.000000015; // * dxMultiplier; //-150
					m_maxDx = .000000025; // * dxMultiplier; //250
				}

				m_xStep = ( m_maxX - m_minX ) / SCOREARRAYSIZE;
				m_dxStep = ( m_maxDx - m_minDx ) / SCOREARRAYSIZE;

				for( unsigned int i = 0; i < SCOREARRAYSIZE; ++i )
					for( unsigned int j = 0; j < SCOREARRAYSIZE; ++j )
						m_targetScore[ i ][ j ] = 0;
				m_targetPointCounter = 0;

				float prevValue;
				float curValue = iTestData.Get( 0 );
				for( unsigned int i = 1; i < iTestData.GetSize(); ++i )
				{
					prevValue = curValue;
					curValue = iTestData.Get( i );
					UpdateScore( prevValue, curValue, iDeltaT * 25, false ); // 25 - workaround
				}
				m_initialized = true;
			}

			this->ResetScore();
		}

		void ResetScore()
		{
			m_scorePointCounter = 0;
			for( unsigned int i = 0; i < SCOREARRAYSIZE; ++i )
				for( unsigned int j = 0; j < SCOREARRAYSIZE; ++j )
					m_score[ i ][ j ] = 0;
		}

		void UpdateScore( float prevValue, float curValue, float iTimeStep, bool iIsScore = true )
		{
			float dx = ( curValue - prevValue ) / iTimeStep;
			int cellXNumber = ( int )( ( curValue - m_minX ) / m_xStep );
			int cellDxNumber = ( int )( ( dx - m_minDx ) / m_dxStep );

			if( cellXNumber >= SCOREARRAYSIZE || cellXNumber < 0 ) {
				std::cerr << "WLEMassonSimulator: X score data is out of range: " << cellXNumber << ".\n";
				return;
			}
			if( cellDxNumber >= SCOREARRAYSIZE || cellDxNumber < 0 ) {
				std::cerr << "WLEMassonSimulator: dX score data is out of range" << cellDxNumber << ".\n";
				return;
			}

			if( iIsScore ) {
				++( m_score[ cellXNumber ][ cellDxNumber ] );
				++m_scorePointCounter;
			}
			else {
				++( m_targetScore[ cellXNumber ][ cellDxNumber ] );
				++m_targetPointCounter;
			}
		}

		float ComputeScore() const
		{
			if( m_scorePointCounter == 0 )
				return 1000;
			float res = 0;
			for( unsigned int i = 0; i < SCOREARRAYSIZE; ++i )
				for( unsigned int j = 0; j < SCOREARRAYSIZE; ++j )
					res += pow( fabs( m_score[ i ][ j ] / m_scorePointCounter - m_targetScore[ i ][ j ] / m_targetPointCounter ), 2 );
			return sqrt( res );
		}

		std::string m_type;
		bool  m_initialized;
		float m_score [SCOREARRAYSIZE][SCOREARRAYSIZE];
		float m_targetScore [SCOREARRAYSIZE][SCOREARRAYSIZE];
		float m_minX;
		float m_maxX;
		float m_minDx;
		float m_maxDx;
		float m_xStep;
		float m_dxStep;
		unsigned int m_scorePointCounter;
		unsigned int m_targetPointCounter;
	};

	float m_deltaT;

	WScoreContainer  m_I100Score;
	WScoreContainer  m_I400Score;
	WScoreContainer  m_I700Score;
	WScoreContainer  m_V0Score;
	WScoreContainer  m_V20Score;
	WScoreContainer  m_V40Score;
};

#endif
