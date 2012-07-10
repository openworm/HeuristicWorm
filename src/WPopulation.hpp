/*
 * File:   WPopulation.hpp
 * Author: Alexander Dibert
 *
 * Created on October 9, 2011, 11:58 AM
 */

#ifndef _WPOPULATION_H
#define	_WPOPULATION_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>

/*
 * WGenomeStruct template struct.
 * Keeps the genome structure
 */

template < typename Type > struct WGenomeStruct // TODO class? Get Set?
{
public:
	WGenomeStruct()
	{
		m_minMax[ 0 ].reserve( 30 );
		m_minMax[ 1 ].reserve( 30 );
		m_minMax[ 2 ].reserve( 30 );
		m_name.reserve( 30 );
		m_length = 0;
    }

	void AddParam( Type iMin, Type iMax, std::string iName )
	{
		m_minMax[ 0 ].push_back( iMin );
		m_minMax[ 1 ].push_back( iMax );
		m_minMax[ 2 ].push_back( iMax - iMin );
		m_name.push_back( iName );
		++m_length;
	}

	void AddTemplateLine( std::string str, std::string symb ) // Gene names should be already known,
	{
		m_templateStrings.push_back( str );
		for( size_t i = 0; i < m_name.size(); ++i ) {
			if( !symb.compare( m_name[ i ] ) ) {
				m_templateSymbols.push_back( i );
				return;
			}
		}
		m_templateSymbols.push_back( -1 );
	}

public:
	std::vector< Type > m_minMax[3];
	std::vector< std::string > m_name;
	std::vector< std::string > m_templateStrings;  // A template used in serialization.
	std::vector< int> m_templateSymbols;
	unsigned int m_length;
};

//---------------------------------------------------------------------------------

/*
 * Genome template class.
 * Represents information about one geneme.
 * Contains genome structure as a static field
 */
template < typename Type > class WGenome
{
public:
	WGenome();
	WGenome( size_t iSize );
	WGenome( const WGenome & );
	~WGenome();

	void Reserve( size_t iSize );
	void Resize( size_t iSize );
	void AddGene( Type iItem );

	WGenome< Type > operator +( const WGenome< Type > & iGenome ) const;
	WGenome< Type > operator -() const;
	WGenome< Type > operator -( const WGenome< Type > & iGenome ) const;
	WGenome< Type > operator *( Type iItem ) const;
	WGenome< Type > operator /( int iScalar ) const;
	WGenome< Type > & operator =( const WGenome< Type > & iGenome );

	bool operator < ( WGenome< Type > const & iGenome ) const { return m_score < iGenome.GetScore(); }

	void Normalize();

	void SetGene( int i, Type iItem );
	Type GetGene( int i ) const { return m_genome[ i ]; }
	Type & operator []( int i ) { return m_genome[ i ]; }

	void SetScore( float iScore );
	void SetTestScore( float iScore );
	float GetScore() const { return m_score; }
	float GetTestScore() const { return m_testScore; }

	size_t GetLength() const { return m_genome.size(); }

	void Serialize( std::ostream & iStream, bool iuseTemplate = false ) const;
	void Deserialize( std::istream & iStream );

public:
	static WGenomeStruct<Type> m_genomeStruct;

private:
	std::vector< Type > m_genome;
	float m_score;
	float m_testScore;

	typedef typename std::vector< Type >::const_iterator TGenomeVecConstIter;
	typedef typename std::vector< Type >::iterator TGenomeVecIter;
};

//-----------------------------------------------------------------------
template < typename Type > WGenomeStruct< Type > WGenome< Type >::m_genomeStruct;

//-----------------------------------------------------------------------
template < typename Type > WGenome< Type >::WGenome()
		: m_genome( 0 )
		, m_score( std::numeric_limits< float >::infinity() )
		, m_testScore( std::numeric_limits< float >::infinity() ) {}

//-----------------------------------------------------------------------
template < typename Type > WGenome< Type >::WGenome( size_t iSize )
		: m_genome( 0 )
		, m_score( std::numeric_limits< float >::infinity() )
		, m_testScore ( std::numeric_limits< float >::infinity() )
{
	m_genome.reserve( iSize );
}

//-----------------------------------------------------------------------
template < typename Type > WGenome< Type >::WGenome( const WGenome & iGenome )
		: m_genome( 0 )
		, m_score( std::numeric_limits< float >::infinity() )
		, m_testScore( std::numeric_limits< float >::infinity() )
{
	for( size_t i = 0; i < iGenome.GetLength(); ++i ) {
		m_genome.push_back( iGenome.GetGene( i ) );
	}
}

//-----------------------------------------------------------------------
template < typename Type > WGenome< Type >::~WGenome() {}

//-----------------------------------------------------------------------
template < typename Type > void WGenome< Type >::Reserve( size_t iSize )
{
	m_genome.reserve( iSize );
}

//-----------------------------------------------------------------------
template < typename Type > void WGenome< Type >::Resize( size_t iSize )
{
	m_genome.resize( iSize );
}

//-----------------------------------------------------------------------
template < typename Type > void WGenome< Type >::AddGene( Type iItem )
{
	m_genome.push_back( iItem );
}

//-----------------------------------------------------------------------
template < typename Type > WGenome< Type > WGenome < Type >::operator +( const WGenome< Type > & iGenome ) const
{
	WGenome< Type > oGenome = WGenome< Type >( m_genome.size() );
	for( size_t i = 0; i < m_genome.size(); ++i ) {
		oGenome.AddGene( m_genome[ i ] + iGenome.GetGene( i ) );
	}
	return oGenome;
}

//-----------------------------------------------------------------------
template < typename Type > WGenome< Type > WGenome < Type >::operator -() const
{
	WGenome< Type > oGenome = WGenome< Type >( m_genome.size() );
	for( size_t i = 0; i < m_genome.size(); ++i ) {
        oGenome.Add( - m_genome[ i ] );
	}
	return oGenome;
}

//-----------------------------------------------------------------------
template < typename Type > WGenome < Type > WGenome < Type >::operator -( const WGenome< Type > & iGenome ) const
{
    WGenome< Type > oGenome( m_genome.size() );
	for( size_t i = 0; i < m_genome.size(); ++i ) {
		oGenome.AddGene( m_genome[ i ] - iGenome.GetGene( i ) );
	}
	return oGenome;
}

//-----------------------------------------------------------------------
template < typename Type > WGenome< Type > WGenome < Type >::operator *( Type iItem ) const
{
	WGenome< Type > oGenome( m_genome.size() );
	for( size_t i = 0; i < m_genome.size(); ++i ) {
		oGenome.AddGene( m_genome[ i ] * iItem );
	}
	return oGenome;
}

//-----------------------------------------------------------------------
template < typename Type > WGenome< Type > WGenome < Type >::operator /( int iScalar ) const
{
	WGenome< Type > oGenome( m_genome.size() );
	for( size_t i = 0; i < m_genome.size(); ++i ) {
		oGenome.AddGene( m_genome[ i ] / iScalar );
	}
	return oGenome;
}

//-----------------------------------------------------------------------
template < typename Type > WGenome< Type > & WGenome < Type >::operator =( const WGenome< Type > & iGenome )
{
	m_genome.clear();
	for( size_t i = 0; i < iGenome.GetLength(); ++i ) {
		m_genome.push_back( iGenome.GetGene( i ) );
	}
	m_score = iGenome.GetScore();
	m_testScore = iGenome.GetTestScore();
	return *this;
}

//-----------------------------------------------------------------------
template < typename Type > void WGenome < Type >::Normalize()
{
	float diff;
	for( size_t i = 0; i < m_genomeStruct.m_length; ++i )
	{
		if( m_genome[ i ] > m_genomeStruct.m_minMax[ 1 ][ i ] ) {
			diff = m_genomeStruct.m_minMax[ 1 ][ i ] - m_genomeStruct.m_minMax[ 0 ][ i ];
			m_genome[ i ] = m_genomeStruct.m_minMax[ 0 ][ i ] + std::fmod( m_genome[ i ] - m_genomeStruct.m_minMax[ 0 ][ i ], diff );
		}
		if( m_genome[ i ] < m_genomeStruct.m_minMax[ 0 ][ i ] ) {
			diff = m_genomeStruct.m_minMax[ 1 ][ i ] - m_genomeStruct.m_minMax[ 0 ][ i ];
			m_genome[ i ] = m_genomeStruct.m_minMax[ 1 ][ i ] - std::fmod( m_genomeStruct.m_minMax[ 0 ][ i ] - m_genome[ i ], diff );
		}
	}
}

//-----------------------------------------------------------------------
template < typename Type > void WGenome< Type >::SetGene( int i, Type iItem )
{
	m_genome[ i ] = iItem;
}

//-----------------------------------------------------------------------
template < typename Type > void WGenome< Type >::SetScore( float iScore )
{
	m_score = iScore;
}

//-----------------------------------------------------------------------
template < typename Type > void WGenome< Type >::SetTestScore( float iScore )
{
	m_testScore = iScore;
}

//-----------------------------------------------------------------------
template < typename Type > void WGenome< Type >::Serialize( std::ostream & iStream, bool iuseTemplate ) const
{
	if( !iuseTemplate )
	{
		for( TGenomeVecConstIter it = m_genome.begin(); it != m_genome.end(); ++it )
			iStream << *it << std::endl;
		iStream << m_score << std::endl;
		iStream << m_testScore << std::endl;
		return;
	}

	std::vector< std::string >::const_iterator strIter = m_genomeStruct.m_templateStrings.cbegin();
	std::vector< int >::const_iterator symbIter = m_genomeStruct.m_templateSymbols.cbegin();
	for( ; strIter != m_genomeStruct.m_templateStrings.end() && symbIter != m_genomeStruct.m_templateSymbols.end(); ++strIter, ++symbIter ) {
		iStream << *strIter;
		if( *symbIter != -1  ) {
			iStream << m_genome[ *symbIter ];
		}
	}
	if( strIter != m_genomeStruct.m_templateStrings.end() )
		iStream << *strIter;
}

//-----------------------------------------------------------------------
template < typename Type > void WGenome< Type >::Deserialize( std::istream & iStream )
{
	for( TGenomeVecIter it = m_genome.begin(); it != m_genome.end(); ++it )
		iStream >> *it;
	iStream >> m_score;
	iStream >> m_testScore;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

/*
 * Population class.
 * Contains a vector of genomes.
 */
template < typename Type > class WPopulation
{
public:
	WPopulation();
	WPopulation( size_t iSize );
	~WPopulation();

	void Reserve( size_t iSize );
	void Resize( size_t iSize );
	void Clear();

	void AddGenome( const WGenome < Type > & iItem );
	void SetGene( int iGenomeNum, int iGenNum, Type iValue );
	WGenome< Type > const & GetGenome( int i ) const;

	float GetBestScore() const;
	void SetBestScore( float iScore );
	int GetBestNum() const;
	void SetBestNum( int iNum );
	WGenome< Type > & GetBestGene();
	void SetGenomeLength( size_t iLength ) { m_genomeLength = iLength; }
	unsigned int GetSize() { return m_genomeVec.size(); }

	WGenome< Type > & operator [](int i);

	void Serialize( std::ostream & iStream, bool iuseTemplate = false ) const;
	void Deserialize( std::istream & iStream );
	void Serialize( std::string iPrefix, bool iuseTemplate = false ) const;
	void Deserialize( std::string iPrefix );

private:
	std::vector< WGenome< Type > > m_genomeVec;

	size_t m_genomeLength;
	float m_bestScore;
	int m_bestNum;

	typedef typename std::vector< WGenome< Type > >::iterator TGenomeVecIter;
	typedef typename std::vector< WGenome< Type > >::const_iterator TGenomeVecConstIter;
};

//-----------------------------------------------------------------------
template < typename Type > WPopulation< Type >::WPopulation()
		: m_genomeVec( 0 )
		, m_bestScore( std::numeric_limits< float >::infinity() )
		, m_bestNum( 0 )
{
    m_genomeVec.reserve( 100 );
}

//-----------------------------------------------------------------------
template < typename Type > WPopulation< Type >::WPopulation( size_t iSize )
		: m_genomeVec( 0 )
		, m_bestScore( std::numeric_limits< float >::infinity() )
		, m_bestNum( 0 )
{
    m_genomeVec.reserve( iSize );
}

//-----------------------------------------------------------------------
template < typename Type > WPopulation< Type >::~WPopulation() {}

//-----------------------------------------------------------------------
template < typename Type > void WPopulation< Type >::Reserve( size_t iSize )
{
	m_genomeVec.reserve( iSize );
}

//-----------------------------------------------------------------------
template < typename Type > void WPopulation< Type >::Resize( size_t iSize )
{
	m_genomeVec.resize( iSize );
}

//-----------------------------------------------------------------------
template < typename Type > void WPopulation< Type >::Clear()
{
	m_genomeVec.clear();
	m_bestScore = 10000000;
	m_bestNum = 0;
}

//-----------------------------------------------------------------------
template < typename Type > void WPopulation< Type >::AddGenome( const WGenome< Type > & iItem )
{
	// if( iItem.GetLength() != m_genomeLength )
		// std::cout << "ERROR: Population: genome add: wrong genome length" << std::endl;
	m_genomeVec.push_back( iItem );
}

//-----------------------------------------------------------------------
template < typename Type > void WPopulation< Type >::SetGene( int iGenomeNum, int iGenNum, Type iValue )
{
	m_genomeVec[ iGenomeNum ].SetGene( iGenNum, iValue );
}

//-----------------------------------------------------------------------
template < typename Type > WGenome< Type > const & WPopulation < Type >::GetGenome( int iNum ) const
{
	return m_genomeVec[ iNum ];
}

//-----------------------------------------------------------------------
template < typename Type > WGenome< Type > & WPopulation < Type >::operator []( int i )
{
	return m_genomeVec[ i ];
}

//-----------------------------------------------------------------------
template <typename Type> float WPopulation< Type >::GetBestScore() const
{
    return m_bestScore;
}

//-----------------------------------------------------------------------
template <typename Type> void WPopulation< Type >::SetBestScore( float iScore )
{
    m_bestScore = iScore;
}

//-----------------------------------------------------------------------
template <typename Type> int WPopulation< Type >::GetBestNum() const
{
    return m_bestNum;
}

//-----------------------------------------------------------------------
template <typename Type> void WPopulation< Type >::SetBestNum( int iNum )
{
    m_bestNum = iNum;
}

//-----------------------------------------------------------------------
template <typename Type> WGenome< Type > & WPopulation< Type >::GetBestGene()
{
    return m_genomeVec[ m_bestNum ];
}

//-----------------------------------------------------------------------
template <typename Type> void WPopulation< Type >::Serialize ( std::ostream& iStream, bool iuseTemplate ) const
{
	for( TGenomeVecConstIter it = m_genomeVec.begin(); it != m_genomeVec.end(); ++it )
		it->Serialite( iStream, iuseTemplate );
}

//-----------------------------------------------------------------------
template <typename Type> void WPopulation< Type >::Deserialize ( std::istream& iStream )
{
	for( TGenomeVecIter it = m_genomeVec.begin(); it != m_genomeVec.end(); ++it )
		it->Derialite( iStream );
}

//-----------------------------------------------------------------------
template <typename Type> void WPopulation< Type >::Serialize ( std::string iPrefix, bool iuseTemplate ) const
{
	unsigned int counter = 0;
	std::ostringstream sstream;

	std::cout << "Going to serialize " << m_genomeVec.size() << " genomes" << std::endl;

	for( TGenomeVecConstIter it = m_genomeVec.begin(); it != m_genomeVec.end(); ++it )
	{
		sstream.str("");
		sstream << iPrefix << "_" << counter << ".txt";
		std::string tmp = sstream.str();
		std::ofstream stream( sstream.str() );

		if( !stream.is_open() ) {
			std::cerr << "Population serializer: Error: cand open output file: " << sstream.str() << std::endl;
			return;
		}

		it->Serialize( stream, iuseTemplate );
		stream.close();
		++counter;
	}
}

//-----------------------------------------------------------------------
template <typename Type> void WPopulation< Type >::Deserialize ( std::string iPrefix )
{
	unsigned int counter = 0;
	std::ostringstream sstream;

	std::cout << "Going to deserialize " << m_genomeVec.size() << " genomes" << std::endl;

	for( TGenomeVecIter it = m_genomeVec.begin(); it != m_genomeVec.end(); ++it )
	{
		sstream.str("");
		sstream << iPrefix << "_" << counter << ".txt";
		std::ifstream stream( sstream.str() );

		if( !stream.is_open() ) {
			std::cout << "Population deserializer: Error: cand open output file: " << sstream.str() << std::endl;
			return;
		}

		it->Resize( m_genomeLength );
		it->Deserialize( stream );
		stream.close();
		++counter;
	}
}

#endif	/* _WPOPULATION_H */

