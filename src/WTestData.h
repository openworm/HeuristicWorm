/*
 * WTestData.h
 *  Created on: Oct 11, 2011
 *  Author: Alexander Dibert
 */

#ifndef WTESTDATA_H_
#define WTESTDATA_H_

#include <string>

/*
 * Test class. Keeps and provides an access to the test/train data.
 */

template <typename Type> class WTest {

public:
	WTest();
	WTest( const WTest & );

	void SetName( std::string const & iName );
	void SetType( std::string const & iType );
	void Add ( Type const & iData );
	Type Get ( int iNum ) const;
	size_t GetSize() const { return m_data.size(); }
	std::string const & GetName() const { return m_name; }
	std::string const & GetType() const { return m_type; }

private:
	std::string m_type;
	std::string m_name;
	std::vector< Type > m_data;
};

template < typename Type > WTest <Type>::WTest() {}

template < typename Type > WTest <Type>::WTest( const WTest & iTest )
{
	m_name = iTest.m_name;
	m_type = iTest.m_type;
	for( unsigned int i = 0; i < iTest.GetSize(); ++i )
		m_data.push_back( iTest.Get( i ) );
}

template < typename Type > void WTest <Type>::SetType( std::string const & iType ) {
	m_type = iType;
}

template < typename Type > void WTest <Type>::SetName( std::string const & iName ) {
	m_name = iName;
}

template < typename Type > void WTest <Type>::Add( Type const & iData ) {
	m_data.push_back( iData );
}

template < typename Type > Type WTest <Type>::Get( int iNum ) const {
	return m_data[ iNum ];
}

//-----------------------------------------------

template < typename Type > class WTestData {

public:
	WTestData();
	void Add( WTest< Type > const & iData );
	Type Get( int iTestNum, int iNum ) const;
	WTest< Type > const & Get( int iTestNum ) const;
	size_t GetSize() const { return m_dataVec.size(); }
	size_t GetSize( int iTestNum ) const { return m_dataVec[ iTestNum ].GetSize(); }

private:
	std::vector< WTest< Type > > m_dataVec;
};

template <typename Type> WTestData< Type >::WTestData() {}

template <typename Type> void WTestData< Type >::Add( WTest< Type > const & iData ) {
	m_dataVec.push_back( iData );
}

template <typename Type> Type WTestData< Type >::Get( int iTestNum, int iNum ) const {
	return m_dataVec[ iTestNum ].Get( iNum );
}

template <typename Type> WTest< Type > const & WTestData< Type >::Get( int iTestNum ) const {
	return m_dataVec[ iTestNum ];
}

#endif /* WTESTDATA_H_ */
