/*
 *  WSolver.h
 *  Created on: Feb 07, 2012
 *  Author: Alexander Dibert
 */
#ifndef WSolver_H_
#define WSolver_H_

#include "WPopulation.hpp"

class WSolver
{
public:
	virtual void Init( const WGenome< float > & iGenome ) = 0;
	virtual void Train( WGenome< float > const & iGenome, std::vector< float > & oResult ) = 0;
	virtual void Test( WGenome< float > const & iGenome, std::vector< float > & oResult ) = 0;
};

#endif