/*
 * WRandom.cpp
 *
 *  Created on: Nov 14, 2011
 *  Author: Alexander Dibert
 */

#include <math.h>
#include "WRandom.h"


float WRand( float iMin, float iMax ) {
	return iMin + ( float )rand() / RAND_MAX * ( iMax - iMin );
}

float WNormRand( float iVariance ) {
	return ( WRand( -0.5, 0.5 ) + WRand( -0.5, 0.5 ) + WRand( -0.5, 0.5 ) + WRand( -0.5, 0.5 ) + WRand( -0.5, 0.5 ) + WRand( -0.5, 0.5 ) + WRand( -0.5, 0.5 ) + WRand( -0.5, 0.5 ) +
			 WRand( -0.5, 0.5 ) + WRand( -0.5, 0.5 ) + WRand( -0.5, 0.5 ) + WRand( -0.5, 0.5 ) ) * sqrt( iVariance );
}

