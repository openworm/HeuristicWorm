/*
 * Helpers.cpp
 * Created on: Nov 14, 2011
 * Author: Alexander Dibert
 */

#include "Helpers.h"

//-----------------------------------------------------------------------
std::string Helpers::ITOA( int i ) {
	std::stringstream out;
	out << i;
	return out.str();
}
