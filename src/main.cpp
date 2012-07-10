/*
 * File:   main.cpp
 * Author: Alexander Dibert
 *
 * Created on October 9, 2011, 11:55 AM
 */
#include <fstream>
#include <sstream>
#include <string.h>
#include <getopt.h>

#include "WDiffEvolution.h"
#include "WDiffEvolutionSA.h"
#include "WEvolutionStrategy.h"
#include "WEvolutionStrategyModified.h"
#include "WPureRandom.h"
#include "WSimpleSimulator.h"
#include "WLeMassonSimulator.h"
#include "WJordanSolver.h"

int main(int argc, char* argv[])
{
	std::ifstream input( "input.txt" );
	if( !input.is_open() ) {
		std::cerr << "Error opening input config file\n";
		return 1;
	}
	std::ifstream trainingData( "trainingData.txt" );
	if( !trainingData.is_open() ) {
		std::cerr << "Error opening training data file\n";
		return 1;
	}
	std::ifstream testData( "testData.txt" );
	if( !testData.is_open() ) {
		std::cerr << "Error opening testsData file\n";
		return 1;
	}

	const char * options = "htgsLSRE";
	const struct option long_options[] = {
		{ "help", no_argument, NULL, 'h' },
		{ "generate", no_argument, NULL, 'g' },
		{ "step", no_argument, NULL, 's' },
		{ "lemasson", no_argument, NULL, 'L' },
		{ "test", no_argument, NULL, 't' },
		{ "self_adaptive", no_argument, NULL, 'S' },
		{ "pure_random", no_argument, NULL, 'R' },
		{ "evolution_strategy", no_argument, NULL, 'E' },
		{ NULL, 0, NULL, 0 }
	};

	int res;
	int optionIndex;
	bool generate = false;
	bool step = false;
	bool printHelp = false;
	bool test = false;
	bool selfAdaptive = false;
	char algType = 'D';
	char simulatorType = 'S';
	char solverType = 'J';

	while( ( res = ::getopt_long( argc, argv, options, long_options, &optionIndex ) ) !=-1 ) {
		switch( res )
		{
			case 'g':
				generate = true;
				break;
			case 's':
				step = true;
				break;
			case 't':
				test = true;
				break;
			case 'L':
				simulatorType = 'L';
				break;
			case 'S':
				selfAdaptive = true;
				break;
			case 'R':
				algType = 'R';
				break;
			case 'E':
				algType = 'E';
				break;
			case 'h':
			default:
				printHelp = true;
				break;
		}
	}

	if( printHelp )
	{
		std::cout << "-h[--help] Print this message." << std::endl;
		std::cout << "-g[--generate] Generate population and write it into files in ./CurrPopulation/." << std::endl; // TODO correct dir
		std::cout << "-s[--step] Read generation from files from ./CurrPopulation, calculate fitness and write it back to files." << std::endl;
		std::cout << "-L[--lemasson] Use Le Masson's type fitness function." << std::endl;
		std::cout << "-S[--self_adaptive] Use self adaptive (or just more sophisticated) version of an algorithm." << std::endl;
		std::cout << "-R[--pure_random] Use a random guess to construct a new population." << std::endl;
		std::cout << "-E[--evolution_strategy] Ude Evolution strategy algorythm." << std::endl;
		return 0;
	}

	// TODO MORE ENCAPSULATION!!!

	if( generate ) // TODO Make it able to create dir if needed.
	{
		WDiffEvolution DE( NULL, false );
		DE.Init( input );
		DE.GeneratePopulation();
		DE.SerializePopulation();
		return 0;
	}

	WSolver * solver;
	switch( solverType )
	{
		case 'j':
		default:
			solver = new WJordanSolver();
			break;
	}

	WSimulator * simulator;
	switch( simulatorType )
	{
		case 'L':
		{
			simulator = new WLeMassonSimulator( solver, trainingData, testData );
			break;
		}
		case 'S':
		default:
			simulator = new WSimpleSimulator( solver, trainingData, testData );
			break;
	}

	if( test )
	{
		WDiffEvolution DE( simulator, true );
		DE.Init( input );
		DE.Start();
		delete simulator;
		delete solver;
		return 0;
	}

	WHeuristicAlg< float, float > * alg;
	switch( algType )
	{
	case 'D':
	{
		if( selfAdaptive )
			alg = new WDiffEvolutionSA( simulator, false );
		else
			alg = new WDiffEvolution( simulator, false );
		break;
	}
	case 'R':
	{
		alg = new WPureRandom( simulator, false );
		break;
	}
	case 'E':
	{
		if( selfAdaptive )
			alg = new WEvolutionStrategyModified( simulator, false );
		else
			alg = new WEvolutionStrategy( simulator, false );
		break;
	}
	default:
		std::cerr << "Unknown alg type;" << std::endl;
	}


	if( step )  // TODO Run generation if can't load population.
				// TODO Run selection.
	{
		alg->Init( input );
		alg->DeserializePopulation();
		alg->Fitness();
		alg->SerializePopulation();
		delete simulator;
		delete solver;
		delete alg;
		return 0;
	}

	// Common launch.
	alg->Init( input );
	alg->GeneratePopulation();
	alg->Start();
	delete simulator;
	delete solver;
	delete alg;
    return 0;
}



