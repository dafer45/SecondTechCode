/* Copyright 2018 Kristofer Björnson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "TBTK/Model.h"
#include "TBTK/Plotter.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Streams.h"

using namespace std;
using namespace TBTK;
using namespace Plot;

///////////////////////
// Model parameters. //
///////////////////////
//Number of sites and states to extract the probability density for.
const int SIZE_X = 500;
const int NUM_STATES = 7;
double t = 1;


/////////////////
// Potentials. //
/////////////////
//Enum type for specifying the type of the potential.
enum PotentialType{
	InfiniteSquareWell,
	SquareWell,
	HarmonicOscillator,
	DoubleWell,
	Step,
	Barrier
};

//The type of the potential.
PotentialType potentialType;

//Infinite square well.
complex<double> infiniteSquareWell(int x){
	return 0;
}

//Square well.
complex<double> squareWell(int x){
	//Parameters.
	const double WELL_DEPTH = -5e-3;
	const int WELL_LEFT_BOUNDARY_SITE = 200;
	const int WELL_RIGHT_BOUNDARY_SITE = 300;

	if(x < WELL_LEFT_BOUNDARY_SITE || x > WELL_RIGHT_BOUNDARY_SITE)
		return 0;
	else
		return WELL_DEPTH;
}

//Harmonic oscillator.
complex<double> harmonicOscillator(int x){
	//Parameters.
	const double CURVATURE = 1e-6;

	return CURVATURE*pow(x - SIZE_X/2, 2);
}

//Double well.
complex<double> doubleWell(int x){
	//Parameters.
	const double QUADRATIC = -1e-6;
	const double QUARTIC = 3e-11;

	return QUADRATIC*pow(x - SIZE_X/2, 2)
		+ QUARTIC*pow(x - SIZE_X/2, 4);
}

//Step.
complex<double> step(int x){
	//Parameters.
	const double STEP_HEIGHT = 3e-3;

	if(x < SIZE_X/2)
		return 0;
	else
		return STEP_HEIGHT;
}

//Barrier.
complex<double> barrier(int x){
	//Parameters.
	const double BARRIER_POTENTIAL_LEFT = 0;
	const double BARRIER_POTENTIAL_BARRIER = 4e-3;
	const double BARRIER_POTENTIAL_RIGHT = 2e-3;
	const double BARRIER_LEFT_BOUNDARY_SITE = 225;
	const double BARRIER_RIGHT_BOUNDARY_SITE = 275;

	if(x < BARRIER_LEFT_BOUNDARY_SITE)
		return BARRIER_POTENTIAL_LEFT;
	else if(x < BARRIER_RIGHT_BOUNDARY_SITE)
		return BARRIER_POTENTIAL_BARRIER;
	else
		return BARRIER_POTENTIAL_RIGHT;
}

//Callback that dynamically returns the current potential on a given
//site.
complex<double> potential(const Index &toIndex, const Index &fromIndex){
	int x = fromIndex[0];

	switch(potentialType){
	case InfiniteSquareWell:
		return infiniteSquareWell(x);
	case SquareWell:
		return squareWell(x);
	case HarmonicOscillator:
		return harmonicOscillator(x);
	case DoubleWell:
		return doubleWell(x);
	case Step:
		return step(x);
	case Barrier:
		return barrier(x);
	default:
		Streams::out << "Error. This should never happen.";
		exit(1);
	}
}

////////////////////
// Visualization. //
////////////////////
//Returns the current potential on the Array format.
Array<double> getPotential(){
	Array<double> p({SIZE_X});
	for(unsigned int x = 0; x < SIZE_X; x++)
		p[{x}] = real(potential({x}, {x}));

	return p;
}

//Returns the minimum value in the potential.
double getMin(const Array<double> &potential){
	double min = potential[{0}];
	for(unsigned int x = 0; x < SIZE_X; x++)
		if(potential[{x}] < min)
			min = potential[{x}];

	return min;
}

//Scales the probability densities such that NUM_STATES evenly spaced states
//can be stacked on top of each other without the probability densities
//overlapping with each other.
void scaleProbabilityDensities(
	Array<double> &probabilityDensities,
	double min,
	double max
){
	double maxProbabilityDensity = 0;
	for(unsigned int state = 0; state < NUM_STATES; state++){
		for(unsigned int x = 0; x < (unsigned int)SIZE_X; x++){
			if(probabilityDensities[{state, x}] > maxProbabilityDensity){
				maxProbabilityDensity
					= probabilityDensities[
						{state, x}
					];
			}
		}
	}
	double scaleFactor = (max - min)/(maxProbabilityDensity*NUM_STATES)/2.;
	for(unsigned int state = 0; state < NUM_STATES; state++)
		for(unsigned int x = 0; x < (unsigned int)SIZE_X; x++)
			probabilityDensities[{state, x}] *= scaleFactor;
}

//Shift the probability densities such that the zero level for a given state is
//at its eigenvalue.
void shiftProbabilityDensities(
	Array<double> &probabilityDensities,
	const Property::EigenValues &eigenValues
){
	for(unsigned int state = 0; state < NUM_STATES; state++)
		for(unsigned int x = 0; x < (unsigned int)SIZE_X; x++)
			probabilityDensities[{state, x}] += eigenValues(state);
}

//Plot the potential and probability densities and save the results to file.
void plot(
	Array<double> probabilityDensities,
	const Property::EigenValues &eigenValues,
	const string &filename
){
	//Get the current potential on array format.
	Array<double> potential = getPotential();

	//Set the minimum bound for the plot to be the minimum of the
	//potential. Set the maximum bound to be the energy of the
	//state that is one higher than the last state for which the
	//probability density is calculated.
	double min = getMin(potential);
	double max = eigenValues(NUM_STATES);

	//Scale the probability densities such that NUM_STATES evenly
	//spaced states can be stacked on top of each other without the
	//probability densities overlapping with each other.
	scaleProbabilityDensities(
		probabilityDensities,
		min,
		max
	);

	//Shift the probability densities such that the zero level for
	//a given state is at its energy.
	shiftProbabilityDensities(probabilityDensities, eigenValues);

	//Plot the potential and probability densities.
	Plotter plotter;
	plotter.setBoundsY(min, max);
	plotter.setHold(true);
	plotter.plot(
		potential,
		Decoration({224, 64, 64}, Decoration::LineStyle::Line, 2)
	);
	for(int state = 0; state < NUM_STATES; state++){
		plotter.plot(
			probabilityDensities.getSlice({state, IDX_ALL})
		);
	}
	plotter.save(filename);
}

///////////
// Main. //
///////////
int main(int argc, char **argv){
	//Create the Model.
	Model model;
	for(unsigned int x = 0; x < SIZE_X; x++){
		//Kinetic terms.
		model << HoppingAmplitude(2*t, {x}, {x});
		if(x + 1 < SIZE_X)
			model << HoppingAmplitude(-t, {x + 1}, {x}) + HC;

		//Potential term.
		model << HoppingAmplitude(potential, {x}, {x});
	}
	model.construct();

	//Setup the Solver and PropertyExtractor.
	Solver::Diagonalizer solver;
	solver.setModel(model);
	PropertyExtractor::Diagonalizer propertyExtractor(solver);

	//List of potentials to run the calculation for.
	vector<PotentialType> potentialTypes = {
		InfiniteSquareWell,
		SquareWell,
		HarmonicOscillator,
		DoubleWell,
		Step,
		Barrier
	};

	//List of filenames to save the results to.
	vector<string> filenames = {
		"figures/InfiniteSquareWell.png",
		"figures/SquareWell.png",
		"figures/HarmonicOscillator.png",
		"figures/DoubleWell.png",
		"figures/Step.png",
		"figures/Barrier.png"
	};

	//Run the calculation and plot the result for each potential.
	for(unsigned int n = 0; n < potentialTypes.size(); n++){
		//Set the current potential.
		potentialType = potentialTypes[n];

		//Run the solver.
		solver.run();

		//Calculate the probability density for the first NUM_STATES.
		Array<double> probabilityDensities({NUM_STATES, SIZE_X});
		for(unsigned int state = 0; state < NUM_STATES; state++){
			for(unsigned int x = 0; x < (unsigned int)SIZE_X; x++){
				complex<double> amplitude
					= propertyExtractor.getAmplitude(
						state,
						{x}
					);
				probabilityDensities[{state, x}]
					= pow(abs(amplitude), 2);
			}
		}

		//Calculate the eigenvalues.
		Property::EigenValues eigenValues
			= propertyExtractor.getEigenValues();

		//Plot the results.
		plot(
			probabilityDensities,
			eigenValues,
			filenames[n]
		);
	}

	return 0;
}
