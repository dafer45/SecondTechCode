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
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/Visualization/MatPlotLib/Plotter.h"

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

int main(int argc, char **argv){
	//Parameters.
	const unsigned int SIZE_X = 20;
	const unsigned int SIZE_Y = 20;
	double t = 1;
	int state = 0;

	//Create the Model.
	Model model;
	for(unsigned int x = 0; x < SIZE_X; x++){
		for(unsigned int y = 0; y < SIZE_Y; y++){
			if(x+1 < SIZE_X){
				model << HoppingAmplitude(
					-t,
					{x + 1,	y},
					{x,	y}
				) + HC;
			}
			if(y+1 < SIZE_Y){
				model << HoppingAmplitude(
					-t,
					{x, y + 1},
					{x, y}
				) + HC;
			}
		}
	}
	model.construct();

	//Setup and run the Solver.
	Solver::Diagonalizer solver;
	solver.setModel(model);
	solver.run();

	//Setup the PropertyExtractor.
	PropertyExtractor::Diagonalizer propertyExtractor(solver);

	//Print the eigenvalue for the given state.
	Streams::out << "The energy of state " << state << " is "
		<< propertyExtractor.getEigenValue(state) << "\n";

	//Calculate the probability density for the given state.
	Array<double> probabilityDensity({SIZE_X, SIZE_Y});
	for(unsigned int x = 0; x < SIZE_X; x++){
		for(unsigned int y = 0; y < SIZE_Y; y++){
			//Get the probability amplitude at site (x, y) for the
			//given state.
			complex<double> amplitude
				= propertyExtractor.getAmplitude(
					state,
					{x, y}
				);

			//Calculate the probability density.
			probabilityDensity[{x, y}] = pow(
				abs(amplitude),
				2
			);
		}
	}

	//Plot the probability density.
	Plotter plotter;
	plotter.plot(probabilityDensity);
	plotter.save("figures/ProbabilityDensity.png");

	return 0;
}
