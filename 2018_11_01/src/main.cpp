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

#include "TBTK/AbstractIndexFilter.h"
#include "TBTK/Model.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/TBTK.h"
#include "TBTK/Visualization/MatPlotLib/Plotter.h"

using namespace std;
using namespace TBTK;
using namespace Visualization::MatPlotLib;

//Parameters.
const unsigned int SIZE = 41;
const unsigned int SIZE_X = SIZE;
const unsigned int SIZE_Y = SIZE;
const double OUTER_RADIUS = SIZE/2;
const double INNER_RADIUS = SIZE/8;
double t = 1;
int state = 0;

//IndexFilter.
class MyIndexFilter : public AbstractIndexFilter{
public:
	//Implements AbstractIndexFilter::clone().
	MyIndexFilter* clone() const{
		return new MyIndexFilter();
	}

	//Implements AbstractIndexFilter::isIncluded().
	bool isIncluded(const Index &index) const{
		//Extract x and y from the Index.
		int x = index[0];
		int y = index[1];

		//Calculate the distance from the center.
		double r = sqrt(
			pow(abs(x - (int)SIZE_X/2), 2)
			+ pow(abs(y - (int)SIZE_Y/2), 2)
		);

		//Return true if the distance is less than the outer radius of
		//the annulus, but larger than the inner radius.
		if(r < OUTER_RADIUS && r > INNER_RADIUS)
			return true;
		else
			return false;
	}
};

int main(int argc, char **argv){
	//Initialize TBTK.
	Initialize();

	//Create filter.
	MyIndexFilter filter;

	//Create the Model.
	Model model;
	model.setFilter(filter);
	for(unsigned int x = 0; x < SIZE_X; x++){
		for(unsigned int y = 0; y < SIZE_Y; y++){
			model << HoppingAmplitude(
				-t,
				{x + 1,	y},
				{x,	y}
			) + HC;
			model << HoppingAmplitude(
				-t,
				{x, y + 1},
				{x, y}
			) + HC;
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
	Array<double> probabilityDensity({SIZE_X, SIZE_Y}, 0);
	for(unsigned int x = 0; x < SIZE_X; x++){
		for(unsigned int y = 0; y < SIZE_Y; y++){
			if(!filter.isIncluded({x, y}))
				continue;
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
