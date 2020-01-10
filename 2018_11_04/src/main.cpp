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

#include "TBTK/Array.h"
#include "TBTK/Model.h"
#include "TBTK/PropertyExtractor/Diagonalizer.h"
#include "TBTK/Solver/Diagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/TBTK.h"

using namespace std;
using namespace TBTK;

int main(int argc, char **argv){
	//Parameters.
	const int SIZE_X = 4;
	const int SIZE_Y = 3;
	double t = 1;

	//Create the Model.
	Model model;
	for(unsigned int x = 0; x < SIZE_X; x++){
		for(unsigned int y = 0; y < SIZE_Y; y++){
			model << HoppingAmplitude(
				4*t,
				{x, y},
				{x, y}
			);
			if(x + 1 < SIZE_X){
				model << HoppingAmplitude(
					-t,
					{x + 1,	y},
					{x,	y}
				) + HC;
			}
			if(y + 1 < SIZE_Y){
				model << HoppingAmplitude(
					-t,
					{x, y + 1},
					{x, y}
				) + HC;
			}
		}
	}
	model.construct();

	//Get the HoppingAmplitudeSet from the Model and extract the basis
	//size.
	const HoppingAmplitudeSet &hoppingAmplitudeSet
		= model.getHoppingAmplitudeSet();
	unsigned int basisSize = hoppingAmplitudeSet.getBasisSize();

	//Initialize the Hamiltonian on a format most suitable for the
	//algorithm at hand.
	Array<complex<double>> hamiltonian({basisSize, basisSize}, 0.);

	//Iterate over the HoppingAmplitudes.
	for(
		HoppingAmplitudeSet::ConstIterator iterator
			= hoppingAmplitudeSet.cbegin();
		iterator != hoppingAmplitudeSet.cend();
		++iterator
	){
		//Extract the amplitude and physical indices from the
		//HoppingAmplitude.
		complex<double> amplitude = (*iterator).getAmplitude();
		const Index &toIndex = (*iterator).getToIndex();
		const Index &fromIndex = (*iterator).getFromIndex();

		//Convert the physical indices to linear indices.
		unsigned int row = hoppingAmplitudeSet.getBasisIndex(toIndex);
		unsigned int column = hoppingAmplitudeSet.getBasisIndex(
			fromIndex
		);

		//Write the amplitude to the Hamiltonian that will be used in
		//this algorithm.
		hamiltonian[{row, column}] += amplitude;
	}

	//Print the Hamiltonian.
	for(unsigned int row = 0; row < basisSize; row++){
		for(unsigned int column = 0; column < basisSize; column++){
			Streams::out << real(hamiltonian[{row, column}])
				<< "\t";
		}
		Streams::out << "\n";
	}

	return 0;
}
