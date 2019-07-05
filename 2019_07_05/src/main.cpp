/* Copyright 2019 Kristofer Björnson
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

#include "TBTK/BrillouinZone.h"
#include "TBTK/Model.h"
#include "TBTK/Plotter.h"
#include "TBTK/Property/DOS.h"
#include "TBTK/PropertyExtractor/BlockDiagonalizer.h"
#include "TBTK/Range.h"
#include "TBTK/Solver/BlockDiagonalizer.h"
#include "TBTK/Streams.h"
#include "TBTK/UnitHandler.h"
#include "TBTK/Vector3d.h"

using namespace std;
using namespace TBTK;
using namespace Plot;

complex<double> i(0, 1);

int main(int argc, char **argv){
	//Set the natural units for this calculation.
	UnitHandler::setScales({"1 C", "1 pcs", "1 eV", "1 Ao", "1 K", "1 s"});

	//Define parameters.
	double t = 3;	//eV
	double a = 2.5;	//Ångström
	unsigned int BRILLOUIN_ZONE_RESOLUTION = 1000;
	vector<unsigned int> numMeshPoints = {
		BRILLOUIN_ZONE_RESOLUTION,
		BRILLOUIN_ZONE_RESOLUTION
	};
	const int K_POINTS_PER_PATH = 100;
	const double ENERGY_LOWER_BOUND = -10;
	const double ENERGY_UPPER_BOUND = 10;
	const int ENERGY_RESOLUTION = 1000;

	//Setup lattice vector. By using three three-dimensional vectors
	//instead of two two-dimensional vectors, the cross product expression
	//for the reciprocal vectors can be expressed in terms of cross
	//products.
	Vector3d r[3];
	r[0] = Vector3d({a,	0,		0});
	r[1] = Vector3d({-a/2,	a*sqrt(3)/2,	0});
	r[2] = Vector3d({0,	0,		a});

	//Define nearest neighbor vectors for A site.
	Vector3d r_AB[3];
	r_AB[0] = (r[0] + 2*r[1])/3.;
	r_AB[1] = -r[1] + r_AB[0];
	r_AB[2] = -r[0] - r[1] + r_AB[0];

	//Calculate the reciprocal lattice vectors.
	Vector3d k[3];
	for(unsigned int n = 0; n < 3; n++){
		k[n] = 2*M_PI*r[(n+1)%3]*r[(n+2)%3]/(
			Vector3d::dotProduct(r[n], r[(n+1)%3]*r[(n+2)%3])
		);
	}

	//Setup the BrillouinZone.
	BrillouinZone brillouinZone(
		{
			{k[0].x, k[0].y},
			{k[1].x, k[1].y}
		},
		SpacePartition::MeshType::Nodal
	);

	//Create mesh.
	vector<vector<double>> mesh = brillouinZone.getMinorMesh(
		numMeshPoints
	);

	//Setup model.
	Model model;
	for(unsigned int n = 0; n < mesh.size(); n++){
		//Get the Index representation of the current k-point.
		Index kIndex = brillouinZone.getMinorCellIndex(
			mesh[n],
			numMeshPoints
		);

		//Calculate the matrix element.
		Vector3d k({mesh[n][0], mesh[n][1], 0});
		complex<double> h_01 = -t*(
			exp(-i*Vector3d::dotProduct(k, r_AB[0]))
			+ exp(-i*Vector3d::dotProduct(k, r_AB[1]))
			+ exp(-i*Vector3d::dotProduct(k, r_AB[2]))
		);

		//Add the matrix element to the model.
		model << HoppingAmplitude(
			h_01,
			{kIndex[0], kIndex[1], 0},
			{kIndex[0], kIndex[1], 1}
		) + HC;
	}
	model.construct();

	//Setup the solver.
	Solver::BlockDiagonalizer solver;
	solver.setModel(model);
	solver.run();

	//Setup the property extractor.
	PropertyExtractor::BlockDiagonalizer propertyExtractor(solver);
	propertyExtractor.setEnergyWindow(
		ENERGY_LOWER_BOUND,
		ENERGY_UPPER_BOUND,
		ENERGY_RESOLUTION
	);

	//Calculate the density of states.
	Property::DOS dos = propertyExtractor.calculateDOS();
	Plotter plotter;
	plotter.setLabelX("Energy");
	plotter.setLabelY("DOS");
	plotter.plot(dos, 0.03);
	plotter.save("figures/DOS.png");

	//Define high symmetry points.
	Vector3d Gamma({0,		0,			0});
	Vector3d M({M_PI/a,		-M_PI/(sqrt(3)*a),	0});
	Vector3d K({4*M_PI/(3*a), 	0,			0});

	//Define paths between high symmetry points.
	vector<vector<Vector3d>> paths = {
		{Gamma, M},
		{M, K},
		{K, Gamma}
	};

	//Calculate the band structure along the path Gamma -> M -> K -> Gamma.
	Array<double> bandStructure({2, 3*K_POINTS_PER_PATH}, 0);
	Range interpolator(0, 1, K_POINTS_PER_PATH);
	for(unsigned int p = 0; p < 3; p++){
		//Select the start and end points for the current path.
		Vector3d startPoint = paths[p][0];
		Vector3d endPoint = paths[p][1];

		//Loop over a single path.
		for(unsigned int n = 0; n < K_POINTS_PER_PATH; n++){
			//Interpolate between the paths start and end point.
			Vector3d k = (
				interpolator[n]*endPoint
				+ (1 - interpolator[n])*startPoint
			);

			//Get the Index representation of the current k-point.
			Index kIndex = brillouinZone.getMinorCellIndex(
				{k.x, k.y},
				numMeshPoints
			);

			//Extract the eigenvalues for the current k-point.
			bandStructure[{0, n+p*K_POINTS_PER_PATH}] = propertyExtractor.getEigenValue(kIndex, 0);
			bandStructure[{1, n+p*K_POINTS_PER_PATH}] = propertyExtractor.getEigenValue(kIndex, 1);
		}
	}

	//Find max and min value for the band structure.
	double min = bandStructure[{0, 0}];
	double max = bandStructure[{1, 0}];
	for(unsigned int n = 0; n < 3*K_POINTS_PER_PATH; n++){
		if(min > bandStructure[{0, n}])
			min = bandStructure[{0, n}];
		if(max < bandStructure[{1, n}])
			max = bandStructure[{1, n}];
	}

	//Plot the band structure.
	plotter.clear();
	plotter.setHold(true);
	plotter.setLabelX("k");
	plotter.setLabelY("Energy");
	plotter.plot(bandStructure.getSlice({0, _a_}));
	plotter.plot(bandStructure.getSlice({1, _a_}));
	plotter.plot({K_POINTS_PER_PATH, K_POINTS_PER_PATH}, {min, max});
	plotter.plot({2*K_POINTS_PER_PATH, 2*K_POINTS_PER_PATH}, {min, max});
	plotter.save("figures/BandStructure.png");

	return 0;
}
