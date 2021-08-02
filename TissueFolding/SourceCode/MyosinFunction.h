/*
 * MyosinFunction.h
 *
 *  Created on: 21 Sep 2015
 *      Author: melda
 */

#ifndef MYOSINFUNCTION_H_
#define MYOSINFUNCTION_H_


#include <stdio.h>
#include <iostream>
//#include <vector>
using namespace std;

class MyosinFunction{
private:

public:
	int Id;							///< The unique identification number of the myosin function
	int initTime;					///< The application time of the myosin response, in time steps. If the input value(sec) does not fit with the time steps, it will be rounded down.
	bool isApical;					///< The boolean stating if the myosin function is for the apical surface of selected tissue layer
	bool isLateral;					///< The boolean stating if the myosin function is applied laterally
	bool isPolarised;				///< The boolean stating if the myosin function is is polarised, defaul value is false, and applies a uniform contraction of the surface
	bool applyToColumnarLayer;		///< Boolean stating if the myosin levels should be applied to columnar layer
	bool applyToPeripodialMembrane; ///< Boolean stating if the myosin levels should be applied to peripodial membrane
	bool manualStripes;				///< Boolean stating if the myosin levels are based on manual stripes
	bool useEllipses;				///< Boolean stating if the myosin levels are based on marker ellipses
	double stripeSize1, stripeSize2;
	double initialPoint, endPoint;
	double manualCMyoEq;
	double manualOrientation[2];
	double ellipseLateralcEq;
	double ellipsecEq;
	int nGridX;						///< The number of grid points that discretise the tissue in x
	int nGridY;						///< The number of grid points that discretise the tissue in y
	double **EquilibriumMyosinMatrix;	///<The matrix of equilibrium myosin levels. It is a matrix of doubles for equilibrium concentration at each grid point. The dimensions of the matrix are equal to (MyosinFunction::nGridX, MyosinFunction::nGridY), and set in constructor of the MyosinFunction.
	double ***MyosinOrientationMatrix;	///<The matrix of equilibrium myosin orientations. It is a matrix of double triplets for equilibrium orientation at each grid point. The dimensions of the matrix are equal to (MyosinFunction::nGridX [3], MyosinFunction::nGridY [3]), and set in constructor of the MyosinFunction.
	//double kMyo = 0.09873; 				///<myosin update rate constant, in units (1/time) time beeing the same units as dt, currently sec;
	MyosinFunction(int id, bool isApical, bool isPolarised, int initTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, int nX, int nY, double** MyoMat, double** angleMat){
		/**
		 * Define the inputs!
		 */
		//this ->kMyo = 0.09873;
		this->manualStripes = false;
		this->useEllipses = false;
		this->isLateral = false;
		this->Id = id;
		this->initTime = initTime;
		this->applyToColumnarLayer = applyToColumnarLayer;
		this->applyToPeripodialMembrane = applyToPeripodialMembrane;
		this->isApical = isApical;
		this->isPolarised = isPolarised;
		this ->nGridX = nX;
		this ->nGridY = nY;
		EquilibriumMyosinMatrix = new double*[(const int) nGridX];
		for (int i=0; i<nGridX; ++i){
			EquilibriumMyosinMatrix[i] = new double[(const int) nGridY];
			for (int j=0; j<nGridY; ++j){
				EquilibriumMyosinMatrix[i][j] = MyoMat[i][j];
			}
		}
		MyosinOrientationMatrix = new double**[(const int) nGridX];
		for (int i=0; i<nGridX; ++i){
			MyosinOrientationMatrix[i] = new double*[(const int) nGridY];
			for (int j=0; j<nGridY; ++j){
				// PI = 3.14159265359, converting angle to degrees
				MyosinOrientationMatrix[i][j] = new double[2];
				double tet = angleMat[i][j] * 3.14159265359/180.0;
				MyosinOrientationMatrix[i][j] =  new double[2];
				MyosinOrientationMatrix[i][j][0] = cos(tet);
				MyosinOrientationMatrix[i][j][1] = sin(tet);
			}
		}
		//Filling up unused parameters
		this ->stripeSize1 = 0;
		this ->stripeSize2 = 0;
		this ->initialPoint = 0;
		this ->endPoint = 0;
		this ->manualCMyoEq = 0;
		this->manualOrientation[0] = 1;
		this->manualOrientation[1] = 0;
		this->ellipseLateralcEq = 0;
		this->ellipsecEq = 0;
	} ///< The constructor of MyosinFunction.
	MyosinFunction(int id, bool isApical, bool isPolarised, int initTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, double stripeSize1, double stripeSize2, double initialPoint, double endPoint, double manualcEq, double tetha){
		/**
		 *  Define the inputs!
		 */
		//this ->kMyo = 0.09873;
		this->manualStripes = true;
		this->useEllipses = false;
		this->isLateral = false;
		this->Id = id;
		this->initTime = initTime;
		this->applyToColumnarLayer = applyToColumnarLayer;
		this->applyToPeripodialMembrane = applyToPeripodialMembrane;
		this->isApical = isApical;
		this->isPolarised = isPolarised;
		this ->nGridX = 0;
		this ->nGridY = 0;
		this ->stripeSize1 = stripeSize1;
		this ->stripeSize2 = stripeSize2;
		this ->initialPoint = initialPoint;
		this ->endPoint = endPoint;
		this ->manualCMyoEq = manualcEq;
		tetha = tetha * 3.14159265359/180.0;
		manualOrientation[0] = cos(tetha);
		manualOrientation[1] = sin(tetha);
		//Filling up unused parameters
		EquilibriumMyosinMatrix = new double*[(const int) nGridX];
		for (int i=0; i<nGridX; ++i){
			EquilibriumMyosinMatrix[i] = new double[(const int) nGridY];
			for (int j=0; j<nGridY; ++j){
				EquilibriumMyosinMatrix[i][j] = 0.0;
			}
		}
		MyosinOrientationMatrix = new double**[(const int) nGridX];
		for (int i=0; i<nGridX; ++i){
			MyosinOrientationMatrix[i] = new double*[(const int) nGridY];
			for (int j=0; j<nGridY; ++j){
				MyosinOrientationMatrix[i][j] = new double[2];
				MyosinOrientationMatrix[i][j] =  new double[2];
				MyosinOrientationMatrix[i][j][0] = 0;
				MyosinOrientationMatrix[i][j][1] = 0;
			}
		}
		this->ellipseLateralcEq = 0;
		this->ellipsecEq = 0;
	} ///< The constructor of MyosinFunction.


	MyosinFunction(int id, bool isApical, bool isLateral, int initTime, bool applyToColumnarLayer, bool applyToPeripodialMembrane, double ellipseLateralcEq, double ellipsecEq){
		this->manualStripes = false;
		this->useEllipses = true;
		this->isLateral = isLateral;
		this->Id = id;
		this->initTime = initTime;
		this->applyToColumnarLayer = applyToColumnarLayer;
		this->applyToPeripodialMembrane = applyToPeripodialMembrane;
		this->isApical = isApical;
		this->isPolarised = isPolarised;
		this->ellipseLateralcEq = ellipseLateralcEq;
		this->ellipsecEq = ellipsecEq;
		//Filling up unused parameters
		this ->nGridX = 0;
		this ->nGridY = 0;
		EquilibriumMyosinMatrix = new double*[(const int) nGridX];
		for (int i=0; i<nGridX; ++i){
			EquilibriumMyosinMatrix[i] = new double[(const int) nGridY];
			for (int j=0; j<nGridY; ++j){
				EquilibriumMyosinMatrix[i][j] = 0.0;
			}
		}
		MyosinOrientationMatrix = new double**[(const int) nGridX];
		for (int i=0; i<nGridX; ++i){
			MyosinOrientationMatrix[i] = new double*[(const int) nGridY];
			for (int j=0; j<nGridY; ++j){
				MyosinOrientationMatrix[i][j] = new double[2];
				MyosinOrientationMatrix[i][j] =  new double[2];
				MyosinOrientationMatrix[i][j][0] = 0;
				MyosinOrientationMatrix[i][j][1] = 0;
			}
		}
		this ->stripeSize1 = 0;
		this ->stripeSize2 = 0;
		this ->initialPoint = 0;
		this ->endPoint = 0;
		this ->manualCMyoEq = 0;
		this->manualOrientation[0] = 1;
		this->manualOrientation[1] = 0;
	}

	~MyosinFunction(){
		for (int i=0; i<nGridX; ++i){
			for (int j=0; j<nGridY; ++j){
				delete[] MyosinOrientationMatrix[i][j];
			}
		}
		for (int i=0; i<nGridX; ++i){
			delete[] MyosinOrientationMatrix[i];
			delete[] EquilibriumMyosinMatrix[i];
		}
		delete[] MyosinOrientationMatrix;
		delete[] EquilibriumMyosinMatrix;
	};


	int getGridX(){
		return nGridX;
	}///< This function returns MyosinFunction#nGridX.

	int getGridY(){
		return nGridY;
	}///< This function returns MyosinFunction#nGridY.

	double** getEquilibriumMyoMatrix(){
		return EquilibriumMyosinMatrix;
	}///< This function returns MyosinFunction#GrowthMatrix.

	double getEquilibriumMyoMatrixElement(int i, int j){
		return EquilibriumMyosinMatrix[i][j];
	}///< This function returns a the equilibrium myosin concentration at grid point [i]\[j\] (in dimensions MyosinFunction#nGridX, MyosinFunction#nGridY).

	double*** getOrientationMatrix(){
		return MyosinOrientationMatrix;
	}///< This function returns MyosinFunction#MyosinOrientationMatrix.

	double getOrientationMatrixElement(int i, int j, int k){
		return MyosinOrientationMatrix[i][j][k];
	}///< This function returns the myosin orientation at grid point [i]\[j\] (in dimensions MyosinFunction#nGridX, MyosinFunction#nGridY as in normalised vector [ DV axis (x), AP axis (y)] ).

	void setOrientationMatrixElement(double ex, double ey, int i, int j){
		MyosinOrientationMatrix[i][j][0] = ex;
		MyosinOrientationMatrix[i][j][1] = ey;
	}///< This function sets the myosin orientation at grid point [i]\[j\] (in dimensions MyosinFunction#nGridX, MyosinFunction#nGridY), to the growth rate [ex, ez] in the format [ DV axis (x), AP axis (y)] The vector MUST be normalised.

	void writeSummary(ofstream &saveFileSimulationSummary){
		/**
		 *  This function will write the MyosinFunction details into the simulation summary file, provided as the first input.
		 *  The output should look like: \n
		 *			Growth Type:  growth From File (3)
		 *			Initial time(sec): GridBasedGrowthFunction#initTime	FinalTime time(sec): GridBasedGrowthFunction#endTime	Growth matrix mesh size: GridBasedGrowthFunction#nGridX GridBasedGrowthFunction#nGridY
		 */
		saveFileSimulationSummary<<"MyosinFunction:"<<endl;
		saveFileSimulationSummary<<"	Initiation time(sec): ";
		saveFileSimulationSummary<<initTime;
		saveFileSimulationSummary<<"	Is it uni-polar myosin: ";
		saveFileSimulationSummary<<isPolarised;
		saveFileSimulationSummary<<"	Applied to apical layer: ";
		saveFileSimulationSummary<<isApical;
		saveFileSimulationSummary<<"	Applied to columnar layer: ";
		saveFileSimulationSummary<<applyToColumnarLayer;
		saveFileSimulationSummary<<"	Applied to peripodial membrane: ";
		saveFileSimulationSummary<<applyToPeripodialMembrane;
		saveFileSimulationSummary<<"	Equilibrium myosin concentration matrix mesh size: ";
		saveFileSimulationSummary<<nGridX <<" "<<nGridY<<endl;
	}///< The
};



#endif /* MYOSINFUNCTION_H_ */
