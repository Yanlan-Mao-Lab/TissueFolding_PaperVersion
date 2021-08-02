/*
 * CellMigration.h
 *
 *  Created on: 14 Sep 2016
 *      Author: melda
 */

#ifndef CELLMIGRATION_H_
#define CELLMIGRATION_H_


#include <stdio.h>
#include <iostream>
#include <vector>
#include "Node.h"
#include "ShapeBase.h"

#include <omp.h>


using namespace std;
/**
 *  The Newton-Raphson solver class
 *  */
class CellMigration{
private:
	int numberOfElements; 		//< Number of elements of the system
	double migrationRate;		//< The rate of cell migration as a fraction of curent volume that leaves the element per hour.
	double migrationOrigin[2];	//< The origin (2D) of cell migration on the peripodial membrane
	double r1, r2; 			//<The size of the tissue ellipsoid, r1 is the longer axis.
	double centreInZAtOriginOfMigration; //< The mid point of the tissue at the point where the mimgration is initiated
	double* 	elementAnglesList;
	double* 	listOfLeavingVolumePerStep;
	double*		listOfGainedVolumePerStep;
	int*		listOfNumberOfNeigsEachElementSendsMaterialTo;
	int*		listOfNumberOfNeigsEachElementGetsMaterialFrom;
	double*		listOfRateFractions; //< The list containing the weighing for the migration rate. Colmunar cells will have dampened migration, cells will stop at the edge.
	vector< vector<int> > elementConnectivityMap; //the array of vectors, that gives the

	//double migrationRadius;		//< The radius of the region on peripodial membrane that will initiate migration
	double calculateAngleOfVector(double dx, double dy, double dz);
	double calculateXYPositionFractionScaledToEllipse(double posX, double posY);

public:
	CellMigration(int nElements, double rate); 	//<Constructer of the cell migration events
	~CellMigration();						//<Desturctor of the cell migration events

	void assignElementConnectivity(vector <Node*>& Nodes, vector <ShapeBase*>& Elements);
	void assignOriginOfMigration(vector <ShapeBase*>& Elements, double originAngle);
	void assignElementRadialVectors(vector <ShapeBase*>& Elements);
	void updateMigrationLists(vector <ShapeBase*>& Elements, double dt);
	void updateVolumesWithMigration(vector <ShapeBase*>& Elements);
	void updateMigratingElements(vector <ShapeBase*>& Elements);
	void generateListOfRateFractions(vector <ShapeBase*>& Elements);
	double getRateFractionForElement(int i);
	double getMigrationAngleForElement(int i);


};


#endif /* CELLMIGRATION_H_ */
