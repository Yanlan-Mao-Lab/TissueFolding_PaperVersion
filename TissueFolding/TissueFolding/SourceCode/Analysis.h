/*
 * Analysis.h
 *
 *  Created on: 25 Jul 2016
 *      Author: melda
 */

#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include <vector>
#include "Node.h"
#include "ShapeBase.h"
using namespace std;

class Analysis{

protected:
   //double getApicalSideLengthAverage();
	int Dim ;		//dimensions of the system, should be 3D.
	int nNodes;		//number of nodes of the system to be analysed
	int nElements; //number of elements of the system to be analysed

	ofstream saveFileVolumeMaps;
	ofstream saveFileLenghtMeasurements;
	ofstream saveFileKinkPositions;
	ofstream saveFileFoldMarkers;
	ofstream saveFileCircumference;
	string saveDirectoryString;
	vector <int> apicalContourLineDVNodeIds;
	vector <int> basalContourLineDVNodeIds;
	//vector<vector<bool> > connectivityMap;
	bool** connectivityMap;
	void setUpConnectivityMap(vector<Node*> &nodes);
public:
	double yPosForSideDVLine;
	double relativeYPosSideDVLine;
	vector<double> apicalContourLineDVSelectedYPositionsX;
	vector<double> apicalContourLineDVSelectedYPositionsZ;
	Analysis(int dim, string saveDirectoryToDisplayString,vector<Node*> &nodes, double boundingBoxWidth);
	~Analysis();

	void calculateBoundingBoxSizeAndAspectRatio(int timeInSec,double boundingBoxLength, double boundingBoxWidth);
	void sortPositionMinToMax(vector<Node*> &nodes, int axisToSortWith, vector <int> &linkToArrayToSort );
	void sortPointsMinToMaxBasedFirstArray(vector<double> &x, vector<double> &z, vector<int> &baseNodeId0, vector <int> &baseNodeId1);
	void sortPointsMinToMaxBasedOnInitialPos(vector<double> &base, vector<double> &x, vector<double> &z);
	void setUpContourLinesDV(vector<Node*> &nodes, double boundingBoxWidth);
	void calculateContourLineLengthsDV(vector<Node*> &nodes);
	void setUpSideApicalDVContour(vector<Node*> &nodes, double boundingBoxWidth);
	void updateSideContourPosition(double boundingBoxSizeY);
	void findApicalKinkPointsDV(int timeInSec, double boundingBoxXMin,  double boundingBoxLength, double boundingBoxWidth, vector<Node*> &nodes);
	void setUpContourLinesAP(vector<Node*> &nodes);
	void calculateContourLineLengthsAP(vector<Node*> &nodes);
	void calculateTissueVolumeMap(vector<ShapeBase*> &elements, int timeInSec, double boundingBoxXMin, double boundingBoxYMin, double boundingBoxLength, double boundingBoxWidth);
	void saveNodesOnFold(int timeInSec, vector<Node*>& Nodes);
	void saveApicalCircumferencePosition(int timeInSec, vector<Node*> &nodes);

};

#endif /* ANALYSIS_H_ */
