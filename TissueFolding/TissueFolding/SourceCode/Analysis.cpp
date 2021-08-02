/*
 * Analysis.cpp
 *
 *  Created on: 25 Jul 2016
 *      Author: melda
 */



#include "Analysis.h"

Analysis::Analysis(int dim, string saveDirectoryToDisplayString, vector<Node*> &nodes, double boundingBoxWidth){
	Dim = dim;
	nNodes = nodes.size();
	yPosForSideDVLine = 0.1;
	relativeYPosSideDVLine = 0.3;
	saveDirectoryString  = saveDirectoryToDisplayString;
	string saveFileString = saveDirectoryString +"/Save_AnalysisVolumeMaps";
	const char* name_saveFileVolumeMaps = saveFileString.c_str();;
	saveFileVolumeMaps.open(name_saveFileVolumeMaps, ofstream::out);
	if (!(saveFileVolumeMaps.good() && saveFileVolumeMaps.is_open())){
		cerr<<"Cannot open the save file to write analysis: "<<name_saveFileVolumeMaps<<endl;
	}
	string saveFileLenghtMeasurementsString = saveDirectoryString +"/Save_AnalysisLengthMeasurements";
	const char* name_saveFileLenghtMeasurements = saveFileLenghtMeasurementsString.c_str();;
	saveFileLenghtMeasurements.open(name_saveFileLenghtMeasurements, ofstream::out);
	if (!(saveFileLenghtMeasurements.good() && saveFileLenghtMeasurements.is_open())){
		cerr<<"Cannot open the save file to write analysis: "<<name_saveFileLenghtMeasurements<<endl;
	}
	string saveFileKinkPositionsString = saveDirectoryString +"/Save_AnalysisKinkPositions";
	const char* name_saveFileKinkPositions = saveFileKinkPositionsString.c_str();;
	saveFileKinkPositions.open(name_saveFileKinkPositions, ofstream::out);
	if (!(saveFileKinkPositions.good() && saveFileKinkPositions.is_open())){
		cerr<<"Cannot open the save file to write analysis: "<<name_saveFileKinkPositions<<endl;
	}
	string saveFileCircumferenceString = saveDirectoryString +"/Save_AnalysisCircumference";
	const char* name_saveFileCircumference = saveFileCircumferenceString.c_str();;
	saveFileCircumference.open(name_saveFileCircumference, ofstream::out);
	if (!(saveFileCircumference.good() && saveFileCircumference.is_open())){
		cerr<<"Cannot open the save file to write analysis: "<<name_saveFileCircumference<<endl;
	}
	string saveFileFoldMarkersString = saveDirectoryString +"/Save_AnalysisNodesOnFold";
	const char* name_saveFileFoldMarkers = saveFileFoldMarkersString.c_str();
	saveFileFoldMarkers.open(name_saveFileFoldMarkers, ofstream::out);
	if (!(saveFileFoldMarkers.good() && saveFileFoldMarkers.is_open())){
		cerr<<"Cannot open the save file to write analysis: "<<name_saveFileFoldMarkers<<endl;
	}
	setUpContourLinesDV(nodes, boundingBoxWidth);
	setUpConnectivityMap(nodes);
}


Analysis::~Analysis(){
	for (int i=0; i<nNodes; ++i){
		delete[] connectivityMap[i];
	}
    delete[] connectivityMap;
}

void Analysis::saveApicalCircumferencePosition(int timeInSec, vector<Node*> &nodes){
	//cout<<"inside save circ" <<endl;
	//[time] [nodeId] [at-midline] [x] [y] [z]
	for (vector<Node*>::iterator itNode = nodes.begin(); itNode<nodes.end(); ++itNode){
		if ((*itNode)->atCircumference && (*itNode)->tissuePlacement == 1){
			saveFileCircumference<<timeInSec<<" ";
			saveFileCircumference<<(*itNode)->Id<<" ";
			if ((*itNode)->Position[1] == 0){
				saveFileCircumference<<1<<" ";
			}
			else{
				saveFileCircumference<<0<<" ";
			}
			saveFileCircumference<<(*itNode)->Position[0]<<" ";
			saveFileCircumference<<(*itNode)->Position[1]<<" ";
			saveFileCircumference<<(*itNode)->Position[2]<<endl;
		}
		if (!(*itNode)->atCircumference && (*itNode)->tissuePlacement == 1 && (*itNode)->Position[1] == 0 ){
			saveFileCircumference<<timeInSec<<" ";
			saveFileCircumference<<(*itNode)->Id<<" ";
			saveFileCircumference<<1<<" ";
			saveFileCircumference<<(*itNode)->Position[0]<<" ";
			saveFileCircumference<<(*itNode)->Position[1]<<" ";
			saveFileCircumference<<(*itNode)->Position[2]<<endl;
		}

	}
}

void Analysis::setUpConnectivityMap(vector<Node*> &nodes){
	connectivityMap = new bool* [(const int) nNodes];
	for (int i=0; i<nNodes; ++i){
		connectivityMap[i] = new bool [(const int) nNodes];
		std::fill(connectivityMap[i],connectivityMap[i]+nNodes,false);
		connectivityMap[i][i] = true;
	}
	for (vector<Node*>::iterator itNode = nodes.begin(); itNode<nodes.end(); ++itNode){
		int n = (*itNode)->immediateNeigs.size();
		int id = (*itNode)->Id;
		for (int i=0; i<n; ++i){
			int currNeigId = nodes[(*itNode)->immediateNeigs[i]]->Id;
			connectivityMap[currNeigId][id] = true;
			connectivityMap[id][currNeigId] = true;
		}
	}
	cerr<<"finished connectivity map"<<endl;
}

void Analysis::calculateBoundingBoxSizeAndAspectRatio(int timeInSec, double boundingBoxLength, double boundingBoxWidth){
	// format is [time in sec] [DV length in microns] [AP length in microns] [aspect ratio DV/AP] [...] then contour lengths will follow in coming functions.
	saveFileLenghtMeasurements<<timeInSec<<" ";
	saveFileLenghtMeasurements.width(6);
	saveFileLenghtMeasurements.precision(3);
	saveFileLenghtMeasurements<<boundingBoxLength<<" ";
	saveFileLenghtMeasurements.width(6);
	saveFileLenghtMeasurements.precision(3);
	saveFileLenghtMeasurements<<boundingBoxWidth<<" ";
	saveFileLenghtMeasurements.width(6);
	saveFileLenghtMeasurements.precision(3);
	saveFileLenghtMeasurements<<boundingBoxLength/boundingBoxWidth<<" ";

}

void Analysis::setUpContourLinesDV(vector<Node*> &nodes, double boundingBoxWidth){
	double posThreshold = 10E-2;
	double negThreshold = (-1.0)*posThreshold;
	for (int i=0 ; i<nNodes; ++i){
		if (nodes[i]->tissuePlacement == 0 ){ //tissue placement 0 is for basal nodes
			if (nodes[i]->Position[1] > negThreshold && nodes[i]->Position[1] < posThreshold){
				basalContourLineDVNodeIds.push_back(i);
			}
		}
		if (nodes[i]->tissuePlacement == 1 ){ //tissue placement 1 is for apical nodes
			if (nodes[i]->Position[1] > negThreshold && nodes[i]->Position[1] < posThreshold){
				apicalContourLineDVNodeIds.push_back(i);
			}
		}
	}
	setUpSideApicalDVContour(nodes, boundingBoxWidth);
	sortPositionMinToMax(nodes, 0, basalContourLineDVNodeIds); //sort basalContourLineDVNodeIds with x position going from minimum to maximum
	sortPositionMinToMax(nodes, 0, apicalContourLineDVNodeIds); //sort apicalContourLineDVNodeIds with x position going from minimum to maximum
	cout<<" Apical contour Ids: ";
	for (int i=0; i<apicalContourLineDVNodeIds.size(); ++i){
		cout<<apicalContourLineDVNodeIds[i]<<" ";
	}
	cout<<endl;
}


void Analysis::calculateContourLineLengthsDV(vector<Node*> &nodes){
	double currContourLength = 0;
	//calculating contour length for basal DV axis:
	int n = basalContourLineDVNodeIds.size();
	int id0 = basalContourLineDVNodeIds[0];
	for (int i=1;i<n;++i){
		int id1 = basalContourLineDVNodeIds[i];
		double dx = nodes[id0]->Position[0] - nodes[id1]->Position[0];
		double dy = nodes[id0]->Position[1] - nodes[id1]->Position[1];
		double dz = nodes[id0]->Position[2] - nodes[id1]->Position[2];
		double d = dx*dx + dy*dy + dz*dz;
		d = pow(d,0.5);
		currContourLength += d;
		id0 = id1;
	}
	//writing the basal contour length:
	saveFileLenghtMeasurements.width(9);
	saveFileLenghtMeasurements.precision(6);
	saveFileLenghtMeasurements<<currContourLength<<" ";
	//calculating apical contour length:
	currContourLength = 0;
	n = apicalContourLineDVNodeIds.size();
	id0 = apicalContourLineDVNodeIds[0];
	for (int i=1;i<n;++i){
		int id1 = apicalContourLineDVNodeIds[i];
		double dx = nodes[id0]->Position[0] - nodes[id1]->Position[0];
		double dy = nodes[id0]->Position[1] - nodes[id1]->Position[1];
		double dz = nodes[id0]->Position[2] - nodes[id1]->Position[2];
		double d = dx*dx + dy*dy + dz*dz;
		d = pow(d,0.5);
		currContourLength += d;
		id0 = id1;
	}
	//writing the apical contour length:
	saveFileLenghtMeasurements.width(6);
	saveFileLenghtMeasurements.precision(3);
	saveFileLenghtMeasurements<<currContourLength<<endl;

}

void Analysis::updateSideContourPosition(double boundingBoxWidth){
	yPosForSideDVLine = relativeYPosSideDVLine * boundingBoxWidth;
}

void Analysis::setUpSideApicalDVContour(vector<Node*> &nodes, double boundingBoxWidth){
	//cout<<"entered setUpSideApicalDVContour"<<endl;
	apicalContourLineDVSelectedYPositionsX.clear();
	apicalContourLineDVSelectedYPositionsZ.clear();
	updateSideContourPosition(boundingBoxWidth);
	vector<double> initialPositionsInX;
	//cout<<" y cut off positions: "<<yPosForSideDVLine<<" relative to tissue: "<<relativeYPosSideDVLine<<endl;
	for (vector<Node*>::iterator itNode = nodes.begin(); itNode<nodes.end(); ++itNode){
		if ((*itNode)->tissuePlacement == 1){
			//cout<<" checking node "<<(*itNode)->Id<<" # of immediate neigs: "<<(*itNode)->immediateNeigs.size()<<endl;
			//checking only apical nodes
			int nNeig = (*itNode)->immediateNeigs.size();
			double diff0 = (*itNode)->Position[1]-yPosForSideDVLine;
			for (int i=0; i<nNeig; ++i){
				if ((*itNode)->Id < (*itNode)->immediateNeigs[i]
				    && nodes[(*itNode)->immediateNeigs[i]]->tissuePlacement == 1 ){
					//cout<<" checking node "<<(*itNode)->Id<<" aginst node "<<nodes[(*itNode)->immediateNeigs[i]]->Id<<endl;
					//checking only neigs with id numbers larger than mine, therefore have not been checked yet:
					//checking only apical neigs
					int id1 = (*itNode)->immediateNeigs[i];
					double diff1 = nodes[id1]->Position[1]-yPosForSideDVLine;
					if (diff0*diff1 < 0){
						//the two nodes are on different sides of the line!:
						//find the intersection point
						double vec0To1[3];
						for (int k=0; k<3; ++k){
							vec0To1[k] = nodes[id1]->Position[k] - (*itNode)->Position[k];
						}
						double scaleFactor = (yPosForSideDVLine-(*itNode)->Position[1])/vec0To1[1];
						double pointX = (*itNode)->Position[0] + vec0To1[0]*scaleFactor;
						double pointZ = (*itNode)->Position[2] + vec0To1[2]*scaleFactor;
						apicalContourLineDVSelectedYPositionsX.push_back(pointX);
						apicalContourLineDVSelectedYPositionsZ.push_back(pointZ);
						//calculating the same point in the junction in original positions:
						double vec0To1InitialPos[3];
						for (int k=0; k<3; ++k){
							vec0To1InitialPos[k] = nodes[id1]->initialPosition[k] - (*itNode)->initialPosition[k];
						}
						double pointXInitialPos =(*itNode)->initialPosition[0] + vec0To1InitialPos[0]*scaleFactor;
						initialPositionsInX.push_back(pointXInitialPos);
						//cout<<"node: "<<(*itNode)->Id<<" initial pos: "<<(*itNode)->initialPosition[0]<<" "<<(*itNode)->initialPosition[1]<<" "<<(*itNode)->initialPosition[2]<<endl;
						//cout<<"node: "<<(*itNode)->Id<<" current pos: "<<(*itNode)->Position[0]<<" "<<(*itNode)->Position[1]<<" "<<(*itNode)->Position[2]<<endl;
						//cout<<"node: "<<id1<<" initial pos: "<<nodes[id1]->initialPosition[0]<<" "<<nodes[id1]->initialPosition[1]<<" "<<nodes[id1]->initialPosition[2]<<endl;
						//cout<<"node: "<<id1<<" current pos: "<<nodes[id1]->Position[0]<<" "<<nodes[id1]->Position[1]<<" "<<nodes[id1]->Position[2]<<endl;
						//cout<<"vec0To1: "<<vec0To1[0]<<" "<<vec0To1[1]<<" "<<vec0To1[2]<<" vec0To1InitialPos: "<<vec0To1InitialPos[0]<<" "<<vec0To1InitialPos[1]<<" "<<vec0To1InitialPos[2]<<endl;
						//cout<<"scaleFactor "<<scaleFactor<<" selected current point "<<pointX<<" "<<yPosForSideDVLine<<" "<<pointZ<<" selected initial x: "<<pointXInitialPos<<endl;
					}
				}
			}
		}
	}
	sortPointsMinToMaxBasedOnInitialPos(initialPositionsInX, apicalContourLineDVSelectedYPositionsX,apicalContourLineDVSelectedYPositionsZ);
	//sortPointsMinToMaxBasedFirstArray(apicalContourLineDVSelectedYPositionsX, apicalContourLineDVSelectedYPositionsZ, apicalContourLineDVSelectedYNodeCouples0,apicalContourLineDVSelectedYNodeCouples1);
	//int nXpos = apicalContourLineDVSelectedYPositionsX.size();
	//for (int i=0;i<nXpos;++i){
	//	cout<<" side contour: "<<i<<" of "<<apicalContourLineDVSelectedYPositionsX.size()<<": "<<apicalContourLineDVSelectedYPositionsX[i]<<" "<<apicalContourLineDVSelectedYPositionsZ[i]<<" base: "<<initialPositionsInX[i]<<endl;
	//}
	//cout<<"finalised setUpSideApicalDVContour"<<endl;
}

void Analysis::findApicalKinkPointsDV(int timeInSec,  double boundingBoxXMin, double boundingBoxLength, double boundingBoxWidth, vector<Node*> &nodes){
	cout<<" inside findApicalKinkPointsDV"<<endl;
	setUpSideApicalDVContour(nodes,boundingBoxWidth);
	for (int DVLineIndex = 0; DVLineIndex<2; DVLineIndex++){
		//by default calculating for the mid line;
		int n = apicalContourLineDVNodeIds.size();
		if (n<3){
			return;
		}
		int id0 = apicalContourLineDVNodeIds[0];
		int id1 = apicalContourLineDVNodeIds[1];
		if (DVLineIndex ==1){
			n = apicalContourLineDVSelectedYPositionsX.size();
			if (n<3){
				cout<<"too few nodes - cannot continue on contour calculation on the side"<<endl;
				return;
			}
			id0 = 0;
			id1 = 1;
		}
		for (int i=0;i<n-2;++i){
			int id2 = apicalContourLineDVNodeIds[i+2];
			if (DVLineIndex ==1){
				id2 =2;
			}
			double dx,dz;
			double posId1[3];
			if (DVLineIndex == 1){
				dx = apicalContourLineDVSelectedYPositionsX[id0] - apicalContourLineDVSelectedYPositionsX[id1];
				dz = apicalContourLineDVSelectedYPositionsZ[id0] - apicalContourLineDVSelectedYPositionsZ[id1];
				posId1[0]=apicalContourLineDVSelectedYPositionsX[id1];
				posId1[1]=yPosForSideDVLine;
				posId1[2]=apicalContourLineDVSelectedYPositionsZ[id1];
			}
			else{
				//calculating the first slope between node id 0 and node id1:
				dx = nodes[id0]->Position[0] - nodes[id1]->Position[0];
				dz = nodes[id0]->Position[2] - nodes[id1]->Position[2];
				posId1[0]=nodes[id1]->Position[0];
				posId1[1]=nodes[id1]->Position[1];
				posId1[2]=nodes[id1]->Position[2];
			}
			double m01 = dz/dx;
			//the slope for node id1 and node id2
			if (DVLineIndex == 1){
				dx = apicalContourLineDVSelectedYPositionsX[id1] - apicalContourLineDVSelectedYPositionsX[id2];
				dz = apicalContourLineDVSelectedYPositionsZ[id1] - apicalContourLineDVSelectedYPositionsZ[id2];
			}
			else{
				dx = nodes[id1]->Position[0] - nodes[id2]->Position[0];
				dz = nodes[id1]->Position[2] - nodes[id2]->Position[2];
			}
			double m12 = dz/dx;
			double meanSlope = ( m01 + m12 ) * 0.5;
			double midPointX01, midPointX12;
			if (DVLineIndex == 1){
				midPointX01 = (apicalContourLineDVSelectedYPositionsX[id0] + apicalContourLineDVSelectedYPositionsX[id1])*0.5;
				midPointX12 = (apicalContourLineDVSelectedYPositionsX[id1] + apicalContourLineDVSelectedYPositionsX[id2])*0.5;

			}
			else{
				midPointX01 = (nodes[id0]->Position[0] + nodes[id1]->Position[0])*0.5;
				midPointX12 = (nodes[id1]->Position[0] + nodes[id2]->Position[0])*0.5;
			}
			double secondDerivative = (m01 - m12)/(midPointX01 - midPointX12);
			double radiusofCurvature = pow((1 + meanSlope*meanSlope ),1.5)/secondDerivative;
			double relativeX;
			if (DVLineIndex == 1){
				relativeX = (apicalContourLineDVSelectedYPositionsX[id1] - boundingBoxXMin)/boundingBoxLength;
			}
			else{
				relativeX = (nodes[id1]->Position[0] - boundingBoxXMin)/boundingBoxLength;
			}
			id0 = id1;
			id1 = id2;
			double yPos = 0;
			if (DVLineIndex ==1){
				yPos = yPosForSideDVLine;
			}
			saveFileKinkPositions<<timeInSec<<" ";
			saveFileKinkPositions<<id1<<" ";
			saveFileKinkPositions.width(6);
			saveFileKinkPositions.precision(3);
			saveFileKinkPositions<<relativeX<<" "<<yPos<<" "<<radiusofCurvature<<" "<<posId1[0]<<" "<<posId1[1]<<" "<<posId1[2]<<endl;
		}
	}
	cout<<" finished findApicalKinkPointsDV "<<endl;
}

void Analysis::sortPositionMinToMax(vector<Node*> &nodes, int axisToSortWith, vector <int> &linkToArrayToSort ){
	int n = linkToArrayToSort.size();
	bool swapped = true;
	while (swapped){
		swapped = false;
		for(int i=1; i<n; ++i){
			double pos1 = nodes[linkToArrayToSort[i]]->Position[axisToSortWith];
			double pos0 = nodes[linkToArrayToSort[i-1]]->Position[axisToSortWith];
			if(pos1<pos0){
				int temp=linkToArrayToSort[i-1];
				linkToArrayToSort[i-1]=linkToArrayToSort[i];
				linkToArrayToSort[i]=temp;
				swapped = true;
			}
		}
	}
}


void Analysis::sortPointsMinToMaxBasedOnInitialPos(vector<double> &base, vector<double> &x, vector<double> &z){
	int n = x.size();
	bool swapped = true;
	while (swapped){
		swapped = false;
		for(int i=1; i<n; ++i){
			double pos1 = base[i];
			double pos0 = base[i-1];
			if(pos1<pos0){
				double tempX=x[i-1];
				double tempZ=z[i-1];
				double tmpBase = base[i-1];
				x[i-1]=x[i];
				x[i]=tempX;
				z[i-1]=z[i];
				z[i]=tempZ;
				base[i-1]=base[i];
				base[i]=tmpBase;
				swapped = true;
			}
		}
	}
}
void Analysis::sortPointsMinToMaxBasedFirstArray(vector<double> &x, vector<double> &z, vector<int> &baseNodeId0, vector <int> &baseNodeId1 ){
	int n = x.size();
	bool swapped = true;
	while (swapped){
		swapped = false;
		for(int i=1; i<n; ++i){
			double pos1 = x[i];
			double pos0 = x[i-1];
			if(pos1<pos0){
				double tempX=x[i-1];
				double tempZ=z[i-1];
				double tempBase0 = baseNodeId0[i-1];
				double tempBase1 = baseNodeId1[i-1];
				x[i-1]=x[i];
				x[i]=tempX;
				z[i-1]=z[i];
				z[i]=tempZ;
				baseNodeId0[i-1]=baseNodeId0[i];
				baseNodeId0[i]=tempBase0;
				baseNodeId1[i-1]=baseNodeId1[i];
				baseNodeId1[i]=tempBase1;
				swapped = true;
			}
		}
	}
}

void Analysis::calculateTissueVolumeMap	(vector<ShapeBase*> &elements, int timeInSec, double boundingBoxXMin, double boundingBoxYMin, double boundingBoxLength, double boundingBoxWidth){
	//Format:
	//[time] [element id] [relative pos in X] [relative pos in Y] [reference shape volume] [current ideal volume] [current emergent volume]
	double totTissueIdealVolume = 0;
	double totTissueEmergentVolume = 0;
	vector<ShapeBase*>::iterator itEle;
	for (itEle=elements.begin(); itEle<elements.end(); ++itEle){
		//Write if explicit ECM
		//double* c = new double[3];
		//c = (*itEle)->getCentre();
    	(*itEle)->calculateRelativePosInBoundingBox(boundingBoxXMin, boundingBoxYMin,boundingBoxLength, boundingBoxWidth);
    	double* ReletivePos = new double[2];
    	(*itEle)->getRelativePosInBoundingBox(ReletivePos);
    	double currEmergentVolume = (*itEle) -> calculateCurrentGrownAndEmergentVolumes();
		double currIdealVolume = (*itEle)->GrownVolume;
		totTissueIdealVolume  += currIdealVolume;
		totTissueEmergentVolume  += currEmergentVolume;
		(*itEle)->calculateEmergentShapeOrientation();
		saveFileVolumeMaps<<timeInSec<<" ";
		saveFileVolumeMaps<<(*itEle) ->Id;
		saveFileVolumeMaps<<" "<<ReletivePos[0]<<" "<<ReletivePos[1]<<" ";
		saveFileVolumeMaps<<(*itEle)->ReferenceShape->Volume<<" ";
		saveFileVolumeMaps<<currIdealVolume<<" ";
		saveFileVolumeMaps<<currEmergentVolume<<" ";
		saveFileVolumeMaps<<(*itEle)->emergentShapeLongAxis[0]<<" "<<(*itEle)->emergentShapeLongAxis[1]<<" ";
		saveFileVolumeMaps<<(*itEle)->emergentShapeShortAxis[0]<<" "<<(*itEle)->emergentShapeShortAxis[1]<<" ";
		saveFileVolumeMaps<<(*itEle)->isECMMimicing<<" ";
		saveFileVolumeMaps<<endl;
		delete[] ReletivePos;
	}
	cout<<" Time: "<<timeInSec<<" total tissue ideal volume: "<<totTissueIdealVolume<<" total tissue emergent volume:" <<totTissueEmergentVolume<<endl;
}

void Analysis::saveNodesOnFold(int timeInSec, vector<Node*>& Nodes){
	int foldIds[(const int)nNodes];
	std::fill(foldIds,foldIds+nNodes,-1);
	int nextFoldId = 0;
	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		int id = (*itNode)->Id;
		if ((*itNode)->onFoldInitiation){
			int activeFoldId = nextFoldId;
			vector<int> idsToChange;
			//now go through connectivity map to see if any neigs have been assigned anything so far
			for (int i=0; i<id; ++i){
				if (connectivityMap[id][i]){
					if (foldIds[i] != -1){
						activeFoldId = foldIds[i];
						idsToChange.push_back(activeFoldId);
					}
				}
			}
			foldIds[id] =activeFoldId;
			int n = idsToChange.size();
			for (int i=0; i<n; ++i){
				for (int j=0; j<nNodes; ++j){
					if (foldIds[j] == idsToChange[i]){
						foldIds[j] = activeFoldId;
					}
				}
			}
			if (activeFoldId == nextFoldId){
				nextFoldId++;
			}
		}
	}

	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->onFoldInitiation){
			saveFileFoldMarkers.width(10);
			saveFileFoldMarkers<<timeInSec;
			saveFileFoldMarkers.width(10);
			saveFileFoldMarkers<<(*itNode)->Id;
			saveFileFoldMarkers.width(10);
			saveFileFoldMarkers.precision(3);
			saveFileFoldMarkers<<(*itNode)->Position[0];
			saveFileFoldMarkers.width(10);
			saveFileFoldMarkers.precision(3);
			saveFileFoldMarkers<<(*itNode)->Position[1];
			saveFileFoldMarkers.width(10);
			saveFileFoldMarkers.precision(3);
			saveFileFoldMarkers<<(*itNode)->Position[2];
			saveFileFoldMarkers.width(10);
			saveFileFoldMarkers.precision(3);
			saveFileFoldMarkers<<foldIds[(*itNode)->Id];
			saveFileFoldMarkers<<endl;

			cout.width(10);
			cout<<timeInSec;
			cout.width(10);
			cout<<(*itNode)->Id;
			cout.width(10);
			cout.precision(3);
			cout<<(*itNode)->Position[0];
			cout.width(10);
			cout.precision(3);
			cout<<(*itNode)->Position[1];
			cout.width(10);
			cout.precision(3);
			cout<<(*itNode)->Position[2];
			cout.width(10);
			cout.precision(3);
			cout<<foldIds[(*itNode)->Id];
			cout<<endl;

		}
	}
}


