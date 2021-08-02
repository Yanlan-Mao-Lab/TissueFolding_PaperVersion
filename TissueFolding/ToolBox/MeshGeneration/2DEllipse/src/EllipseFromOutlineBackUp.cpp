using namespace std;
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdio.h>
#include "math.h"
#include <vector>

class EllipseLayoutGenerator{
public:
	EllipseLayoutGenerator(double r1_1, double r1_2, double r2_1, double r2_2, double sideLen, bool ysymmetry, double condensation);
	double pi;
	double r1[2];
	double r2[2];
	double r1Init[2];
	double r2Init[2];
	double LinkerR1[4];	//Inner and outer radia of the loop, 4 needed as they are ellipses
	double LinkerR2[2];	//The width of the loop, from mid point to apical side and to basal side
	double condensation;
	double sideLen;
	double sideLenInit; 
	vector <double> posx,posy;//,posz;
	vector <int> tissueType; //0: columnar layer, 1: peripodium
	vector <int> atBorder;	//1 id the node is at the circumference of the setup, 0 otherwise
	double tetha;
	double dtet[4];
	double Linkerdtet[4];
	double Circumference[4];
    	double LinkerCircumference[4];
	double linkerSideLen;
	int Vborder[4];
	int LinkerVborder[4];
	vector <double> LinkerOuterRingPosz,LinkerOuterRingPosx,LinkerInnerRingPosz,LinkerInnerRingPosx;
	bool tethaAtZero;
	bool r1failed;	
	vector <int> Links0,Links1;	
	vector <int*>triangles;
	vector <int> linkerLinks0,linkerLinks1;	
	vector <int*>linkerTriangles;
	vector <double> trianglePeripodialness;
	vector <double> LinkerPosz, LinkerPosx;
	vector <int>   LinkerAtBorder;
	int 		nNonDuplicateLinkers;
	vector <bool>   duplicateLinkerNode;
	vector <int> 	LinkerTissuePlacement;
	vector <double> LinkerVectorsOnCircumferenceX;
	vector <double> LinkerVectorsOnCircumferenceY;
	vector <int> SortedCircumferenceForLinkers;
	double dForNodeGeneration;
	bool addPeripodial;
	int peripodialLayers;
	double lumenHeight;
	double peripodialHeight;
	bool symmetricY;

	void calculateCircumference();
	void calculateCurrentBorderNumber();
	void calculatedtet();
	bool updateRadia();
	void addMidLine();
	void Tesselate2D();
	void readInTesselation2D();
	void linkerTesselate2D();
	void linkerReadInTesselation2D();
	void writeVectors2D(ofstream &vectorsForGnuplot);
	void writeNodes2D(ofstream &nodesForGnuplot);
	void writeMeshFileForSimulation(double zHeight, int zLayers);
	void writeTissueWeights(ofstream& MeshFile, vector <double> &recordedPeripodialness);
	void writeTriangleMeshFileForSimulation(double zHeight);
	void addEquidistantRingMin();
	void minimiseCurrentPositions(bool addingPeripodial, int n , double d, int iterLen, int iterWid, vector<double> &CurrPosX, vector <double> &CurrPosY);
	void cleanUpCurrentPositions(double d, vector<double> &CurrPosX, vector <double> &CurrPosY);
	void cleanUpAgainsAlreadyExisting(double d, vector<double> &CurrPosX, vector <double> &CurrPosY);
	void addRectangle();
	void calculateAverageSideLength();
	void calculatePeripodialMembraneParameters(double ABHeight, int ABLayers, double peripodialHeightFrac, double lumenHeightFrac, double peripodialSideCurveFrac);
	void calculatePeripodialAttachVectors(bool symmetricY);	
	void extractBorder(bool symmetricY, vector<double>& x, vector <double>& y, vector <int> &ids);
	void sortPositions(bool symmetricY, vector<double>& x, vector <double>& y, vector <int> &ids);
	void calculateVectors(bool symmetricY, vector<double>& vecx, vector <double>& vecy, vector <int>& ids );

	void calculateLinkerCircumference();
	void calculateLinkerBorderNumber();
	void calculateLinkerdtet();
	void calculateLinkerRings();
	void translateLinkers(double dz, double dx);
	void calculateTrianglePeripodialness(double midZ);
	void scaleInputOutline(vector <float>& x, vector <float>&y);
	void calculateBoundingBox(double* boundingBox, vector <float>&  x, vector <float>& y);
	void bringThicknessCentreToOrigin(vector <float>&  x, vector <float>& y);
	void smooth(int pointAvr, vector <float>&  x, vector <float>& y);
	void scaleToTissueSize(double* boundingBox, vector <float>&  x, vector <float>& y);
	void reOrderToFixNodesAtZeroTips(vector <float>&  x, vector <float>& y);
	void arrangeSideLength(vector <float>&  x, vector <float>& y);
	void clearSymmetric(vector <float>&  x, vector <float>& y);
	void bluntTip(vector <float>&  x, vector <float>& y);
	void addLayersToOutline(int nLayer, double* bounidngBox, vector <float>&  x, vector <float>& y);
	void writeInPosVectors(vector <float>&  x, vector <float>& y);
};

EllipseLayoutGenerator::EllipseLayoutGenerator(double r1_1, double r1_2, double r2_1, double r2_2, double sideLen, bool ysymmetry, double condensation){
	pi = 3.14;
	this->r1[0] = r1_1;
	this->r1[1] = r1_2;
	this->r2[0] = r2_1;
	this->r2[1] = r2_2;
	this->r1Init[0] = r1_1;
	this->r1Init[1] = r1_2;
	this->r2Init[0] = r2_1;
	this->r2Init[1] = r2_2;

	this->condensation = condensation;
	this->sideLen = sideLen*condensation;
	this->sideLenInit =this->sideLen;
	this->symmetricY = ysymmetry;
	tethaAtZero=true;
	r1failed = false;
	//peripodial related parameter setting
	addPeripodial = false;
	peripodialLayers = 0;
	lumenHeight = 0.0;
	peripodialHeight = 0.0;
	LinkerR1[0] = 0.0;
	LinkerR1[1] = 0.0;
	LinkerR2[0] = 0.0;
	LinkerR2[1] = 0.0;
	nNonDuplicateLinkers= 0;
}

void EllipseLayoutGenerator::calculateCircumference(){
	for (int iterSide = 0 ;iterSide <4; iterSide++){
		int iterLen,iterWid;
		if (iterSide == 0){iterLen = 0; iterWid = 0;}		//iterSide == 0 is putting in (+)ve x - (+)ve y 
		else if (iterSide == 1){iterLen = 1; iterWid = 0;}	//iterSide == 1 is putting in (-)ve x - (+)ve y 
		else if (iterSide == 2){iterLen = 1; iterWid = 1;}	//iterSide == 2 is putting in (-)ve x - (-)ve y 
		else if (iterSide == 3){iterLen = 0; iterWid = 1;}	//iterSide == 3 is putting in (+)ve x - (-)ve y 
		double h = ((r1[iterLen]-r2[iterWid])*(r1[iterLen]-r2[iterWid]))/((r1[iterLen]+r2[iterWid])*(r1[iterLen]+r2[iterWid]));
		Circumference[iterSide] = pi * (r1[iterLen] + r2[iterWid]) * (1.0 + (3*h / (10.0 + pow(4.0 - 3*h,0.5) ))); 
	}
}

void EllipseLayoutGenerator::calculateCurrentBorderNumber(){
	for (int i =0; i<4; ++i){
		Vborder[i] = Circumference[i] / (sideLen*1.8);
	}
}

void EllipseLayoutGenerator::calculatedtet(){
	for (int i =0; i<4; ++i){
		int n  = Vborder[i]/2;	
		dtet[i] = pi/n;
	}
}

void EllipseLayoutGenerator::linkerTesselate2D(){
	ofstream pointsForTesselation;
	pointsForTesselation.open("./LinkerPoints.poly",ofstream::trunc);
	//writing nodes:	
	int nLinker = LinkerOuterRingPosz.size() + LinkerInnerRingPosz.size();
	pointsForTesselation<< nLinker<<" 2 0 1"<<endl; //dim, attribute number , border markers on or of
	int counter = 0;	
	for (int i =0; i< LinkerOuterRingPosz.size(); ++i){
		pointsForTesselation<<counter<<" "<<LinkerOuterRingPosz[i]<<" "<<LinkerOuterRingPosx[i]<<endl;
		counter++;
		cout<<"OuterRing: "<<i<<" : "<<LinkerOuterRingPosz[i]<<" "<<LinkerOuterRingPosx[i]<<endl;
	}
	for (int i =0; i< LinkerInnerRingPosz.size(); ++i){
		pointsForTesselation<<counter<<" "<<LinkerInnerRingPosz[i]<<" "<<LinkerInnerRingPosx[i]<<endl;
		counter++;
		cout<<"InnerRing: "<<i<<" : "<<LinkerInnerRingPosz[i]<<" "<<LinkerInnerRingPosx[i]<<endl;
	}
	//writing segments:
	//the first set of segments will link all the outer nodes, that is n-1 segments, 
	//the second set will be the same for inner layer
	//then there is the segments coundint the two layers
	//There will be one hole, looping over the inner side
	int nSeg = nLinker-2+2;
	pointsForTesselation<< nSeg<<" 1"<<endl; //nSeg, border markers on or of
	counter = 0;
	//writing the segments of outer loop
	int nOuter = LinkerOuterRingPosz.size();
	for (int i =0; i< nOuter-1; ++i){
		pointsForTesselation<<counter<<" "<<i<<" "<<i+1<<endl;
		counter++;
	}
	//writing the segments of the inner loop
	int nInner = LinkerInnerRingPosz.size();
	for (int i =0; i< nInner-1; ++i){
		pointsForTesselation<<counter<<" "<<nOuter+i<<" "<<nOuter+i+1<<endl;
		counter++;
	}
	//writing the connections between layers:
	pointsForTesselation<<counter<<" "<<nOuter-1<<" "<<nOuter+nInner-1<<endl;
	counter++;
	pointsForTesselation<<counter<<" "<<nOuter<<" "<<0<<endl;
	//writing holes:
	pointsForTesselation<<"1"<<endl; //1 hole
	//hole loops over the inner side:
	for (int i =0; i< nInner-1; ++i){
		pointsForTesselation<<0<<" "<<nOuter+i<<" "<<nOuter+i+1<<endl;
		counter++;
	}
	pointsForTesselation<<0<<" "<<nOuter+nInner-1<<" "<<nOuter<<endl;
	//finished writing the holes
	double maxArea = (linkerSideLen/condensation) * (linkerSideLen/condensation)*pow(3.0,0.5)/2.0 * 1/2.0;
	maxArea = maxArea*1.5;
	cerr<<"Max Area: "<<maxArea<<endl;

	ostringstream Convert;
	Convert << maxArea; // Use some manipulators
	string maxAreaStr = Convert.str(); // Give the result to the string
	string sysCommand = "/home/melda/Documents/TissueFolding/ToolBox/MeshGeneration/triangle/triangle -pq25a"+maxAreaStr+" ./LinkerPoints.poly  ";
	cerr<<"Running triangulation with: "<<sysCommand<<endl;
	system(sysCommand.c_str());
}

void EllipseLayoutGenerator::linkerReadInTesselation2D(){
	cout<<" start of read in linker tesselation, size of points: "<<posx.size()<<" size of tissue shape: "<<tissueType.size()<<endl;	
	ifstream MeshFile;	
	MeshFile.open("./LinkerPoints.1.ele", ifstream::in);
	//readHeader:
	int ntri;
	MeshFile>>ntri;
	int nodesPerTri;
	MeshFile>>nodesPerTri;
	int nAttributes;
	MeshFile>>nAttributes;
	for (int i=0; i<ntri; ++i){
		int triId;
		MeshFile>>triId;
		int* pnts;
		pnts = new int[3];
		for (int j=0;j<3;++j){
			MeshFile >> pnts[j];
		}
		linkerTriangles.push_back(pnts);
		linkerLinks0.push_back(pnts[0]);
		linkerLinks1.push_back(pnts[1]);
		linkerLinks0.push_back(pnts[1]);
		linkerLinks1.push_back(pnts[2]);
		linkerLinks0.push_back(pnts[2]);
		linkerLinks1.push_back(pnts[0]);
	}
	MeshFile.close();
	ifstream NodeFile;	
	NodeFile.open("./LinkerPoints.1.node", ifstream::in);
	int nNode;
	int bordersMarked;
	NodeFile>>nNode;
	NodeFile>>nodesPerTri;	//dimensions
	NodeFile>>nAttributes;	//attribute number
	NodeFile>>bordersMarked;	//border markers on or off
	for (int i=0; i<nNode; ++i){
		int nodeId;
		NodeFile>>nodeId;
		//cerr<<"Reading node : "<<i<<"("<<nodeId<<") of "<<nNode<<endl;
		double pnts[2];
		for (int j=0;j<2;++j){
			NodeFile >> pnts[j];	
		}
		//reading the flag for border nodes, and generating the vector for it
		NodeFile>>bordersMarked;
		LinkerAtBorder.push_back(bordersMarked);
		//cerr<<"		read pos: "<< pnts[0]<<" "<<pnts[1]<<" borders? "<<bordersMarked<<endl;
		LinkerPosz.push_back(pnts[0]);	
		LinkerPosx.push_back(pnts[1]);
		
	}
	cout<<" end of read in tesselation, size of points: "<<posx.size()<<" size of tissue shape: "<<tissueType.size()<<endl;
	NodeFile.close();
}


void EllipseLayoutGenerator::Tesselate2D(){
	ofstream vectorsForGnuplot,nodesForGnuplot;
	nodesForGnuplot.open("./NodesPreTesselation.out",ofstream::trunc);	
	writeNodes2D(nodesForGnuplot);
	ofstream pointsForTesselation;
	pointsForTesselation.open("./Points.node",ofstream::trunc);
	pointsForTesselation<<posx.size()<<" 2 0 1"<<endl; //dim, attribute number , border markers on or of
	for (int i =0; i< posx.size(); ++i){
		pointsForTesselation<<i<<" "<<posx[i]<<" "<<posy[i]<<endl;
	}
	double maxArea = (sideLen/condensation) * (sideLen/condensation)*pow(3.0,0.5)/2.0 * 1/2.0;
	maxArea = maxArea*1.5;
	cerr<<"Max Area: "<<maxArea<<endl;


	ostringstream Convert;
	Convert << maxArea; // Use some manipulators
	string maxAreaStr = Convert.str(); // Give the result to the string
	string sysCommand = "/home/melda/Documents/TissueFolding/ToolBox/MeshGeneration/triangle/triangle -q33a"+maxAreaStr+" ./Points.node  ";
	cerr<<"Running triangulation with: "<<sysCommand<<endl;
	system(sysCommand.c_str());
}

void EllipseLayoutGenerator::readInTesselation2D(){
	cout<<" start of read in tesselation, size of points: "<<posx.size()<<" size of tissue shape: "<<tissueType.size()<<endl;
		
	ifstream MeshFile;	
	MeshFile.open("./Points.1.ele", ifstream::in);
	//readHeader:
	int ntri;
	MeshFile>>ntri;
	int nodesPerTri;
	MeshFile>>nodesPerTri;
	int nAttributes;
	MeshFile>>nAttributes;
	for (int i=0; i<ntri; ++i){
		int triId;
		MeshFile>>triId;
		int* pnts;
		pnts = new int[3];
		for (int j=0;j<3;++j){
			MeshFile >> pnts[j];
		}
		triangles.push_back(pnts);
		Links0.push_back(pnts[0]);
		Links1.push_back(pnts[1]);
		Links0.push_back(pnts[1]);
		Links1.push_back(pnts[2]);
		Links0.push_back(pnts[2]);
		Links1.push_back(pnts[0]);
	}
	MeshFile.close();
	ifstream NodeFile;	
	NodeFile.open("./Points.1.node", ifstream::in);
	int nNode;
	int bordersMarked;
	NodeFile>>nNode;
	NodeFile>>nodesPerTri;	//dimensions
	NodeFile>>nAttributes;	//attribute number
	NodeFile>>bordersMarked;	//border markers on or off
	for (int i=0; i<nNode; ++i){
		int nodeId;
		NodeFile>>nodeId;
		//cerr<<"Reading node : "<<i<<"("<<nodeId<<") of "<<nNode<<endl;
		double pnts[2];
		for (int j=0;j<2;++j){
			NodeFile >> pnts[j];	
		}
		//reading the flag for border nodes, and generating the vector for it
		NodeFile>>bordersMarked;
		atBorder.push_back(bordersMarked);
		//cerr<<"		read pos: "<< pnts[0]<<" "<<pnts[1]<<" borders? "<<bordersMarked<<endl;
		if (nodeId<posx.size()){
			//cerr<<"		overwriting existing placement on vector"<<endl;
			posx[nodeId]=pnts[0];	
			posy[nodeId]=pnts[1];
			tissueType[nodeId] = 0; //all nodes are columnar at this stage	
		}
		else{
			//cerr<<"		adding points to system"<<endl;
			posx.push_back(pnts[0]);	
			posy.push_back(pnts[1]);
			tissueType.push_back(0); //all nodes are columnar at this stage		
		}
	}
	cout<<" end of read in tesselation, size of points: "<<posx.size()<<" size of tissue shape: "<<tissueType.size()<<endl;
	
	NodeFile.close();
	//writing for gnuplot vectors:
	ofstream vectorsForGnuplot,nodesForGnuplot;;
	vectorsForGnuplot.open("./VectorsPostTesselation.out",ofstream::trunc);
	nodesForGnuplot.open("./NodesPostTesselation.out",ofstream::trunc);
	writeVectors2D(vectorsForGnuplot);	
	writeNodes2D(nodesForGnuplot);
}

void EllipseLayoutGenerator::writeVectors2D(ofstream &vectorsForGnuplot){
	for (int i =0; i< Links1.size(); ++i){
		double x = posx[Links0[i]];
		double y = posy[Links0[i]];
		double dx =  posx[Links1[i]] - x;
		double dy =  posy[Links1[i]] - y;
		vectorsForGnuplot<<x<<"  "<<y<<" "<<dx<<" "<<dy<<endl;
	}	
}
void EllipseLayoutGenerator::writeNodes2D(ofstream &nodesForGnuplot){
	for (int i =0; i< posx.size(); ++i){
		double x = posx[i];
		double y = posy[i];
		nodesForGnuplot<<x<<"  "<<y<<endl;
	}	
}

void EllipseLayoutGenerator::calculateAverageSideLength(){
	int nTri = triangles.size();
	double SumSide =0.0;	
	for (int i =0; i<nTri; ++i){
		int IDcorner0, IDcorner1;
		int node0, node1;
		double dx, dy, side;		
		for (int j=0; j<3; ++j){		
			if (j==0){IDcorner0 = 0, IDcorner1 = 1;}
			else if (j==1){IDcorner0 = 1, IDcorner1 = 2;}
			else if (j==2){IDcorner0 = 2, IDcorner1 = 0;}
			node0 = triangles[i][IDcorner0];
			node1 = triangles[i][IDcorner1];
			dx = posx[node0] - posx[node1];
			dy = posy[node0] - posy[node1];
			side = pow(dx*dx + dy*dy , 0.5);
			SumSide +=side;
		}
	}
	SumSide /= (3.0 * nTri);
	cerr<<"Average side length of a triangle: "<<SumSide<<" -- desired length was: "<<sideLen/condensation<<endl;
}

void EllipseLayoutGenerator::calculatePeripodialAttachVectors(bool symmetricY){
	//generate copies;
	vector<double> x, y;
	vector <int> ids;
	extractBorder(symmetricY, x,y,ids);
	sortPositions(symmetricY,x,y,ids);
	calculateVectors(symmetricY,x,y,ids);
	for(int i=0; i<x.size(); ++i){
		LinkerVectorsOnCircumferenceX.push_back(x[i]);
		LinkerVectorsOnCircumferenceY.push_back(y[i]);
		SortedCircumferenceForLinkers.push_back(ids[i]);
		cout<<"pos: "<<posx[ids[i]]<<" "<<posy[ids[i]]<<" newbase: "<<posx[ids[i]]+x[i]<<" "<<posy[ids[i]]+y[i]<<endl;
	}
}

void EllipseLayoutGenerator::extractBorder(bool symmetricY, vector<double>& x, vector <double>& y, vector <int> &ids){
	int n = posx.size();
	for(int i=0; i<n; ++i){
		if(atBorder[i] == 1){
			x.push_back(posx[i]);
			y.push_back(posy[i]);
			ids.push_back(i);
		}	
	}
	if (symmetricY){
		//find maximum and minimum in x:
		n = x.size();
		int maxId, minId;
		double maxX= -1000, minX = 1000;
		for(int i=0; i<n; ++i){
			if (x[i] <minX){
				minX  = x[i];
				minId = ids[i];				
			}
			if (x[i] >maxX){
				maxX  = x[i];
				maxId = ids[i];				
			}	
		}
		//I need to keep these two points, and remove all remaining nodes at y = 0;
		int k = 0;		
		while(k<n){
			if (y[k]<1E-4 && y[k]>-1E-4){
				//this node is at the y=0 border, check if it is one of hte end points:
				if ( ids[k] == minId || ids[k] == maxId){
					//no need to do anything it is one of the borders				
				} 		
				else{
					x.erase (x.begin()+k);
					y.erase (y.begin()+k);
					ids.erase (ids.begin()+k);
					k--;
					n--;					
				}
			}
			k++;
		}
	}
}

void EllipseLayoutGenerator::sortPositions(bool symmetricY, vector<double>& x, vector <double>& y, vector <int> &ids){
	int n=x.size();	
	//calculate centre:
	double c[2] = {0.0,0.0};
	for (int j =0 ; j<n; ++j){
		c[0] += x[j];	
		c[1] += y[j];
	}
	c[0] /= n;
	c[1] /= n;
	if (symmetricY){
		c[1] = 0;
	};
	//ordering the position list clockwise rotation
	//calculate the angle of each circumference point	
	vector <double> angles;
	for (int j =0 ; j<n; ++j){
		double currX =x[j]-c[0];
		double currY =y[j]-c[1];
		double tet = atan2(currY,currX);
		if (tet<0){tet += 2.0*3.14;}
		angles.push_back(tet);
	}
	//sort id list:
	bool swapped = true;
	while (swapped){
		swapped = false;
		for(int i=1; i<n; ++i){
			if(angles[i]<angles[i-1]){
				int temp=ids[i-1];
				ids[i-1]=ids[i];
				ids[i]=temp;
				double td = angles[i-1];
				angles[i-1]=angles[i];
				angles[i]=td;
				swapped = true;
				//cout<<"swapped "<<i <<" with "<<i-1<<endl;
			}
		}
	}
}

void EllipseLayoutGenerator::calculateVectors(bool symmetricY, vector<double>& vecx, vector <double>& vecy, vector <int>& ids ){
	//calculatecentre:
	int n=ids.size();	
	//calculate centre:
	double c[2] = {0.0,0.0};
	for (int j =0 ; j<n; ++j){
		c[0] += posx[ids[j]];	
		c[1] += posy[ids[j]];
	}
	c[0] /= n;
	c[1] /= n;
	if (symmetricY){
		c[1] = 0;
	};
	int nodeId0;
	int nodeId1;
	int nodeId2;
	for (int i =0 ; i<n; ++i){
		if (i == 0){
			if (!symmetricY){
				nodeId0 = ids[i];
				nodeId1 = ids[i+1];
				nodeId2 = ids[n-1];
			}
			else{
				vecx[i] = 1.0;
				vecy[i] = 0.0;
				continue;				
			}
		}else if (i == n-1){
			if (!symmetricY){
				nodeId0 = ids[i];
				nodeId1 = ids[0];
				nodeId2 = ids[i-1];
			}
			else{
				vecx[i] = -1.0;
				vecy[i] = 0.0;
				continue;			
			}
		}
		else{
			nodeId0 = ids[i];
			nodeId1 = ids[i+1];
			nodeId2 = ids[i-1];
		}
		double vec1[2];
		double vec2[2];
		vec1[0] = (posx[nodeId0] - posx[nodeId1]);
		vec1[1] = (posy[nodeId0] - posy[nodeId1]);
		vec2[0] = (posx[nodeId0] - posx[nodeId2]);
		vec2[1] = (posy[nodeId0] - posy[nodeId2]);
		//normalise vectors: 
		double mag = pow(vec1[0]*vec1[0] + vec1[1]*vec1[1],0.5);
		vec1[0] /= mag; vec1[1] /= mag;
		mag = pow(vec2[0]*vec2[0] + vec2[1]*vec2[1],0.5);
		vec2[0] /= mag; vec2[1] /= mag;
		//cout<<"nodeids: "<< nodeId0<<" "<<nodeId1<<" "<<nodeId2<<" vec1: "<<vec1[0]<<" "<<vec1[1]<<" vec2: "<<vec2[0]<<" "<<vec2[1]<<endl;
		//sum:
		vec1[0] += vec2[0];
		vec1[1] += vec2[1];
		mag = vec1[0]*vec1[0] + vec1[1]*vec1[1];
		//cout<<"sum: "<<vec1[0]<<" "<<vec1[1]<<" mag: "<<mag<<endl;
		if (mag < 1E-5){
			//the two nodes are linear, the resulting vector is of zero length.
			//I will rotate one of the vectors 90 degrees and it will be pointing in the correct orientation, the actual direction can be fixed in the below loop:
			//this vector is already normalised, I am skipping the normalisation step
			vec1[0] = -1.0*vec2[1];
			vec1[1] = vec2[0];
		}
		else{
			mag = pow(mag,0.5);
			vec1[0] /= mag; vec1[1] /= mag;
		}
		//now I have the vector to point out from the base node 0. BUT, this will point outwards only if the tissue curvature is convex at all points
		//I need to check if it actually is pointing out, as the experimentally driven tissues can be concave at points.
		//the cross produc of vec[2] and the vector to the cell centre should have the opposite sign with the corss product of my orientation vector and vector 2.
		double vecCentre[2];
		vecCentre[0] = c[0] - posx[nodeId0] ;
		vecCentre[1] = c[1] - posy[nodeId0] ;
		mag = pow(vecCentre[0]*vecCentre[0] + vecCentre[1]*vecCentre[1],0.5);
		vecCentre[0] /= mag; vecCentre[1] /= mag;
		double crossZ1 = vec2[0]*vecCentre[1] - vec2[1]*vecCentre[0];
		double crossZ2 = vec2[0]*vec1[1] - vec2[1]*vec1[0];
		if (crossZ1*crossZ2 > 0 ){
			//the vectors are pointing to the same direction! Need to rotate vec1 180 degrees:
			vec1[0] *= -1.0;
			vec1[1] *= -1.0;
		}
		vecx[i] = vec1[0];
		vecy[i] = vec1[1];
	}
}

void EllipseLayoutGenerator::calculatePeripodialMembraneParameters(double ABHeight, int ABLayers, double peripodialHeightFrac, double lumenHeightFrac, double peripodialSideCurveFrac){
	lumenHeight = lumenHeightFrac*ABHeight;
	peripodialHeight = peripodialHeightFrac*ABHeight;
    	linkerSideLen = ABHeight / ABLayers;	
	peripodialLayers = ceil(peripodialHeight/linkerSideLen );
	double dzPer = peripodialHeight/peripodialLayers;
	double dzCol = linkerSideLen;
	double tissueTop = ABHeight*(1+ lumenHeightFrac + peripodialHeightFrac);
	double midPoint = ABHeight* (1+ 0.5*lumenHeightFrac);
	LinkerR1[0] = tissueTop-midPoint;
	LinkerR1[1] = midPoint;
	LinkerR1[2] = ABHeight * 0.5 * lumenHeightFrac;
	LinkerR1[3] = ABHeight * 0.5 * lumenHeightFrac;
	LinkerR2[0] = ABHeight * peripodialSideCurveFrac;
	LinkerR2[1] = ABHeight * 0.5 * lumenHeightFrac;
	calculateLinkerCircumference();
	cout<<"LinkerSideLen: "<<linkerSideLen<<"LinkerR1: "<<LinkerR1[0]<<" "<<LinkerR1[1]<<" "<<LinkerR1[2]<<" "<<LinkerR1[3]<<" LinkerR2 "<<LinkerR2[0]<<" "<<LinkerR2[1]<<endl;
	cout<<"circumferences: "<<LinkerCircumference[0]<<" "<<LinkerCircumference[1]<<" "<<LinkerCircumference[2]<<" "<<LinkerCircumference[3]<<endl;
	// add the nodes on peripodial side to outer ring, go from apical to basal (do not add the borders, only the mid range)
	// then we will ad the outer loop, then we will add the nodes corresponding to columnar side, starting from basal side.
	//the purpose is to draw a loop.
	for (int i=0; i<peripodialLayers-1; ++i){
		LinkerOuterRingPosz.push_back(LinkerR1[2]+(i+1)*dzPer);
		LinkerOuterRingPosx.push_back(0);
	}
	calculateLinkerBorderNumber();
	calculateLinkerdtet();
	calculateLinkerRings();	
	for (int i=0; i<ABLayers-1; ++i){
		LinkerOuterRingPosz.push_back(-1.0*midPoint+(i+1)*dzCol);
		LinkerOuterRingPosx.push_back(0);
	}
	translateLinkers(midPoint,0.0); //dz, dx
	linkerTesselate2D();
	linkerReadInTesselation2D();
	calculateTrianglePeripodialness(midPoint);
	for (int i=0; i<LinkerPosz.size(); ++i){
		cout<<"linker pos after triangulation: "<<LinkerPosz[i]<<" "<<LinkerPosx[i]<<endl;
	}
	for (int i=0; i<LinkerPosz.size(); ++i){
		if (LinkerAtBorder[i] == 1){
			if (LinkerPosx[i]>-1E-4 && LinkerPosx[i]<1E-4){
				//this is a node that will correspond to a columnar node or peripodial cap node.
				cout<<" LinkerPos: "<<LinkerPosx[i]<<" "<<LinkerPosz[i]<<" declared duplicate"<<endl;
				duplicateLinkerNode.push_back(true);
			}
			else{
				//this is a node that is unique
				duplicateLinkerNode.push_back(false);
				cout<<" LinkerPos: "<<LinkerPosx[i]<<" "<<LinkerPosz[i]<<" border but not  duplicate"<<endl;
				nNonDuplicateLinkers++;		
			}
			//this is a border node on the linker list, see which border it belongs to:
			//if the distance to midline is within smallRadius (+ -) 10%SideLen, then it is apical 
			double dx = (LinkerPosx[i]-0.0);
			double dz = (LinkerPosz[i]-midPoint);
			double d = pow (dx*dx + dz*dz,0.5);
			if (d<LinkerR1[2]+0.1*linkerSideLen){
				LinkerTissuePlacement.push_back(1);	//apical	
			}
			else{
				LinkerTissuePlacement.push_back(0);	//basal
			}

		}
		else{
		//node is not at border, it is not duplicate, and it is midline:
			duplicateLinkerNode.push_back(false);
			LinkerTissuePlacement.push_back(2);
			nNonDuplicateLinkers++;
			cout<<" LinkerPos: "<<LinkerPosx[i]<<" "<<LinkerPosz[i]<<" not border not  duplicate"<<endl;
		}		
	}
}

void EllipseLayoutGenerator::calculateTrianglePeripodialness(double midZ){
	int n = linkerTriangles.size();
	for (int i=0; i<n; ++i){
		//cout<<" calculating peripodialness for triangle: "<<i<<" of "<<n<<endl;
		int id0= linkerTriangles[i][0];
		int id1= linkerTriangles[i][1];
		int id2= linkerTriangles[i][2];
		double centre[2] = {1.0/3.0*(LinkerPosx[id0]+LinkerPosx[id1]+LinkerPosx[id2]), 1.0/3.0*(LinkerPosz[id0]+LinkerPosz[id1]+LinkerPosz[id2])};
		//cout<<" midpoint: "<<midZ<<endl;
		//cout<<" pos: "<<endl;
		//cout<<LinkerPosx[id0]<<" "<<LinkerPosz[id0]<<endl;
		//cout<<LinkerPosx[id1]<<" "<<LinkerPosz[id1]<<endl;
		//cout<<LinkerPosx[id2]<<" "<<LinkerPosz[id2]<<endl;
		//the vector from the centre of the linker loops ({0, midZ}) to the triangle centre (centre[2])		
		double vec[2] = {centre[0], centre[1]-midZ};
		//cout<<"vec: "<<vec[0]<<" "<<vec[1]<<endl;
		double mag = pow(vec[0] * vec[0] + vec[1]*vec[1],0.5);
		vec[0] /= mag;
		vec[1] /= mag;
		//cout<<"after normalisation, vec: "<<vec[0]<<" "<<vec[1]<<" (mag: "<<mag<<")"<<endl;
		//dot product of this vector with the vector (0, -1):		
		double dot = -1.0*vec[1];
		//the corresponding angle:
		double tet = acos(dot);
		//cout<<"dotp: "<<dot<<" angle: "<<tet<<endl;
		//from 0 to 180, the angle will determine the peripodial character of the linker zone.
		//For angle = 0. tissue is columanr, for angle 90, tissue is 50% periipodial, for angle 180, it is peripodial.
		double peripodialness = tet/3.14159265359;
		trianglePeripodialness.push_back(peripodialness);

		cout<<"centre: "<<centre[0]<<" "<<centre[1]<<" tet: "<<tet<<" peripodialness: "<<peripodialness<<endl;
	} 
}

void EllipseLayoutGenerator::calculateLinkerCircumference(){
    for (int iterSide = 0 ;iterSide <4; iterSide++){
        int iterLen,iterWid;
        if (iterSide == 0){iterLen = 0; iterWid = 0;}		//iterSide == 0 is putting in (+)ve x - (+)ve y outer layer (peripodial side)
        else if (iterSide == 1){iterLen = 1; iterWid = 0;}	//iterSide == 1 is putting in (-)ve x - (+)ve y outer layer (columnar side)
        else if (iterSide == 2){iterLen = 2; iterWid = 1;}	//iterSide == 2 is putting in (+)ve x - (+)ve y inner layer (peripodial side)
        else if (iterSide == 3){iterLen = 3; iterWid = 1;}	//iterSide == 3 is putting in (-)ve x - (+)ve y inner layer (columnar side)
        double h = ((LinkerR1[iterLen]-LinkerR2[iterWid])*(LinkerR1[iterLen]-LinkerR2[iterWid]))/((LinkerR1[iterLen]+LinkerR2[iterWid])*(LinkerR1[iterLen]+LinkerR2[iterWid]));
        LinkerCircumference[iterSide] = pi * (LinkerR1[iterLen] + LinkerR2[iterWid]) * (1.0 + (3*h / (10.0 + pow(4.0 - 3*h,0.5) )));
    }

}

void EllipseLayoutGenerator::calculateLinkerBorderNumber(){
    for (int i =0; i<4; ++i){
        LinkerVborder[i] = LinkerCircumference[i] / (linkerSideLen*1.8);
    }
}

void EllipseLayoutGenerator::calculateLinkerdtet(){
	for (int i =0; i<4; ++i){
		int n  = LinkerVborder[i]/2;
		Linkerdtet[i] = pi/n;
	}
}

void EllipseLayoutGenerator::calculateLinkerRings(){
	int iterLen, iterWid;
	for (int i = 0; i<4; i++){ //this is for left and right sides of outer(i=0,1) or inner (i = 2,3) loop
		if (i == 0){iterLen = 0; iterWid = 0;}		//iterSide == 0 is putting in (+)ve x - (+)ve y (outer loop)
		else if (i == 1){iterLen = 1; iterWid = 0;}	//iterSide == 1 is putting in (-)ve x - (+)ve y (outer loop)
		else if (i == 2){iterLen = 2; iterWid = 1;} //iterSide == 0 is putting in (+)ve x - (+)ve y (inner loop)
		else if (i == 3){iterLen = 3; iterWid = 1;}	//iterSide == 1 is putting in (-)ve x - (+)ve y (inner loop)
		double CQ = LinkerCircumference[i]/4.0;
		int n= CQ / linkerSideLen;
		Linkerdtet[i]= pi/2.0/n;
		double d = CQ / n;
		cout<<" circumference of ring: "<<CQ<<" number of nodes: "<<n<<" dtet: "<<Linkerdtet[i]<<" d: "<<d<<endl;

		dForNodeGeneration  = d;
		double x0 = LinkerR1[i];
		double y0 = 0;
		vector <double> CurrPosX, CurrPosY;
		CurrPosX.push_back(x0);
		CurrPosY.push_back(y0);
		tetha = Linkerdtet[i];
		for (int j= 0; j<n-1; ++j){
			double x = LinkerR1[iterLen]*cos(tetha);
			double y = LinkerR1[iterWid]*sin(tetha);
			CurrPosX.push_back(x);
			CurrPosY.push_back(y);
			tetha += Linkerdtet[i];
		}
		CurrPosX.push_back(0);
		CurrPosY.push_back(LinkerR2[iterWid]);
		minimiseCurrentPositions(true , n , d, iterLen, iterWid, CurrPosX,CurrPosY);
		cleanUpCurrentPositions(d, CurrPosX,CurrPosY);
		for (int j = 0; j<CurrPosX.size(); ++j){
			cout<<"node "<<j<<" pos: "<<CurrPosX[j]<<" "<<CurrPosY[j]<<endl;
		}
		n = CurrPosX.size();
		int startindex = 0, endIndex = n;
		//rotate the values for nodes towards the columnar side.
		if (i == 1 || i == 3 ){endIndex = n-1;for (int j = 0; j<n; ++j){CurrPosX[j] *= (-1.0);}}
		
		if (i ==0){
			//record outer ring, positive side:
			for (int j = startindex; j<endIndex; ++j){
				LinkerOuterRingPosz.push_back(CurrPosX[j]);
				LinkerOuterRingPosx.push_back(CurrPosY[j]);
			}
		}else if (i ==1){
			//record the outer ring negative side, I will add them in reverse order to keep nodes sorted:
			for (int j = endIndex-1; j>=startindex; --j){
				LinkerOuterRingPosz.push_back(CurrPosX[j]);
				LinkerOuterRingPosx.push_back(CurrPosY[j]);
			}	
		}
		if (i ==2){
			//record inner ring, positive side:
			for (int j = startindex; j<endIndex; ++j){
				LinkerInnerRingPosz.push_back(CurrPosX[j]);
				LinkerInnerRingPosx.push_back(CurrPosY[j]);
			}
		}else if (i ==3){
			//record the inner ring negative side, I will add them in reverse order to keep nodes sorted:
			for (int j = endIndex-1; j>=startindex; --j){
				//record outer ring, positive side:
				LinkerInnerRingPosz.push_back(CurrPosX[j]);
				LinkerInnerRingPosx.push_back(CurrPosY[j]);
			}	
		}
	}
}

void EllipseLayoutGenerator::minimiseCurrentPositions(bool addingPeripodial, int n , double d, int iterLen, int iterWid, vector<double> &CurrPosX, vector <double> &CurrPosY){
	double ksp = 0.2;
	int niteration = 100000;
	for (int a = 0; a< niteration; ++a){
		if (a < 100) {ksp = 0.2;}
		else if (a<500) {ksp = 0.1;}
		else if (a<5000)  {ksp =0.05;}
		else if (a<10000)  {ksp =0.01;}
		else if (a<20000)  {ksp =0.005;}
		else  {ksp =0.001;}
		for (int i= 1; i<n; ++i){
			double vec[2]={0.0,0.0};
			double v[2]={0.0,0.0};
			//pre neig:
			vec[0] = CurrPosX[i-1] - CurrPosX[i];
			vec[1] = CurrPosY[i-1] - CurrPosY[i];	
			double dsq = vec[0]*vec[0] + vec[1]*vec[1];
			double dmag = pow(dsq,0.5);
			double F =  (dmag-d)*ksp;
			v[0] = F*vec[0]/dmag;
			v[1] = F*vec[1]/dmag;
			//cerr<<"vel of node : "<<i<<" "<<v[0]<<" "<<v[1]<<endl;
			//post neig:
			vec[0] = CurrPosX[i+1] - CurrPosX[i];
			vec[1] = CurrPosY[i+1] - CurrPosY[i];	
			dsq = vec[0]*vec[0] + vec[1]*vec[1];
			dmag = pow(dsq,0.5);
			//cerr<<"dmag of node : "<<i<<" "<<dmag<<endl;	
			F =  (dmag-d)*ksp;
			v[0] += F*vec[0]/dmag;
			v[1] += F*vec[1]/dmag;
			//cerr<<"vel of node : "<<i<<" "<<v[0]<<" "<<v[1]<<endl;
			//add energy for deviating from ellipse:
			dsq = CurrPosX[i]*CurrPosX[i] + CurrPosY[i]*CurrPosY[i];
			dmag = pow(dsq,0.5);
			double ellipseSum;
			if (addingPeripodial){
				ellipseSum = (CurrPosX[i]/LinkerR1[iterLen]) * (CurrPosX[i]/LinkerR1[iterLen]) +  (CurrPosY[i]/LinkerR2[iterWid]) * (CurrPosY[i]/LinkerR2[iterWid]);
			}
			else{
				ellipseSum = (CurrPosX[i]/r1[iterLen]) * (CurrPosX[i]/r1[iterLen]) +  (CurrPosY[i]/r2[iterWid]) * (CurrPosY[i]/r2[iterWid]);
			}
			if (ellipseSum > 1){
				ellipseSum = (-1.0)*pow(( ellipseSum-1),0.5);		
			}
			else{
				ellipseSum = pow(( 1.0 - ellipseSum),0.5);	
			}
			F = ellipseSum*ksp;
			v[0] += F*CurrPosX[i]/dmag;
			v[1] += F*CurrPosY[i]/dmag;
			//cerr<<"vel of node : "<<i<<" "<<v[0]<<" "<<v[1]<<endl;	
			CurrPosX[i] += v[0];
			CurrPosY[i] += v[1];
		}
	}
}

void EllipseLayoutGenerator::cleanUpCurrentPositions(double d, vector<double> &CurrPosX, vector <double> &CurrPosY){
	//now cheking iof any of the non-neighbouring nodes are too close
	// I will remove them from the list, therefore I am nt checking first or last ones:
	double threshold = d*0.7;
	threshold *= threshold;	
	int n = CurrPosX.size();
	int i=1;	
	while (i < n-1){
		bool deleteNode = false;
		for (int j=0; j< n; ++j){
			if (CurrPosX[i]<0 || CurrPosY[i]<0){deleteNode = true;}
			if (i != j && deleteNode ==false){
				double dx = CurrPosX[i] - CurrPosX[j];
				double dy = CurrPosY[i] - CurrPosY[j];
				double d2 = dx*dx + dy*dy;
				//cerr<<"i: "<<i<<" j "<<j<<" d: "<<d2 <<" threshold : "<<threshold <<" pos: "<<CurrPosX[i]<<" "<<CurrPosY[i]<<endl;
				if (d2 < threshold){deleteNode = true;}
			}
			if(deleteNode){
				vector<double>::iterator it;
  				it = CurrPosX.begin();
				it +=i;
				CurrPosX.erase(it);
				it = CurrPosY.begin();	
				it +=i;
				CurrPosY.erase(it);
				n = CurrPosX.size();
				i = 0;		
				break;				
			}			
		}
		i++;	
	}
}

void EllipseLayoutGenerator::translateLinkers(double dz, double dx){
	int n = LinkerOuterRingPosz.size();
	for (int i=0; i<n; ++i){
		LinkerOuterRingPosz[i] += dz;
		LinkerOuterRingPosx[i] += dx;
	}
	n = LinkerInnerRingPosz.size();
	for (int i=0; i<n; ++i){
		LinkerInnerRingPosz[i] += dz;
		LinkerInnerRingPosx[i] += dx;
	}
}

void EllipseLayoutGenerator::writeMeshFileForSimulation(double zHeight, int zLayers){	
	vector <double> recordedNodesX, recordedNodesY, recordedNodesZ;
	vector <double> recordedPeripodialness;
	cout<<" Writing triangle file, size of points: "<<posx.size()<<" size of tissue shape: "<<tissueType.size()<<endl;
	ofstream MeshFile;
	MeshFile.open("./MeshFile.out",ofstream::trunc);
	cout<<" nNonDuplicateLinkers: "<<nNonDuplicateLinkers<<endl;
	int nNodes = posx.size()*(zLayers+1);
	if (addPeripodial){
		nNodes +=  posx.size()*(peripodialLayers+1) + SortedCircumferenceForLinkers.size()*nNonDuplicateLinkers;
	}
	MeshFile<<nNodes<<endl;
	cerr<<"posx.size(): "<<posx.size()<<" "<<posx.size()*(zLayers+1)<<endl;
	//Adding the basal layer:
	int n = posx.size();
	for (int i=0;i<n;++i){
		MeshFile<<posx[i]<<"	"<<posy[i]<<"	"<<0<<"    0	"<<tissueType[i]<<"	"<<atBorder[i]<<endl;	//x-coord   y-coord   z-coord   basal-node identifier(0) tissueType(columnar-0, peripodium-2) flag for border nodes
		recordedNodesX.push_back(posx[i]);
		recordedNodesY.push_back(posy[i]);
		recordedNodesZ.push_back(0);
	}
	//Adding columnar nodes, midline and apical surface:
	double  dzHeight = zHeight/zLayers;
	double currzHeight = dzHeight;
	for (int layers=0; layers<zLayers; ++layers){
		int nodePositionIdentifier = 2; //0=basal. 1=apical, 2=midline
		if (layers == zLayers-1){nodePositionIdentifier=1;}
		for (int i=0;i<n;++i){
			MeshFile<<posx[i]<<"	"<<posy[i]<<"	"<<currzHeight<<"    "<<nodePositionIdentifier<<"	"<<tissueType[i]<<"	"<<atBorder[i]<<endl;	
			recordedNodesX.push_back(posx[i]);
			recordedNodesY.push_back(posy[i]);
			recordedNodesZ.push_back(currzHeight);	
		}
		currzHeight += dzHeight;
	}
	cout<<" number of recorded nodes (columnar finished) : "<<recordedNodesX.size()<<" calculated node size: "<<nNodes<<endl;
	if(addPeripodial){
		currzHeight -= dzHeight;
		//Adding peripodial cap nodes:
		currzHeight  += lumenHeight;
		dzHeight = peripodialHeight/peripodialLayers;
		for (int layers=0; layers<peripodialLayers+1; ++layers){
			int nodePositionIdentifier = 2; //0=basal, 1=apical, 2=midline
			if (layers == 0){nodePositionIdentifier=1;} //apical for the bottom layer looking into lumen
			else if (layers == peripodialLayers){nodePositionIdentifier=0;} //basal for the top layer
			cout<<" adding layer : "<<layers<<" height: "<<	currzHeight<<" peripodialLayers: "<<peripodialLayers<<" nodePositionIdentifier: "<<nodePositionIdentifier<<endl;		
			for (int i=0;i<n;++i){
				MeshFile<<posx[i]<<"	"<<posy[i]<<"	"<<currzHeight<<"    "<<nodePositionIdentifier<<"	"<<1<<"	"<<0<<endl;	//1 0 : tissue type peripodial, at border false.
				recordedNodesX.push_back(posx[i]);
				recordedNodesY.push_back(posy[i]);
				recordedNodesZ.push_back(currzHeight);
			}
			currzHeight += dzHeight;
		}
		cout<<" number of recorded nodes (peri cap finished) : "<<recordedNodesX.size()<<" calculated node size: "<<nNodes<<endl;
		//Adding peripodial loop nodes:
		//trial for one:
		//calculate centre from borders:
		double c[2] = {0.0,0.0};
		for (int j =0 ; j<SortedCircumferenceForLinkers.size(); ++j){
			int id = SortedCircumferenceForLinkers[j];			
			c[0] += posx[id];	
			c[1] += posy[id];
		}
		c[0] /= SortedCircumferenceForLinkers.size();
		c[1] /= SortedCircumferenceForLinkers.size();
		if (symmetricY){
			c[1] = 0;
		};
		int nCirc = SortedCircumferenceForLinkers.size();
		for (int i =0 ; i<nCirc;++i){
			int id  = SortedCircumferenceForLinkers[i];
			//the normalised vector from centre to node position
			//double vec[2] = {posx[id] - c[0], posy[id] - c[1]};
			//double mag = pow(vec[0]*vec[0]+ vec[1]*vec[1],0.5);
			//vec[0] /= mag; vec[1] /= mag;
			double vec[2] = {1, 0};
			//calculate the angle between the orientation of the linker vector and the vector from centre to node:
			double c = vec[0]*LinkerVectorsOnCircumferenceX[i] +vec[1]*LinkerVectorsOnCircumferenceY[i];
			double s = sin(acos(c));
			//double rotMat[3][3] = {{c,-s,0},{s,c,0},{0,0,1}};
			//rotate the node list of the linkers with rotmat
			//translate the origin onto the base node xy
			//add the nodes to node list
			for (int j=0; j<LinkerPosx.size(); ++j){
				if (!duplicateLinkerNode[j]){
					//rotate:			
					double pos[3] = {c*LinkerPosx[j], s*LinkerPosx[j], LinkerPosz[j]};
					//translate:
					pos[0] += posx[id];
					pos[1] += posy[id];
					//add:
					MeshFile<<pos[0]<<"	"<<pos[1]<<"	"<<pos[2]<<"    "<<LinkerTissuePlacement[j]<<"	"<<2<<"	"<<0<<endl;	//2 0 : tissue type linker, at border false.
					//book keeping					
					recordedNodesX.push_back(pos[0]);
					recordedNodesY.push_back(pos[1]);
					recordedNodesZ.push_back(pos[2]);
				}
						
			}
			cout<<" number of recorded nodes (after i-"<<i<<") : "<<recordedNodesX.size()<<" calculated node size: "<<nNodes<<endl;
		}

	}
	
	//Adding elements:
	int nTri = triangles.size();
	int nLinkerTri = linkerTriangles.size()*(SortedCircumferenceForLinkers.size()-1);
	//int nLinkerTri = linkerTriangles.size();
	int nElements = nTri*(zLayers +peripodialLayers)+nLinkerTri;
	
	MeshFile<<nElements<<endl;
	int currOffset =0;
	//loop for columnar elements:
	for (int layers=0; layers<zLayers; ++layers){
		for (int i =0; i<nTri; ++i){
			double refPos[6][3];
			int nodes[6];
			for (int j=0;j<3;++j){			
				nodes[j]=triangles[i][j]+currOffset;
			}
			for (int j=3;j<6;++j){			
				nodes[j]=triangles[i][j-3]+currOffset+n;
			}	


			for (int j=0;j<6;++j){
				refPos[j][0] = recordedNodesX[nodes[j]];
				refPos[j][1] = recordedNodesY[nodes[j]];
				refPos[j][2] = recordedNodesZ[nodes[j]];
			}
			//writing Shapetype
			MeshFile<<"1	";
			//writing nodes of prism
			for (int j=0;j<6;++j){
				MeshFile<<nodes[j]<<"     ";
			}
			//writing positions for reference prism:
			for (int j=0;j<6;++j){
				MeshFile<<refPos[j][0]<<"   "<<refPos[j][1]<<"   "<<refPos[j][2]<<"    ";		
			}
			/*//basal nodes:
			double currzHeight = dzHeight*layers;
			for (int j=0;j<3;++j){
				MeshFile<<posx[triangles[i][j]]<<"   "<<posy[triangles[i][j]]<<"   "<<currzHeight<<"    ";		
			}
			currzHeight = dzHeight*(layers+1);
			//apical nodes:
			for (int j=0;j<3;++j){
				MeshFile<<posx[triangles[i][j]]<<"   "<<posy[triangles[i][j]]<<"   "<<currzHeight<<"   ";		
			}*/
			MeshFile<<endl;
			recordedPeripodialness.push_back(0.0); //This is a columnar element, the peripodialness is zero
		}
		currOffset += n;
	}
	//loop for peripodial elements:
	if(addPeripodial){
		//move one more layer up on nodes, I do not want to connect top layer of columnar to bottom of peripodial, I want to continue on peroipodial:
		currOffset += n;
		for (int layers=0; layers<peripodialLayers; ++layers){
			for (int i =0; i<nTri; ++i){
				double refPos[6][3];
				int nodes[6];
				for (int j=0;j<3;++j){			
					nodes[j]=triangles[i][j]+currOffset;
				}
				for (int j=3;j<6;++j){			
					nodes[j]=triangles[i][j-3]+currOffset+n;
				}
				for (int j=0;j<6;++j){
					refPos[j][0] = recordedNodesX[nodes[j]];
					refPos[j][1] = recordedNodesY[nodes[j]];
					refPos[j][2] = recordedNodesZ[nodes[j]];
				}
				//writing Shapetype
				MeshFile<<"1	";
				//writing nodes of prism
				for (int j=0;j<6;++j){
					MeshFile<<nodes[j]<<"     ";
				}
				//writing positions for reference prism:
				for (int j=0;j<6;++j){
					MeshFile<<refPos[j][0]<<"   "<<refPos[j][1]<<"   "<<refPos[j][2]<<"    ";		
				}
				//writing positions for reference prism:
				//basal nodes:
				/*double currzHeight = dzHeight*layers;
				for (int j=0;j<3;++j){
					MeshFile<<posx[triangles[i][j]]<<"   "<<posy[triangles[i][j]]<<"   "<<currzHeight<<"    ";		
				}
				currzHeight = dzHeight*(layers+1);
				//apical nodes:
				for (int j=0;j<3;++j){
					MeshFile<<posx[triangles[i][j]]<<"   "<<posy[triangles[i][j]]<<"   "<<currzHeight<<"   ";		
				}*/
				MeshFile<<endl;
				recordedPeripodialness.push_back(1.0); //This is a peripodial element, the peripodialness is zero
			}
			currOffset += n;
		}
		//loop for adding linker zone elements:
		//displaying the linker nodes:
		for (int i=0; i<LinkerPosz.size(); i++){
			cout<<"Linker: "<<i<<" pos: "<<LinkerPosx[i]<<" 0 "<<LinkerPosz[i]<<" duplicate? "<<duplicateLinkerNode[i]<<" position: "<<LinkerTissuePlacement[i]<<endl;
		}
		cout<<" number of recorded nodes: "<<recordedNodesX.size()<<" calculated node size: "<<nNodes<<endl;
		int nCirc = SortedCircumferenceForLinkers.size();
		for (int i =0 ; i<nCirc-1;++i){
			int circBaseId0 = SortedCircumferenceForLinkers[i];
			int circBaseId1 = SortedCircumferenceForLinkers[i+1];
			cout<<"circBaseIds: "<<circBaseId0<<" "<<circBaseId1<<endl;
			double refPos[6][3];
			int nodes[6];
			for (int k=0;k<linkerTriangles.size();++k){
				//cout<<" working on trianlge "<<k<<" of "<<linkerTriangles.size()<<endl;
				for (int j=0;j<3;++j){
					int idOnLinkerOrder = linkerTriangles[k][j];
					//cout<<" corner: "<<j<<" idOnLinkerOrder: "<<idOnLinkerOrder<<endl;
					if (duplicateLinkerNode[idOnLinkerOrder]){
						//cout<<" the node "<<j<<" is duplicate, pos: "<<LinkerPosx[idOnLinkerOrder]<<" 0 "<<LinkerPosz[idOnLinkerOrder]<<endl;
						//this is a duplicate node, I will go up on the circumference point, untill
						//I find the closest point in z (should be overlaping)
						//Then add that node to the node list.
						double dzCol = zHeight/zLayers;
						double dzPer = peripodialHeight/peripodialLayers;
						int currLayer = 0;					
						for (int k=0; k<zLayers+peripodialLayers+2; k++){
							//searching the columnar list:
							//cout<<" checking heights for k:"<<k<<" 	zLayers: "<<zLayers<<" 	dzCol: "<<dzCol<< " dzPer "<<dzPer<<endl;				
							double hToCompare = k*dzCol;							
							if (k<zLayers+1 && LinkerPosz[idOnLinkerOrder] > hToCompare-0.01 && LinkerPosz[idOnLinkerOrder] < hToCompare+0.01){
								//cout<<"equivalent layer is on columnar, k: "<<k<<endl;									
								currLayer = k;
								break;
							}
							else{
								double hToCompare = zLayers*dzCol + lumenHeight+ (k-zLayers-1)*dzPer;
								if( LinkerPosz[idOnLinkerOrder] > hToCompare-0.01 && LinkerPosz[idOnLinkerOrder] < hToCompare+0.01){
									//cout<<"equivalent layer is on peripodial, k: "<<k<<endl;									
									currLayer = k;
									break;										
								}						
							}		
						}
						nodes[j]= circBaseId0 + currLayer*n;
						nodes[j+3]= circBaseId1 + currLayer*n;
					}
					else{
						int currOffsetForThisCircumferenceNode = i*nNonDuplicateLinkers;
						int currOffsetForMainTissueNodes = n*(peripodialLayers+1+zLayers+1);
						int idAfterDuplicateClean = idOnLinkerOrder;
						for (int aa = 0; aa<idOnLinkerOrder; ++aa){
							if (duplicateLinkerNode[aa]){
								idAfterDuplicateClean--;							
							}	
						}
						//cout<<"currOffsetForThisCircumferenceNode: "<<currOffsetForThisCircumferenceNode<<" currOffsetForMainTissueNodes: "<<currOffsetForMainTissueNodes<<endl;					
						nodes[j]=idAfterDuplicateClean+currOffsetForThisCircumferenceNode+currOffsetForMainTissueNodes;
						nodes[j+3] = nodes[j]+nNonDuplicateLinkers;
					}	
					refPos[j][0] = recordedNodesX[nodes[j]];
					refPos[j][1] = recordedNodesY[nodes[j]];
					refPos[j][2] = recordedNodesZ[nodes[j]];
					refPos[j+3][0] = recordedNodesX[nodes[j+3]];
					refPos[j+3][1] = recordedNodesY[nodes[j+3]];
					refPos[j+3][2] = recordedNodesZ[nodes[j+3]];
				}
				//cout<<"Nodes: "<<nodes[0]<<"  "<<nodes[1]<<" "<<nodes[2]<<" "<<nodes[3]<<" " <<nodes[4]<<" " <<nodes[5]<<endl;	
				//writing Shapetype
				MeshFile<<"1	";
				//writing nodes of prism
				for (int j=0;j<6;++j){
					MeshFile<<nodes[j]<<"     ";
				}
				//writing positions for reference prism:
				for (int j=0;j<6;++j){
					MeshFile<<refPos[j][0]<<"   "<<refPos[j][1]<<"   "<<refPos[j][2]<<"    ";		
				}
				MeshFile<<endl;
				recordedPeripodialness.push_back(trianglePeripodialness[k]); //This is a linker element, the peripodialness is calculated arter linker triangulation
			}
		}
	}
	writeTissueWeights(MeshFile,recordedPeripodialness);
}

void EllipseLayoutGenerator::writeTissueWeights(ofstream& MeshFile, vector <double> &recordedPeripodialness){
	if (!addPeripodial){
		MeshFile<<0<<endl;
	}else{
		MeshFile<<1<<endl;
		int n = recordedPeripodialness.size();	
		for (int i=0; i<n; ++i){
			MeshFile<<recordedPeripodialness[i]<<endl;	//recording peripodialness		
		}	
	}
}

void EllipseLayoutGenerator::scaleInputOutline(vector <float>& x, vector <float>&y){
	double* boundingBox; //[minX , minY, maxX, maxY]
	boundingBox = new double[4];
	bringThicknessCentreToOrigin(x,y);
	smooth(100,x,y);
	calculateBoundingBox(boundingBox,x,y);
	scaleToTissueSize(boundingBox,x,y);
	reOrderToFixNodesAtZeroTips(x,y);
	arrangeSideLength(x,y);
	if (symmetricY){
		clearSymmetric(x,y);
	}
	bluntTip(x,y);
	calculateBoundingBox(boundingBox,x,y);
	//addLayersToOutline(2,boundingBox,x,y);
	writeInPosVectors(x,y);
	int n=x.size();
	for (int i=0; i<n; ++i){
		cout<<" node: "<<i<<" "<<x[i]<<" "<<y[i]<<endl;	
	}
}

void EllipseLayoutGenerator::bluntTip(vector <float>&  x, vector <float>& y){
	//I will cut the tip until the thinnest point becomes sligthly larger than the side length
	//this is necessary for stability, and avoiding tiny triangles.
	//find the y tip index
	//I know first node is at y=0. I need the second one:
	int n = x.size();
	int yTipIndex = 0;	
	for (int i=1; i<n; ++i){
		if (y[i] ==0){
			yTipIndex = i;
			break;
		}
	}
	cout<<"yTipIndex: "<<yTipIndex<<endl;
	for (int i=0; i<n; ++i){
		cout<<" currnodelist: "<<i<<" "<<x[i]<<" "<<y[i]<<endl;	
	}
	//now go back till you reach the height of sideLength
	int firstBreak = 0;	
	for (int i=yTipIndex; i<n; --i){
		if (y[i] > sideLen*0.9){
			firstBreak = i-1;
			break;		
		}	
	}
	//now go forward till you reach the height of side Length
	int secondBreak = n;
	for (int i=yTipIndex; i<n; i++){
		if (y[i] > sideLen*0.9){
			secondBreak = i;
			break;		
		}	
	}
	//take backup
	vector <float> tmpX, tmpY;
	for (int i=0; i<n; ++i){
		tmpX.push_back(x[i]);
		tmpY.push_back(y[i]);
	}
	x.clear();
	y.clear();
	//now remove those items;
	for (int i=0; i<firstBreak; ++i){
		x.push_back(tmpX[i]);
		y.push_back(tmpY[i]);
	}
	x.push_back(tmpX[yTipIndex]);
	y.push_back(tmpY[yTipIndex]);
	for (int i=secondBreak; i<n; ++i){
		x.push_back(tmpX[i]);
		y.push_back(tmpY[i]);
	}		
}

void EllipseLayoutGenerator::writeInPosVectors(vector <float>&  x, vector <float>& y){
	int n=x.size();
	for (int i=0; i<n; ++i){
		posx.push_back(x[i]);
		posy.push_back(y[i]);
		tissueType.push_back(0); //0: columnar layer, 1: peripodium
	}	
}

void EllipseLayoutGenerator::addLayersToOutline(int nLayer, double* boundingBox, vector <float>&  x, vector <float>& y){
	//scale down to have about 0.75 sideLen gap on each side:
	float dx = boundingBox[2] - boundingBox[0];
	float dy = boundingBox[3] - boundingBox[1];
	if(symmetricY){	
		dy *=2.0;	
	}
	float scaleX = (dx - 2.0*0.75*sideLen)/dx;
	float scaleY = (dy - 2.0*0.75*sideLen)/dy;
	//cout<<"scaleX: "<<scaleX<<" scaleY: "<<scaleY<<endl;
	int n = x.size();
	for (int j=0; j<nLayer; ++j){
		for (int i=n*j; i<n*(j+1); ++i){
			x.push_back(scaleX*x[i]);
			y.push_back(scaleY*y[i]);		
		}	
	}
}

void EllipseLayoutGenerator::clearSymmetric(vector <float>&  x, vector <float>& y){		
	if (symmetricY){
		int n=x.size();
		vector<float> tmpX, tmpY;;
		for (int i=0; i<n; ++i){
			if (y[i]>=0){
				tmpX.push_back(x[i]);
				tmpY.push_back(y[i]);		
			}		
		}
		x.clear();
		y.clear();
		n = tmpX.size();
		for (int i=0; i<n; ++i){
			x.push_back(tmpX[i]);		
			y.push_back(tmpY[i]);		
		}
	}
}

void EllipseLayoutGenerator::calculateBoundingBox(double* boundingBox, vector <float>&  x, vector <float>& y){
	double minX = 1000, maxX = -1000, minY = 1000, maxY = -1000;
	int n = x.size();
	for (int i=0; i<n; ++i){
		if (x[i]<minX){minX = x[i];};	
		if (x[i]>maxX){maxX = x[i];};	
		if (y[i]<minY){minY = y[i];};	
		if (y[i]>maxY){maxY = y[i];};	
	}
	boundingBox[0] = minX;
	boundingBox[1] = minY;
	boundingBox[2] = maxX;
	boundingBox[3] = maxY;
}

void EllipseLayoutGenerator::bringThicknessCentreToOrigin(vector <float>&  x, vector <float>& y){
	int n = x.size();
	//bring the thickest point in y to y=0:
	//find the x min:
	float minX=1000;	
	int minXId;
	for (int i=0; i<n; ++i){
		if (x[i]<minX){minX = x[i];minXId = i;};	
	}	
	float dy = -y[minXId];
	for (int i=0; i<n; ++i){
		y[i] +=	dy;
			
	}
	//bring x =0 to align with y peak:
	//find the y peak:
	float maxY=-1000;	
	int maxYId;
	for (int i=0; i<n; ++i){
		if (y[i]>maxY){maxY = y[i];maxYId = i;};	
	}
	float dx = -x[maxYId];
	for (int i=0; i<n; ++i){
		x[i] +=	dx;
		
	}
}

void EllipseLayoutGenerator::smooth(int pointAvr, vector <float>&  x, vector <float>& y){
	int n = x.size();
	float tmpCoord[(const int) n][2];
	float range = (pointAvr-1.0)/2.0;
	/*for (int i=range; i<n-range; ++i){
		float currSum[2] = {0.0,0.0};
		for (int j=-range; j<=range; ++j){
			currSum[0] += x[i+j];
			currSum[1] += y[i+j];
		}
		currSum[0] /=pointAvr;
		currSum[1] /=pointAvr;
		tmpCoord[i][0] = currSum[0];
		tmpCoord[i][1] = currSum[1];
	}
	for (int i=range; i<n-range; ++i){
		x[i] = tmpCoord[i][0];
		y[i] = tmpCoord[i][1];	
	}*/
	for (int i=0; i<n; ++i){
		float currSum[2] = {0.0,0.0};
		for (int j=-range; j<=range; ++j){
			int index = i+j;
			if (index <0){
				index += n;
			}
			if (index >=n){
				index -=n;			
			}
			currSum[0] += x[index];
			currSum[1] += y[index];
		}
		currSum[0] /=pointAvr;
		currSum[1] /=pointAvr;
		tmpCoord[i][0] = currSum[0];
		tmpCoord[i][1] = currSum[1];
		//x[i] = currSum[0];
		//y[i] = currSum[1];	
	}
	for (int i=range; i<n-range; ++i){
		x[i] = tmpCoord[i][0];
		y[i] = tmpCoord[i][1];	
	}	
}

void EllipseLayoutGenerator::scaleToTissueSize(double* boundingBox, vector <float>&  x, vector <float>& y){
	//bounding box format: [minX , minY, maxX, maxY]	
	float sizeOutline[2] = {boundingBox[2] - boundingBox[0],  boundingBox[3] - boundingBox[1]};
	float sizeTissue[2] = {r1[0]+r1[1],  r2[0]+r2[1]};
	float xMultiplier = sizeTissue[0] / sizeOutline[0];
	float yMultiplier = sizeTissue[1] / sizeOutline[1];
	int n = x.size();	
	for (int i=0; i<n; ++i){
		x[i] *= xMultiplier;
		y[i] *= yMultiplier;
	}
}

void EllipseLayoutGenerator::reOrderToFixNodesAtZeroTips(vector <float>&  x, vector <float>& y){
	//find first point to pass from (-) y to (+) y
	int index = 0;
	bool foundTransition = false;
	bool foundPositive = false;
	bool foundNegative = false;
	int posIndex = 0;
	int negIndex = 0;
	int n=x.size();
	bool flippedDirection = false;
	int counter =0;
	while (!foundTransition){
		if(y[index] > 0){
			foundPositive = true;
			posIndex = index;
			index --;
			if (foundPositive && foundNegative){
				foundTransition = true;
			}
		}
		else{
			foundNegative = true;
			negIndex = index;
			index ++;
			if (foundPositive && foundNegative){
				foundTransition = true;
			}
		}
		if (index >= n){
			index -=n;
		}
		if (index <0){
			index +=n;
		}
	}
	//I have found the point that is the transition in y at the positive x tip.
	//Now I need to intrapolate and put a point to y=0;
	float dx = x[posIndex] - x[negIndex];
	float dy = y[posIndex] - y[negIndex];
	float dy0 = -1.0*y[negIndex];
	float multiplier = dy0/dy;
	float newX = multiplier * dx + x[negIndex];
	vector<float> tmpX, tmpY;	
	tmpX.push_back(newX);
	tmpY.push_back(0.0);	
	n=x.size();
	if (posIndex>negIndex){	
		for (int i=posIndex; i<n; ++i){	
			tmpX.push_back(x[i]);
			tmpY.push_back(y[i]);	
		}
		for (int i=0; i<negIndex+1; ++i){	
			tmpX.push_back(x[i]);
			tmpY.push_back(y[i]);	
		}
	}
	else{
		for (int i=posIndex; i<n; ++i){	
			tmpX.push_back(x[i]);
			tmpY.push_back(y[i]);	
		}			
	}
	//Now I need to repear the procedure again, to add another point to the x tip:
	// I know first transition is at point 0 now, so I can move on and start from 1
	//find first negative:
	index = 1;
	foundTransition = false;
	foundPositive = false;
	foundNegative = false;
	posIndex = 0;
	negIndex = 0;
	n=tmpX.size();
	while (!foundNegative){	
		if(tmpY[index] < 0){
			foundNegative = true;
			negIndex = index;
			index --;
			if (index <0){
				index +=n;		
			}
			posIndex = index;
		}
		index++;
	}
	//interpolate:
	dx = tmpX[posIndex] - tmpX[negIndex];
	dy = tmpY[posIndex] - tmpY[negIndex];
	dy0 = -1.0*tmpY[negIndex];
	multiplier = dy0/dy;
	newX = multiplier * dx + tmpX[negIndex];
	x.clear();
	y.clear();	
	n=tmpX.size();
	for (int i =0; i<=posIndex; ++i){
		x.push_back(tmpX[i]);	
		y.push_back(tmpY[i]);			
	}
	x.push_back(newX);
	y.push_back(0);
	for (int i =posIndex+1; i<n; ++i){
		x.push_back(tmpX[i]);	
		y.push_back(tmpY[i]);			
	}
}

void EllipseLayoutGenerator::arrangeSideLength(vector <float>&  x, vector <float>& y){
	//I am scaling the nodes to be double the side length, the triangulation can fill in between	
	float scaledSideLen = sideLen*2.0;	
	vector<int> tmpX, tmpY;	
	int n=x.size();
	int yTipIndex =0;
	for (int i=0; i<n; ++i){
		//find the tip in Y (that is, the transition from (+)ve to (-)ve y).
		//I want that node in position:		
		if (y[i] < 0 && i>0){
			yTipIndex = i-1;
			break;	
		}			
	}
	for (int i=0; i<n; ++i){
		tmpX.push_back(x[i]);			
		tmpY.push_back(y[i]);			
	}
	x.clear();
	y.clear();
	n=tmpX.size();
	x.push_back(tmpX[0]);
	y.push_back(tmpY[0]);
	float lastX=x[0];
	float lastY=y[0];
	float s2 = scaledSideLen*scaledSideLen;
	for (int i=1; i<yTipIndex; ++i){
		float dx= lastX - tmpX[i]; 
		float dy= lastY - tmpY[i]; 
		float d2 = dx*dx + dy*dy;
		if (d2 > s2){
			x.push_back(tmpX[i]);
			y.push_back(tmpY[i]);
			lastX = tmpX[i];
			lastY = tmpY[i];	
		}
	}
	//remove all node up until the distance to the yTip is above scaledSideLen:
	float d = 0.1*scaledSideLen;
	while (d<scaledSideLen){
		n=x.size();
		float dx= x[n-1] - tmpX[yTipIndex]; 
		float dy= y[n-1] - tmpY[yTipIndex]; 
		float d2 = dx*dx + dy*dy;
		d = pow (d2,0.5);
		if (d<scaledSideLen){
			x.pop_back();
			y.pop_back();		
		}
	}
	//add tip:
	x.push_back(tmpX[yTipIndex]);
	y.push_back(tmpY[yTipIndex]);
	lastX = tmpX[yTipIndex];
	lastY = tmpY[yTipIndex];
	//now calculate last distance again, if it is more than 1.5 scaledSideLen, then add one in between:
	n=x.size();
	float dx= x[n-2] - x[n-1]; 
	float dy= y[n-2] - y[n-1]; 
	float d2 = dx*dx + dy*dy;
	d = pow (d2,0.5);
	if (d>1.5*scaledSideLen){
		float midX = 0.5* (x[n-2] + x[n-1]);
		float midY = 0.5* (y[n-2] + y[n-1]);
		x.pop_back();	//remove tip point
		y.pop_back();		
		x.push_back(midX);	//add new point
		y.push_back(midY);
		x.push_back(tmpX[yTipIndex]);	//add tip point
		y.push_back(tmpY[yTipIndex]);
	}
	//I also need the tip index in x indice recording, not tmpX (which is the full set).
	int yTipIndexInXCoord = x.size()-1;
	//now continue loop:
	n=tmpX.size();
	for (int i=yTipIndex+1; i<n; ++i){
		float dx= lastX - tmpX[i]; 
		float dy= lastY - tmpY[i]; 
		float d2 = dx*dx + dy*dy;
		if (d2 > s2){
			x.push_back(tmpX[i]);
			y.push_back(tmpY[i]);
			lastX = tmpX[i];
			lastY = tmpY[i];	
		}
	}
	//check loop connection between the last node and the first node:
	//remove all nodes up until the distance is above scaledSideLen
	d = 0.1*scaledSideLen;
	while (d<scaledSideLen){
		n=x.size();
		float dx= x[0] - x[n-1]; 
		float dy= y[0] - y[n-1]; 
		float d2 = dx*dx + dy*dy;
		d = pow (d2,0.5);
		if (d<scaledSideLen){
			x.pop_back();
			y.pop_back();		
		}
	}
	//now calculate last distance again, if it is more than 1.5 scaledSideLen, then add one in between:
	n=x.size();
	dx= x[0] - x[n-1]; 
	dy= y[0] - y[n-1]; 
	d2 = dx*dx + dy*dy;
	d = pow (d2,0.5);
	if (d>1.5*scaledSideLen){
		n=x.size();
		float midX = 0.5* (x[0] + x[n-1]);
		float midY = 0.5* (y[0] + y[n-1]);
		x.push_back(midX);
		y.push_back(midY);
	}
	//Now I will smooth the data once more, then fix initial poonts again, anf it will be ready:
	float fixedPos[2][2] = {{x[0], y[0]},{x[yTipIndexInXCoord],y[yTipIndexInXCoord]}};
	smooth(3,x,y);
	if (x[0] < fixedPos[0][0]){
		x[0] = fixedPos[0][0];	//due to smoothing, the node above first node may have moved further down x, I do not want to bring this closer to centre
	}
	y[0] = fixedPos[0][1];
	if (x[yTipIndexInXCoord] > fixedPos[1][0]){
		x[yTipIndexInXCoord] = fixedPos[1][0];	//same in negative direction
	}	
	y[yTipIndexInXCoord] = fixedPos[1][1];
}


bool readinputs (int argc, char **argv, double* parameters, ifstream& inputOutline){
	if (argc<2){
		cerr<<"Please give shape type: 1  = ellipse, 2 = rectangle"<<endl;
		return 0;	
	}
	else {
		const char* inpstring = argv[1];
		parameters[0] = atof(inpstring);
	}
	int offset = 0;
	if (parameters[0] == 1){
		if(argc<6){
			cerr<<"Please give lenght1, length2, width1 and width2 (in this order)"<<endl;
			return 0;
		}
		else{
			const char* inpstring1 = argv[2];
			parameters[1] = atof(inpstring1);
			const char* inpstring2 = argv[3];
			parameters[2] = atof(inpstring2);		
			const char* inpstring3 = argv[4];
			parameters[3] = atof(inpstring3);
			const char* inpstring4 = argv[5];
			parameters[4] = atof(inpstring4);			
		}
		offset = 6;
	}
	if(argc<offset+1){
		cerr<<"Please give the ABHeight of the tissue"<<endl;
		return 0;	
	}
	else{
		const char* inpstring = argv[offset];
		parameters[5] = atof(inpstring);		
	}
	if(argc<offset+2){
		cerr<<"Please give the side length of prisms"<<endl;
		return 0;
	}
	else{
		const char* inpstring = argv[offset+1];
		parameters[6] = atof(inpstring);		
	}
	if(argc<offset+3){
		cerr<<"Please give the number of layers in z"<<endl;
		return 0;	
	}
	else{
		const char* inpstring = argv[offset+2];
		parameters[7] = atof(inpstring);		
	}
	if(argc<offset+4){
		cerr<<"Please give if the mesh will be symmetric in y(bool)"<<endl;
		return 0;	
	}
	else{
		const char* inpstring = argv[offset+3];
		parameters[8] = atof(inpstring);		
	}
	if(argc<offset+5){
		cerr<<"Please give the input file path"<<endl;
		return 0;	
	}
	else{
		const char* inpstring = argv[offset+4];
		inputOutline.open(inpstring,ifstream::in);
		if (!inputOutline.good() || !inputOutline.is_open()){
			cerr<<"input outline cannot be opened: "<<inpstring<<endl;
			return 0;
		}
	}
	return 1;
}

void printOutSystemSetup(double * parameters){
	string ShapeType = "NULL";
	if (parameters[0] == 1){ShapeType = "ellipse";}
	if (parameters[0] == 2){ShapeType = "rectangle";}
	cerr<<"Generating mesh for "<<ShapeType<<" :"<<endl;
	cerr<<"	DVRadia    : "<<parameters[1]<<" "<<parameters[2]<<endl;
	cerr<<"	APRadia    : "<<parameters[3]<<" "<<parameters[4]<<endl;
	cerr<<"	ABHeight   : "<<parameters[5]<<endl;
	cerr<<"	sideLength : "<<parameters[6]<<endl;
	cerr<<"	ABLayers   : "<<parameters[7]<<endl;
	cerr<<"	symmetricY : "<<parameters[8]<<endl;
}

void readInOutline(vector <float>& x, vector <float>&y, ifstream& inputOutline){
	int n, id,currX, currY;
	inputOutline>>n;	
	for(int i=0;i<n;++i){
		inputOutline>>id;	
		inputOutline>>currX;
		inputOutline>>currY;
		x.push_back((float) currX);
		y.push_back((float) currY);
	}
}

int main(int argc, char **argv)
{	
	double * parameters;
	parameters = new double[9];
	ifstream inputOutline;
	bool success = readinputs (argc, argv, parameters,inputOutline);
	if (!success) {
		cerr<<"Error in parameter inputs, quitting"<<endl;	
		return 1;
	};
	int 	GlobalShape;
	double 	DVRadius[2];
	double 	APRadius[2];
	double 	ABHeight;
	double 	sideLength;
	int    	ABLayers;
	bool 	symmetricY = false;
	
	if (success) {
		printOutSystemSetup(parameters);	
		GlobalShape = parameters[0];
		DVRadius[0] = parameters[1];
		DVRadius[1] = parameters[2];
		APRadius[0] = parameters[3];
		APRadius[1] = parameters[4];
		ABHeight = parameters[5];
		sideLength = parameters[6];
		ABLayers = parameters[7];
		if (parameters[8] == 1 ){ 		
			symmetricY = true;
		}
	}
	EllipseLayoutGenerator Lay01(DVRadius[0],DVRadius[1], APRadius[0], APRadius[1], sideLength, symmetricY, 0.9);
	bool addPeripodial = true;
	double peripodialHeightFrac = 0.45;
	double lumenHeightFrac = 0.25;
	double peripodialSideCurveFrac = 5.52 /ABHeight ; //5.52 is the side thickness Maria Measured for 48 hr discs
	//double peripodialSideCurveFrac = 7.1 /ABHeight ; //7.1 is the side thickness Maria Measured for 72 hr discs
	if(addPeripodial){
		Lay01.addPeripodial = true;
		Lay01.calculatePeripodialMembraneParameters(ABHeight, ABLayers, peripodialHeightFrac, lumenHeightFrac,peripodialSideCurveFrac);
	}
	vector <float> x, y;	
	readInOutline(x,y,inputOutline);
	Lay01.scaleInputOutline(x,y);
	
	//cout<<"before Tesselate2D, n =  "<<	Lay01.posx.size()<<endl;
	Lay01.Tesselate2D();
	Lay01.readInTesselation2D();
	//now I have the triangulated mesh. If I have peripodial, I need to get the circumference.
	//Then I will calculate the vectors pointing out from each of the circumference nodes.
	//Then I will use those values to add the elements in mesh generation.
	if(addPeripodial){
		Lay01.calculatePeripodialAttachVectors(symmetricY);
	}
	Lay01.calculateAverageSideLength();
	Lay01.writeMeshFileForSimulation(ABHeight,ABLayers);	
	//output the points for plotting:
	int n=Lay01.posx.size();	
	//cerr<<"r1: "<<Lay01.r1[0]<<" "<<Lay01.r1[1]<<" r2: "<<Lay01.r2[0]<<" "<<Lay01.r2[1]<<endl;
	Lay01.calculateAverageSideLength();
}


