#ifndef Prism_H
#define Prism_H


#include "ShapeBase.h"
#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

class Prism : public ShapeBase{

protected:
    void setTissueCoordsRotationsBuffers();
    void getCurrRelaxedShape(gsl_matrix * CurrRelaxedShape);
    void setShapeFunctionDerivatives(gsl_matrix * ShapeFuncDer,double eta, double zeta, double nu);
    void setShapeFunctionDerivativeStack(gsl_matrix* ShapeFuncDer, gsl_matrix* ShapeFuncDerStack);
	void setCoeffMat();
	void calculateDVector();
	void calculateD81Tensor();
	void calculateCurrk(boost::numeric::ublas::matrix<double> &currk, boost::numeric::ublas::matrix<double> &currB, boost::numeric::ublas::matrix<double>& currBE, boost::numeric::ublas::matrix<double> &currBo, double eta, double zeta, double nu);
    void calculateCurrNodalForces(gsl_matrix *gslcurrge, gsl_matrix *gslcurrgv, gsl_matrix *gslcurrF, gsl_matrix* displacementPerDt, int pointNo);
    void calculateCurrTriPointFForRotation(gsl_matrix *currF,int pointNo);
    void calculateNormalToBottom();
	void calculateReferenceNormalToBottom();
	void calculateNormalToTop();
	void calculateReferenceNormalToTop();
	void getCurrentAlignmentSides(double*, double*);
	void getCurrentAlignmentFaces(double* RefSide, double* ShapeSide, double* RefFace, double* ShapeFace);
	void updateAlignmentTurn();
	//void updateReferenceShapeBaseFromBuffer();
	//void resetBuffersAfterGrowth();
	//void calculateZVecForTissueCoordAlignment(double* u);
	//void calculateXVecForTissueCoordAlignment(double* u);
	void calculateReferenceVolume();
	void calculatePlaneNormals(double** normals);
	void assignNodalVector(double* vec, int id0, int id1);
	bool checkNodePlaneConsistency(double** normals);
	void setInitialEdgeLenghts();
	void checkEdgeLenghtsForBinding(vector<int>& masterIds, vector<int>& slaveIds);
	double getApicalSideLengthAverage();
	double getBasalSideLengthAverage();
	void assignExposedSurfaceAreaIndices(vector <Node*>& Nodes);
	double initialApilcalEdgeLengthsSq[3];
	double initialBasalEdgeLengthsSq[3];
public:
	Prism(int* NodeIds,vector<Node*>& Nodes, int CurrId, bool thereIsPlasticDeformation);
	~Prism();
	void  setElasticProperties(double EApical, double EBasal, double EMid,  double EECM, double v);
	void  fillLateralNeighbours(vector<Node*>& Nodes, vector<int>& lateralNeigbours );
	void  calculateBasalNormal(double * normal);
	void  calculateApicalNormalCurrentShape();
	void  AlignReferenceBaseNormalToZ();
    void  calculateElementShapeFunctionDerivatives();

	//void  calculateForces(int RKid, double** SystemForces);
	void  checkHealth();
	void getApicalTriangles(vector <int> &ApicalTriangles);
	int getCorrecpondingApical(int currNodeId);
	bool IsThisNodeMyBasal(int currNodeId);
	bool IsThisNodeMyApical(int currNodeId);
	double getElementHeight();
	void AddPackingToSurface(int tissueplacementOfPackingNode, double Fx, double Fy,double Fz,  double **PackingForces,vector<Node*> &Nodes, bool& allCornersFixedX, bool& allCornersFixedY, bool& allCornersFixedZ);
	void getRelevantNodesForPacking(int TissuePlacementOfPackingNode, int& id1, int& id2, int& id3);
	bool IsPointCloseEnoughForPacking(double* Pos,  float threshold, int TissuePlacementOfPackingNode);
	void calculateNormalForPacking(int tissuePlacementOfNormal);
	void calculateApicalArea();
	void calculateBasalArea();
	void calculateMyosinForcesAreaBased(double forcePerMyoMolecule);
	void calculateMyosinForcesTotalSizeBased(double forcePerMyoMolecule);
	void distributeMyosinForcesAreaBased(bool isIsotropic, bool apical, double forcePerMyoMolecule);
	void distributeMyosinForcesTotalSizeBased(bool isIsotropic, bool apical, double forcePerMyoMolecule);
	//void fillLateralNeighbours();
	void updateElasticProperties();
	void setBasalNeigElementId(vector<ShapeBase*>& elementsList);
	void constructElementStackList(const int discretisationLayers, vector<ShapeBase*>& elementsList);
	void getApicalNodePos(double* posCorner);
	void getBasalNodePos(double* posCorner);
	void getApicalCentre(double* centre);
	void getBasalCentre(double* centre);
	void getLumenFacingNodeIds(int* nodeIds,int& numberOfTriangles);
	void getReferenceApicalCentre(double* centre);
	void getReferenceBasalCentre(double* centre);
	double* getApicalMinViscosity(vector<Node*> Nodes);
	double* getBasalMinViscosity(vector<Node*> Nodes);
	bool IspointInsideTriangle(int tissueplacement,double x, double y,double z);
	void checkRotationConsistency3D();
	void copyElementInformationAfterRefinement(ShapeBase* baseElement, int layers, bool thereIsPlasticDeformation);
	bool areNodesDirectlyConnected(int node0, int node1);
};

#endif
