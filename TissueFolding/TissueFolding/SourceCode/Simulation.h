#ifndef Simulation_H
#define Simulation_H

#include <stdio.h>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>

#include "ModelInputObject.h"
#include "ShapeBase.h"
#include "Node.h"
#include "GrowthFunctionBase.h"
#include "GrowthFunctionTypes.h"
#include "MyosinFunction.h"
#include "NewtonRaphsonSolver.h"
#include "CellMigration.h"
#include "MuscleFibre.h"


#include <omp.h>
//test for rici pull

class ModelInputObject;
using namespace std;

class Simulation {
private:
	int currElementId;
	ModelInputObject* ModInp;
	ofstream saveFileMesh;
	ofstream saveFileTensionCompression;
    ofstream saveFileGrowth;
    ofstream saveFileGrowthRate;
	ofstream saveFileForces;
	ofstream saveFileProteins;
	ofstream saveFilePhysicalProp;
	ofstream saveFileSpecificType;
	ofstream saveFileGrowthRedistribution;
	ofstream saveFileNodeBinding;
	ofstream saveFileCollapseAndAdhesion;
	ofstream saveFilePacking;
	ofstream saveFileSimulationSummary;
	ifstream saveFileToDisplayMesh;
	ifstream inputMeshFile;
	ifstream saveFileToDisplayTenComp;
    ifstream saveFileToDisplayGrowth;
    ifstream saveFileToDisplayGrowthRate;
	ifstream saveFileToDisplayForce;
	ifstream saveFileToDisplayProteins;
	ifstream saveFileToDisplayPhysicalProp;
	ifstream saveFileToDisplaySpecificType;
	ifstream saveFileToDisplayPacking;
	ifstream saveFileToDisplayVel;
	ifstream saveFileToDisplaySimSum;
	ifstream saveFileToDisplaySpecificNodeTypes;
	ifstream saveFileToDisplayGrowthRedistribution;
	ifstream saveFileToDisplayNodeBinding;
	ifstream saveFileToDisplayCollapseAndAdhesion;

	bool TensionCompressionSaved;
    bool GrowthSaved;
    bool GrowthRateSaved;
	bool ForcesSaved;
	bool ProteinsSaved;
	bool physicalPropertiesSaved;
	bool PackingSaved;
	bool growthRedistributionSaved;
	bool nodeBindingSaved;
	bool specificElementTypesRecorded;
	bool collapseAndAdhesionSaved;
	int	 nCircumferencialNodes;
	int dorsalTipIndex,ventralTipIndex,anteriorTipIndex,posteriorTipIndex;
	double StretchDistanceStep;
	bool recordForcesOnFixedNodes;
	double boundingBoxSize[3];
	//double columnarBoundingBoxSize[3];
	//double peripodialBoundingBoxSize[3];
	bool ContinueFromSave;
    int growthRotationUpdateFrequency;
    //vector <Node*> symmetricYBoundaryNodes;
    //vector <Node*> symmetricXBoundaryNodes;
    vector <int> AblatedNodes;
    bool checkedForCollapsedNodesOnFoldingOnce;
    bool boundLateralElements;

	bool readModeOfSim(int& i, int argc, char **argv);
	bool readParameters(int& i, int argc, char **argv);
	bool readOutputDirectory(int& i, int argc, char **argv);
	bool readSaveDirectoryToDisplay(int& i, int argc, char **argv);
	bool openFilesToDisplay();
	bool readSystemSummaryFromSave();
	bool readSpecificNodeTypesFromSave();
	void initiateNodesFromSave();
	bool readNodeDataToContinueFromSave();
	void initiateNodesFromMeshInput();
	void initiateElementsFromSave();
	bool readElementDataToContinueFromSave();
	void initiateElementsFromMeshInput();
	void initiatePrismFromSave();
	bool readShapeData(int i);
	void initiateTriangleFromSave(double height);
	void initiatePrismFromMeshInput();
	void initiateTetrahedraFromMeshInput();
	void initiateTriangleFromMeshInput();
	void initiatePrismFromSaveForUpdate(int k);
	void initiateTriangleFromSaveForUpdate(int k,double height);
	void removeElementFromEndOfList();
	void updateNodeNumberFromSave();
	void updateNodePositionsFromSave();
	void updateElementStatesFromSave();
	void updateForcesFromSave();
	void updateTensionCompressionFromSave();
    void updateGrowthFromSave();
    void updateGrowthRateFromSave();
    void updateProteinsFromSave();
    void updatePhysicalPropFromSave();
    void updatePackingFromSave();
    void updateNodeBindingFromSave();
    void  updateCollapseAndAdhesionFromSave();
    void updateGrowthRedistributionFromSave();
	void readTensionCompressionToContinueFromSave();
    void readGrowthToContinueFromSave();
    void readGrowthRateToContinueFromSave();
    void readProteinsToContinueFromSave();
    void readPhysicalPropToContinueFromSave();
    void readGrowthRedistributionToContinueFromSave();
    void readNodeBindingToContinueFromSave();
    void readCollapseAndAdhesionToContinueFromSave();
	bool readFinalSimulationStep();
	void reInitiateSystemForces(int oldSize);
	bool checkInputConsistency();
	void setDefaultParameters();
	bool openFiles();
	bool reOpenOutputFile();
	void initiateSystemForces();
	bool initiateMesh(int MeshType,float zHeight);
	bool initiateMesh(int MeshType, int Row, int Column, float SideLength, float zHeight);
	//bool initiateMesh(int MeshType, string inputtype, float SideLength, float zHeight );
	bool initiateMesh(int MeshType);
	void readInTissueWeights();
	bool checkIfTissueWeightsRecorded();
	bool generateColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList);
	void sortColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList);
	void clearCircumferenceDataFromSymmetricityLine();
	//void removeSymmetryBorderFromColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList);
	void getAverageSideLength(double& periAverageSideLength, double& colAverageSideLength);
	int getElementArrayIndexFromId(int currId);
	int getNodeArrayIndexFromId(int currId);
	bool isColumnarLayer3D();
	bool checkIfThereIsPeripodialMembrane();
	void setLinkerCircumference();
	bool calculateTissueHeight();
	void assignInitialZPositions();
	bool addStraightPeripodialMembraneToTissue();
	bool addCurvedPeripodialMembraneToTissue();
	void calculateDiscretisationLayers(double &hColumnar, int& LumenHeightDiscretisationLayers, double &hLumen, double &peripodialHeight, int& peripodialHeightDiscretisationLayers, double& hPeripodial);
	void fillColumnarBasedNodeList(vector< vector<int> > &ColumnarBasedNodeArray, vector <int> &ColumnarCircumferencialNodeList);
	void calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, double* pos, double sideThickness);
	void calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, int nodeId2, double* pos, double sideThickness);
	void addNodesForPeripodialOnOuterCircumference (vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray, double hColumnar, int LumenHeightDiscretisationLayers, double hLumen, int peripodialHeightDiscretisationLayers, double hPeripodial);
	void addNodesForPeripodialOnColumnarCircumference (vector< vector<int> > &ColumnarBasedNodeArray, int LumenHeightDiscretisationLayers, double hLumen, int peripodialHeightDiscretisationLayers, double hPeripodial);
	void addLateralPeripodialElements(int LumenHeightDiscretisationLayers, int peripodialHeightDiscretisationLayers, vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray);
	void addNodesForPeripodialOnCap(vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &PeripodialCapNodeArray, int TissueHeightDiscretisationLayers, int LumenHeightDiscretisationLayers, int peripodialHeightDiscretisationLayers, double hPeripodial);
	void constructTriangleCornerListOfApicalSurface( vector< vector<int> > &TriangleList);
	void addCapPeripodialElements( vector< vector<int> > &TriangleList, vector< vector<int> > &PeripodialCapNodeArray, int peripodialHeightDiscretisationLayers);
	void correctCircumferentialNodeAssignment(vector< vector<int> > OuterNodeArray);


	void initiateSinglePrismNodes(float zHeight);
	void initiateSinglePrismElement();
	void initiateNodesByRowAndColumn(int Row, int Column,  float SideLength, float zHeight);
	void fixApicalBasalNodes(vector<int> &NodesToFix);
	void checkForNodeFixing();
	void checkForNodeBinding();
	void clearUpRigidFixedNodesFromSlaves();
	bool bindEllipseAxes();
	bool bindCircumferenceXY();
	bool areNodesToCollapseOnLateralECM(int slaveNodeId, int masterNodeId);
	bool checkEdgeLenghtsForBindingPotentiallyUnstableElements();
	bool bindPeripodialToColumnar();
	bool areNodesOnNeighbouingElements(int masterNoeId, int slaveNodeId);
	void manualAdhesion(int masterNodeId,int slaveNodeId);
	double distanceSqBetweenNodes(int id0, int id1);
	bool isAdhesionAllowed(int masterNodeId, int slaveNodeId);
	bool checkForElementFlippingUponNodeCollapse(vector<int> &newCollapseList, double* avrPos);
	bool adhereNodes();
	bool updatePositionsOfNodesCollapsingInStages();
	void induceClones();
	void initiateElementsByRowAndColumn(int Row, int Column);
	void assignPhysicalParameters();
	void checkForZeroExternalViscosity();
	void calculateStiffnessMatrices();
    void calculateShapeFunctionDerivatives();
    void assignNodeMasses();
    void updateNodeMasses();
    void assignElementalSurfaceAreaIndices();
    void updateNodeViscositySurfaces();
    void updateElementToConnectedNodes(vector <Node*>& Nodes);
	void assignConnectedElementsAndWeightsToNodes();
	void fixAllD(int i, bool fixWithViscosity);
	void fixAllD(Node* currNode, bool fixWithViscosity);
	void fixX(int i, bool fixWithViscosity);
	void fixX(Node* currNode, bool fixWithViscosity);
	void fixY(int i, bool fixWithViscosity);
	void fixY(Node* currNode, bool fixWithViscosity);
	void fixZ(int i, bool fixWithViscosity);
	void fixZ(Node* currNode, bool fixWithViscosity);
	void zeroForcesOnNode(int i);
    void processDisplayDataAndSave();
	//void updateDisplaySaveValuesFromRK();
	void saveStep();
	void writeSimulationSummary();
	void writeSpecificNodeTypes();
	void writeMeshFileSummary();
	void writeGrowthRatesSummary();
	void writeMyosinSummary();
	void writeECMSummary();
	void writeECMProperties();
	void writeActinSummary();
	void writeExperimentalSummary();
	void writePipetteSumary();
	void writeSaveFileStepHeader();
	void writeNodes();
	void writeElements();
	void writeSaveFileStepFooter();
	void writeTensionCompression();
    void writeGrowth();
	void writeForces();
	void writePacking();
	void writeProteins();
	void writePhysicalProp();
	void writeNodeBinding();
	void writeCollapseAndAdhesion();
	void writeGrowthRedistribution();
	void calculateMyosinForces();
	void cleanUpMyosinForces();
	void checkForMyosinUpdates();
	void updateEquilibriumMyosinsFromInputSignal(MyosinFunction* currMF);
	void calculateGrowth();
	void calculateShapeChange();
	void cleanUpGrowthRates();
	void assignIfElementsAreInsideEllipseBands();
	void checkForPinningPositionsUpdate();
	void updateRelativePositionsToApicalPositioning();
	void updatePinningPositions();
    void updateGrowthRotationMatrices();
	void cleanUpShapeChangeRates();
	void calculateGrowthUniform(GrowthFunctionBase* currGF);
	void calculateGrowthRing(GrowthFunctionBase* currGF);
	void calculateGrowthGridBased(GrowthFunctionBase* currGF);
	//void calculatePeripodialGrowthGridBased(GrowthFunctionBase* currGF);
	void calculateShapeChangeUniform (GrowthFunctionBase* currSCF);
	void calculateShapeChangeMarkerEllipseBased (GrowthFunctionBase* currSCF);
	void changeCellShapesInSystem();
	void changeCellShapeRing(int currIndexForParameters);
	void setStretch();
	void setUpClampBorders(vector<int>& clampedNodeIds);
	void moveClampedNodesForStretcher();
	void moveAFMBead();
	void recordForcesOnClampBorders();
	void setupPipetteExperiment();
	void setUpAFM();
    void addPipetteForces(gsl_matrix *gExt);
    void addMyosinForces(gsl_matrix* gExt);
	void laserAblate(double OriginX, double OriginY, double Radius);
	void laserAblateTissueType(int ablationType);
	void fillInNodeNeighbourhood();
	void setBasalNeighboursForApicalElements();
	void conserveColumnVolume();
	void fillInElementColumnLists();
	void updateElementVolumesAndTissuePlacements();
	void updateElasticPropertiesForAllNodes(); //When a simulation is read from save, the modified physical properties will be read.  The elasticity calculation parameters (lambda, mu, and related matrices) should be updated after reading the files. This function is called upon finishing reading a save file, to ensure the update.

	void clearNodeMassLists();
	void clearLaserAblatedSites();
    void manualPerturbationToInitialSetup(bool deform, bool rotate);
    void pokeElement(int elementId, double dx, double dy, double dz);
    void addCurvatureToColumnar(double h);
    void thickecECM();
    void addSoftPeriphery(double* fractions);
    void setupYsymmetricity();
    void setupXsymmetricity();
    void ablateSpcific();
    void setUpECMMimicingElements();
    void assigneElementsAtTheBorderOfECM();
    void assigneElementsAtTheBorderOfActin();
    void setUpActinMimicingElements();
    void clearScaleingDueToApikobasalRedistribution();
    void checkForVolumeRedistributionInTissue();
    void assignCompartment();
    //void setSymmetricNode(Node* currNode, double yLimPos);


public:
    bool savedStepwise;
	ofstream outputFile;
	bool displayIsOn;
	bool DisplaySave;
	bool reachedEndOfSaveFile;
	float dt;
	int timestep;
	double currSimTimeSec;
	double SimLength;	//in seconds
	string saveDirectory;
	string saveDirectoryToDisplayString;
	string inputMeshFileName;
	string outputFileString;
	bool saveImages;
	bool saveData;
	string name_saveFile;
	int imageSaveInterval;
	int dataSaveInterval;
	double EApical,EBasal,EMid;
	double EColumnarECM, EPeripodialECM;
	double poisson;
	double discProperApicalViscosity;
	double discProperBasalViscosity;
	double discProperMidlineViscosity;
	int noiseOnPysProp[4];
	bool zeroExternalViscosity[3]; //The boolean stating if there is zero external viscosity on any of the 3 dimensions
	bool extendExternalViscosityToInnerTissue;
	vector <bool> changedECM;
	double externalViscosityDPApical;
	double externalViscosityDPBasal;
	double externalViscosityPMApical;
	double externalViscosityPMBasal;
	double externalViscosityLZApical;
	double externalViscosityLZBasal;
	int MeshType;
	int Row;
	int Column;
	float SideLength;
	float zHeight;
	//bool fixWithExternalViscosity;
	bool ApicalNodeFixWithExternalViscosity;
	bool BasalNodeFixWithExternalViscosity;
	bool CircumferentialNodeFixWithHighExternalViscosity[5];
	bool NotumNodeFixWithExternalViscosity;
	double fixingExternalViscosity[3];
	bool ApicalNodeFix[3];
	bool BasalNodeFix[3];
	bool NotumNodeFix[3];
	double notumFixingRange[2];
	bool CircumferentialNodeFix[5][3]; //row 0: apical circumferece x,y,z ; row 1: basal circumference x,y,z; row 2: linker apical circumference x,y,z, row 3: linker basal circumference x,y,z, row 4: all circumference x,y,z
	double PeripodialElasticity;
	double peripodialApicalViscosity;
	double peripodialBasalViscosity;
	double peripodialMidlineViscosity;
	double currViscMidline;
	double PeripodialThicnessScale;
	double PeripodialLateralThicnessScale; //the thickness of the side region linking two layers, as a fraction of tissue thickness
	double lumenHeight;
	double lumenHeightScale;
	//linker zone parameters
	bool BaseLinkerZoneParametersOnPeripodialness;
	double LinkerZoneApicalElasticity;
	double LinkerZoneBasalYoungsModulus;
	double linkerZoneApicalViscosity;
	double linkerZoneBasalViscosity;
	double linkerZoneMidlineViscosity;

	int nGrowthFunctions;
	bool GridGrowthsPinnedOnInitialMesh;
	int nGrowthPinning; ///number of times growth pinning positions will be updated, the relative positions of elements in tissue.
	int* growthPinUpdateTime;
	bool* growthPinUpdateBools;
	int gridGrowthsInterpolationType; //0 = no interpolation, step function, 1 = linear interpolation (default = 1).
	vector <double***> GrowthMatrices;
	vector<GrowthFunctionBase*> GrowthFunctions;

	int nShapeChangeFunctions;
	vector<GrowthFunctionBase*> ShapeChangeFunctions;
	double shapeChangeECMLimit;

	bool thereIsPlasticDeformation;
	bool plasticDeformationAppliedToPeripodial;
	bool plasticDeformationAppliedToColumnar;
	bool volumeConservedInPlasticDeformation;
	double plasticDeformationHalfLife;
	double zRemodellingLowerThreshold;
	double zRemodellingUpperThreshold;

	int nMyosinFunctions;
	vector<MyosinFunction*> myosinFunctions;
	double kMyo;
	double forcePerMyoMolecule;
	bool thereIsMyosinFeedback;
	double MyosinFeedbackCap;
	bool addCurvatureToTissue;
	double tissueCurvatureDepth;
	vector <Node*> Nodes;
	vector <ShapeBase*> Elements;
	int nElements;
	int nNodes;
	double** SystemForces;
	double** PackingForces;
	//double** PackingForcesPreviousStep;
	//double** PackingForcesTwoStepsAgoStep;
	double** FixedNodeForces;
	vector <double> randomForces;
	double randomForceMean;
	double randomForceVar;
	double SystemCentre[3];
	bool needPeripodialforInputConsistency;
	bool thereIsPeripodialMembrane;
	bool AddPeripodialMembrane;
    bool symmetricY;
    bool symmetricX;
    bool addingRandomForces;
    bool conservingColumnVolumes;
	bool stretcherAttached;
	bool thereIsArtificaialRelaxation;
	bool relaxECMInArtificialRelaxation;
	double artificialRelaxationTime;
	vector <int> leftClampBorder;
	vector <int> rightClampBorder;
	double leftClampForces[3];
	double rightClampForces[3];
	bool DVClamp;
	int distanceIndex;	//the index of dimension for stretcher clamp position,
	bool PipetteSuction;
	bool ApicalSuction;
	bool TissueStuckOnGlassDuringPipetteAspiration;
	vector <int> TransientZFixListForPipette;
	int StretchInitialTime, StretchEndTime;
	int PipetteInitialStep;
	int nPipetteSuctionSteps;
	vector<double> pipetteSuctionTimes;
	vector<double> pipetteSuctionPressures;
	double pipetteCentre[3];
	double pipetteDepth;
	double pipetteInnerRadius;
    double pipetteThickness;
	double pipetteInnerRadiusSq;
	double effectLimitsInZ[2];
	double SuctionPressure[3];
	double StretchMin, StretchMax, StretchStrain;
	vector <int*> TrianglesToDraw;
	vector <double*> NodesToDraw;
	double TissueHeight;
	int TissueHeightDiscretisationLayers;
	double boundingBox[2][3];
	vector <int> pacingNodeCouples0;
	vector <int> pacingNodeCouples1;
	vector <bool> pacingNodeCouplesHaveAdhered;
	vector <int> pacingNodeSurfaceList0; // this is the id of the base node that is packing
	vector <int> pacingNodeSurfaceList1; // lists 1 to 3 are the edges of the triangle
	vector <int> pacingNodeSurfaceList2;
	vector <int> pacingNodeSurfaceList3;
	vector <int> pacingNodeEdgeList0; // this is the id of the base node that is packing
	vector <int> pacingNodeEdgeList1; // lists 1 to 2 are the nodes of the packing edge
	vector <int> pacingNodeEdgeList2;
	vector <int> pacingNodePointList0; // this is the id of the base node that is packing
	vector <int> pacingNodePointList1; // list 1 is the packing node

	vector <int> initialSignsSurfacex;
	vector <int> initialSignsSurfacey;
	vector <int> initialSignsSurfacez;
	vector <int> initialSignsEdgex;
	vector <int> initialSignsEdgey;
	vector <int> initialSignsEdgez;
	vector <int> initialSignsPointx;
	vector <int> initialSignsPointy;
	vector <int> initialSignsPointz;

	vector <double> initialWeightSurfacex;
	vector <double> initialWeightSurfacey;
	vector <double> initialWeightSurfacez;
	vector <double> initialWeightEdgex;
	vector <double> initialWeightEdgey;
	vector <double> initialWeightEdgez;
	vector <double> initialWeightPointx;
	vector <double> initialWeightPointy;
	vector <double> initialWeightPointz;

	int	nMarkerEllipseRanges;
	vector<double> markerEllipseBandXCentres;
	vector<double> markerEllipseBandR1Ranges;
	vector<double> markerEllipseBandR2Ranges;
    bool 	thereIsECMChange;
	vector <int> numberOfECMChangeEllipseBands;
	vector< vector<int> > ECMChangeEllipseBandIds;
	vector <double> ECMChangeBeginTimeInSec;
	vector <double> ECMChangeEndTimeInSec;
	vector <double>	ECMStiffnessChangeFraction;
	vector <double> ECMRenewalHalfLifeTargetFraction;
	vector <double> ECMViscosityChangeFraction;
	vector <bool> 	changeApicalECM;
	vector <bool> 	changeBasalECM;
	vector <bool> 	ECMChangeTypeIsEmergent;
	double notumECMChangeInitTime,notumECMChangeEndTime;
	double notumECMChangeFraction;
	double hingeECMChangeInitTime,hingeECMChangeEndTime;
	double hingeECMChangeFraction;
	double pouchECMChangeInitTime,pouchECMChangeEndTime;
	double pouchECMChangeFraction;

	int numberOfMyosinAppliedEllipseBands;
	vector <int> myosinEllipseBandIds;
	
	vector <bool> startedStiffnessPerturbation;
	bool ThereIsStiffnessPerturbation;
	vector <bool> ThereIsApicalStiffnessPerturbation;
	vector <bool> ThereIsBasalStiffnessPerturbation;
	vector <bool> ThereIsWholeTissueStiffnessPerturbation;
	vector <bool> ThereIsBasolateralStiffnessPerturbation;
	vector <bool> ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation;
	vector <double> stiffnessChangedToFractionOfOriginal;
	vector <double> stiffnessPerturbationBeginTimeInSec;
	vector <double> stiffnessPerturbationEndTimeInSec;
	vector <int> numberOfStiffnessPerturbationAppliesEllipseBands;
	vector< vector<int> > stiffnessPerturbationEllipseBandIds;

	double packingDetectionThreshold;
	double packingDetectionThresholdGrid[10][5];
	double packingThreshold;
	double packingMultiplier;
	double sigmoidSaturationForPacking;

	//AFM bead
	vector <int> nodesPackingToBead;
	double packingToBeadThreshold;
	double beadR;
	double 	beadPos[3];
	vector <double > initialWeightPackingToBeadx;
	vector <double > initialWeightPackingToBeady;
	vector <double > initialWeightPackingToBeadz;
	vector <double > distanceToBead;

	//soft periphery parameters:
	bool 	softPeriphery;
	double 	softDepth;
	double 	softnessFraction;
	bool 	softPeripheryBooleans[4]; //  [applyToApical]  [applyToBasal]  [applyToColumnar]  [applyToPeripodial]
    bool 	thereIsAdhesion;
    bool    collapseNodesOnAdhesion;
    bool 	adherePeripodialToColumnar;
    bool	thereIsEmergentEllipseMarking;
    bool    thereNodeCollapsing;

	bool implicitPacking;
	bool thereIsCellMigration;
	CellMigration* cellMigrationTool;
	vector <double> drawingPointsX, drawingPointsY, drawingPointsZ;

	bool thereIsExplicitECM;
	bool addLateralECMManually;
	double lateralECMThickness;
	bool thereIsExplicitActin;
	double ECMRenawalHalfLife; //The half life for ECM renewal inside plastic deformation

	int numberOfClones;
	vector<double> cloneInformationX;
	vector<double> cloneInformationY;
	vector<double> cloneInformationR;
	vector<bool> cloneInformationUsingAbsolueGrowth;
	vector<double> cloneInformationGrowth;

	bool encloseTissueBetweenSurfaces;
	double initialZEnclosementBoundaries[2];
	double finalZEnclosementBoundaries[2];
	double initialTimeToEncloseTissueBetweenSurfacesSec;
	double finalTimeToEncloseTissueBetweenSurfacesSec;
	double zEnclosementBoundaries[2];
	double packingToEnclosingSurfacesThreshold;
	double packingDetectionToEnclosingSurfacesThreshold;
	vector <int> nodesPackingToPositiveSurface;
	vector <int> nodesPackingToNegativeSurface;
	vector <double> initialWeightPackingToPositiveSurface;
	vector <double> initialWeightPackingToNegativeSurface;


	vector< vector<int> > ellipseIdsForBaseAxisBinding;
	vector< vector<bool> > ellipseBasesAreBoundOnAxis;
	bool thereIsCircumferenceXYBinding;

	int nApikobasalVolumeRedistributionFunctions;

	vector<double> apikobasalVolumeRedistributionBeginTimeInSec;
	vector<double> apikobasalVolumeRedistributionEndTimeInSec;
	vector <bool> apikobasalVolumeRedistributionFunctionShrinksApical;
	vector<int> apikobasalVolumeRedistributionFunctionEllipseNumbers;
	vector<vector <int> > apikobasalVolumeRedistributionFunctionEllipseBandIds;
	vector<double> apikobasalVolumeRedistributionScales;

	NewtonRaphsonSolver *NRSolver;

	//packi
	Simulation();
	~Simulation();
	void assignTips();
	bool readExecutableInputs(int argc, char **argv);
	bool initiateSystem();
	void calculateSystemCentre();
	void cleanMatrixUpdateData();
	void resetForces(bool resetPacking);
	void calculateApicalSize();
	void calculateBoundingBox();
    void calculateZProjectedAreas();
    void correctzProjectedAreaForMidNodes();
    void clearProjectedAreas();
    void checkForExperimentalSetupsBeforeIteration();
    void checkForExperimentalSetupsWithinIteration();
    void checkForExperimentalSetupsAfterIteration();
    void setLateralElementsRemodellingPlaneRotationMatrices();

    void updateChangeForExplicitECM(int idOfCurrentECMPerturbation);
    void updateECMRenewalHalflifeMultiplier(int idOfCurrentECMPerturbation);
    void updateChangeForViscosityBasedECMDefinition(int idOfCurrentECMPerturbation);
    void calculateChangeRatesForECM(int idOfCurrentECMPerturbation);
    void checkECMChange();
    void updateStiffnessChangeForActin(int idOfCurrentStiffnessPerturbation);
    void calculateStiffnessChangeRatesForActin(int idOfCurrentStiffnessPerturbation);
    void checkStiffnessPerturbation();
    void checkForEllipseIdUpdateWithECMDegradation();
    void checkForEmergentEllipseFormation();

    void updateOnFoldNodesFromCollapse();
    void artificialRelax();
    void checkEllipseAllocationWithCurvingNodes();
    void updateEllipseWithCollapse();
    void checkForLeftOutElementsInEllipseAssignment();
    bool runOneStep();
    void updatePlasticDeformation();
    void updateStepNR();
    void calculateNumericalJacobian(bool displayMatricesDuringNumericalCalculation);
    void updateElementPositionsinNR(gsl_matrix* uk);
    void updateNodePositionsNR(gsl_matrix* uk);
    void calculateRandomForces();
    void addRandomForces(gsl_matrix* gExt);
    void smallStrainrunOneStep();
    void packToPipetteWall();
    void calculatePacking2(double PeriThreshold, double ColThreshold);
    void addToEdgeList(Node* nodePointer, ShapeBase* elementPointer, vector <int> & edgeNodeData0, vector <int> & edgeNodeData1);
    bool checkIfPointIsOnEdge(int node0, int node1, double x, double y, double z, double& vx, double& vy, double& vz);
   	void detectPacingNodes();
   	void assignFoldRegionAndReleasePeripodial(Node* NodeMAster, Node* NodeSlave);
   	void detectPacingCombinations();
   	void cleanUpPacingCombinations();

   	void calculatePackingToEnclosingSurfacesJacobian3D(gsl_matrix* K);
	void calculatePackingToAFMBeadJacobian3D(gsl_matrix* K);
   	void calculatePackingForcesToEnclosingSurfacesImplicit3D();
   	void detectPacingToEnclosingSurfacesNodes();
   	void detectPackingToAFMBead();

    void calculatePacking();
    void calculatePackingK(gsl_matrix* K);
    void calculatePackingNumerical(gsl_matrix* K);
    void calculatePackingForcesImplicit3D();
    void calculatePackingForcesExplicit3D();
    void calculatePackingJacobian3D(gsl_matrix* K);
    void addValueToMatrix(gsl_matrix* K, int i, int j, double value);
    void addPackingForces(gsl_matrix* gExt);
	void checkPackingToPipette(bool& packsToPip, double* pos, double* pipF,double mass);
	void getNormalAndCornerPosForPacking(Node* NodePointer, ShapeBase* ElementPointer, double* normalForPacking,double* posCorner);
	void getApicalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner);
	void getBasalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner);
	inline void CapPackingForce(double& Fmag);
	void bringMyosinStimuliUpToDate();
	void redistributePeripodialMembraneForces(int RKId);
	void updateElementPositions();
	void updateMasterSlaveNodesInBinding();
	void updateElementPositionsSingle(int i );
	bool initiateSavedSystem();
	void updateOneStepFromSave();
	void alignTissueDVToXPositive();
	void alignTissueAPToXYPlane();
	bool checkFlip();
	void flagElementsThatNeedRefinement();
	void refineElements();
	void addNodesForRefinement(ShapeBase* currElement, int* newNodeIdList);
	void addElementsForRefinement(int* elementsIdsOnColumn, int* newNodeIdList);
	void wrapUpAtTheEndOfSimulation();
	void writeRelaxedMeshFromCurrentState();
	void writeMeshRemovingAblatedRegions();
	void calculateDVDistance();
	void TissueAxisPositionDisplay();
	void coordinateDisplay();
	void correctTiltedGrowthForGrowthFunction(GrowthFunctionBase* currGF);
	void calculateTiltedElementPositionOnBase(ShapeBase* currElement);
	void calculateBaseElementsFinalPosition(int Id,double DVGrowth, double APGrowth, double ABGrowth);
	void calculateCurrentElementsFinalPosition(ShapeBase* currElement);
	void fixNode0InPosition(double x, double y, double z);

    void addNodesForSideECMOnOuterCircumference (vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray , double hColumnar);
    void addSideECMElements(vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray);
    bool addSideECMLayer();

};

#endif
