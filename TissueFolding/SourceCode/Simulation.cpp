
#include "Simulation.h"
#include "Prism.h"
#include "RandomGenerator.h"
#include <string.h>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>


using namespace std;

Simulation::Simulation(){
	savedStepwise = true;
	currElementId = 0;
	ModInp = new ModelInputObject();
	SystemCentre[0]=0.0; SystemCentre[1]=0.0; SystemCentre[2]=0.0;
	TissueHeight = 0.0;
	TissueHeightDiscretisationLayers= 1;
	timestep = 0;
	currSimTimeSec = 0;
	reachedEndOfSaveFile = false;
	AddPeripodialMembrane = false;
	thereIsPeripodialMembrane = false;
	needPeripodialforInputConsistency = false;
	conservingColumnVolumes = false;
	lumenHeight = -20;
	boundingBoxSize[0]=1000.0; boundingBoxSize[1]=1000.0; boundingBoxSize[2]=1000.0;
	//columnarBoundingBoxSize[0]=1000.0; columnarBoundingBoxSize[1]=1000.0; columnarBoundingBoxSize[2]=1000.0;
	//peripodialBoundingBoxSize[0]=1000.0; peripodialBoundingBoxSize[1]=1000.0; peripodialBoundingBoxSize[2]=1000.0;
	ContinueFromSave = false;
    growthRotationUpdateFrequency = 60.0/dt;
    nElements = 0;
    nNodes = 0;
    if (growthRotationUpdateFrequency<1) {growthRotationUpdateFrequency =1;}
	setDefaultParameters();
	implicitPacking = true;
	thereIsAdhesion = false;
	collapseNodesOnAdhesion = false;
	thereNodeCollapsing = false;
	adherePeripodialToColumnar = false;
	thereIsEmergentEllipseMarking = false;
	thereIsArtificaialRelaxation = false;
	artificialRelaxationTime = -1;
	relaxECMInArtificialRelaxation = false;
	checkedForCollapsedNodesOnFoldingOnce = false;
}

Simulation::~Simulation(){
    cerr<<"destructor for simulation called"<<endl;
	delete ModInp;
	if (nGrowthPinning>0){
		delete [] growthPinUpdateTime;
		delete [] growthPinUpdateBools;
	}
	//n nodes
	for (int j=0;j<nNodes;++j){
		delete[] SystemForces[j];
		delete[] PackingForces[j];
		//delete[] PackingForcesPreviousStep[j];
		//delete[] PackingForcesTwoStepsAgoStep[j];
		delete[] FixedNodeForces[j];
	}
	delete[] SystemForces;
    delete[] PackingForces;
    //delete[] PackingForcesPreviousStep;
    //delete[] PackingForcesTwoStepsAgoStep;
    delete[] FixedNodeForces;
    cout<<"deleting elements"<<endl;
	while(!Elements.empty()){
		ShapeBase* tmp_pt;
		tmp_pt = Elements.back();
		Elements.pop_back();
		delete tmp_pt;
        //cerr<<"Element list size: "<<Elements.size()<<endl;
	}
    cout<<"deleting nodes"<<endl;
	while(!Nodes.empty()){
		Node* tmp_pt;
		tmp_pt = Nodes.back();
		Nodes.pop_back();
		delete tmp_pt;
	}
    cout<<"deleting GrowthFunctions"<<endl;
	while(!GrowthFunctions.empty()){
		GrowthFunctionBase* tmp_GF;
		tmp_GF = GrowthFunctions.back();
		GrowthFunctions.pop_back();
		delete tmp_GF;
        //cerr<<"GrowtFuncitons list size: "<<GrowthFunctions.size()<<endl;
	}
    cout<<"deleting NRSolver"<<endl;
    delete NRSolver;
    cout<<"deletion complete"<<endl;

}

void Simulation::setDefaultParameters(){
	dt = 0.01;				//sec
	SimLength = 10.0; 		//10 sec of simulation
	saveImages = false;		//do not save simulation images by default
	saveData = false;		//do not save simulation data by default
	imageSaveInterval = 60.0/dt;	//save images every minute
    dataSaveInterval  = 60.0/dt;	//save data every minute
	saveDirectory = "Not-Set";	//the directory to save the images and data points
	saveDirectoryToDisplayString  = "Not-Set"; //the file whcih will be read and displayed - no simulation
	EApical = 10.0;
	EBasal = 10.0;
	EMid = 10.0;
	PeripodialElasticity = 10.0;
	EColumnarECM = 10.0;
	EPeripodialECM = 10.0;
	poisson = 0.3;
	discProperApicalViscosity = 0.0;
	discProperBasalViscosity = 0.0;
	for (int i=0; i<3; ++i){
		zeroExternalViscosity[i] = true;
	}
	memset(noiseOnPysProp,0,4*sizeof(int));
	// The default input is a calculated mesh of width 4 elements, each element being 2.0 unit high
	// and having 1.0 unit sides of the triangles.
	MeshType = 2;
	Row = 4;
	Column = Row-2;
	SideLength=1.0;
	zHeight = 2.0;
	//fixWithExternalViscosity = false;
	ApicalNodeFixWithExternalViscosity = false;
	BasalNodeFixWithExternalViscosity = false;
	NotumNodeFixWithExternalViscosity = false;
	for (int i=0; i<5; ++i){
		CircumferentialNodeFixWithHighExternalViscosity[i] = false;
		for (int j=0; j<3; j++){
			CircumferentialNodeFix[i][j] = false;
			ApicalNodeFix[j] = false;
			BasalNodeFix[j] = false;
			NotumNodeFix[j] = false;
			fixingExternalViscosity[j] = 0;
		}
	}
	notumFixingRange[0] = 1.1;
	notumFixingRange[1] = -0.1;
	nGrowthFunctions = 0;
	GridGrowthsPinnedOnInitialMesh = false;
	nGrowthPinning = 0;
	gridGrowthsInterpolationType = 1;
	nShapeChangeFunctions = 0;
	TensionCompressionSaved = true;
    GrowthSaved = true;
    GrowthRateSaved = true;
	ForcesSaved = true;
	ProteinsSaved = true;
	PackingSaved = true;
	growthRedistributionSaved = true;
	nodeBindingSaved = true;
	physicalPropertiesSaved = true;
	collapseAndAdhesionSaved = true;
	PeripodialElasticity = 0.0;
	peripodialApicalViscosity = discProperApicalViscosity;
	peripodialBasalViscosity  = discProperBasalViscosity;
	PeripodialThicnessScale = 1.0;
	PeripodialLateralThicnessScale = 0.3;
	lumenHeightScale = 0.3;
	dorsalTipIndex = 0;
	ventralTipIndex = 1;
	anteriorTipIndex = 0;
	posteriorTipIndex = 0;
	stretcherAttached = false;
	recordForcesOnFixedNodes = false;
	distanceIndex = false;
	StretchDistanceStep = 0.0;
	DVClamp = true;
	StretchInitialTime = -100;
	StretchEndTime = -100;
	PipetteSuction = false;
	PipetteInitialStep= -100;
	nPipetteSuctionSteps = 0;
	ApicalSuction = true;
	TissueStuckOnGlassDuringPipetteAspiration = true;
	pipetteCentre[0] = 0.0;
	pipetteCentre[1] = 0.0;
	pipetteCentre[2] = 0.0;
	pipetteDepth = 0.0;
	pipetteInnerRadius =0.0;
	SuctionPressure[0] = 0.0;
	SuctionPressure[1] = 0.0;
	SuctionPressure[2] = 0.0;
	pipetteInnerRadiusSq = pipetteInnerRadius*pipetteInnerRadius;
    pipetteThickness = 11.7;
	effectLimitsInZ[0] = pipetteCentre[2] - pipetteDepth;
	effectLimitsInZ[1] = pipetteCentre[2] + pipetteDepth;

	boundingBox[0][0] =  1000.0;	//left x
	boundingBox[0][1] =  1000.0;	//low y
	boundingBox[0][2] =  1000.0;	//bottom z
	boundingBox[1][0] = -1000.0;	//right x
	boundingBox[1][1] = -1000.0;	//high y
	boundingBox[1][2] = -1000.0;	//top z

	addCurvatureToTissue = false;
	tissueCurvatureDepth = 0.0;
	symmetricY = false;
	symmetricX = false;
	addingRandomForces = false;
	randomForceMean = 0.0;
	randomForceVar = 0.0;
	packingDetectionThreshold = 5.5;
	packingThreshold = 6.0;

	softPeriphery = false;
	softDepth = 0.0;
	softnessFraction = 0.0;
	softPeripheryBooleans[0] = false; //apical
	softPeripheryBooleans[1] = false; //basal
	softPeripheryBooleans[2] = false; //columnar
	softPeripheryBooleans[3] = false; //peripodial

	thereIsPlasticDeformation = false;
	plasticDeformationAppliedToPeripodial = false;
	plasticDeformationAppliedToColumnar = false;
	volumeConservedInPlasticDeformation = false;
	plasticDeformationHalfLife = 0.0;
	zRemodellingLowerThreshold = 0.5;
	zRemodellingUpperThreshold = 2.0;

	kMyo = 0.09873;
	forcePerMyoMolecule = 1.0;
	thereIsMyosinFeedback = false;
	MyosinFeedbackCap = 0.0;

	BaseLinkerZoneParametersOnPeripodialness = true;
	LinkerZoneApicalElasticity = 0.0;
	LinkerZoneBasalYoungsModulus = 0.0;
	linkerZoneApicalViscosity = discProperApicalViscosity;
	linkerZoneBasalViscosity = discProperBasalViscosity;

	extendExternalViscosityToInnerTissue = false;
	externalViscosityDPApical = 0.0;
	externalViscosityDPBasal  = 0.0;
	externalViscosityPMApical = 0.0;
	externalViscosityPMBasal  = 0.0;
	externalViscosityLZApical = 0.0;
	externalViscosityLZBasal  = 0.0;

	numberOfMyosinAppliedEllipseBands = 0;


	thereIsECMChange = false;

	thereIsCellMigration = false;
	thereIsExplicitECM = false;
	addLateralECMManually = false;
	lateralECMThickness = 0.0;
	ECMRenawalHalfLife = 0.0;
	thereIsExplicitActin = false;

	nMarkerEllipseRanges = 0;
	ThereIsStiffnessPerturbation = false;

	encloseTissueBetweenSurfaces =  false;
	zEnclosementBoundaries[0] = -1000;
	zEnclosementBoundaries[1] = 1000;
	initialZEnclosementBoundaries[0] = -1000;
	initialZEnclosementBoundaries[1] = 1000;
	finalZEnclosementBoundaries[0] = -1000;
	finalZEnclosementBoundaries[1] = 1000;

	thereIsCircumferenceXYBinding= false;

	notumECMChangeInitTime = 10000000.0;
	notumECMChangeEndTime = 0.0;
	notumECMChangeFraction =1.0;
	hingeECMChangeInitTime = 10000000.0;
	hingeECMChangeEndTime = 0.0;
	hingeECMChangeFraction =1.0;
	pouchECMChangeInitTime = 10000000.0;
	pouchECMChangeEndTime = 0.0;
	pouchECMChangeFraction =1.0;

	boundLateralElements = false;

	beadR = 0;
	beadPos[0] = -1000.0;
	beadPos[1] = -1000.0;
	beadPos[2] = -1000;
}

bool Simulation::readExecutableInputs(int argc, char **argv){
	int i = 1;
	bool Success = true;
	while(i<argc){
		const char *inptype = argv[i];
		if (string(inptype) == "-mode"){
			Success = readModeOfSim(i, argc, argv);
		}
		else if (string(inptype) == "-i"){
			Success = readParameters(i, argc, argv);
		}
		else if (string(inptype) == "-od"){
			Success = readOutputDirectory(i, argc, argv);
		}
		else if (string(inptype) == "-dInput"){
			Success = readSaveDirectoryToDisplay(i, argc, argv);
		}
		else {
			cerr<<"Please enter a valid option key: {-mode,-i, -od, -dInput}, current string: "<<inptype<<endl;
			return false;
		}
		i++;
		if (!Success){
			return Success;
		}
	}
	Success = checkInputConsistency();
	return Success;
}

bool Simulation::readModeOfSim(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
		cerr<<" input the mode of simulation: {DisplaySave, SimulationOnTheGo, ContinueFromSave, Default}"<<endl;
		return false;
	}
	const char* inpstring = argv[i];
	if (string(inpstring) == "DisplaySave"){
		DisplaySave = true;
		return true;
	}
	else if (string(inpstring) == "SimulationOnTheGo" || string(inpstring) == "Default"){
		DisplaySave = false;
		return true;
	}
	else if (string(inpstring) == "ContinueFromSave"){
		ContinueFromSave = true;
		DisplaySave = false;
		//DisplaySave = true;
		return true;
	}
	else{
		cerr<<"Please provide input mode: -mode {DisplaySave, SimulationOnTheGo, ContinueFromSave, Default}";
		return false;
	}
}

bool Simulation::readParameters(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
		cerr<<" input the model input file"<<endl;
		return false;
	}
	//cerr<<"Reading parameter input file"<<endl;
	const char* inpstring = argv[i];
	ModInp->Sim=this;
	ModInp->parameterFileName =  inpstring;
	//cerr<<" Reading parameters from file: "<<ModInp->parameterFileName<<endl;
	bool Success = ModInp->readParameters();
	if (!Success){
		return Success;
	}
	return true;
}

bool Simulation::readSaveDirectoryToDisplay(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
		cerr<<" input the save directory, contents of which will be displayed"<<endl;
		return false;
	}
	const char* inpstring = argv[i];
	saveDirectoryToDisplayString = string(inpstring);
	return true;
}

bool Simulation::readOutputDirectory(int& i, int argc, char **argv){
	i++;
	if (i >= argc){
		cerr<<" input the save directory"<<endl;
		return false;
	}
	const char* inpstring = argv[i];
	//This will set the save directory, but will not change the safe file boolean.
	//If your model input file states no saving, then the error and output files
	//will be directed into this directory, but the frame saving will not be toggled
	saveDirectory= string(inpstring);
	return true;
}

bool Simulation::readFinalSimulationStep(){
	bool success  = openFilesToDisplay();
	if (!success){
		return false;
	}
	//backing up data save interval and time step before reading the system summary:
	double timeStepCurrentSim = dt;
	int dataSaveIntervalCurrentSim = dataSaveInterval;
	//bool ZerothFrame = true;
	//reading system properties:
	success  = readSystemSummaryFromSave();
	if (!success){
		return false;
	}
	string currline;
	while(reachedEndOfSaveFile == false){
		//skipping the header:
		getline(saveFileToDisplayMesh,currline);
		cerr<<" currline in read last step: "<<currline<<endl;
		if(saveFileToDisplayMesh.eof()){
			reachedEndOfSaveFile = true;
			break;
		}
		success = readNodeDataToContinueFromSave();
		cout<<"dt after  readNodeDataToContinueFromSave "<<dt<<" timeStepCurrentSim: "<<timeStepCurrentSim<<" dataSaveInterval: "<<dataSaveInterval<<" dataSaveIntervalCurrentSim: "<<dataSaveIntervalCurrentSim<<endl;
		if (!success){
			return false;
		}
		readElementDataToContinueFromSave();
		cout<<"dt after  readElementDataToContinueFromSave "<<dt<<" timeStepCurrentSim: "<<timeStepCurrentSim<<" dataSaveInterval: "<<dataSaveInterval<<" dataSaveIntervalCurrentSim: "<<dataSaveIntervalCurrentSim<<endl;
		//assignPhysicalParameters();
		//cout<<" reading tension"<<endl;
		if (TensionCompressionSaved){
			readTensionCompressionToContinueFromSave();
		}
		//cout<<" reading growth"<<endl;
        if (GrowthSaved){
            readGrowthToContinueFromSave();
        }
		//cout<<" reading growth rates"<<endl;
        if (GrowthRateSaved){
            readGrowthRateToContinueFromSave();
        }
		//cout<<" reading proteins"<<endl;
        if (ProteinsSaved){
        	readProteinsToContinueFromSave();
		}
		//cout<<" reading physical properties"<<endl;
        if (physicalPropertiesSaved){
        	readPhysicalPropToContinueFromSave();
        }
		//cout<<" reading forces"<<endl;
		if (ForcesSaved){
    		updateForcesFromSave();
    	}
		//cout<<" reading packing"<<endl;
    	if (PackingSaved){
    		updatePackingFromSave();
    	}
		//cout<<" reading growth redistribution"<<endl;
    	if (growthRedistributionSaved){
    		readGrowthRedistributionToContinueFromSave();
    	}
		//cout<<" reading node binding"<<endl;
    	if (nodeBindingSaved){
    		readNodeBindingToContinueFromSave();
    	}
		//cout<<" reading collapse and adhesion"<<endl;
    	if(collapseAndAdhesionSaved){
    		readCollapseAndAdhesionToContinueFromSave();
    	}
		//if (ZerothFrame){
		//	ZerothFrame = false;
		//	cout<<"dt after  ZerothFrame if clause "<<dt<<" timeStepCurrentSim: "<<timeStepCurrentSim<<" dataSaveInterval: "<<dataSaveInterval<<endl;
		//}
		//else{
			timestep = timestep + dataSaveInterval;
			currSimTimeSec += dt*dataSaveInterval;
		//}
		cout<<"current time step "<<timestep<<" currSimTimeSec: "<<currSimTimeSec<<endl;

		//skipping the footer:
		getline(saveFileToDisplayMesh,currline);
		while (currline.empty() && !saveFileToDisplayMesh.eof()){
			//skipping empty line
			getline(saveFileToDisplayMesh,currline);
		}
		//Now I will save it again:
		cout<<"saving"<<endl;
		saveStep();
	}
	cout<<"carry out updates"<<endl;
	updateElementVolumesAndTissuePlacements();
	updateElasticPropertiesForAllNodes();
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	clearLaserAblatedSites();
    calculateShapeFunctionDerivatives();
	updateElementPositions();
	calculateBoundingBox();
	bringMyosinStimuliUpToDate();
	//During a simulation, the data is saved after the step is run, and the time step is incremented after the save. Now I have read in the final step,
	//I need to increment my time step to continue the next time step from here.
	//currSimTimeSec += dt*dataSaveInterval;
	//timestep++;

	//bring the time step and data save time steps to the main modelinput:
	dataSaveInterval = dataSaveIntervalCurrentSim;
	dt = timeStepCurrentSim;

	return true;
}

void Simulation::updateMasterSlaveNodesInBinding(){
	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		for (int dim = 0; dim<3; ++dim){
			if ((*itNode)->slaveTo[dim] > -1){
				int slaveDof = (*itNode)->Id*3+dim;
				int masterDof = (*itNode)->slaveTo[dim]*3+dim;
				vector <int> fix;
				fix.push_back(slaveDof);
				fix.push_back(masterDof);
				NRSolver->slaveMasterList.push_back(fix);
				NRSolver->boundNodesWithSlaveMasterDefinition = true;
			}
		}
	}
}

bool Simulation::checkInputConsistency(){
	if (ContinueFromSave &&  saveDirectoryToDisplayString == "Not-Set"){
		cerr <<"The mode is set to continue from saved simulation, please provide an input directory containing the saved profile, using -dInput"<<endl;
		return false;
	}
	if (saveData || saveImages){
		if (saveDirectory == "Not-Set"){
			cerr <<"Modelinput file requires saving, please provide output directory, using -od tag"<<endl;
			return false;
		}
	}
	if (DisplaySave){
		if (saveDirectoryToDisplayString == "Not-Set"){
			cerr <<"The mode is set to display from save, please provide an input directory, using -dInput"<<endl;
			return false;
		}
	}
	if (AddPeripodialMembrane == false){
		for (int i=0; i<nGrowthFunctions; ++i){
			if(GrowthFunctions[i]->applyToPeripodialMembrane){
				cerr<<"There is no peripodial membrane, while growth function "<<i<<" is applicable to peropodial membrane, further checks needed"<<endl;
				needPeripodialforInputConsistency = true;
			}
		}
	}
	return true;
}

bool Simulation::initiateSystem(){
	bool Success = openFiles();
	if (!Success){
		return Success;
	}
	if (MeshType == 1){
		Success = initiateMesh(MeshType, zHeight); //zHeight
	}
	else if (MeshType == 2){
		Success = initiateMesh(MeshType, Row, Column,  SideLength,  zHeight);
	}
	else if(MeshType == 4){
		cout<<"mesh type : 4"<<endl;
		Success = initiateMesh(MeshType);
	}
	if (!Success){
		return Success;
	}
	if (symmetricY || symmetricX){
		 clearCircumferenceDataFromSymmetricityLine();
	}
	cout<<"cleared symmetricity"<<endl;

	Success = checkIfThereIsPeripodialMembrane();
	cout<<"calculating tissue height"<<endl;
	Success = calculateTissueHeight(); //Calculating how many layers the columnar layer has, and what the actual height is.
	cout<<"calculating bounding box"<<endl;
	calculateBoundingBox();
	cout<<"assigning Initial Z Positions"<<endl;

	assignInitialZPositions();
	if (!Success){
		return Success;
	}
	if (AddPeripodialMembrane){
		if (thereIsPeripodialMembrane){
			Success = false;
			cerr<<"Error-there is already peripodial membrane added to the mesh, but modelinput file asks to add another"<<endl;
		}
		else{
			//Success = addCurvedPeripodialMembraneToTissue();
			Success = addStraightPeripodialMembraneToTissue();
		    calculateBoundingBox();
			if (Success){
				thereIsPeripodialMembrane = true;
			}
		}
	}
	if (needPeripodialforInputConsistency){
		if (!thereIsPeripodialMembrane){
			cerr<<"There is no peripodial membrane but at least one growth function desires one"<<endl;
			Success = false;
		}
	}
	if (thereIsExplicitECM){
		if (addLateralECMManually){
			addSideECMLayer();
		}
		setUpECMMimicingElements();
		//TO DO: NEED TO UPDATE THE PHYSICAL PROPERTIES AFTER THIS!
	}
	if (addCurvatureToTissue){
		addCurvatureToColumnar(tissueCurvatureDepth);
	}
	if (thereIsExplicitActin){
		setUpActinMimicingElements();
	}
	setBasalNeighboursForApicalElements();
	fillInNodeNeighbourhood();
	fillInElementColumnLists();
	checkForNodeFixing();
	assignTips();
	if (!Success){
		return Success;
	}
	initiateSystemForces();
	calculateSystemCentre();
	assignPhysicalParameters();
	checkForZeroExternalViscosity();
	calculateShapeFunctionDerivatives();
	assignNodeMasses();
	assignElementalSurfaceAreaIndices();
	assignConnectedElementsAndWeightsToNodes();
    alignTissueDVToXPositive();
    //alignTissueAPToXYPlane();
    calculateBoundingBox();
    calculateDVDistance();
    assignCompartment();
	vector<ShapeBase*>::iterator itElement;
    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	(*itElement)->calculateRelativePosInBoundingBox(boundingBox[0][0],boundingBox[0][1],boundingBoxSize[0],boundingBoxSize[1]);
    }
    cout<<" positions in relative bounding box calculated"<<endl;

    updateRelativePositionsToApicalPositioning();

    cout<<" positions are update to apical"<<endl;

    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        (*itElement)->setInitialRelativePosInBoundingBox();
    }
    cout<<" initial positions in relative bounding box set"<<endl;
    induceClones();
    //bool softHinge = true;
    //double hingeLimits[2] ={0.40, 0.65};
    //double softnessLevel = 0.5;
    //if (softHinge){
    //	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    //		(*itElement)->assignSoftHinge(hingeLimits[0], hingeLimits[1],softnessLevel);
    //	}
    //}
	if (stretcherAttached){
		setStretch();
	}
    cout<<"setting the pipette"<<endl;
    //setUpAFM();
	if(PipetteSuction){
		setupPipetteExperiment();
	}

	if (symmetricY){
		setupYsymmetricity();
	}
	if (symmetricX){
		setupXsymmetricity();
	}
	assignIfElementsAreInsideEllipseBands();
	if (ContinueFromSave){
		cout<<"Reading Final SimulationStep: "<<endl;
		Success = readFinalSimulationStep();
		if (!Success){
			return Success;
		}
	}
	if (saveData){
		cout<<"writing the summary current simulation parameters"<<endl;
		writeSimulationSummary();
		writeSpecificNodeTypes();
	}
	nNodes = Nodes.size();
	nElements = Elements.size();
	//initiating the NR solver object, that generates the necessary matrices
	//for solving the tissue dynamics
    NRSolver = new NewtonRaphsonSolver(Nodes[0]->nDim,nNodes);
    if (ContinueFromSave){
    	updateMasterSlaveNodesInBinding();
    }
    checkForNodeBinding();
	//peripodial z binding
    if (adherePeripodialToColumnar){
		bool thereIsBinding = bindPeripodialToColumnar();
		if (thereIsBinding){
			NRSolver->boundNodesWithSlaveMasterDefinition = true;
		}
    }
	if (thereIsCellMigration) {
    	cout<<"initiation of cell migration"<<endl;
    	double cellMigrationOriginAngle = M_PI/2.0;
    	cellMigrationTool = new CellMigration(nElements,0.5); //50% leaves the tissue region per hour
        cellMigrationTool->assignElementRadialVectors(Elements);
        cellMigrationTool->assignOriginOfMigration(Elements, cellMigrationOriginAngle);
        cellMigrationTool->assignElementConnectivity(Nodes,Elements);
        cellMigrationTool->generateListOfRateFractions(Elements);
    }
    if (thereIsPlasticDeformation){
    	setLateralElementsRemodellingPlaneRotationMatrices();
    }
    cout<<" system initiated"<<endl;
	return Success;
}

void Simulation::checkForNodeBinding(){
	if (thereIsCircumferenceXYBinding){
		bool thereIsBinding = bindCircumferenceXY();
		if (thereIsBinding){
			NRSolver->boundNodesWithSlaveMasterDefinition = true;
		}
	}
	if (ellipseIdsForBaseAxisBinding.size()>0){
		bool thereIsBinding = bindEllipseAxes();
		if (thereIsBinding){
			//I will not equate the values, if the parameter of NRSolver
			//has been modified to be true before, it should stay true.
			//If my decisin from abpve function is false(I have not bound any ellipses)
			//then I should not alter this parameter;
			NRSolver->boundNodesWithSlaveMasterDefinition = true;
		}
	}
}



bool Simulation::bindEllipseAxes(){
	bool thereIsBinding = false;
	int dim = 3;
	int nEllipseFunctions = ellipseIdsForBaseAxisBinding.size();
	for (int ellipseFunctionIterator=0; ellipseFunctionIterator<nEllipseFunctions; ++ellipseFunctionIterator){
		//checking one ellipse base binding rule:
		vector <int> nodeIds;
		int nBoundEllipses = ellipseIdsForBaseAxisBinding[ellipseFunctionIterator].size();
		int masterNodeId =  nNodes*2;
		cout<<"checking binding rule :"<<ellipseFunctionIterator<<" - axis: ";
		cout<<ellipseBasesAreBoundOnAxis[ellipseFunctionIterator][0]<<" ";
		cout<<ellipseBasesAreBoundOnAxis[ellipseFunctionIterator][1]<<" ";
		cout<<ellipseBasesAreBoundOnAxis[ellipseFunctionIterator][2]<<" ";
		cout<<" - ellipses: ";
		for (int ellipseIdIterator =0; ellipseIdIterator<nBoundEllipses;++ellipseIdIterator){
			cout<<ellipseIdsForBaseAxisBinding[ellipseFunctionIterator][ellipseIdIterator]<<"  ";
		}
		cout<<endl;

		//loop over all nodes, to see if
		//* node is inside an ellipse band
		//* it is basal
		//* it is inside an ellipse of interest
		for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ((*itNode)->insideEllipseBand && (*itNode)->tissuePlacement == 0){//inside ellipse band and basal
				//checking one ellipse base binding rule:
				for (int ellipseIdIterator =0; ellipseIdIterator<nBoundEllipses;++ellipseIdIterator){
					if ((*itNode)->coveringEllipseBandId == ellipseIdsForBaseAxisBinding[ellipseFunctionIterator][ellipseIdIterator]){
						//cout<<"Node : "<<(*itNode)->Id<<" is basal and inside ellipse band: "<<(*itNode)->coveringEllipseBandId<<endl;
						nodeIds.push_back((*itNode)->Id);
						if ((*itNode)->Id<masterNodeId){
							//I will avoid making a fixed node my master:
							bool skipThisNode = false;
							for (int dimIterator = 0; dimIterator <3; ++dimIterator){
								if (ellipseBasesAreBoundOnAxis[ellipseFunctionIterator][dimIterator] &&
										(*itNode)->FixedPos[dimIterator]){
									skipThisNode = true;
									break;
								}
							}
							if (!skipThisNode){
								masterNodeId = (*itNode)->Id;
							}
						}
						break;
					}
				}
			}
		}
		if (masterNodeId ==  nNodes*2){
			//I could not pick a master node, as all the nodes within the marker ellipse region
			//are already bound in at least on of their qualifying axes. Remove this fixed axis
			//from the requirements of the ellipse and restart simulation!
			cout<<" skipping ellipse binding function "<<ellipseFunctionIterator<<", due to clash with node fixing!"<<endl;
			cerr<<" skipping ellipse binding function "<<ellipseFunctionIterator<<", due to clash with node fixing!"<<endl;
			return thereIsBinding;
		}
		//Now I have a list of all nodes that are to be bound;
		//I have the minimum node id, that is to be the master of this selection
		int n= nodeIds.size();
		int dofXmaster  = masterNodeId*dim;
		int dofYmaster  = dofXmaster+1;
		int dofZmaster  = dofXmaster+2;
		for (int i=0;i<n;++i){
			if (nodeIds[i] != masterNodeId){
				int slaveNodeId = nodeIds[i];
				int dofXslave  = slaveNodeId*dim;
				int dofYslave  = dofXslave+1;
				int dofZslave  = dofXslave+2;
				if (ellipseBasesAreBoundOnAxis[ellipseFunctionIterator][0]){
					//x axis is bound:
					//see explanation under bindCircumferenceXY function for this check
					if (!Nodes[slaveNodeId]->FixedPos[0]){
						vector <int> fixX;
						fixX.push_back(dofXslave);
						fixX.push_back(dofXmaster);
						NRSolver->slaveMasterList.push_back(fixX);
						Nodes[slaveNodeId]->slaveTo[0] = masterNodeId;
						Nodes[masterNodeId]->isMaster[0] = true;
						thereIsBinding = true;
					}
				}
				if (ellipseBasesAreBoundOnAxis[ellipseFunctionIterator][1]){
					//y axis is bound:
					if (!Nodes[slaveNodeId]->FixedPos[1]){
						vector <int> fixY;
						fixY.push_back(dofYslave);
						fixY.push_back(dofYmaster);
						NRSolver->slaveMasterList.push_back(fixY);
						Nodes[slaveNodeId]->slaveTo[1] = masterNodeId;
						Nodes[masterNodeId]->isMaster[1] = true;
						thereIsBinding = true;
					}
				}
				if (ellipseBasesAreBoundOnAxis[ellipseFunctionIterator][2]){
					//z axis is bound:
					if (!Nodes[slaveNodeId]->FixedPos[2]){
						vector <int> fixZ;
						fixZ.push_back(dofZslave);
						fixZ.push_back(dofZmaster);
						NRSolver->slaveMasterList.push_back(fixZ);
						Nodes[slaveNodeId]->slaveTo[2] = masterNodeId;
						Nodes[masterNodeId]->isMaster[2] = true;
						thereIsBinding = true;
					}
				}
			}
		}
	}
	return thereIsBinding;
}

bool Simulation::bindPeripodialToColumnar(){
	bool thereIsBinding = false;
	int dim = 3;
	vector<int> masterIds,slaveIds;
	//make a list of master slave couples between columnar and peripodial:
	for (vector<Node*>::iterator itPeripodialNode = Nodes.begin();itPeripodialNode<Nodes.end(); ++itPeripodialNode){
		if ((*itPeripodialNode)->tissueType != 1 || (*itPeripodialNode)->tissuePlacement != 1 ){
			//not peripodial node or is not apical
			continue;
		}
		if ((*itPeripodialNode)->hasLateralElementOwner){
			continue;
		}
		for (vector<Node*>::iterator itColumnarNode = Nodes.begin();itColumnarNode<Nodes.end(); ++itColumnarNode){
			if ((*itColumnarNode)->tissueType != 0 || (*itColumnarNode)->tissuePlacement != 1){
				//not columnar node or is not apical
				continue;
			}
			if ((*itColumnarNode)->hasLateralElementOwner){
				continue;
			}
			double threshold = 1E-5;
			if(    (*itColumnarNode)->Position[0] < (*itPeripodialNode)->Position[0]+threshold
				&& (*itColumnarNode)->Position[0] > (*itPeripodialNode)->Position[0]-threshold
				&& (*itColumnarNode)->Position[1] < (*itPeripodialNode)->Position[1]+threshold
				&& (*itColumnarNode)->Position[1] > (*itPeripodialNode)->Position[1]-threshold
					){
				masterIds.push_back((*itColumnarNode)->Id);
				slaveIds.push_back((*itPeripodialNode)->Id);
				(*itColumnarNode)->attachedToPeripodial=true;
				cout<<" master: "<<(*itColumnarNode)->Id<<" slave: "<<(*itPeripodialNode)->Id<<endl;
			}
		}
	}
	//Now go through the list, check if the binding is feasible:
	int n = masterIds.size();
	for (int idMasterSlaveCouple=0;idMasterSlaveCouple<n;++idMasterSlaveCouple){
		int masterNodeId = masterIds[idMasterSlaveCouple];
		int slaveNodeId = slaveIds[idMasterSlaveCouple];

		for (int i=0; i<3; ++i){ //int i = 2;// the z dimension!
			if (Nodes[masterNodeId]->FixedPos[i]){
				Nodes[slaveNodeId]->FixedPos[i]=true;
			}
			//not using an else, as the slave could be fixed in given dimension independent of the master
			if (!Nodes[slaveNodeId]->FixedPos[i]){
				int dofmaster = masterNodeId*dim+i;
				int dofslave  = slaveNodeId*dim+i;
				//if slave is already slave of another node
				//make the master of the slave the new slave
				//   algorithm will take care of the rest.
				if (Nodes[slaveNodeId]->slaveTo[i] > -1){
					//Current slave has a master. If this master is the current master, or the master of current maste, than dont do anything, all is fine:
					if(Nodes[slaveNodeId]->slaveTo[i] == masterNodeId || Nodes[slaveNodeId]->slaveTo[i] == Nodes[masterNodeId]->slaveTo[i] ){
						continue;
					}
					dofslave = Nodes[slaveNodeId]->slaveTo[i]*dim+i;
					slaveNodeId = Nodes[slaveNodeId]->slaveTo[i];
				}
				//cout<<"DOF not fixed,  initial              : "<<dofmaster<<" "<<dofslave<<endl;
				//check if the master dof is already bound to something:
				NRSolver->checkMasterUpdate(dofmaster,masterNodeId);
				//It may have been that the slave was the master of the master node.
				//Now I have updated the master to the slave, and they are equal.
				//I do not need to add anything, as the master-slave relation is already implemented.
				//cout<<"DOF not fixed, after master update   : "<<dofmaster<<" "<<dofslave<<endl;
				if (dofmaster != dofslave){
					bool continueAddition =  NRSolver->checkIfCombinationExists(dofslave,dofmaster);
					//cout<<"DOF not fixed, continueAddition? : "<<continueAddition<<endl;
					if (continueAddition){
						bool madeChange = NRSolver->checkIfSlaveIsAlreadyMasterOfOthers(dofslave,dofmaster);
						if (madeChange){
							for (int nodeIt = 0 ; nodeIt<Nodes.size(); ++nodeIt){
								if(Nodes[nodeIt]->slaveTo[i]==slaveNodeId){
									Nodes[nodeIt]->slaveTo[i]=masterNodeId;
								}
							}
						}
						vector <int> fixDOF;
						fixDOF.push_back(dofslave);
						fixDOF.push_back(dofmaster);
						NRSolver->slaveMasterList.push_back(fixDOF);
						Nodes[slaveNodeId]->slaveTo[i] = masterNodeId;
						Nodes[masterNodeId]->isMaster[i] = true;
						thereIsBinding = true;
					}
				}
			}
		}
	}
	return thereIsBinding;
}

bool Simulation::areNodesToCollapseOnLateralECM(int slaveNodeId, int masterNodeId){
	if(Nodes[slaveNodeId]->hasLateralElementOwner || Nodes[masterNodeId]->hasLateralElementOwner){
		cout<<"binding nodes: "<<slaveNodeId<<" "<<masterNodeId<<" on single element collapse, but will not collapse nodes, as they are too close on lateral element"<<endl;
		return true;
	}
	int nElement = Nodes[slaveNodeId]->connectedElementIds.size();
	for (int elementCounter = 0; elementCounter<nElement; elementCounter++){
		int elementId = Nodes[slaveNodeId]->connectedElementIds[elementCounter];
		if (Elements[elementId]->isECMMimimcingAtCircumference){
			cout<<"binding nodes: "<<slaveNodeId<<" "<<masterNodeId<<" on single element collapse, but will not collapse nodes, as they are too close on circumferential element - slave"<<endl;
			return true;
		}
	}
	nElement = Nodes[masterNodeId]->connectedElementIds.size();
	for (int elementCounter = 0; elementCounter<nElement; elementCounter++){
		int elementId = Nodes[masterNodeId]->connectedElementIds[elementCounter];
		if (Elements[elementId]->isECMMimimcingAtCircumference){
			cout<<"binding nodes: "<<slaveNodeId<<" "<<masterNodeId<<" on single element collapse, but will not collapse nodes, as they are too close on circumferential element - master"<<endl;
			return true;
		}
	}
	return false;
}

bool Simulation::checkEdgeLenghtsForBindingPotentiallyUnstableElements(){
	bool thereIsBinding = false;
	int dim = 3;
	vector<int> masterIdsBulk,slaveIdsBulk;
	if (boundLateralElements == false){
		for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			boundLateralElements = true;
			//cout<<"checking element "<<(*itElement)->Id<<endl;
			(*itElement)->checkEdgeLenghtsForBinding(masterIdsBulk,slaveIdsBulk);
			//if ((*itElement)->Id == 17851){
			//	cout<<" diaplaying node positions for element 17851"<<endl;
			//	(*itElement)->displayPositions();
			//}
			int selectedPair[2] = {0,0};
			if((*itElement)->isECMMimimcingAtCircumference && (*itElement)->tissuePlacement == 0){//check lateral elements and bind them:
				int* nodeIds = (*itElement)->getNodeIds();
				double* vec= new double[3];
				for (int idim =0; idim<3; idim++){
					vec[idim] = Nodes[nodeIds[0]]->Position[idim]-Nodes[nodeIds[1]]->Position[idim];
				}
				double L = (*itElement)->normaliseVector3D(vec);
				if (L < 0.5){
					selectedPair[0] = 0;
					selectedPair[1] = 1;
				}
				else{
					for (int idim =0; idim<3; idim++){
						vec[idim] = Nodes[nodeIds[0]]->Position[idim]-Nodes[nodeIds[2]]->Position[idim];
					}
					double L = (*itElement)->normaliseVector3D(vec);
					if (L < 0.5){
						selectedPair[0] = 0;
						selectedPair[1] = 2;
					}
					else{
						selectedPair[0] = 1;
						selectedPair[1] = 2;
					}
				}
				//cout<<"selected "<<selectedPair[0]<<" "<<selectedPair[1]<<endl;
				masterIdsBulk.push_back(nodeIds[selectedPair[0]]);
				slaveIdsBulk.push_back(nodeIds[selectedPair[1]]);
				masterIdsBulk.push_back(nodeIds[selectedPair[0]+3]);
				slaveIdsBulk.push_back(nodeIds[selectedPair[1]+3]);
				for (int i=1; i < TissueHeightDiscretisationLayers;++i){
					int idOfElementOnSameColumn = (*itElement)->elementsIdsOnSameColumn[i];
					nodeIds = Elements[idOfElementOnSameColumn]->getNodeIds();
					masterIdsBulk.push_back(nodeIds[selectedPair[0]]);
					slaveIdsBulk.push_back(nodeIds[selectedPair[1]]);
					masterIdsBulk.push_back(nodeIds[selectedPair[0]+3]);
					slaveIdsBulk.push_back(nodeIds[selectedPair[1]+3]);
				}
			}
		}
	}

	//masterIdsBulk.push_back(11124);
	//masterIdsBulk.push_back(11130);

	//slaveIdsBulk.push_back(5);
	//slaveIdsBulk.push_back(543);

	//clean duplicates:
	vector<int> masterIds,slaveIds;
	int n = masterIdsBulk.size();
	//n=0;
	//masterIds.push_back(5467);
	//slaveIds.push_back(5659);
	for (int i=0;i<n;++i){
		bool duplicate = false;
		for (unsigned int j=0; j<masterIds.size(); ++j){
			if (masterIdsBulk[i] == masterIds[j] && slaveIdsBulk[i] == slaveIds[j]){
				duplicate = true;
			}
			if (masterIdsBulk[i] == slaveIds[j] && slaveIdsBulk[i] == masterIds[j]){
				duplicate = true;
			}
		}
		if (!duplicate){
			masterIds.push_back(masterIdsBulk[i]);
			slaveIds.push_back(slaveIdsBulk[i]);
		}
	}
	//Now go through the list, check if the binding is feasible:
	n = masterIds.size();
	for (int idMasterSlaveCouple=0;idMasterSlaveCouple<n;++idMasterSlaveCouple){
		int masterNodeId = masterIds[idMasterSlaveCouple];
		int slaveNodeId = slaveIds[idMasterSlaveCouple];
		if (binary_search(Nodes[masterNodeId]->collapsedWith.begin(), Nodes[masterNodeId]->collapsedWith.end(),slaveNodeId)){
			//the couple is already collapsed;
			//cout<<" the couple is already collapsed, slave master list size should be non-zero, size: "<<NRSolver->slaveMasterList.size()<<endl;
			//cout<<" is there binding? "<<NRSolver->boundNodesWithSlaveMasterDefinition<<endl;
			continue;
		}
		bool  nodesHaveLateralOwners = areNodesToCollapseOnLateralECM(slaveNodeId,masterNodeId);
		if(!nodesHaveLateralOwners){
			cout<<"the nodes do not have lateral owners, continue to collapse"<<endl;
			Nodes[slaveNodeId]->collapseOnNode(Nodes, masterNodeId);
		}
		//cout<<" arranging positions"<<endl;
		for(int i=0; i<dim; ++i){
			if (Nodes[masterNodeId]->FixedPos[i]){
				Nodes[slaveNodeId]->FixedPos[i]=true;
			}
			//not using an else, as the slave could be fixed in given dimension independent of the master
			if (!Nodes[slaveNodeId]->FixedPos[i]){
				int dofmaster = masterNodeId*dim+i;
				int dofslave  = slaveNodeId*dim+i;
				//if slave is already slave of another node
				//make the master of the slave the new slave
				//   algorithm will take care of the rest.
				if (Nodes[slaveNodeId]->slaveTo[i] > -1){
					//Current slave has a master. If this master is the current master, or the master of current maste, than dont do anything, all is fine:
					if(Nodes[slaveNodeId]->slaveTo[i] == masterNodeId || Nodes[slaveNodeId]->slaveTo[i] == Nodes[masterNodeId]->slaveTo[i] ){
						cout<<"slave "<<slaveNodeId<<" is already bound to master "<<masterNodeId<<endl;
						continue;
					}
					dofslave = Nodes[slaveNodeId]->slaveTo[i]*dim+i;
					slaveNodeId = Nodes[slaveNodeId]->slaveTo[i];
				}
				//cout<<"DOF not fixed,  initial              : "<<dofmaster<<" "<<dofslave<<endl;
				//check if the master dof is already bound to something:
				NRSolver->checkMasterUpdate(dofmaster,masterNodeId);
				//It may have been that the slave was the master of the master node.
				//Now I have updated the master to the slave, and they are equal.
				//I do not need to add anything, as the master-slave relation is already implemented.
				//cout<<"DOF not fixed, after master update   : "<<dofmaster<<" "<<dofslave<<endl;
				if (dofmaster != dofslave){
					bool continueAddition =  NRSolver->checkIfCombinationExists(dofslave,dofmaster);
					//cout<<"DOF not fixed, continueAddition? : "<<continueAddition<<endl;
					if (continueAddition){
						bool madeChange = NRSolver->checkIfSlaveIsAlreadyMasterOfOthers(dofslave,dofmaster);
						if (madeChange){
							cout<<"slave "<<slaveNodeId<<" is already master of others"<<endl;
							for (int nodeIt = 0 ; nodeIt<Nodes.size(); ++nodeIt){
								if(Nodes[nodeIt]->slaveTo[i]==slaveNodeId){
									Nodes[nodeIt]->slaveTo[i]=masterNodeId;
								}
							}
						}
						vector <int> fixDOF;
						fixDOF.push_back(dofslave);
						fixDOF.push_back(dofmaster);
						NRSolver->slaveMasterList.push_back(fixDOF);
						Nodes[slaveNodeId]->slaveTo[i] = masterNodeId;
						Nodes[masterNodeId]->isMaster[i] = true;
						thereIsBinding = true;
					}
				}
			}
		}
	}
	if(thereIsBinding){
		updateElementPositions();
	}
	return thereIsBinding;
}

bool Simulation::bindCircumferenceXY(){
	bool thereIsBinding = false;
	int dim = 3;
		vector <int> nodeIds;
		for (int i=0; i<nNodes; ++i){
			if (Nodes[i]->tissuePlacement == 0){
				//basal node
				if (Nodes[i]->atCircumference){
					//at basal circumference:
					nodeIds.push_back(i);
				}
			}
		}
		int n = nodeIds.size();
		for (int i=0; i<n; ++i){
			//cout<<" in circumference binding, checking nodeId: "<<nodeIds[i]<<" ";
			int masterNodeId = nodeIds[i];
			int dofXmaster = masterNodeId*dim;
			int dofYmaster = masterNodeId*dim+1;
			int currNodeId = masterNodeId;
			bool reachedTop = false;
			while (!reachedTop){
					//Nodes[currNodeId]->tissuePlacement != 1 || (thereIsPeripodialMembrane && )){ //while node is not apical
				int nConnectedElements = Nodes[currNodeId]->connectedElementIds.size();
				for (int j=0;j<nConnectedElements;++j){
					int elementId = Nodes[currNodeId]->connectedElementIds[j];
					bool IsBasalOwner = Elements[elementId]->IsThisNodeMyBasal(currNodeId);
					if (IsBasalOwner){
						int slaveNodeId = Elements[elementId]->getCorrecpondingApical(currNodeId);
						if (Nodes[masterNodeId]->FixedPos[0]){
							//master node is fixed in X, slave needs to be fixed as well, and no need for binding calculations.
							Nodes[slaveNodeId]->FixedPos[0]=true;
						}
						if (Nodes[masterNodeId]->FixedPos[1]){
							Nodes[slaveNodeId]->FixedPos[1]=true;
						}
						//If the DoF is fixed rigidly by node fixing functions, then I
						//will not make it a slave of any node. Making it a slave leads to displacement
						//I would need to check for fixed nodes and correct my forces/jacobian twice.
						//First fix jacobian and forces, then bind nodes, then fix jacobian again to
						//clear up all the fixed nodes that were slaves and bound to other nodes.
						//If I do the binding first,without the first fixing I will carry the load of the fixed node onto the master node unnecessarily
						//Simpler and cleaner to not make them slaves at all.
						if (!Nodes[slaveNodeId]->FixedPos[0]){
							int dofXslave  = slaveNodeId*dim;
							vector <int> fixX;
							fixX.push_back(dofXslave);
							fixX.push_back(dofXmaster);
							NRSolver->slaveMasterList.push_back(fixX);
							Nodes[slaveNodeId]->slaveTo[0] = masterNodeId;
							Nodes[masterNodeId]->isMaster[0] = true;
							thereIsBinding = true;
						}
						if (!Nodes[slaveNodeId]->FixedPos[1]){
							int dofYslave  = slaveNodeId*dim+1;
							vector <int> fixY;
							fixY.push_back(dofYslave);
							fixY.push_back(dofYmaster);
							NRSolver->slaveMasterList.push_back(fixY);
							Nodes[slaveNodeId]->slaveTo[1] = masterNodeId;
							Nodes[masterNodeId]->isMaster[1] = true;
							thereIsBinding = true;
						}
						currNodeId = slaveNodeId;
						//cout<<" found slave: "<<slaveNodeId<<endl;
						break;
					}
				}
				if (!thereIsPeripodialMembrane && Nodes[currNodeId]->tissuePlacement == 1){
					//there is no peripodial and I have reached apical nodes
					reachedTop = true;
				}
				if (thereIsPeripodialMembrane && Nodes[currNodeId]->tissuePlacement == 0){
					//there is peripodial and I have reached basal nodes again
					reachedTop = true;
				}
				//cout<<" tissue placement of next node: "<<Nodes[currNodeId]->tissuePlacement<<endl;
			}
		}
		return thereIsBinding;
}

void Simulation::setLateralElementsRemodellingPlaneRotationMatrices(){
	calculateSystemCentre();
	double cx = SystemCentre[0];
	double cy = SystemCentre[1];
	if (symmetricX){
		cx = 0;
	}
	if (symmetricY){
		cy = 0;
	}
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->tissueType == 2){
			(*itElement)->setLateralElementsRemodellingPlaneRotationMatrix(cx,cy);
		}
	}
}

void Simulation::assignNodeMasses(){
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
	    (*itElement)->assignVolumesToNodes(Nodes);
		//(*itElement)->assignSurfaceAreaToNodes(Nodes);
	}
}

void Simulation::assignElementalSurfaceAreaIndices(){
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->assignExposedSurfaceAreaIndices(Nodes);
	}
}

void Simulation::assignConnectedElementsAndWeightsToNodes(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->assignElementToConnectedNodes(Nodes);
	}
}

bool Simulation::openFiles(){
	bool Success;
	if (saveData){
		//Mesh information at each step:
		string saveFileString = saveDirectory +"/Save_Frame";
		const char* name_saveFileMesh = saveFileString.c_str();
		saveFileMesh.open(name_saveFileMesh, ofstream::out);
		if (saveFileMesh.good() && saveFileMesh.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileMesh<<endl;
			Success = false;
		}

		saveFileString = saveDirectory +"Save_Summary";
		const char* name_saveFileSimulationSummary = saveFileString.c_str();
		saveFileSimulationSummary.open(name_saveFileSimulationSummary, ofstream::out);
		if (saveFileSimulationSummary.good() && saveFileSimulationSummary.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileSimulationSummary<<endl;
			Success = false;
		}

		//tension compression information at each step:
		saveFileString = saveDirectory +"/Save_TensionCompression";
		const char* name_saveFileTenComp = saveFileString.c_str();
		cout<<"opening the file" <<name_saveFileTenComp<<endl;
		saveFileTensionCompression.open(name_saveFileTenComp, ofstream::binary);
		if (saveFileTensionCompression.good() && saveFileTensionCompression.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileTenComp<<endl;
			Success = false;
		}

        //Growth information at each step (Fg):
        saveFileString = saveDirectory +"/Save_Growth";
        const char* name_saveFileGrowth = saveFileString.c_str();
        cout<<"opening the file" <<name_saveFileGrowth<<endl;
        saveFileGrowth.open(name_saveFileGrowth, ofstream::binary);
        if (saveFileGrowth.good() && saveFileGrowth.is_open()){
            Success = true;
        }
        else{
            cerr<<"could not open file: "<<name_saveFileGrowth<<endl;
            Success = false;
        }

        //Growth rate information at each step (rx,ry,rz, (as in exp(rx*dt) )for display purposes only):
        saveFileString = saveDirectory +"/Save_GrowthRate";
        const char* name_saveFileGrowthRate = saveFileString.c_str();
        cout<<"opening the file" <<name_saveFileGrowthRate<<endl;
        saveFileGrowthRate.open(name_saveFileGrowthRate, ofstream::binary);
        if (saveFileGrowthRate.good() && saveFileGrowthRate.is_open()){
            Success = true;
        }
        else{
            cerr<<"could not open file: "<<name_saveFileGrowthRate<<endl;
            Success = false;
        }

		//Force information at each step:
		saveFileString = saveDirectory +"/Save_Force";
		const char* name_saveFileForces = saveFileString.c_str();
		saveFileForces.open(name_saveFileForces, ofstream::binary);
		if (saveFileForces.good() && saveFileForces.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileForces<<endl;
			Success = false;
		}
		//opening the packing information file:
		saveFileString = saveDirectory +"/Save_Packing";
		const char* name_saveFilePacking = saveFileString.c_str();
		saveFilePacking.open(name_saveFilePacking, ofstream::binary);
		if (saveFilePacking.good() && saveFilePacking.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFilePacking<<endl;
			Success = false;
		}
		//opening the protein information file:
		saveFileString = saveDirectory +"/Save_Proteins";
		const char* name_saveFileProteins = saveFileString.c_str();
		saveFileProteins.open(name_saveFileProteins, ofstream::binary);
		if (saveFileProteins.good() && saveFileProteins.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileProteins<<endl;
			Success = false;
		}
		//opening the physical property information file:
		saveFileString = saveDirectory +"/Save_PhysicalProp";
		const char* name_saveFilePhysicalProp = saveFileString.c_str();
		saveFilePhysicalProp.open(name_saveFilePhysicalProp, ofstream::binary);
		if (saveFilePhysicalProp.good() && saveFilePhysicalProp.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFilePhysicalProp<<endl;
			Success = false;
		}

		//opening the specific Node/Element Type file:
		saveFileString = saveDirectory +"/Save_SpecificElementAndNodeTypes";
		const char* name_saveFileSpecificType = saveFileString.c_str();
		saveFileSpecificType.open(name_saveFileSpecificType, ofstream::binary);
		if (saveFileSpecificType.good() && saveFileSpecificType.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileSpecificType<<endl;
			Success = false;
		}



		//opeining the growth redistribution information file:
		saveFileString = saveDirectory +"/Save_GrowthRedistribution";
		const char* name_saveFileGrowRedist = saveFileString.c_str();
		cout<<"opening the file" <<name_saveFileGrowRedist<<endl;
		saveFileGrowthRedistribution.open(name_saveFileGrowRedist, ofstream::binary);
		if (saveFileGrowthRedistribution.good() && saveFileGrowthRedistribution.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileGrowRedist<<endl;
			Success = false;
		}

		//node binding information at each step:
		saveFileString = saveDirectory +"/Save_NodeBinding";
		const char* name_saveFileNodeBind = saveFileString.c_str();
		cout<<"opening the file" <<name_saveFileNodeBind<<endl;
		//saveFileNodeBinding.open(name_saveFileNodeBind, ofstream::binary);
		saveFileNodeBinding.open(name_saveFileNodeBind, ofstream::out);
		if (saveFileNodeBinding.good() && saveFileNodeBinding.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileNodeBind<<endl;
			Success = false;
		}
		//node collapse and adhesion information at each step:
		saveFileString = saveDirectory +"/Save_CollapseAndAdhesion";
		const char* name_saveFileCollapseAndAdhesion = saveFileString.c_str();
		cout<<"opening the file" <<name_saveFileCollapseAndAdhesion<<endl;
		saveFileCollapseAndAdhesion.open(name_saveFileCollapseAndAdhesion, ofstream::binary);
		if (saveFileCollapseAndAdhesion.good() && saveFileCollapseAndAdhesion.is_open()){
			Success = true;
		}
		else{
			cerr<<"could not open file: "<<name_saveFileCollapseAndAdhesion<<endl;
			Success = false;
		}
	}
	if (saveDirectory == "Not-Set"){
		cerr<<"Output directory is not set, outputting on Out file in current directory"<<endl;
		outputFileString = "./Out";
	}
	else {
		outputFileString = saveDirectory +"/Out";
	}
	const char* name_outputFile = outputFileString.c_str();
	outputFile.open(name_outputFile, ofstream::out);
	if (outputFile.good() && outputFile.is_open()){
		Success = true;
	}
	else{
		cerr<<"could not open file: "<<name_outputFile<<endl;
		Success = false;
	}
	return Success;
}

bool Simulation::reOpenOutputFile(){
	outputFile.close();
	const char* name_outputFile = outputFileString.c_str();
	outputFile.open(name_outputFile, ofstream::out);
	if (!(outputFile.good() && outputFile.is_open())){
		cerr<<"at step: "<<currSimTimeSec<<" could not open file: "<<name_outputFile<<endl;
		return false;
	}
	return true;
}

void Simulation::writeSimulationSummary(){
	saveFileSimulationSummary<<"TimeStep(sec):  ";
	saveFileSimulationSummary<<dt<<endl;
	saveFileSimulationSummary<<"DataSaveInterval(sec):  ";
	saveFileSimulationSummary<<dataSaveInterval*dt<<endl;
	saveFileSimulationSummary<<"ModelinputName:  ";
	saveFileSimulationSummary<<ModInp->parameterFileName<<endl<<endl;
	saveFileSimulationSummary<<"Mesh_Type:  ";
	saveFileSimulationSummary<<MeshType<<endl;
	saveFileSimulationSummary<<"	Symmetricity-x: "<<symmetricX<<" Symmetricity-y: "<<symmetricY<<endl;
	writeMeshFileSummary();
	writeGrowthRatesSummary();
	writeMyosinSummary();
	writeECMSummary();
	writeActinSummary();
	writeExperimentalSummary();
}

void Simulation::writeSpecificNodeTypes(){
	//WCounting the number of elements and nodes for, actin elements, ECM mimicing elements
	//Elements marked by ellipses and nodes marked by ellipses.
	int counterForActinMimicingElements = 0;
	int counterForECMMimicingElements = 0;
	int counterForMarkerEllipsesOnElements = 0;
	int counterForMarkerEllipsesOnNodes = 0;
	int numberOfCounterEllipses = 0;

	for (vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->isActinMimicing){
			counterForActinMimicingElements++;
		}
		if ((*itElement)->isECMMimicing){
			counterForECMMimicingElements++;
		}
		if ((*itElement)->insideEllipseBand){
			counterForMarkerEllipsesOnElements++;
			if (numberOfCounterEllipses<(*itElement)->coveringEllipseBandId){
				numberOfCounterEllipses=(*itElement)->coveringEllipseBandId;
			}
		}
	}
	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->insideEllipseBand){
			counterForMarkerEllipsesOnNodes++;
		}
	}
	//Writing explicit actin layer:
	saveFileSpecificType.write((char*) &counterForActinMimicingElements, sizeof counterForActinMimicingElements);
	for (vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->isActinMimicing){
			saveFileSpecificType.write((char*) &(*itElement)->Id, sizeof (*itElement)->Id);
		}
	}

	//Writing explicit ECM layer:
	saveFileSpecificType.write((char*) &counterForECMMimicingElements, sizeof counterForECMMimicingElements);
	for (vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->isECMMimicing){
			saveFileSpecificType.write((char*) &(*itElement)->Id, sizeof (*itElement)->Id);
		}
	}

	//Writing marker ellipses for elements:
	saveFileSpecificType.write((char*) &counterForMarkerEllipsesOnElements, sizeof counterForMarkerEllipsesOnElements);
	for (vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->insideEllipseBand){
			saveFileSpecificType.write((char*) &(*itElement)->Id, sizeof (*itElement)->Id);
			saveFileSpecificType.write((char*) &(*itElement)->coveringEllipseBandId, sizeof (*itElement)->coveringEllipseBandId);
		}
	}

	//Writing marker ellipses for nodes:
	saveFileSpecificType.write((char*) &counterForMarkerEllipsesOnNodes, sizeof counterForMarkerEllipsesOnNodes);
	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->insideEllipseBand){
			saveFileSpecificType.write((char*) &(*itNode)->Id, sizeof (*itNode)->Id);
			saveFileSpecificType.write((char*) &(*itNode)->coveringEllipseBandId, sizeof (*itNode)->coveringEllipseBandId);
		}
	}
	cout<<"wrote specific element types: "<<counterForMarkerEllipsesOnNodes<<" "<<counterForMarkerEllipsesOnElements<<" "<<counterForECMMimicingElements<<" "<<counterForActinMimicingElements<<endl;
}

void Simulation::writeMeshFileSummary(){
	if ( MeshType == 2){
		saveFileSimulationSummary<<"	Row: ";
		saveFileSimulationSummary<<Row;
		saveFileSimulationSummary<<"	Column: ";
		saveFileSimulationSummary<<Column;
		saveFileSimulationSummary<<"	SideLength: ";
		saveFileSimulationSummary<<SideLength;
		saveFileSimulationSummary<<"	zHeight: ";
		saveFileSimulationSummary<<zHeight;
		saveFileSimulationSummary<<endl;
	}
	if ( MeshType == 4){
		saveFileSimulationSummary<<"	MeshFileName: ";
		saveFileSimulationSummary<<inputMeshFileName;
		saveFileSimulationSummary<<endl;
	}
}

void Simulation::writeGrowthRatesSummary(){
	saveFileSimulationSummary<<endl;
	for (int i=0; i<nGrowthFunctions; ++i){
		if(GrowthFunctions[i]->Type == 3){ //grid based function does not need dt in summary
			GrowthFunctions[i]->writeSummary(saveFileSimulationSummary);
		}
		else{
			GrowthFunctions[i]->writeSummary(saveFileSimulationSummary,dt);
		}
	}
}

void Simulation::writeMyosinSummary(){
	saveFileSimulationSummary<<endl;
	for (int i=0; i<nMyosinFunctions; ++i){
		myosinFunctions[i]->writeSummary(saveFileSimulationSummary);
	}
}


void Simulation::writeECMSummary(){
	saveFileSimulationSummary<<endl;
	saveFileSimulationSummary<<"There is explicit ECM:	"<<thereIsExplicitECM<<endl;
	if (thereIsExplicitECM){
		writeECMProperties();
	}
}

void Simulation::writeECMProperties(){
	saveFileSimulationSummary<<"	Columnar ECM Stiffness(Pa): "<<EColumnarECM<<endl;
	saveFileSimulationSummary<<"	Peripodial ECM Stiffness(Pa): "<<EPeripodialECM<<endl;
	saveFileSimulationSummary<<"	ECM remodelling half life (hour): "<<ECMRenawalHalfLife/3600.0<<endl;
	saveFileSimulationSummary<<"	Is there perturbation to the ECM: "<<thereIsECMChange<<endl;
	if (thereIsECMChange){
		int n = ECMChangeBeginTimeInSec.size();
		saveFileSimulationSummary<<"there are "<<n<<" ECM perturbations"<<endl;
		for (int ECMperturbationIndex =0;ECMperturbationIndex<n; ++ECMperturbationIndex ){
			saveFileSimulationSummary<<"		stiffness alteration begins at  "<<ECMChangeBeginTimeInSec[ECMperturbationIndex]/3600.0<<" hrs and ends at "<<ECMChangeEndTimeInSec[ECMperturbationIndex]/3600.0<<endl;
			saveFileSimulationSummary<<"		final fraction of ECM stiffness  "<<ECMStiffnessChangeFraction[ECMperturbationIndex]<<" times original values."<<endl;
			saveFileSimulationSummary<<"		final fraction of ECM renewal time  "<<ECMRenewalHalfLifeTargetFraction[ECMperturbationIndex]<<" times original values."<<endl;
			saveFileSimulationSummary<<"		final fraction of ECM viscosity  "<<ECMViscosityChangeFraction[ECMperturbationIndex]<<" times original values."<<endl;
			saveFileSimulationSummary<<"		stiffness alteration applied to apical ECM: "	<<changeApicalECM[ECMperturbationIndex]<<endl;
			saveFileSimulationSummary<<"		stiffness alteration applied to basal  ECM: "	<<changeBasalECM[ECMperturbationIndex]<<endl;
			saveFileSimulationSummary<<"		stiffness alteration applied to ellipse bands: ";
			for (int i=0; i<numberOfECMChangeEllipseBands[ECMperturbationIndex]; ++i){
				saveFileSimulationSummary<<" "<<ECMChangeEllipseBandIds[ECMperturbationIndex][i];
			}
			saveFileSimulationSummary<<endl;
		}
	}
}

void Simulation::writeActinSummary(){
	saveFileSimulationSummary<<endl;
	saveFileSimulationSummary<<"There is explicit actin:  "<<thereIsExplicitActin<<endl;
}

void Simulation::writeExperimentalSummary(){
	writePipetteSumary();
}

void Simulation::writePipetteSumary(){
	saveFileSimulationSummary<<endl;
	saveFileSimulationSummary<<"There is pipette aspiration:  "<<PipetteSuction<<endl;
	if (PipetteSuction){
		saveFileSimulationSummary<<"	is the tissue stuck on the glass: "<<TissueStuckOnGlassDuringPipetteAspiration<<endl;
		saveFileSimulationSummary<<"	suction on apical side: "<<ApicalSuction<<endl;
		saveFileSimulationSummary<<"	pipette position [centre(x,y,z), effective pippette suction depth, microns]: "<<pipetteCentre[0]<<" "<<pipetteCentre[1]<<" "<<pipetteCentre[2]<<" "<<pipetteDepth<<endl;
		saveFileSimulationSummary<<"	pipette radia [inner & outer(microns)]: "<<pipetteInnerRadius<<" "<<pipetteInnerRadius+pipetteThickness<<endl;
		saveFileSimulationSummary<<"	Steps for pipette suction "<<nPipetteSuctionSteps<<endl;
		saveFileSimulationSummary<<"	suction times (sec): ";
		for (int i=0; i<nPipetteSuctionSteps; ++i){
			saveFileSimulationSummary<<pipetteSuctionTimes[i]<<" ";
		}
		saveFileSimulationSummary<<endl;
		saveFileSimulationSummary<<"	suction pressures (Pa): ";
		for (int i=0; i<nPipetteSuctionSteps; ++i){
			saveFileSimulationSummary<<pipetteSuctionPressures[i]<<" ";
		}
	}

}

void Simulation::writeRelaxedMeshFromCurrentState(){
	string meshSaveString = saveDirectory +"/MeshFromEndPoint.mesh";
	const char* name_meshSaveString = meshSaveString.c_str();;
	ofstream file;
	file.open(name_meshSaveString, ofstream::out);
	//count non ablated nodes:
	file<<nNodes;
	file<<endl;
	for (int i=0; i<nNodes; ++i){
		file << Nodes[i]->Position[0];
		file<<" \t";
		file << Nodes[i]->Position[1];
		file<<" \t";
		file << Nodes[i]->Position[2];
		file<<" \t";
		file << Nodes[i]->tissuePlacement;
		file<<" \t";
		file << Nodes[i]->tissueType;
		file<<" \t";
		file << Nodes[i]->atCircumference;
		file << endl;
	}
	file<<nElements;
	file<<endl;
	for (int i=0; i<nElements; ++i){
		file<< Elements[i]->getShapeType();
		file<<" \t";
		const int n = Elements[i]->getNodeNumber();
		int* NodeIds;
		NodeIds = new int[n];
		NodeIds = Elements[i]->getNodeIds();
		for (int j=0; j<n; ++j){
			file<< NodeIds[j];
			file<<" \t";
		}
		for (int j=0; j<n; ++j){
			file << Nodes[NodeIds[j]]->Position[0];
			file<<" \t";
			file << Nodes[NodeIds[j]]->Position[1];
			file<<" \t";
			file << Nodes[NodeIds[j]]->Position[2];
			file<<" \t";
		}
		file << endl;
	}
	file.close();
}

bool Simulation::openFilesToDisplay(){
	string saveFileString = saveDirectoryToDisplayString +"/Save_Frame";
	const char* name_saveFileToDisplayMesh = saveFileString.c_str();
	saveFileToDisplayMesh.open(name_saveFileToDisplayMesh, ifstream::in);
	if (!(saveFileToDisplayMesh.good() && saveFileToDisplayMesh.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayMesh<<endl;
		return false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_Summary";
	const char* name_saveFileToDisplaySimSum = saveFileString.c_str();;
	saveFileToDisplaySimSum.open(name_saveFileToDisplaySimSum, ifstream::in);
	if (!(saveFileToDisplaySimSum.good() && saveFileToDisplaySimSum.is_open())){
		cerr<<"Cannot open the simulation summary: "<<name_saveFileToDisplaySimSum<<endl;
		return false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_SpecificElementAndNodeTypes";
	const char* name_saveSpecificElementAndNodeTypes = saveFileString.c_str();;
	saveFileToDisplaySpecificNodeTypes.open(name_saveSpecificElementAndNodeTypes, ifstream::in);
	specificElementTypesRecorded = true;
	if (!(saveFileToDisplaySpecificNodeTypes.good() && saveFileToDisplaySpecificNodeTypes.is_open())){
		cerr<<"Cannot open the specific node types: "<<name_saveSpecificElementAndNodeTypes<<endl;
		specificElementTypesRecorded = false;
		//return false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_TensionCompression";
	const char* name_saveFileToDisplayTenComp = saveFileString.c_str();
	saveFileToDisplayTenComp.open(name_saveFileToDisplayTenComp, ifstream::in);
	if (!(saveFileToDisplayTenComp.good() && saveFileToDisplayTenComp.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayTenComp<<endl;
		TensionCompressionSaved = false;
	}

    saveFileString = saveDirectoryToDisplayString +"/Save_Growth";
    const char* name_saveFileToDisplayGrowth = saveFileString.c_str();
    saveFileToDisplayGrowth.open(name_saveFileToDisplayGrowth, ifstream::in);
    if (!(saveFileToDisplayGrowth.good() && saveFileToDisplayGrowth.is_open())){
        cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayGrowth<<endl;
        GrowthSaved = false;
    }

    saveFileString = saveDirectoryToDisplayString +"/Save_GrowthRate";
	const char* name_saveFileToDisplayGrowthRate = saveFileString.c_str();
	saveFileToDisplayGrowthRate.open(name_saveFileToDisplayGrowthRate, ifstream::in);
	if (!(saveFileToDisplayGrowthRate.good() && saveFileToDisplayGrowthRate.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayGrowthRate<<endl;
		GrowthRateSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_Force";
	const char* name_saveFileToDisplayForce = saveFileString.c_str();;
	saveFileToDisplayForce.open(name_saveFileToDisplayForce, ifstream::in);
	if (!(saveFileToDisplayForce.good() && saveFileToDisplayForce.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayForce<<endl;
		ForcesSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_Proteins";
	const char* name_saveFileToDisplayProteins = saveFileString.c_str();;
	saveFileToDisplayProteins.open(name_saveFileToDisplayProteins, ifstream::in);
	if (!(saveFileToDisplayProteins.good() && saveFileToDisplayProteins.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayProteins<<endl;
		ProteinsSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_Packing";
	const char* name_saveFileToDisplayPacking = saveFileString.c_str();;
	saveFileToDisplayPacking.open(name_saveFileToDisplayPacking, ifstream::in);
	if (!(saveFileToDisplayPacking.good() && saveFileToDisplayPacking.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayPacking<<endl;
		PackingSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_PhysicalProp";
	const char* name_saveFileToDisplayPhysicalProp = saveFileString.c_str();;
	saveFileToDisplayPhysicalProp.open(name_saveFileToDisplayPhysicalProp, ifstream::in);
	if (!(saveFileToDisplayPhysicalProp.good() && saveFileToDisplayPhysicalProp.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayPhysicalProp<<endl;
		physicalPropertiesSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_GrowthRedistribution";
	const char* name_saveFileToDisplayGrowRedist = saveFileString.c_str();
	saveFileToDisplayGrowthRedistribution.open(name_saveFileToDisplayGrowRedist, ifstream::in);
	if (!(saveFileToDisplayTenComp.good() && saveFileToDisplayTenComp.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayGrowRedist<<endl;
		growthRedistributionSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_NodeBinding";
	const char* name_saveFileToDisplayNodeBinding = saveFileString.c_str();
	saveFileToDisplayNodeBinding.open(name_saveFileToDisplayNodeBinding, ifstream::in);
	if (!(saveFileToDisplayNodeBinding.good() && saveFileToDisplayNodeBinding.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayNodeBinding<<endl;
		nodeBindingSaved = false;
	}
	saveFileString = saveDirectoryToDisplayString +"/Save_CollapseAndAdhesion";
	const char* name_saveFileToDisplayCollapseAndAdhesion = saveFileString.c_str();
	saveFileToDisplayCollapseAndAdhesion.open(name_saveFileToDisplayCollapseAndAdhesion, ifstream::in);
	if (!(saveFileToDisplayCollapseAndAdhesion.good() && saveFileToDisplayCollapseAndAdhesion.is_open())){
		cerr<<"Cannot open the save file to display: "<<name_saveFileToDisplayCollapseAndAdhesion<<endl;
		collapseAndAdhesionSaved = false;
	}
	return true;
}

bool Simulation::initiateSavedSystem(){
	bool success  = openFilesToDisplay();
	if (!success){
		return false;
	}
	//reading system properties:
	success  = readSystemSummaryFromSave();
	if (!success){
		return false;
	}
	string currline;
	//skipping the header:
	getline(saveFileToDisplayMesh,currline);
	initiateNodesFromSave();
	initiateElementsFromSave();
	fillInNodeNeighbourhood();
	assignPhysicalParameters();
	initiateSystemForces();
	if (specificElementTypesRecorded){
		success  = readSpecificNodeTypesFromSave();
	}
	if (!success){
		return false;
	}
	cout<<" reading Ten comp"<<endl;
	if (TensionCompressionSaved){
		updateTensionCompressionFromSave();
	}
	cout<<" reading growth"<<endl;

    if (GrowthSaved){
        updateGrowthFromSave();
    }
	cout<<" reading growth rate"<<endl;

    if (GrowthRateSaved){
        updateGrowthRateFromSave();
    }
	cout<<" reading forces"<<endl;

	if (ForcesSaved){
		updateForcesFromSave();
	}
	if (growthRedistributionSaved){
		updateGrowthRedistributionFromSave();
	}
	if (nodeBindingSaved){
		updateNodeBindingFromSave();
	}
	if(collapseAndAdhesionSaved){
		updateCollapseAndAdhesionFromSave();
	}
	updateElementVolumesAndTissuePlacements();
    //cleanMatrixUpdateData();
	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	clearLaserAblatedSites();
    //calculateStiffnessMatrices();
    calculateShapeFunctionDerivatives();
	updateElementPositions();
	//skipping the footer:
	getline(saveFileToDisplayMesh,currline);
	while (currline.empty() && !saveFileToDisplayMesh.eof()){
		//skipping empty line
		getline(saveFileToDisplayMesh,currline);
	}
	if(saveFileToDisplayMesh.eof()){
		reachedEndOfSaveFile = true;
		return true;
	}

	//cout<<"skipped footer: "<<currline<<endl;
	return true;
}

bool Simulation::readSystemSummaryFromSave(){
	string dummystring;
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: TimeStep(sec):"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "TimeStep(sec):"
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: dt value"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dt;

	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: DataSaveInterval(sec):"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "DataSaveInterval(sec):"
	double dummydouble;
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: save interval value"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummydouble;
	dataSaveInterval =  (int) ceil(dummydouble/dt);


	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: ModelinputName:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "ModelinputName: "
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting name of the model input file"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading the model input file
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: Mesh_Type:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "Mesh_Type: "
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: the mesh type(int)"<<endl;
		return false;
	}
	int dummyint;
	saveFileToDisplaySimSum >> dummyint; //reading the mesh type
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: Symmetricity-x:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "Symmetricity-x: "
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: symmetricitX boolean:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> symmetricX;
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: Symmetricity-y:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> dummystring; //reading "Symmetricity-y: "
	if(saveFileToDisplaySimSum.eof()){
		cerr<<"reached the end of summary file, expecting: symmetricitY boolean:"<<endl;
		return false;
	}
	saveFileToDisplaySimSum >> symmetricY;
	cout<<"read the summary, symmetricity data: "<<symmetricX<<" "<<symmetricY<<endl;
	return true;
}

bool Simulation::readSpecificNodeTypesFromSave(){
	int currElementId;
	//read actin layer to display:
	int counterForActinMimicingElements;
	saveFileToDisplaySpecificNodeTypes.read((char*) &counterForActinMimicingElements, sizeof counterForActinMimicingElements);
	if (counterForActinMimicingElements>0){
		thereIsExplicitActin = true;
	}
	for (int i=0; i<counterForActinMimicingElements; i++){
		saveFileToDisplaySpecificNodeTypes.read((char*) &currElementId, sizeof currElementId);
		int currIndice = getElementArrayIndexFromId(currElementId);
		if (currIndice < 0){
			cout<<" error in reading actin mimicking elements "<<endl;
			cerr<<" error in reading actin mimicking elements "<<endl;
			return false;
		}
		Elements[currIndice]->isActinMimicing = true;
	}
	//read in ECM layer to display:
	int counterForECMMimicingElements;
	saveFileToDisplaySpecificNodeTypes.read((char*) &counterForECMMimicingElements, sizeof counterForECMMimicingElements);
	if (counterForECMMimicingElements>0){
		thereIsExplicitECM = true;
	}
	for (int i=0; i<counterForECMMimicingElements; i++){
		saveFileToDisplaySpecificNodeTypes.read((char*) &currElementId, sizeof currElementId);
		int currIndice = getElementArrayIndexFromId(currElementId);
		if (currIndice < 0){
			cout<<" error in reading ECM mimicking elements "<<endl;
			cerr<<" error in reading ECM mimicking elements "<<endl;
			return false;
		}
		Elements[currIndice]->isECMMimicing = true;
		//cout<<" ECM indice: "<<currIndice<<" "<<currElementId<<endl;
	}
	assigneElementsAtTheBorderOfECM();
	assigneElementsAtTheBorderOfActin();
	//read marker ellipses to display for elements:
	int counterForMarkerEllipsesOnElements;
	int currEllipseBandId;
	saveFileToDisplaySpecificNodeTypes.read((char*) &counterForMarkerEllipsesOnElements, sizeof counterForMarkerEllipsesOnElements);
	cout<<"counterForMarkerEllipsesOnElements "<<counterForMarkerEllipsesOnElements<<endl;
	if (counterForMarkerEllipsesOnElements<0 || counterForMarkerEllipsesOnElements> 200000){
		counterForMarkerEllipsesOnElements = 0;
	}
	for (int i=0; i<counterForMarkerEllipsesOnElements; i++){
		saveFileToDisplaySpecificNodeTypes.read((char*) &currElementId, sizeof currElementId);
		saveFileToDisplaySpecificNodeTypes.read((char*) &currEllipseBandId, sizeof currEllipseBandId);
		int currIndice = getElementArrayIndexFromId(currElementId);
		if (currIndice < 0){
			cout<<" error in reading marker ellipses for elements "<<endl;
			cerr<<" error in reading marker ellipses for elements "<<endl;
			return false;
		}
		Elements[currIndice]->insideEllipseBand = true;
		Elements[currIndice]->coveringEllipseBandId = currEllipseBandId;
	}
	//read marker ellipses to display for nodes:
	int counterForMarkerEllipsesOnNodes;
	int currNodeId;
	saveFileToDisplaySpecificNodeTypes.read((char*) &counterForMarkerEllipsesOnNodes, sizeof counterForMarkerEllipsesOnNodes);

	cout<<" to display counterForMarkerEllipsesOnNodes "<<counterForMarkerEllipsesOnNodes<<endl;
	if (counterForMarkerEllipsesOnNodes<0 || counterForMarkerEllipsesOnNodes> 200000){
		counterForMarkerEllipsesOnNodes = 0;
	}
	for (int i=0; i<counterForMarkerEllipsesOnNodes; i++){
		saveFileToDisplaySpecificNodeTypes.read((char*) &currNodeId, sizeof currNodeId);
		saveFileToDisplaySpecificNodeTypes.read((char*) &currEllipseBandId, sizeof currEllipseBandId);
		int currIndice = getNodeArrayIndexFromId(currNodeId);
		if (currIndice < 0){
			cout<<" error in reading marker ellipses for nodes "<<endl;
			cerr<<" error in reading marker ellipses for nodes "<<endl;
			return false;
		}
		Nodes[currIndice]->insideEllipseBand = true;
		Nodes[currIndice]->coveringEllipseBandId = currEllipseBandId;
	}
	cout<<"read specific element types: "<<counterForMarkerEllipsesOnNodes<<" "<<counterForMarkerEllipsesOnElements<<" "<<counterForECMMimicingElements<<" "<<counterForActinMimicingElements<<endl;
	return true;
}

int Simulation::getElementArrayIndexFromId(int currId){
	int currIndice = 0;
	for (vector<ShapeBase*>::iterator itElement = Elements.begin(); itElement < Elements.end(); itElement++){
		if ((*itElement)->getId() == currId){
			return currIndice;
		}
		currIndice++;
	}
	cout<<"Error in getElementArrayIndexFromId, required element Id : "<<currId<<" returnin -10, inducing a crash"<<endl;
	cerr<<"Error in getElementArrayIndexFromId, required element Id : "<<currId<<" returnin -10, inducing a crash"<<endl;
	return -10;
}

int Simulation::getNodeArrayIndexFromId(int currId){
	int currIndice = 0;
	for (vector<Node*>::iterator itNode = Nodes.begin(); itNode < Nodes.end(); itNode++){
		if ((*itNode)->getId() == currId){
			return currIndice;
		}
		currIndice++;
	}
	cout<<"Error in getElementArrayIndexFromId, required node Id : "<<currId<<" returnin -10, inducing a crash"<<endl;
	cerr<<"Error in getElementArrayIndexFromId, required node Id : "<<currId<<" returnin -10, inducing a crash"<<endl;
	return -10;
}

bool Simulation::readNodeDataToContinueFromSave(){
	int n;
	saveFileToDisplayMesh >> n;
	cout<<"number of nodes: "<<n<<endl;
	if(nNodes != n){
		cerr<<"The node number from save file("<<n<<") and model input("<<Nodes.size()<<") are not consistent - cannot continue simulation from save"<<endl;
		return false;
	}
	for (int i=0; i<n; ++i){
		int tissuePlacement, tissueType;
		bool atCircumference;
		saveFileToDisplayMesh >> Nodes[i]->Position[0];
		saveFileToDisplayMesh >> Nodes[i]->Position[1];
		saveFileToDisplayMesh >> Nodes[i]->Position[2];
		saveFileToDisplayMesh >> tissuePlacement;
		saveFileToDisplayMesh >> tissueType;
		saveFileToDisplayMesh >> atCircumference;
		if (Nodes[i]->tissuePlacement != tissuePlacement || Nodes[i]->atCircumference != atCircumference || Nodes[i]->tissueType != tissueType ){
			cerr<<"Node "<<i<<" properties are not consistent  - cannot continue simulation from save "<<endl;
			return false;
		}
	}
	return true;
}

void Simulation::initiateNodesFromSave(){
	saveFileToDisplayMesh >> nNodes;
	cout<<"number of nodes: "<<nNodes<<endl;
	Node* tmp_nd;
	for (int i=0; i<nNodes; ++i){
		double* pos = new double[3];
		int tissuePlacement, tissueType;
		bool atCircumference;
		saveFileToDisplayMesh >> pos[0];
		saveFileToDisplayMesh >> pos[1];
		saveFileToDisplayMesh >> pos[2];
		saveFileToDisplayMesh >> tissuePlacement;
		saveFileToDisplayMesh >> tissueType;
		saveFileToDisplayMesh >> atCircumference;
		tmp_nd = new Node(i, 3, pos,tissuePlacement, tissueType);
		tmp_nd-> atCircumference = atCircumference;
		Nodes.push_back(tmp_nd);
		delete[] pos;
	}
	cout<<"number of nodes: "<<nNodes<<endl;
}

void Simulation::initiateNodesFromMeshInput(){
	int n;
	inputMeshFile >> n;
	Node* tmp_nd;
	for (int i=0; i<n; ++i){
		double* pos = new double[3];
		int tissuePos = -2;
		int tissueType = -2;
		int atCircumference;
		inputMeshFile >> pos[0];
		inputMeshFile >> pos[1];
		inputMeshFile >> pos[2];
		inputMeshFile >> tissuePos;
		inputMeshFile >> tissueType;
		inputMeshFile >> atCircumference;
		tmp_nd = new Node(i, 3, pos,tissuePos, tissueType);
		tmp_nd->atCircumference = atCircumference;
		Nodes.push_back(tmp_nd);
		nNodes = Nodes.size();
		delete[] pos;
	}
}

void Simulation::initiateElementsFromMeshInput(){
	int n;
	inputMeshFile >> n;
	for (int i=0; i<n; ++i){
		int shapeType;
		inputMeshFile >> shapeType;
		if (shapeType == 1){
			initiatePrismFromMeshInput();
		}
		else{
			cerr<<"Element "<<i<<" of "<<n<<": Error in shape type, corrupt save file! - currShapeType: "<<shapeType<<endl;
			break;
		}
	}
}


void Simulation::initiateElementsFromSave(){
	int n;
	saveFileToDisplayMesh >> n;
	cout<<"number of elements: "<<n<<endl;
	for (int i=0; i<n; ++i){
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (shapeType == 1){
			initiatePrismFromSave();
		}
		else{
			cerr<<"Error in shape type, corrupt save file! - currShapeType: "<<shapeType<<endl;
		}
	}
	cout<<"number of elements: "<<nElements<<endl;
}

bool Simulation::readElementDataToContinueFromSave(){
	int n;
	saveFileToDisplayMesh >> n;
	if (nElements != n){
		cerr<<"The element number from save file and model input are not consistent - cannot continue simulation from save"<<endl;
		cerr<<"n: "<<n<<" nElements: "<<endl;
		return false;
	}
	for (int i=0; i<nElements; ++i){
		//string line;
		//getline(saveFileToDisplayMesh, line);

		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (Elements[i]->getShapeType() != shapeType){
			cerr<<"The element type from save file and model input are not consistent - cannot continue simulation from save"<<endl;
			return false;
		}
		if (shapeType == 1){
			bool success = readShapeData(i);
			if (!success){
				cerr<<"Error reading shape data, element: "<<i<<endl;
				return false;
			}
		}
		else if (shapeType == 4){
			double height = 0.0;
			saveFileToDisplayMesh >> height;
			if (Elements[i]->ReferenceShape->height != height){
				cerr<<"The element height from save file and model input are not consistent - cannot continue simulation from save"<<endl;
				return false;
			}
			bool success = readShapeData(i);
			if (!success){
				cerr<<"Error reading shape data, element: "<<i<<endl;
				return false;
			}
		}
		else{
			cerr<<"Error in shape type, corrupt save file! - currShapeType: "<<shapeType<<endl;
		}

	}
	//string line;
	//getline(saveFileToDisplayMesh, line);

	return true;
}

void Simulation::initiatePrismFromSave(){
	//inserts a new prism at order k into elements vector
	//the node ids and reference shape positions
	//will be updated in function: updateShapeFromSave
	int* NodeIds;
	NodeIds = new int[6];
	for (int i =0 ;i<6; ++i){
		NodeIds[i] = 0;
	}
	Prism* PrismPnt01;
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
	PrismPnt01->updateShapeFromSave(saveFileToDisplayMesh);
	Elements.push_back(PrismPnt01);
	nElements = Elements.size();
	currElementId++;
	delete[] NodeIds;
	//cout<<"Element: "<<PrismPnt01->Id<<endl;
	//Elements[Elements.size()-1]->displayPositions();
	//Elements[Elements.size()-1]->displayReferencePositions();
}

bool Simulation::readShapeData(int i){
	bool IsAblated;
	saveFileToDisplayMesh >> IsAblated;
	if ( Elements[i]->IsAblated != IsAblated){
		cerr<<"The element "<<i<<" ablation from save file and model input are not consistent - cannot continue simulation from save"<<endl;
		return false;
	}
	bool success = Elements[i]->readNodeIdData(saveFileToDisplayMesh);
	if (!success){
		cerr<<"The element "<<i<<" node ids from save file and model input are not consistent - cannot continue simulation from save"<<endl;
		return false;
	}
	success = Elements[i]->readReferencePositionData(saveFileToDisplayMesh);
	if (!success){
		cerr<<"The element "<<i<<" reference shape from save file and model input are not consistent - cannot continue simulation from save"<<endl;
		return false;
	}
	return true;
}


void Simulation::initiatePrismFromMeshInput(){
	int* NodeIds;
	NodeIds = new int[6];
	for (int i =0 ;i<6; ++i){
		int savedId;
		inputMeshFile >> savedId;
		NodeIds[i] = savedId;
	}
	Prism* PrismPnt01;
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
	PrismPnt01->updateReferencePositionMatrixFromMeshInput(inputMeshFile);
	PrismPnt01->checkRotationConsistency3D();
	Elements.push_back(PrismPnt01);
	nElements = Elements.size();
	currElementId++;
	delete[] NodeIds;
	//Elements[Elements.size()-1]->displayReferencePositions();
}



void Simulation::reInitiateSystemForces(int oldSize){
	//deleting the old system forces:
	for (int j=0;j<oldSize;++j){
		delete[] SystemForces[j];
		delete[] PackingForces[j];
		//delete[] PackingForcesPreviousStep[j];
		//delete[] PackingForcesTwoStepsAgoStep[j];
		delete[] FixedNodeForces[j];
	}
	delete[] SystemForces;
	delete[] PackingForces;
	//delete[] PackingForcesPreviousStep;
	//delete[] PackingForcesTwoStepsAgoStep;
	delete[] FixedNodeForces;
	//reinitiating with the new size:
	const int n = nNodes;
	SystemForces = new double*[n];
	PackingForces = new double*[n];
	//PackingForcesPreviousStep = new double*[n];
	//PackingForcesTwoStepsAgoStep = new double*[n];
	FixedNodeForces = new double*[n];
	for (int j=0;j<n;++j){
		SystemForces[j] = new double[3];
		PackingForces[j] = new double[3];
		//PackingForcesPreviousStep[j] = new double[3];
		//PackingForcesTwoStepsAgoStep[j] = new double[3];
		FixedNodeForces[j] = new double[3];
		SystemForces[j][0]=0.0;
		SystemForces[j][1]=0.0;
		SystemForces[j][2]=0.0;
		PackingForces[j][0]=0.0;
		PackingForces[j][1]=0.0;
		PackingForces[j][2]=0.0;
		//PackingForcesPreviousStep[j][0] = 0.0;
		//PackingForcesPreviousStep[j][1] = 0.0;
		//PackingForcesPreviousStep[j][2] = 0.0;
		//PackingForcesTwoStepsAgoStep[j][0] = 0.0;
		//PackingForcesTwoStepsAgoStep[j][1] = 0.0;
		//PackingForcesTwoStepsAgoStep[j][2] = 0.0;
		FixedNodeForces[j][0] = 0.0;
		FixedNodeForces[j][1] = 0.0;
		FixedNodeForces[j][2] = 0.0;
	}
}

void Simulation::updateForcesFromSave(){
	for (int i=0;i<nNodes;++i){
		saveFileToDisplayForce.read((char*) &SystemForces[i][0], sizeof SystemForces[i][0]);
		saveFileToDisplayForce.read((char*) &SystemForces[i][1], sizeof SystemForces[i][1]);
		saveFileToDisplayForce.read((char*) &SystemForces[i][2], sizeof SystemForces[i][2]);
	}
	for (int i=0;i<nElements;++i){
		int n = Elements[i]->getNodeNumber();
		for (int j=0; j<n; ++j){
			saveFileToDisplayForce.read((char*) &Elements[i]->MyoForce[j][0], sizeof &Elements[i]->MyoForce[j][0]);
			saveFileToDisplayForce.read((char*) &Elements[i]->MyoForce[j][1], sizeof &Elements[i]->MyoForce[j][1]);
			saveFileToDisplayForce.read((char*) &Elements[i]->MyoForce[j][2], sizeof &Elements[i]->MyoForce[j][2]);
		}
	}
}

void Simulation::updateTensionCompressionFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		for (int j=0; j<6; ++j){
            double S = gsl_matrix_get((*itElement)->Strain,j,0);
            saveFileToDisplayTenComp.read((char*) &S, sizeof S);
            gsl_matrix_set((*itElement)->Strain,j,0,S);
        }
	}
}

void Simulation::updateGrowthRedistributionFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		bool thereIsDistribution = false;
		bool shrinksElement = false;
		saveFileToDisplayGrowthRedistribution.read((char*) &thereIsDistribution, sizeof thereIsDistribution);
		(*itElement)->thereIsGrowthRedistribution = thereIsDistribution;
		saveFileToDisplayGrowthRedistribution.read((char*) &shrinksElement, sizeof shrinksElement);
		(*itElement)->growthRedistributionShrinksElement = shrinksElement;
	}
}

void Simulation::updateNodeBindingFromSave(){
	/*int n = 0; //size Of Master-Slave List
	saveFileToDisplayNodeBinding.read((char*) &n, sizeof n);
	cout<<" n: "<<n<<endl;
	for (int i=0; i<n; ++i){
		int dofSlave, dofMaster;
		saveFileToDisplayNodeBinding.read((char*) &dofSlave, sizeof dofSlave);
		saveFileToDisplayNodeBinding.read((char*) &dofMaster, sizeof dofMaster);
		int dim = dofSlave % 3;
		int nodeSlave = (dofSlave - dim)/3;
		int nodeMaster = (dofMaster - dim)/3;
		Nodes[nodeSlave]->slaveTo[dim] = nodeMaster;
		Nodes[nodeMaster]->isMaster[dim] = true;
	}*/


	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		(*itNode)->slaveTo[0] = -1;
		(*itNode)->slaveTo[1] = -1;
		(*itNode)->slaveTo[2] = -1;
		(*itNode)->isMaster[0] = false;
		(*itNode)->isMaster[1] = false;
		(*itNode)->isMaster[2] = false;
		(*itNode)->attachedToPeripodial = false;
	}
	int n = 0; //size Of Master-Slave List
	saveFileToDisplayNodeBinding>>n;
	for (int i=0; i<n; ++i){
		int dofSlave, dofMaster;
		saveFileToDisplayNodeBinding>>dofSlave;
		saveFileToDisplayNodeBinding>>dofMaster;
		int dim = dofSlave % 3;
		int nodeSlave = (dofSlave - dim)/3;
		int nodeMaster = (dofMaster - dim)/3;
		Nodes[nodeSlave]->slaveTo[dim] = nodeMaster;
		Nodes[nodeMaster]->isMaster[dim] = true;
		if(Nodes[nodeSlave]->tissueType ==1 ){
			Nodes[nodeMaster]->attachedToPeripodial = true;
		}
		//cout<<"read binding: "<<i<<" of "<<n<<" slave :"<<nodeSlave<<" master  "<<nodeMaster<<endl;
	}
}

void Simulation::updateCollapseAndAdhesionFromSave(){
	/*for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		int adheredTo;
		saveFileToDisplayCollapseAndAdhesion.read((char*) &adheredTo, sizeof adheredTo);
		(*itNode)->adheredTo = adheredTo;
		int nCollapsedNodes;
		saveFileToDisplayCollapseAndAdhesion.read((char*) &nCollapsedNodes, sizeof nCollapsedNodes);
		vector<int> collapseList;
		for (int i=0; i<nCollapsedNodes; ++i){
			int nodeId;
			saveFileToDisplayCollapseAndAdhesion.read((char*) &nodeId, sizeof nodeId);
			collapseList.push_back(nodeId);
		}
		(*itNode)->collapsedWith = collapseList;
	}*/
	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		int adheredTo;
		saveFileToDisplayCollapseAndAdhesion.read((char*) &adheredTo, sizeof adheredTo);
		if ((*itNode)->adheredTo<0){
			(*itNode)->adheredTo = adheredTo;
		}
		if (savedStepwise){
		saveFileToDisplayCollapseAndAdhesion.read((char*) &(*itNode)->positionUpdateOngoing, sizeof (*itNode)->positionUpdateOngoing);
		saveFileToDisplayCollapseAndAdhesion.read((char*) &(*itNode)->positionUpdateCounter, sizeof (*itNode)->positionUpdateCounter);
		}
		int nCollapsedNodes;

		saveFileToDisplayCollapseAndAdhesion.read((char*) &nCollapsedNodes, sizeof nCollapsedNodes);
		//vector<int> collapseList;
		for (int i=0; i<nCollapsedNodes; ++i){
			int nodeId;
			saveFileToDisplayCollapseAndAdhesion.read((char*) &nodeId, sizeof nodeId);
			bool nodeOnList = binary_search((*itNode)->collapsedWith.begin(), (*itNode)->collapsedWith.end(),nodeId);
			if (!nodeOnList){
				(*itNode)->collapsedWith.push_back(nodeId);
			}
		}
	}
}

void Simulation::updateGrowthFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		gsl_matrix* currFg = gsl_matrix_calloc(3,3);
        for (int j=0; j<3; ++j){
            for(int k=0; k<3; ++k){
                double Fgjk;
                saveFileToDisplayGrowth.read((char*) &Fgjk, sizeof Fgjk);
                gsl_matrix_set(currFg,j,k,Fgjk);
            }
        }
        (*itElement)->setFg(currFg);
        gsl_matrix_free(currFg);
    }
}

void Simulation::updateGrowthRateFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double rx=0.0,ry=0.0,rz=0.0;
		saveFileToDisplayGrowthRate.read((char*) &rx, sizeof rx);
		saveFileToDisplayGrowthRate.read((char*) &ry, sizeof ry);
		saveFileToDisplayGrowthRate.read((char*) &rz, sizeof rz);
		(*itElement)->setGrowthRateExpFromInput(rx,ry,rz);
	}
}

void Simulation::updateProteinsFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double c[4];
		double cEq[4];
		for (int j = 0; j<4; j++){
			saveFileToDisplayProteins.read((char*) &c[j], sizeof c[j]);
		}
		for (int j = 0; j<4; j++){
			saveFileToDisplayProteins.read((char*) &cEq[j], sizeof cEq[j]);
		}
		(*itElement)->setMyosinLevels(c[0], c[1], c[2], c[3]);
		(*itElement)->setEquilibriumMyosinLevels(cEq[0], cEq[1], cEq[2], cEq[3]);
	}
}

void Simulation::updatePhysicalPropFromSave(){
	readPhysicalPropToContinueFromSave();
}

void Simulation::updatePackingFromSave(){
	//reading the number of packing nodes:
	int n;
	saveFileToDisplayPacking.read((char*) &n, sizeof n);
	//emptying the packing node vectors:
	pacingNodeCouples0.clear();
	pacingNodeCouples1.clear();
	pacingNodeCouplesHaveAdhered.clear();
	//filling in the packing node vectors for step
	for(int i=0; i<n; ++i){
		int elementId;
		saveFileToDisplayPacking.read((char*) &elementId, sizeof elementId);
		pacingNodeCouples0.push_back(elementId);
		saveFileToDisplayPacking.read((char*) &elementId, sizeof elementId);
		pacingNodeCouples1.push_back(elementId);
		pacingNodeCouplesHaveAdhered.push_back(false);
	}
	//reading the forces:
	for(int i=0; i<n; ++i){
		double Fx,Fy,Fz;
		saveFileToDisplayPacking.read((char*) &Fx, sizeof Fx);
		saveFileToDisplayPacking.read((char*) &Fy, sizeof Fy);
		saveFileToDisplayPacking.read((char*) &Fz, sizeof Fz);
		PackingForces[pacingNodeCouples0[i]][0] = Fx;
		PackingForces[pacingNodeCouples0[i]][1] = Fy;
		PackingForces[pacingNodeCouples0[i]][2] = Fz;
		saveFileToDisplayPacking.read((char*) &Fx, sizeof Fx);
		saveFileToDisplayPacking.read((char*) &Fy, sizeof Fy);
		saveFileToDisplayPacking.read((char*) &Fz, sizeof Fz);
		PackingForces[pacingNodeCouples1[i]][0] = Fx;
		PackingForces[pacingNodeCouples1[i]][1] = Fy;
		PackingForces[pacingNodeCouples1[i]][2] = Fz;
	}
}

void Simulation::readTensionCompressionToContinueFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		for (int j=0; j<6; ++j){
            double S;
            saveFileToDisplayTenComp.read((char*) &S, sizeof S);
            gsl_matrix_set((*itElement)->Strain,j,0,S);
		}
	}
}

void Simulation::readGrowthRedistributionToContinueFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		bool thereIsDistribution = false;
		bool shrinksElement = false;
		saveFileToDisplayGrowthRedistribution.read((char*) &thereIsDistribution, sizeof thereIsDistribution);
        (*itElement)->thereIsGrowthRedistribution = thereIsDistribution;
        saveFileToDisplayGrowthRedistribution.read((char*) &shrinksElement, sizeof shrinksElement);
        (*itElement)->growthRedistributionShrinksElement = shrinksElement;
	}
}

void Simulation::readNodeBindingToContinueFromSave(){
	/*int n = 0; //size Of Master-Slave List
	saveFileToDisplayNodeBinding.read((char*) &n, sizeof n);
	for (int i=0; i<n; ++i){
		int dofSlave, dofMaster;
		saveFileToDisplayNodeBinding.read((char*) &dofSlave, sizeof dofSlave);
		saveFileToDisplayNodeBinding.read((char*) &dofMaster, sizeof dofMaster);
		int dim = dofSlave % 3;
		int nodeSlave = (dofSlave - dim)/3;
		int nodeMaster = (dofMaster - dim)/3;
		Nodes[nodeSlave]->slaveTo[dim] = nodeMaster;
		Nodes[nodeMaster]->isMaster[dim] = true;
	}*/
	//clear all first:
	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		(*itNode)->slaveTo[0] = -1;
		(*itNode)->slaveTo[1] = -1;
		(*itNode)->slaveTo[2] = -1;
		(*itNode)->isMaster[0] = false;
		(*itNode)->isMaster[1] = false;
		(*itNode)->isMaster[2] = false;
		(*itNode)->attachedToPeripodial = false;
	}
	int n = 0; //size Of Master-Slave List
	saveFileToDisplayNodeBinding>>n;
	for (int i=0; i<n; ++i){
		int dofSlave, dofMaster;
		saveFileToDisplayNodeBinding>>dofSlave;
		saveFileToDisplayNodeBinding>>dofMaster;
		int dim = dofSlave % 3;
		int nodeSlave = (dofSlave - dim)/3;
		int nodeMaster = (dofMaster - dim)/3;
		Nodes[nodeSlave]->slaveTo[dim] = nodeMaster;
		Nodes[nodeMaster]->isMaster[dim] = true;
		if(Nodes[nodeSlave]->tissueType ==1 ){
			Nodes[nodeMaster]->attachedToPeripodial = true;
		}
		//cout<<"read binding: "<<i<<" of "<<n<<" slave :"<<nodeSlave<<" master  "<<nodeMaster<<endl;
	}
}


void Simulation::readCollapseAndAdhesionToContinueFromSave(){
	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		int adheredTo;
		saveFileToDisplayCollapseAndAdhesion.read((char*) &adheredTo, sizeof adheredTo);
		if ((*itNode)->adheredTo<0){
			(*itNode)->adheredTo = adheredTo;
		}
		if (savedStepwise){
			saveFileToDisplayCollapseAndAdhesion.read((char*) &(*itNode)->positionUpdateOngoing, sizeof (*itNode)->positionUpdateOngoing);
			saveFileToDisplayCollapseAndAdhesion.read((char*) &(*itNode)->positionUpdateCounter, sizeof (*itNode)->positionUpdateCounter);
		}
		int nCollapsedNodes;
		saveFileToDisplayCollapseAndAdhesion.read((char*) &nCollapsedNodes, sizeof nCollapsedNodes);
		//cout<<" node "<<(*itNode)->Id<<" adheredTo " <<(*itNode)->adheredTo<<" positionUpdateOngoing? "<<(*itNode)->positionUpdateOngoing<<" counter: "<<(*itNode)->positionUpdateCounter<<" nCollapsedNodes? "<<nCollapsedNodes<<endl;
		//vector<int> collapseList;
		for (int i=0; i<nCollapsedNodes; ++i){
			int nodeId;
			saveFileToDisplayCollapseAndAdhesion.read((char*) &nodeId, sizeof nodeId);
			//cout<<"read collapse: "<<" node :"<<(*itNode)->Id<<" adheredTo  "<<(*itNode)->adheredTo<<" reading collapse: "<<i<<" of "<<nCollapsedNodes<<" id: "<<nodeId<<endl;
			bool nodeOnList = binary_search((*itNode)->collapsedWith.begin(), (*itNode)->collapsedWith.end(),nodeId);
			if (!nodeOnList){
				(*itNode)->collapsedWith.push_back(nodeId);
			}
			//collapseList.push_back(nodeId);
		}
		(*itNode)->clearDuplicatesFromCollapseList();
		//(*itNode)->collapsedWith = collapseList;
		/*if ((*itNode)->adheredTo >  0|| (*itNode)->collapsedWith.size()>0){
			cout<<" node "<<(*itNode)->Id<<" is collapsed with ";
			for (int i=0; i<(*itNode)->collapsedWith.size(); ++i){
				cout<<(*itNode)->collapsedWith[i]<<" ";
			}
			cout<<" adhered to: "<<(*itNode)->adheredTo<<" is master: "<<(*itNode)->isMaster[0]<<" "<<(*itNode)->isMaster[1]<<" "<<(*itNode)->isMaster[2]<<" slave to: "<<(*itNode)->slaveTo[0]<<" "<<(*itNode)->slaveTo[1]<<" "<<(*itNode)->slaveTo[2]<<endl;
		}*/
	}
}

void Simulation::readGrowthToContinueFromSave(){
	updateGrowthFromSave();
}

void Simulation::readGrowthRateToContinueFromSave(){
	updateGrowthRateFromSave();
}


void Simulation::readProteinsToContinueFromSave(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double c[4];
		double cEq[4];
		for (int j = 0; j<4; j++){
			saveFileToDisplayProteins.read((char*) &c[j], sizeof c[j]);
		}
		for (int j = 0; j<4; j++){
			saveFileToDisplayProteins.read((char*) &cEq[j], sizeof cEq[j]);
		}
		(*itElement)->setMyosinLevels(c[0], c[1], c[2], c[3]);
		(*itElement)->setEquilibriumMyosinLevels(cEq[0], cEq[1], cEq[2], cEq[3]);
	}
}

void Simulation::readPhysicalPropToContinueFromSave(){
	vector<ShapeBase*>::iterator itElement;
	double E;
	double internalViscposity;
	double externalViscosity[3];
	double zRemodelling;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		saveFileToDisplayPhysicalProp.read((char*) &E, sizeof E);
		saveFileToDisplayPhysicalProp.read((char*) &internalViscposity, sizeof internalViscposity);
		saveFileToDisplayPhysicalProp.read((char*) &zRemodelling, sizeof zRemodelling);
		(*itElement)->setYoungsModulus(E);
		(*itElement)->setViscosity(internalViscposity);
		(*itElement)->setZRemodellingSoFar(zRemodelling);
	}
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		for (int i=0; i<3; ++i){
			saveFileToDisplayPhysicalProp.read((char*) &externalViscosity[i], sizeof externalViscosity[i]);
			(*itNode)->externalViscosity[i] = externalViscosity[i];
		}
	}
}

void Simulation::initiatePrismFromSaveForUpdate(int k){
	//inserts a new prism at order k into elements vector
	//the node ids and reference shape positions
	//will be updated in function: updateShapeFromSave
	int* NodeIds;
	NodeIds = new int[6];
	for (int i =0 ;i<6; ++i){
		//the node ids will be updated in function: updateShapeFromSave
		NodeIds[i] = 0;
	}
	Prism* PrismPnt01;
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
	PrismPnt01->updateShapeFromSave(saveFileToDisplayMesh);
	vector<ShapeBase*>::iterator it = Elements.begin();
	it += k;
	Elements.insert(it,PrismPnt01);
	nElements = Elements.size();
	currElementId++;
	delete[] NodeIds;
}


void Simulation::updateOneStepFromSave(){
	cout<<"updating step from save"<<endl;
	string currline;
	//skipping the header:
	getline(saveFileToDisplayMesh,currline);
	if(saveFileToDisplayMesh.eof()){
		reachedEndOfSaveFile = true;
		return;
	}
	//cout<<"skipped header: "<<currline<<endl;
    if (implicitPacking){
        resetForces(true);	// reset packing forces
    }
    else{
        resetForces(false);	// do not reset packing forces
    }
    updateNodeNumberFromSave();
	updateNodePositionsFromSave();
	updateElementStatesFromSave();
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		//This is updating positions from save.
		(*itElement)->updatePositions(Nodes);
	}
	if (TensionCompressionSaved){
        //cout<<"updating tension compression: "<<endl;
		updateTensionCompressionFromSave();
	}
    if (GrowthSaved){
        updateGrowthFromSave();
    }
    if (GrowthRateSaved){
    	updateGrowthRateFromSave();
    }
	if (ForcesSaved){
		updateForcesFromSave();
	}
	if(ProteinsSaved){
		updateProteinsFromSave();
	}
	if(physicalPropertiesSaved){
		updatePhysicalPropFromSave();
	}
	//cout<<"trying to update packing from save, PackingSaved: "<<PackingSaved<<endl;
	if (PackingSaved){
		updatePackingFromSave();
		//cout<<"updated packing, the size of packing couples: "<<pacingNodeCouples0.size()<<endl;
	}
	if (growthRedistributionSaved){
		updateGrowthRedistributionFromSave();
	}
	if (nodeBindingSaved){
		updateNodeBindingFromSave();
	}
	if (collapseAndAdhesionSaved){
		updateCollapseAndAdhesionFromSave();
	}

	clearNodeMassLists();
	assignNodeMasses();
	assignConnectedElementsAndWeightsToNodes();
	clearLaserAblatedSites();
	calculateBoundingBox();
	//calculateColumnarLayerBoundingBox();
	//if (thereIsPeripodialMembrane){
	//	calculatePeripodialBoundingBox();
	//}
	//skipping the footer:
	string currline2;
	getline(saveFileToDisplayMesh,currline2);
	//cout<<"currline 1st reading: "<<currline2<<endl;
	getline(saveFileToDisplayMesh,currline2);
	//cout<<"currline 2nd reading: "<<currline2<<endl;
	while (currline.empty() && !saveFileToDisplayMesh.eof()){
		//skipping empty line
		//cout<<"skipping empty line"<<endl;
		getline(saveFileToDisplayMesh,currline2);
	}
	//if(saveFileToDisplayMesh.eof()){
	//	reachedEndOfSaveFile = true;
	//	return;
	//}
	//cout<<"in step update, skipped footer: "<<currline2<<endl;
	timestep = timestep + dataSaveInterval;
	currSimTimeSec += dt*dataSaveInterval;
}

void  Simulation::updateNodeNumberFromSave(){
	//cout<<"Updating number of nodes from save"<<endl;
	//cout<<"Is save file open: "<<saveFileToDisplay.is_open()<<" is file good? "<<saveFileToDisplay.good()<<endl;
	int n;
	saveFileToDisplayMesh>> n;
	int currNodeNumber = Nodes.size();
	//cout<<"number of nodes from save file: "<<n <<" number of nodes on the vector: "<<currNodeNumber<<endl;
	if (n>currNodeNumber){
		Node* tmp_nd;
		for (int i = 0; i<(n-currNodeNumber); ++i){
			double* pos = new double[3];
			pos[0]=0.0;
			pos[1]=0.0;
			pos[2]=0.0;
			//the positions will be read and updated in function updateNodePositionsFromSave
			tmp_nd = new Node(i, 3, pos,-1, -1);
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			delete[] pos;
		}
	}
	else{
		for (int i = 0; i<(currNodeNumber-n); ++i){
			Node* tmp_nd;
			tmp_nd = Nodes.back();
			Nodes.pop_back();
			nNodes = Nodes.size();
			delete tmp_nd;
		}
	}
	n = Nodes.size();
	if ( n != nNodes){
		//the node number is change, I updated the node list, now I need to fix system forces:
		reInitiateSystemForces(nElements);
	}
	//cout<<"end of funciton, number of nodes from save file: "<<n <<" number of nodes on the vector: "<<Nodes.size()<<endl;
}

void  Simulation::updateElementStatesFromSave(){
	//cout<<"Updating element states from save"<<endl;
	int n;
	saveFileToDisplayMesh >> n;
	int currElementNumber = Elements.size();
	//The elements list is bigger than necessary, I am deleting elements form the end of the list
	//cout<<"number of elements from save file: "<<nElements <<" number of element on the Elements vector: "<<currElementNumber<<endl;
	while(currElementNumber>n){
		//cout<<"removing elements from list, current number vector length: "<<currElementNumber<<endl;
		removeElementFromEndOfList();
		currElementNumber = Elements.size();
	}
	for (int i=0; i<currElementNumber; ++i){
		//cout<<"reading element: "<<i<<endl;
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		int currShapeType = Elements[i]->getShapeType();
		if (shapeType == currShapeType || shapeType == 2 || shapeType ==4){
			//cout<<"shape type is correct, moving on to update"<<endl;
			//The current shape on the list, and the shape I am reading are of same type, I can read it:
			if (shapeType ==4){
				double height;
				saveFileToDisplayMesh >> height;
			}
			Elements[i]->updateShapeFromSave(saveFileToDisplayMesh);
		}
		else{
			//cout<<"shape type is wrong"<<endl;
			//the current element is a different type, I need to insert generate a new element, and insert it here:
			if (shapeType == 1){
				//the new shape is a prism, I will add it now;
				initiatePrismFromSaveForUpdate(i);
			}
			removeElementFromEndOfList();
			currElementNumber = Elements.size();
		}
	}
	while(n>currElementNumber){
		int shapeType;
		saveFileToDisplayMesh >> shapeType;
		if (shapeType == 1){
			int i = Elements.size()-1;
			//this will initiate a prism at the current point in the Elements vector
			//then it will read the node ids and reference positions from the save file
			//it will update the node ids, and reference positions.
			//the normal positions will be updated using function updatePositions, called by updateOneStepFromSave
			initiatePrismFromSaveForUpdate(i);
			currElementNumber = Elements.size();
		}
	}
}

void Simulation::removeElementFromEndOfList(){
	ShapeBase* tmp_element;
	tmp_element = Elements.back();
	//int shapetype = tmp_element->getShapeType();
	Elements.pop_back();
	nElements = Elements.size();
	delete tmp_element;
}

void Simulation::updateNodePositionsFromSave(){
	//cout<<"Updating node positions from save"<<endl;
	//I have already read the current number of nodes in function "updateNodeNumberFromSave",
	//and I updated the number of nodes where necessary
	//the "cursor" in the file progressed, and is at the beginning of positions now
	for (int i=0; i<nNodes; ++i){
		saveFileToDisplayMesh >> Nodes[i]->Position[0];
		saveFileToDisplayMesh >> Nodes[i]->Position[1];
		saveFileToDisplayMesh >> Nodes[i]->Position[2];
		saveFileToDisplayMesh >> Nodes[i]->tissuePlacement;
		saveFileToDisplayMesh >> Nodes[i]->tissueType;
		saveFileToDisplayMesh >> Nodes[i]->atCircumference;
	}
}

bool Simulation::initiateMesh(int MeshType, float zHeight){
	if (MeshType == 1 ){
		initiateSinglePrismNodes(zHeight);
		initiateSinglePrismElement();
	}
	else if (MeshType == 2){
		cerr<<"Error: Too few arguments for mesh by dimensions"<<endl;
		return false;
	}
	else if ( MeshType == 3 || MeshType ==4 ){
		cerr<<"Error: Wrong set of arguments  for mesh triangulation"<<endl;
		return false;
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
		return false;
	}
	return true;
}

bool Simulation::initiateMesh(int MeshType, int Row, int Column, float SideLength, float zHeight){
	if ( MeshType == 1 ){
		cerr<<"Error: Too many arguments for a single element system"<<endl;
		return false;
	}
	if ( MeshType == 2){
		//The necessary parameters:
		//InputMeshParameters:
		//  MeshInputMode(int-seeDocumentation): 2
		//  MeshRow(int): 3
		//  MeshColumn(int): 1
		//  SideLength: 1.0
		//  zHeight: 1.0
		//  ApicalCircumferenceFixZ: 0
		//  ApicalCircumferenceFixXY: 0
		//  BasalCircumferenceFixZ: 0
		//  BasalCircumferenceFixXY: 0
		initiateNodesByRowAndColumn(Row,Column,SideLength,zHeight);
		initiateElementsByRowAndColumn(Row,Column);
	}
	else if ( MeshType == 3 ){
		cerr<<"Error: Wrong set of arguments for mesh triangulation"<<endl;
		return false;
	}
	else if ( MeshType ==4 ){
		cerr<<"Error: Too many arguments for reading the mesh from file"<<endl;
		return false;
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
		return false;
	}
	return true;
}

/*bool Simulation::initiateMesh(int MeshType, string inputtype, float SideLength, float zHeight ){
	if ( MeshType == 1 ){
		cerr<<"Error: Too many arguments for a single element system"<<endl;
		return false;
	}
	if ( MeshType == 2){
		cerr<<"Error: Too few arguments for mesh by dimensions"<<endl;
		return false;
	}
	else if ( MeshType == 3 ){
		//this will be inputting circumference of the tissue, and the sidelength and z-height
		//generate mesh by triangulation
		return false;
	}
	else if ( MeshType == 4 ){
		cerr<<"Error: Too many arguments for reading the mesh from file"<<endl;
		return false;
		//this will be reading full mesh data
		//read mesh data from file
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
		return false;
	}
	return true;
}*/

bool Simulation::initiateMesh(int MeshType){
	if ( MeshType == 1 ){
		cerr<<"Error: Too many arguments for a single element system"<<endl;
		return false;
	}
	if ( MeshType == 2){
		cerr<<"Error: Too few arguments for mesh by dimensions"<<endl;
		return false;
	}
	else if ( MeshType == 3 ){
		cerr<<"Error: Wrong set of arguments  for mesh triangulation"<<endl;
		return false;
	}
	else if ( MeshType ==4 ){
		//this will be reading full mesh data
		//read mesh data from file
		const char* name_inputMeshFile = inputMeshFileName.c_str();;
		inputMeshFile.open(name_inputMeshFile, ifstream::in);
		if (!(inputMeshFile.good() && inputMeshFile.is_open())){
			cerr<<"Cannot open the input mesh file: "<<name_inputMeshFile<<endl;
			return false;
		}
		cout<<"initiating nodes"<<endl;
		initiateNodesFromMeshInput();
		initiateElementsFromMeshInput();
		cout<<" nNodes: "<<nNodes<<" nEle: "<<nElements<<endl;
		bool areTissueWeightsRecorded = checkIfTissueWeightsRecorded();
		if (areTissueWeightsRecorded){
			readInTissueWeights();
			cout<<" read in tissue weights"<<endl;
		}
		inputMeshFile.close();
		//Alter the mesh:
		//delete all nodes and elements with nodes and elements above a certain z:
		vector<int> deletedNodes, deletedElements;
		for (int i=0; i<nNodes; ++i){
			if (Nodes[i]->Position[2] > 14.25){//I do not want any node above 14.25;
				deletedNodes.push_back(i);
			}
		}
		int nAN = deletedNodes.size();
		vector<ShapeBase*>::iterator itElement;
		for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			if(!(*itElement)->IsAblated){
				for (int j =0; j<nAN; ++j){
					bool IsAblatedNow = (*itElement)->DoesPointBelogToMe(deletedNodes[j]);
					if (IsAblatedNow){
							deletedElements.push_back((*itElement)->Id);
					}
				}
			}
		}
		//generate new node list
		//generate new elemetn list
		//update element ids
		//update node ids
	}
	else {
		cerr<<"Error: Mesh Type not recognised"<<endl;
		return false;
	}
	return true;
}

bool Simulation::checkIfTissueWeightsRecorded(){
	bool tissueWeightsRecorded;
	saveFileToDisplayMesh >> tissueWeightsRecorded;
	return tissueWeightsRecorded;
}

void Simulation::readInTissueWeights(){
	vector<ShapeBase*>::iterator itEle;
	double wPeri;
	for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
		inputMeshFile >> wPeri;
		(*itEle)->setGrowthWeightsViaTissuePlacement(wPeri); //peripodialness weight is recorded
	}
}

bool Simulation::generateColumnarCircumferenceNodeList(	vector <int> &ColumnarCircumferencialNodeList){
	//generating a list of nodes that are at the circumference and at the basal surface
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->atCircumference && Nodes[i]->tissuePlacement == 0){ // tissuePlacement = 0 -> basal node
			ColumnarCircumferencialNodeList.push_back(i);

		}
	}
	int n = ColumnarCircumferencialNodeList.size();
	if (n<=0){
		cerr<<"No circumferncial nodes indicated! Cannot generate PeripodialMembrane"<<endl;
		AddPeripodialMembrane = false;
		thereIsPeripodialMembrane = false;
		return false;
	}
	return true;
}

void Simulation::clearCircumferenceDataFromSymmetricityLine(){
	//I will take out anything that was at the border of symmetry, but I need the nodes at the tips. So find the x tips first:
	if (symmetricY){
		double xTipPos = -1000, xTipNeg = 1000;
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ((*itNode)->Position[0] > xTipPos ){
				xTipPos = (*itNode)->Position[0];
			}
			if ((*itNode)->Position[0] < xTipNeg ){
				xTipNeg = (*itNode)->Position[0];
			}
		}
		double yLimPos = 0.1;
		double yLimNeg = (-1.0)*yLimPos;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			double x = (*itNode)->Position[0];
			double y = (*itNode)->Position[1];
			if ( y < yLimPos){
				if ( y  > yLimNeg){
					(*itNode)->atSymmetricityBorder = true;
					fixY((*itNode),false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
					//the node is indeed at the border, BUT, I will remove it only if it is not at the tip:
					bool atTip = false;
					if ( x < xTipPos+yLimPos && x >xTipPos+yLimNeg){
						atTip = true;
					}
					if (!symmetricX){
						//This if clause is checking the x-tip at the negative end.
						//As described above, the node should still be a circumference node,
						//if it is at the line of symmetry, but at the tip as well.
						//BUT, if there is also xSymmetry, this "negative end tip" will
						//be the tip at the x symmetry line. (If I am modelling a quarter of a circle,
						//with x and y symmetry, this will point be the centre of the circle.)
						//Under these conditions, it is not actually at the tip. It should be
						//removed from circumference node list.
						//Therefore, I am carrying the atTip? check for the negative end only if
						//there is no x-symmetry.
						if ( x > xTipNeg+yLimNeg && x < xTipNeg+yLimPos){
							atTip = true;
						}
					}
					if (!atTip){
						//removing node from list:
						(*itNode)->atCircumference = false;
					}
				}
			}
		}
	}
	if (symmetricX){
		double yTipPos = -1000, yTipNeg = 1000;
		for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ((*itNode)->Position[0] > yTipPos ){
				yTipPos = (*itNode)->Position[1];
			}
			if ((*itNode)->Position[0] < yTipNeg ){
				yTipNeg = (*itNode)->Position[1];
			}
		}
		double xLimPos = 0.1;
		double xLimNeg = (-1.0)*xLimPos;
		for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			double x = (*itNode)->Position[0];
			double y = (*itNode)->Position[1];
			if ( x < xLimPos){
				if ( x  > xLimNeg){
					(*itNode)->atSymmetricityBorder = true;
					fixX((*itNode),false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
					//the node is indeed at the border, BUT, I will remove it only if it is not at the tip:
					bool atTip = false;
					if ( y < yTipPos+xLimPos && y >yTipPos+xLimNeg){
						atTip = true;
					}
					if ( y > yTipNeg+xLimNeg && y < yTipNeg+xLimPos){
						atTip = true;
					}
					if (!atTip){
						//removing node from list:
						(*itNode)->atCircumference = false;
					}
				}
			}
		}
	}
}
/*
void Simulation::removeSymmetryBorderFromColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList){
	//I will take out anything that was at the border of symmetry, but I need the nodes at the tips. So find the x tips first:
	double xTipPos = -1000, xTipNeg = 1000;
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->Position[0] > xTipPos ){
			xTipPos = (*itNode)->Position[0];
		}
		if ((*itNode)->Position[0] < xTipNeg ){
			xTipNeg = (*itNode)->Position[0];
		}
	}
	double yLimPos = 0.1;
	double yLimNeg = (-1.0)*yLimPos;
	int n = ColumnarCircumferencialNodeList.size();
	int i=0;
	while (i<n){
		double x = Nodes[ColumnarCircumferencialNodeList[i]]->Position[0];
		double y = Nodes[ColumnarCircumferencialNodeList[i]]->Position[1];
		if ( y < yLimPos){
			if ( y  > yLimNeg){
				Nodes[ColumnarCircumferencialNodeList[i]]->atSymmetricityBorder = true;
				fixY(Nodes[ColumnarCircumferencialNodeList[i]]);
				//the node is indeed at the border, BUT, I will remove it only if it is not at the tip:
				bool atTip = false;
				if ( x < xTipPos+yLimPos && x >xTipPos+yLimNeg){
					atTip = true;
				}
				if ( x > xTipNeg+yLimNeg && x < xTipNeg+yLimPos){
					atTip = true;
				}
				if (!atTip){
					//removing node from list:
					ColumnarCircumferencialNodeList.erase (ColumnarCircumferencialNodeList.begin()+i);
					Nodes[ColumnarCircumferencialNodeList[i]]->atCircumference = false;
					n--;
					i--;
				}
			}
		}
		i++;
	}
}
*/
void Simulation::sortColumnarCircumferenceNodeList(vector <int> &ColumnarCircumferencialNodeList){
	//ordering the circumferencial nodes of the basal surface in clockwise rotation
	int n = ColumnarCircumferencialNodeList.size();
	vector <double> angles;
	for (int j =0 ; j<n; ++j){
		double x = Nodes[ColumnarCircumferencialNodeList[j]]->Position[0];
		double y = Nodes[ColumnarCircumferencialNodeList[j]]->Position[1];
		double tet = atan2(y,x);
		if (tet<0){tet += 2.0*3.14;}
		angles.push_back(tet);
	}

	bool swapped = true;
	while (swapped){
		swapped = false;
		for(int i=1; i<n; ++i){
			if(angles[i]<angles[i-1]){
				int temp=ColumnarCircumferencialNodeList[i-1];
				ColumnarCircumferencialNodeList[i-1]=ColumnarCircumferencialNodeList[i];
				ColumnarCircumferencialNodeList[i]=temp;
				double td = angles[i-1];
				angles[i-1]=angles[i];
				angles[i]=td;
				swapped = true;
				//cout<<"swapped "<<i <<" with "<<i-1<<endl;
			}
		}
	}
}

void Simulation::getAverageSideLength(double& periAverageSideLength, double& colAverageSideLength){
	double dsumPeri =0.0, dsumCol = 0.0;
	int colCounter =0, periCounter=0;
	double packingDetectionThresholdGridCounter[10][4];
	for (int i=0; i<10;++i){
		for (int j=0;j<5;++j){
			packingDetectionThresholdGrid[i][j] = 0.0;
			packingDetectionThresholdGridCounter[i][j] = 0.0;
		}
	}
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->IsAblated){
			//do not count ablated elements
			if ((*itElement)->tissueType==0){ //element belongs to columnar layer
				double currValue = (*itElement)->getApicalSideLengthAverage();
				if(currValue>0){
					//There may be alements with fully collapsed apical areas, where currValue will be zero
					//I do not wish to count them in averaging
					dsumCol += currValue;
					colCounter++;
				}
				//if the element is at the basal surface, also get the basal lengths into consideration:
				if((*itElement)->tissuePlacement == 0 || (*itElement)->spansWholeTissue){
					currValue = (*itElement)->getBasalSideLengthAverage();
					if(currValue>0){
						//There may be alements with fully collapsed apical areas, where currValue will be zero
						//I do not wish to count them in averaging
						dsumCol += currValue;
						colCounter++;
					}
				}
				double* reletivePos = new double[2];
				(*itElement)->getRelativePosInBoundingBox(reletivePos);
				int relX = floor(reletivePos[0]);
				int relY = floor(reletivePos[1]/2.0);
				if (relX < 0) {relX = 0;}
				if (relX > 9) {relX = 9;}
				if (relY < 0) {relY = 0;}
				if (relY > 4) {relY = 4;}
				delete[] reletivePos;
				packingDetectionThresholdGrid[relX][relY]+=currValue;
				packingDetectionThresholdGridCounter[relX][relY]++;

			}
			else{
				dsumPeri += (*itElement)->getApicalSideLengthAverage();
				periCounter++;
			}
		}
	}
	colAverageSideLength = dsumCol / (double)colCounter;
	for (int i=0; i<10;++i){
		for (int j=0;j<5;++j){
			if(packingDetectionThresholdGridCounter[i][j]>0){
				packingDetectionThresholdGrid[i][j] /= packingDetectionThresholdGridCounter[i][j];
			}
		}
	}
	if (periCounter>0){
		periAverageSideLength = dsumPeri / (double)periCounter;
	}
}

bool  Simulation::isColumnarLayer3D(){
	bool ColumnarLayer3D = false;
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->ShapeDim == 3){
			ColumnarLayer3D  = true;
			break;
		}
	}
	return ColumnarLayer3D;
}

bool Simulation::calculateTissueHeight(){
	//check if the columnar layer is made of 3D elements:
	bool ColumnarLayer3D = isColumnarLayer3D();
	TissueHeight = 0;
	//cout<<"inside calculateTissueHeight"<<endl;
	if (ColumnarLayer3D){
		//cout<<"tissue is 3d"<<endl;
		//Find the first basal node on the Nodes List
		//Find the first element that has this node on the elements list
		//Move apically through the elements until you reach the apical surface - an apical node
		//Find a basal node:
		vector<Node*>::iterator itNode;
		bool foundNode = false;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if((*itNode)->tissueType == 0 && (*itNode)->tissuePlacement == 0){ //Columnar node is basal
				foundNode = true;
				break;
			}
		}
		if (!foundNode){
			return false;
		}
		//Find an element using the basal node, and move on the elements apically, until you reach the apical surface:
		int currNodeId = (*itNode)->Id;
		cout<<"found node: "<<foundNode<<" currNodeId: "<<currNodeId<<endl;
		vector<ShapeBase*>::iterator itElement;
		bool foundElement = true;
		TissueHeightDiscretisationLayers = 0;
		while(Nodes[currNodeId]->tissuePlacement != 1 && foundElement){ //while the node is not apical, and I could find the next element
			foundElement = false;
			for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
				bool IsBasalOwner = (*itElement)->IsThisNodeMyBasal(currNodeId);
				if (IsBasalOwner){
					foundElement = true;
					break;
				}
			}
			//cout<<"found element: "<<foundElement<<" currNodeId: "<<currNodeId<<endl;
			double currentH = (*itElement)->getElementHeight();
			TissueHeight += currentH;
			TissueHeightDiscretisationLayers++;
			currNodeId = (*itElement)->getCorrecpondingApical(currNodeId); //have the next node
			//cout<<"corresponding apical node: "<<currNodeId<<" foundElement bool at this point: "<<foundElement<<" tissueplacement: "<<Nodes[currNodeId]->tissuePlacement<<endl;
		}
		if (!foundElement){
			return false;
		}
	}
	else{
		TissueHeight = Elements[0]->ReferenceShape->height;
		if (TissueHeight == 0){
			cout<<"The coulmanr layer is 2D, but the tissue height of the elements is not assigned properly, cannot obtain TissueHeight"<<endl;
			return false;
		}
	}
	//cout<<"checking for peripodial & lumen in calculateTissueHeight"<<endl;

	if (thereIsPeripodialMembrane){
		double columnarTop = -1000.0, peripodialbottom = 1000.0;
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ((*itNode)->tissueType == 0 ){
				//columnar node
				if ((*itNode)->tissuePlacement == 1){
					//apical node:
					if((*itNode)->Position[2]> columnarTop){
						columnarTop = (*itNode)->Position[2];
					}
				}
			}
			else if ((*itNode)->tissueType == 1 ){
				//peripodial node
				if ((*itNode)->tissuePlacement == 0){
					//basal node:
					if((*itNode)->Position[2]< peripodialbottom){
						peripodialbottom = (*itNode)->Position[2];
					}
				}
			}
		}
		lumenHeight = columnarTop - peripodialbottom;
	}
	//cout<<"finalised calculateTissueHeight"<<endl;
	return true;
}

void Simulation::assignInitialZPositions(){
	for(vector<ShapeBase*>::iterator itElement = Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->setInitialZPosition(boundingBox[0][2], TissueHeight);  //minimum z of the tissue and the tissue height given as input
	}
}

void Simulation::calculateStiffnessMatrices(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		cout<<" setting up element :  "<<(*itElement)->Id<<" of "<<nElements<<endl;
		(*itElement)->calculateReferenceStiffnessMatrix();
	}
}

void Simulation::calculateShapeFunctionDerivatives(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		cout<<" setting up element :  "<<(*itElement)->Id<<" of "<<nElements<<endl;
		(*itElement)->calculateElementShapeFunctionDerivatives();
    }
}

void Simulation::fixAllD(Node* currNode, bool fixWithViscosity){
	for (int j =0 ; j<currNode->nDim; ++j){
		if(fixWithViscosity){
			currNode->externalViscosity[j] = fixingExternalViscosity[j];
			//currNode->baseExternalViscosity[j] = currNode->externalViscosity[j];
			currNode->externalViscositySetInFixing[j] = true;
		}
		else{
			currNode->FixedPos[j]=true;
		}
	}
	//cout<<" fixAllD called on node: "<<currNode->Id<<endl;
}

void Simulation::fixAllD(int i, bool fixWithViscosity){
	for (int j =0 ; j<Nodes[i]->nDim; ++j){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[j] = fixingExternalViscosity[j];
			//Nodes[i]->baseExternalViscosity[j] = Nodes[i]->externalViscosity[j];
			Nodes[i]->externalViscositySetInFixing[j] = true;
		}
		else{
			Nodes[i]->FixedPos[j]=true;
		}
	}
}

void Simulation::fixX(Node* currNode, bool fixWithViscosity){
	if(currNode->nDim>0){
		if(fixWithViscosity){
			currNode->externalViscosity[0] = fixingExternalViscosity[0];
			//currNode->baseExternalViscosity[0] = currNode->externalViscosity[0];
			currNode->externalViscositySetInFixing[0] = true;
		}
		else{
			currNode->FixedPos[0]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<currNode->Id<<" does not have x-dimension"<<endl;
	}
}
void Simulation::fixX(int i, bool fixWithViscosity){
	if(Nodes[i]->nDim>0){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[0] = fixingExternalViscosity[0];
			//Nodes[i]->baseExternalViscosity[0] = Nodes[i]->externalViscosity[0];
			Nodes[i]->externalViscositySetInFixing[0] = true;
		}
		else{
			Nodes[i]->FixedPos[0]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<Nodes[i]->Id<<" does not have x-dimension"<<endl;
	}

}

void Simulation::fixY(Node* currNode, bool fixWithViscosity){
	if(currNode->nDim>1){
		if(fixWithViscosity){
			currNode->externalViscosity[1] = fixingExternalViscosity[1];
			//currNode->baseExternalViscosity[1] = currNode->externalViscosity[1];
			currNode->externalViscositySetInFixing[1] = true;
		}
		else{
			currNode->FixedPos[1]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<currNode->Id<<" does not have y-dimension"<<endl;
	}
}

void Simulation::fixY(int i, bool fixWithViscosity){
	if(Nodes[i]->nDim>1){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[1] = fixingExternalViscosity[1];
			//Nodes[i]->baseExternalViscosity[1] = Nodes[i]->externalViscosity[1];
			Nodes[i]->externalViscositySetInFixing[1] = true;
		}
		else{
			Nodes[i]->FixedPos[1]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<Nodes[i]->Id<<" does not have y-dimension"<<endl;
	}
}

void Simulation::fixZ(Node* currNode, bool fixWithViscosity){
	if(currNode->nDim>2){
		if(fixWithViscosity){
			currNode->externalViscosity[2] = fixingExternalViscosity[2];
			//currNode->baseExternalViscosity[2] = currNode->externalViscosity[2];
			currNode->externalViscositySetInFixing[2] = true;
		}
		else{
			currNode->FixedPos[2]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<currNode->Id<<" does not have z-dimension"<<endl;
	}
}

void Simulation::fixZ(int i, bool fixWithViscosity){
	if(Nodes[i]->nDim>2){
		if(fixWithViscosity){
			Nodes[i]->externalViscosity[2] = fixingExternalViscosity[2];
			//Nodes[i]->baseExternalViscosity[2] = Nodes[i]->externalViscosity[2];
			Nodes[i]->externalViscositySetInFixing[2] = true;
		}
		else{
			Nodes[i]->FixedPos[2]=true;
		}
	}
	else{
		cerr<<"ERROR: Node : "<<Nodes[i]->Id<<" does not have z-dimension"<<endl;
	}
}

void Simulation::zeroForcesOnNode(int i){
	double ForceBalance[3];
	ForceBalance[0] = SystemForces[i][0];
	ForceBalance[1] = SystemForces[i][1];
	ForceBalance[2] = SystemForces[i][2];
	for (int i=0;i<nNodes;++i){
		SystemForces[i][0]-=ForceBalance[0];
		SystemForces[i][1]-=ForceBalance[1];
		SystemForces[i][2]-=ForceBalance[2];
	}
}

void Simulation::initiateSystemForces(){
	const int n = nNodes;
	//n nodes
	SystemForces = new double*[n];
	PackingForces = new double*[n];
	//PackingForcesPreviousStep = new double*[n];
	//PackingForcesTwoStepsAgoStep = new double*[n];
	FixedNodeForces = new double*[n];
	for (int j=0;j<n;++j){
		//3 dimensions
		SystemForces[j] = new double[3];
		PackingForces[j] = new double[3];
		//PackingForcesPreviousStep[j] = new double[3];
		//PackingForcesTwoStepsAgoStep[j] = new double[3];
		FixedNodeForces[j] = new double[3];
		SystemForces[j][0]=0.0;
		SystemForces[j][1]=0.0;
		SystemForces[j][2]=0.0;
		PackingForces[j][0]=0.0;
		PackingForces[j][1]=0.0;
		PackingForces[j][2]=0.0;
		//PackingForcesPreviousStep[j][0]=0.0;
		//PackingForcesPreviousStep[j][1]=0.0;
		//PackingForcesPreviousStep[j][2]=0.0;
		//PackingForcesTwoStepsAgoStep[j][0]=0.0;
		//PackingForcesTwoStepsAgoStep[j][1]=0.0;
		//PackingForcesTwoStepsAgoStep[j][2]=0.0;
		FixedNodeForces[j][0] = 0.0;
		FixedNodeForces[j][1] = 0.0;
		FixedNodeForces[j][2] = 0.0;
		//cout<<"systemforces[i][j]: "<<SystemForces[i][0]<<" "<<SystemForces[i][0]<<" "<<SystemForces[i][0]<<endl;
	}
}

void Simulation::initiateSinglePrismNodes(float zHeight){
	double *pos = new double[3];
	Node* tmp_nd;
	pos[0]=0;pos[1]=1;pos[2]=0;
	tmp_nd = new Node(0, 3, pos,0,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	pos[0]=1;pos[1]=0;pos[2]=0;
	tmp_nd = new Node(1, 3, pos,0,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	pos[0]=0;pos[1]=0;pos[2]=0;
	tmp_nd = new Node(2, 3, pos,0,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	pos[0]=0;pos[1]=1;pos[2]=zHeight;
	tmp_nd = new Node(3, 3, pos,1,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	pos[0]=1;pos[1]=0;pos[2]=zHeight;
	tmp_nd = new Node(4, 3, pos,1,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	pos[0]=0;pos[1]=0;pos[2]=zHeight;
	tmp_nd = new Node(5, 3, pos,1,0);
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	delete[] pos;
}

void Simulation::initiateSinglePrismElement(){
	int* NodeIds;
	NodeIds = new int[6];
	for (int i = 0; i < 6 ; i++){
		NodeIds[i]=i;
	}
	Prism* PrismPnt01;
	PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
	Elements.push_back(PrismPnt01);
	nElements = Elements.size();
	currElementId++;
	fixZ(0, BasalNodeFixWithExternalViscosity);
	fixZ(1, BasalNodeFixWithExternalViscosity);
	fixZ(2, BasalNodeFixWithExternalViscosity);
}


void Simulation::initiateNodesByRowAndColumn(int Row, int Column, float SideLength, float zHeight){
	dorsalTipIndex = Row;
	ventralTipIndex = 0;
	//The height of the equilateral triangle with side length: SideLength
	double sqrt3 = 1.7321;
	float h = sqrt3/2*SideLength;
	vector <double> xPos, yPos;
	vector <int> NodesToFix;
	int toprowcounter = 0;	//number of nodes at the terminating end, the top row. I will need this to add the sideway prisms;
	for (int ColCount = 0; ColCount < Column+1; ++ColCount){
		double CurrY = ColCount*h;
		int CurrRowNum = Row + 1 - ColCount;
		double RowOffset = 0.5*SideLength*ColCount;
		for ( int RowCount = 1; RowCount<CurrRowNum+1; ++RowCount){
			double CurrX = RowOffset + RowCount * SideLength;
			xPos.push_back(CurrX);
			yPos.push_back(CurrY);
			if (RowCount ==1 || RowCount == CurrRowNum || ColCount == Column){
				NodesToFix.push_back(xPos.size()-1);
				if(ColCount == Column){
					toprowcounter++;
				}
			}
		}
		if (ColCount>0){
			CurrY = (-1.0)*CurrY;
			for ( int RowCount = 1; RowCount<CurrRowNum+1; ++RowCount){
				double CurrX = RowOffset + RowCount * SideLength;
				xPos.push_back(CurrX);
				yPos.push_back(CurrY);
				if (RowCount ==1 || RowCount == CurrRowNum || ColCount == Column){
					NodesToFix.push_back(xPos.size()-1);
				}
			}
		}
	}
	int n =  xPos.size();
	Node* tmp_nd;
	double* pos = new double[3];
	//Adding the basal level of nodes, all will form columnar elements:
	for (int i =0; i< n; ++i){
		pos[0] = xPos[i];
		pos[1] = yPos[i];
		pos[2] = 0.0;
		tmp_nd = new Node(i, 3, pos,0,0);
		Nodes.push_back(tmp_nd);
		nNodes = Nodes.size();
	}
	//Adding the apical level, all will form columnar elements:
	for (int i =0; i< n; ++i){
		pos[0] = xPos[i];
		pos[1] = yPos[i];
		pos[2] = zHeight;
		tmp_nd = new Node(n+i, 3, pos,1,0);
		Nodes.push_back(tmp_nd);
		nNodes = Nodes.size();
		nNodes = Nodes.size();
	}
	delete[] pos;
}

void Simulation::setLinkerCircumference(){
	//First find the node that is further back on
	vector<Node*>::iterator itNode;
	double maxX = -100, ZofTip = -100;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType == 2){
			//The node is linker, is the x position higher than already recorded?
			if ((*itNode)->Position[0]> maxX){
				maxX  = (*itNode)->Position[0];
				ZofTip = (*itNode)->Position[2];
			}
		}
	}
	double thres = 0.2;
	//cout<<"Setting circumference, maxX: "<<maxX<<" Z of tip: "<<ZofTip<<endl;
	//Now I have the maxX, and the corresponding z height.
	//Declare all linkers at the basal/or apical surface are the circumferential nodes:
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType == 2){
			//The node is linker, if it is apical or basal, it can be circumferntial:
			if ( (*itNode)->tissuePlacement == 0 || (*itNode)->tissuePlacement == 1 ){
				//The node is linker, if it is in the range of the height of the tip then it is circumferential
				if ( (*itNode)->Position[2] < ZofTip+thres && (*itNode)->Position[2] > ZofTip-thres ){
					(*itNode)->atCircumference = true;
				}
			}

		}
	}
}

void Simulation::checkForNodeFixing(){
	//Are there any circumferential node fixing options enabled:
	bool thereIsCircumFix = false;
	if (!thereIsCircumFix){
		for (int i=0;i<5; ++i){
			for (int j=0;j<3; ++j){
				if (CircumferentialNodeFix[i][j] == true){
					thereIsCircumFix = true;
					break;
				}
			}
			if (thereIsCircumFix){
				break;
			}
		}
	}
	//If there is any circumferential Node fixing:
	if (thereIsCircumFix){
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ( (*itNode)->atCircumference){
				//cout<<"Node "<<(*itNode)->Id<<" at circumference"<<endl;
				//The node is at circumference, I would like to fix only the columnar side,
				bool doNotFix = false;
				if ((*itNode)->tissueType == 1){ //the node is peripodial, no not fix this node
					doNotFix = true;
				}
				/*else if ((*itNode)->tissueType == 2){ //the node is linker (as will be the case for all circumference)
					//I will check if any of its neigs are peripodial:
					int nNeig = (*itNode)->immediateNeigs.size();
					for (int k = 0; k<nNeig; ++k){
						if (Nodes[(*itNode)->immediateNeigs[k]]->tissueType == 1){
							doNotFix = true;
							break;
						}
					}
				}*/
				if (!doNotFix){
					//cout<<"Node "<<(*itNode)->Id<<" in fix options"<<endl;
					//Node is at the circumference, now checking for all possibilities:
					// if i == 0 , I am checking for apical circumference
					// if i == 1 , I am checking for basal  circumference
					// if i == 2 , I am checking for the linker apical circumference
					// if i == 3 , I am checking for the linker basal circumference
					// if i == 4 , I am checking for all    circumference
					for (int i=0;i<5; ++i){
						if ( (i == 0 && (*itNode)->tissuePlacement == 1 ) ||  //tissuePlacement == 1 is apical
							 (i == 1 && (*itNode)->tissuePlacement == 0 ) ||  //tissuePlacement == 0 is basal
							 (i == 2 && (*itNode)->tissueType == 2 && (*itNode)->tissuePlacement == 1 ) ||  //tissuePlacement == 1 is apical
							 (i == 3 && (*itNode)->tissueType == 2 && (*itNode)->tissuePlacement == 0 ) ||  //tissuePlacement == 0 is basal
							 (i == 4)){										  //tissuePlacement is irrelevant, fixing all
							//The node is at circumference; if
							if (CircumferentialNodeFix[i][0]){
								fixX((*itNode),CircumferentialNodeFixWithHighExternalViscosity[i]);
							}
							if (CircumferentialNodeFix[i][1]){
								fixY((*itNode),CircumferentialNodeFixWithHighExternalViscosity[i]);

							}
							if (CircumferentialNodeFix[i][2]){
								fixZ((*itNode),CircumferentialNodeFixWithHighExternalViscosity[i]);
							}
						}
					}
				}
			}
		}
	}
	if (BasalNodeFix[0] || BasalNodeFix[1] || BasalNodeFix[2] ){
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ( (*itNode)->tissuePlacement == 0){
				if (BasalNodeFix[0]){
					fixX((*itNode), BasalNodeFixWithExternalViscosity);
				}
				if (BasalNodeFix[1]){
					fixY((*itNode), BasalNodeFixWithExternalViscosity);
				}
				if (BasalNodeFix[2]){
					fixZ((*itNode), BasalNodeFixWithExternalViscosity);
				}
			}
		}
	}
	if (ApicalNodeFix[0] || ApicalNodeFix[1]  || ApicalNodeFix[2]){
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ( (*itNode)->tissuePlacement == 1){
				if (ApicalNodeFix[0]){
					fixX((*itNode), ApicalNodeFixWithExternalViscosity);
				}
				if (ApicalNodeFix[1]){
					fixY((*itNode), ApicalNodeFixWithExternalViscosity);
				}
				if (ApicalNodeFix[2]){
					fixZ((*itNode), ApicalNodeFixWithExternalViscosity);
				}
			}
		}
	}
	if (NotumNodeFix[0] || NotumNodeFix[1]  || NotumNodeFix[2]){
		vector<Node*>::iterator itNode;
		for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ( (*itNode)->tissuePlacement == 0){
				double boundingBoxXMin = boundingBox[0][0];
				double boundingBoxLength = boundingBoxSize[0];
				double relativeXPos = ((*itNode)->Position[0] - boundingBoxXMin)/boundingBoxLength;
				if ( relativeXPos >= notumFixingRange[0] && relativeXPos <= notumFixingRange[1]){
					if (NotumNodeFix[0]){
						fixX((*itNode), NotumNodeFixWithExternalViscosity);
					}
					if (NotumNodeFix[1]){
						fixY((*itNode), NotumNodeFixWithExternalViscosity);
					}
					if (NotumNodeFix[2]){
						fixZ((*itNode), NotumNodeFixWithExternalViscosity);
					}
				}
			}
		}
	}
}

void Simulation::induceClones(){
	//cout<<"inside induce clones"<<endl;
	for (int i=0;i<numberOfClones; ++i){
		double inMicronsX = cloneInformationX[i]*boundingBoxSize[0] + boundingBox[0][0];
		double inMicronsY = cloneInformationY[i]*boundingBoxSize[1] + boundingBox[0][1];
		double inMicronsRadius = cloneInformationR[i];
		//cout<<" clone: "<<i<<" position: "<<inMicronsX<<" "<<inMicronsY<<" radius: "<<inMicronsRadius<<endl;
		double r2 = inMicronsRadius*inMicronsRadius;
		double growthRateORFold = cloneInformationGrowth[i];
		for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			if (!(*itElement)->isECMMimicing && !(*itElement)->IsAblated){
				double* c = new double[3];
				c = (*itElement)->getCentre();
				double dx = inMicronsX - c[0];
				double dy = inMicronsY - c[1];
				double d2 = dx*dx + dy*dy;
				if (d2 < r2){
					//cout<<"mutating element: "<<(*itElement)->Id<<endl;
					if (cloneInformationUsingAbsolueGrowth[i]){
						(*itElement)->mutateElement(0,growthRateORFold); //the mutation is absolute, using an absolute value
					}
					else{
						(*itElement)->mutateElement(growthRateORFold,0); //the mutation is not absolute, using relative values
					}
				}
				delete[] c;
			}
		}
	}
}

void Simulation::initiateElementsByRowAndColumn(int Row, int Column){
	int xinit1 = 0;
	int xinit2 = xinit1+Row+1;
	int xinit3 = 0;
	int xinit4 = xinit1+2*(Row+1)-1;
	int n = Nodes.size() /2.0;
    //initialising the tissue elements:
	for (int ColCount = 0; ColCount < Column; ++ColCount){
		int CurrRowNum = Row + 1 - ColCount;
		for (int RowCount = 0; RowCount<CurrRowNum-1; ++RowCount ){
			int* NodeIds;
			NodeIds = new int[6];

			NodeIds[0] = xinit1+RowCount;
			NodeIds[1] = xinit1+RowCount+1;
			NodeIds[2] = xinit2+RowCount;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			Prism* PrismPnt01;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;

			NodeIds[0] = xinit3+RowCount;
			NodeIds[1] = xinit4+RowCount;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;
		}
		for (int RowCount = 0; RowCount<CurrRowNum-2; ++RowCount ){
			int* NodeIds;
			NodeIds = new int[6];

			NodeIds[0] = xinit2+RowCount;
			NodeIds[1] = xinit1+RowCount+1;
			NodeIds[2] = xinit2+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			Prism* PrismPnt01;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;

			NodeIds[0] = xinit4+RowCount;
			NodeIds[1] = xinit4+RowCount+1;
			NodeIds[2] = xinit3+RowCount+1;
			NodeIds[3] = NodeIds[0] + n;
			NodeIds[4] = NodeIds[1] + n;
			NodeIds[5] = NodeIds[2] + n;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;

		}
		xinit1 = xinit2;
		xinit2 = xinit4 + CurrRowNum-1;
		xinit3 = xinit4;
		xinit4 = xinit2 + CurrRowNum-2;
	}
	cout<<"finalised element initiation"<<endl;
}

void Simulation::calculateSystemCentre(){
	for (int i = 0; i< nNodes; ++i){
		for (int j =0; j<Nodes[i]->nDim; ++j){
			SystemCentre[j] += Nodes[i]->Position[j];
		}
	}
	SystemCentre[0]= SystemCentre[0]/nNodes;
	SystemCentre[1]= SystemCentre[1]/nNodes;
	SystemCentre[2]= SystemCentre[2]/nNodes;
}

void Simulation::addSoftPeriphery(double* fractions){
	double t2 = softDepth*softDepth;
	double midlineFraction = 0.0;
	double currSoftnessFraction = softnessFraction;
	double tissueTypeMultiplier = 1;
	if (softPeripheryBooleans[0] && softPeripheryBooleans[1]){
		midlineFraction = softnessFraction;
	}
	else if (softPeripheryBooleans[0] || softPeripheryBooleans[1]){
		midlineFraction = 1.0 - (1.0-softnessFraction)*0.5;
	}


	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ( (*itNode)->atCircumference ){
			//this node is at circumference, I will calculate the distance of all elements to this node
			//if an element is close enough, update the fraction matrix set above:
			for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
				bool applyToThisElement = true;
				if ( (*itElement)->tissuePlacement == 1 && !softPeripheryBooleans[0]){
					//element is apical, the softness is NOT applied to apical
					applyToThisElement = false;
				}
				if ((*itElement)->tissuePlacement == 0 && !softPeripheryBooleans[1]){
					//element is basal, the softness is NOT applied to apical
					applyToThisElement = false;
				}
				if ((*itElement)->tissueType == 0 && !softPeripheryBooleans[2]){
					//element is columnar, the softness is NOT applied to columnar
					applyToThisElement = false;
				}
				if ((*itElement)->tissueType == 1 && !softPeripheryBooleans[3]){
					//element is peripodial, the softness is NOT applied to peripodial
					applyToThisElement = false;
				}
				if ( applyToThisElement ){
					//scaling the softness fraction if necessary:
					currSoftnessFraction  = softnessFraction;
					if ((*itElement)->tissuePlacement == 2 ){ //element is on the midline
						currSoftnessFraction= midlineFraction;
					}
					//scaling the tissue type multiplier if necessary (important for linker elements:
					if ((*itElement)->tissueType == 2 ){ // element is on the linker zone
						if (softPeripheryBooleans[2] && softPeripheryBooleans[3]){
							//softness is applied to both peripodial and columnar zones, the linkers should not be scaling anything
							tissueTypeMultiplier = 1.0;
						}else if (softPeripheryBooleans[2] ){
							//softness only applied to columnar layer:
							tissueTypeMultiplier = (*itElement)->getColumnarness();
						}
						else if (softPeripheryBooleans[3] ){
							//softness only applied to peripodial layer:
							tissueTypeMultiplier = (*itElement)->getPeripodialness();
						}
					}
					else{ //element is not linker
						tissueTypeMultiplier = 1;
					}
					currSoftnessFraction  = 1.0 - (1.0-currSoftnessFraction)*tissueTypeMultiplier;
					double *c = (*itElement)->getCentre();
					double dx = c[0]-(*itNode)->Position[0];
					double dy = c[1]-(*itNode)->Position[1];
					//double dz = c[2]-(*itNode1)->Position[2];
					double d = dx*dx + dy*dy;
					if (d<t2){
						d = pow(d,0.5);
						double f = currSoftnessFraction + (1.0 - currSoftnessFraction)*d/softDepth;
						if ((softnessFraction< 1.0 && f < fractions[(*itElement)->Id]) ||(softnessFraction> 1.0 && f > fractions[(*itElement)->Id])){
							fractions[(*itElement)->Id] = f;
						}
					}
					delete[] c;
				}
			}
		}
	}
}

void Simulation::assignPhysicalParameters(){

	double* fractions;
	fractions = new double[(const int) nElements];
	for (int i=0; i<nElements; ++i){
		fractions[i] = 1.0;
	}
	if(softPeriphery){
		addSoftPeriphery(fractions);
	}
	//now I have softness fraction table, I can scale the element properties:
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double r = (rand() % 200) / 100.0;	//random number between 0.00 and 2.00
		r = r - 1.0; 						//random number between -1.00 and 1.00
		float noise1 = r*noiseOnPysProp[0];	//percent noise on current element
		r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise2 = r*noiseOnPysProp[1];
		if ((*itElement)->tissueType == 0){ //Element is on the columnar layer
			double currEApical 	= fractions[(*itElement)->Id] * EApical*(1 + noise1/100.0);
			double currEBasal	= fractions[(*itElement)->Id] * EBasal*(1 + noise1/100.0);
			double currEMid		= fractions[(*itElement)->Id] * EMid*(1 + noise1/100.0);
			double currEECM		= fractions[(*itElement)->Id] * EColumnarECM*(1 + noise1/100.0);
			double currPoisson = poisson*(1 + noise2/100);
			if((*itElement)->isECMMimicing){
				//the Poisson ratio is zero so that the ECM layer will not thin!
				currPoisson = 0;
			}
			(*itElement)->setElasticProperties(currEApical,currEBasal,currEMid,currEECM,currPoisson);
			(*itElement)->setViscosity(discProperApicalViscosity,discProperBasalViscosity,discProperMidlineViscosity);
		}
		else if ((*itElement)->tissueType == 1){ //Element is on the peripodial membrane
			double currEApical = fractions[(*itElement)->Id] * PeripodialElasticity*(1 + noise1/100.0);
			double currEBasal = currEApical;
			double currEMid = currEApical;
			double currEECM = fractions[(*itElement)->Id] * EPeripodialECM*(1 + noise1/100.0);
			double currPoisson = poisson*(1 + noise2/100);
			/*if(thereIsExplicitECM){
				//The input file does not define apical and basal stiffness separately
				//for each element. If I have explicit ECM, I will change the basal ECM
				//stiffness such that it will be the same as basal of the columnar layer,
				//therefore the ECM.
				currEBasal = fractions[(*itElement)->Id] * EBasal*(1 + noise1/100.0);
				currEMid = fractions[(*itElement)->Id] * EMid*(1 + noise1/100.0);
			}*/
			if((*itElement)->isECMMimicing){
				//the Poisson ratio is zero so that the ECM layer will not thin!
				currPoisson = 0;
			}
			(*itElement)->setElasticProperties(currEApical,currEBasal,currEMid,currEECM, currPoisson);
			(*itElement)->setViscosity(peripodialApicalViscosity,peripodialBasalViscosity,peripodialMidlineViscosity);
		}
		else if ((*itElement)->tissueType == 2 ){ //Element is on the linker Zone,
			//The elastic properties of the linker zone ECM are always based on the
			//peripodialness factor.
			double currPeripodialECME = fractions[(*itElement)->Id] * EPeripodialECM*(1 + noise1/100.0);
			double currColumnarECME = fractions[(*itElement)->Id] * EColumnarECM*(1 + noise1/100.0);
			double periWeight 		= (*itElement)->getPeripodialness();
			double colWeight = (*itElement)->getColumnarness();
			double currEECM		= colWeight * currColumnarECME + periWeight * currPeripodialECME;
			if (BaseLinkerZoneParametersOnPeripodialness){
				//I will weight the values:
				double currPeripodialE 	= fractions[(*itElement)->Id] * PeripodialElasticity * (1 + noise1/100.0);
				double currEApical 		= fractions[(*itElement)->Id] * EApical * (1 + noise1/100.0);
				double currEBasal 		= fractions[(*itElement)->Id] * EBasal * (1 + noise1/100.0);
				double currEMid 		= fractions[(*itElement)->Id] * EMid * (1 + noise1/100.0);
				double currPoisson      = poisson*(1 + noise2/100);
				currEApical = colWeight * currEApical + periWeight * currPeripodialE;
				currEBasal  = colWeight * currEBasal  + periWeight * currPeripodialE;
				currEMid    = colWeight * currEMid    + periWeight * currPeripodialE;
				(*itElement)->setElasticProperties(currEApical,currEBasal,currEMid,currEECM, currPoisson);
				double currViscApical  = colWeight * discProperApicalViscosity  + periWeight * peripodialApicalViscosity;
				double currViscBasal   = colWeight * discProperBasalViscosity   + periWeight * peripodialBasalViscosity;
				(*itElement)->setViscosity(currViscApical,currViscBasal); //There is no midline in the linker zone.
			}
			else{
				//I have the inputs provided:
				double currEApical 	= fractions[(*itElement)->Id] * LinkerZoneApicalElasticity*(1 + noise1/100.0);
				double currEBasal	= fractions[(*itElement)->Id] * LinkerZoneBasalYoungsModulus*(1 + noise1/100.0);
				double currEMid		= fractions[(*itElement)->Id] * 0.5 * (LinkerZoneApicalElasticity+LinkerZoneBasalYoungsModulus)*(1 + noise1/100.0);
				double currPoisson  = poisson*(1 + noise2/100);


				if((*itElement)->isECMMimicing){
					//the Poisson ratio is zero so that the ECM layer will not thin!
					currPoisson = 0;
				}
				(*itElement)->setElasticProperties(currEApical,currEBasal,currEMid,currEECM, currPoisson);
				(*itElement)->setViscosity(linkerZoneApicalViscosity,linkerZoneBasalViscosity,linkerZoneMidlineViscosity);
			}
		}
	}
	vector<Node*>::iterator itNode;
	for(itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		double r = (rand() % 200) / 100.0;
		r = r - 1.0;
		float noise3 = r*noiseOnPysProp[1];
		noise3 = (1 + noise3/100.0);
		if ((*itNode)->tissueType == 2 ){ //linker, the viscosity can be averaged, or based on indivudual set of parameters:
			if (BaseLinkerZoneParametersOnPeripodialness){
				//take average of the two:
				double currExternalViscApical = 0.5*(externalViscosityDPApical + externalViscosityPMApical);
				double currExternalViscBasal  = 0.5*(externalViscosityDPBasal  + externalViscosityPMBasal);
				(*itNode)->setExternalViscosity(currExternalViscApical,currExternalViscBasal, extendExternalViscosityToInnerTissue);
				//(*itNode)->setExternalViscosity(externalViscosityDPApical*noise3, externalViscosityDPBasal*noise3, externalViscosityPMApical*noise3, externalViscosityPMBasal*noise3);
			}
			else{
				(*itNode)->setExternalViscosity(externalViscosityLZApical*noise3, externalViscosityLZBasal*noise3, extendExternalViscosityToInnerTissue);
				//(*itNode)->setExternalViscosity(externalViscosityLZApical*noise3, externalViscosityLZBasal*noise3, externalViscosityLZApical*noise3, externalViscosityLZBasal*noise3);
			}
		}else if ((*itNode)->tissueType == 0){ //disc proper
			(*itNode)->setExternalViscosity(externalViscosityDPApical*noise3, externalViscosityDPBasal*noise3, extendExternalViscosityToInnerTissue);
			//(*itNode)->setExternalViscosity(externalViscosityDPApical*noise3, externalViscosityDPBasal*noise3, externalViscosityPMApical*noise3, externalViscosityPMBasal*noise3);
		}else if ((*itNode)->tissueType == 1){
			(*itNode)->setExternalViscosity(externalViscosityPMApical*noise3, externalViscosityPMBasal*noise3, extendExternalViscosityToInnerTissue);
		}
	}
	delete[] fractions;
}

void Simulation::checkForZeroExternalViscosity(){
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		for (int i=0; i<3; ++i){
			if ((*itNode)->externalViscosity[i] > 0){
				zeroExternalViscosity[i] = false;
			}
		}
	}
    if (zeroExternalViscosity[0] || zeroExternalViscosity[1] || zeroExternalViscosity[2]){
    	//At least one of the dimensions have zero viscosity
    	if(nNodes < 3) {
    		cerr<<"Cannot run zero viscosity simulation with less than 3 nodes, run simulation at your own risk!"<<endl;
    	}
        //fix x,y,z for 0
        /*Nodes[0]->FixedPos[0] = true;
        Nodes[0]->FixedPos[1] = true;
        Nodes[0]->FixedPos[2] = true;
        //fix x and z for 1
        Nodes[1]->FixedPos[0] = true;
        Nodes[1]->FixedPos[2] = true;
        //fix z for 2
        Nodes[2]->FixedPos[2] = true;

        //Nodes[4]->FixedPos[0] = true;
        //Nodes[4]->FixedPos[1] = true;
        //Nodes[4]->FixedPos[2] = true;

        */
        if (!symmetricY){
        	//TO DO!! In the node fixing options, there should be at least 2 different axes fixed.
        	//If there are any node fixes, on the surfaces in z, then I do not need to do fix z:
        	bool sufficientXFix = false;
        	bool sufficientYFix = false;
        	bool sufficientZFix = false;
        	for (int a=0; a<3; ++a){
        		//if zeroExternalViscosity[0] is false, it means there is already  sufficient x fix
        		//OR, if I spotted sufficientXFix in previous nodes, it persists,
        		//OR, there is circumferential fix on x dimension
        		sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || CircumferentialNodeFix[a][0];
        		sufficientYFix = sufficientYFix || !zeroExternalViscosity[1] || CircumferentialNodeFix[a][1];
        		sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || CircumferentialNodeFix[a][2];
        	}
			//if there are no z fixe on the circumference, check the whole surfaces:
        	if (!sufficientXFix){
        		sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || BasalNodeFix[0] ;
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || ApicalNodeFix[0] ;
			}
        	if (!sufficientYFix){
        		sufficientYFix = sufficientYFix || !zeroExternalViscosity[1] || BasalNodeFix[1] ;
        		sufficientYFix = sufficientYFix || !zeroExternalViscosity[1] || ApicalNodeFix[1] ;
        	}
        	if (!sufficientZFix){
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || BasalNodeFix[2];
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || ApicalNodeFix[2];
			}

        	if (!sufficientXFix){
        		Nodes[ventralTipIndex]->FixedPos[0] = true;
        	}
        	if (!sufficientYFix){
        		Nodes[ventralTipIndex]->FixedPos[1] = true;
        	}
        	if (!sufficientZFix){
        		Nodes[ventralTipIndex]->FixedPos[2] = true;
        	}
        	if (!sufficientYFix){
        		Nodes[dorsalTipIndex]->FixedPos[1] = true;
        	}
			if (!sufficientZFix){
				Nodes[dorsalTipIndex]->FixedPos[2] = true;
			}
			if (!sufficientZFix){
				//if there is symmetricity, then the mid-line nodes will be fixed in y, and I do not need to fix the third node.
				// in fact, fixing the position of the third node in z will cause problems.
				if (dorsalTipIndex != 1 && 	ventralTipIndex!= 1 ){
					Nodes[1]->FixedPos[2] = true;
				}
				else if (dorsalTipIndex != 2 && ventralTipIndex!= 2 ){
					Nodes[2]->FixedPos[2] = true;
				}
				else if (dorsalTipIndex != 3 && ventralTipIndex!= 3 ){
					Nodes[3]->FixedPos[2] = true;
				}
			}
        }
        else{
        	//There is symmetricity, then there is y fix, if there are any node fixes, on the surfaces in z, then I do not need to do fix z:
        	bool sufficientZFix = false;
        	for (int a=0; a<3; ++a){
        		sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || CircumferentialNodeFix[a][2];
			}
        	//if there are no z fixe on the circumferenece, check the whole surfaces:
			if (!sufficientZFix){
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || BasalNodeFix[2];
				sufficientZFix = sufficientZFix || !zeroExternalViscosity[2] || ApicalNodeFix[2];
			}
        	bool sufficientXFix = false;
			for (int a=0; a<3; ++a){
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || CircumferentialNodeFix[a][0] ;
			}
			//if there are no x fixes on the circumferenece, check the whole surfaces:
			if (!sufficientXFix){
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || BasalNodeFix[0] ;
				sufficientXFix = sufficientXFix || !zeroExternalViscosity[0] || ApicalNodeFix[0] ;
			}
			//there is no circumference or surface fixing of suffieint nature, then I should fix the nodes:
			if (!sufficientXFix){
              	Nodes[ventralTipIndex]->FixedPos[0] = true;
			}
			Nodes[ventralTipIndex]->FixedPos[1] = true;
			if (!sufficientZFix){
				Nodes[ventralTipIndex]->FixedPos[2] = true;
				Nodes[dorsalTipIndex]->FixedPos[2] = true;
			}
			//cout<<"node fixing sufficiency, sufficientXFix: "<<sufficientXFix<<" sufficientZFix: "<<sufficientZFix<<endl;
/*
				vector<ShapeBase*>::iterator itElement;
				bool foundElement = false;
				int apicalId = 0;
				//find the apical node corresponding to the basal ventral tip:
				for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
					bool IsBasalOwner = (*itElement)->IsThisNodeMyBasal(ventralTipIndex);
					if (IsBasalOwner){
						foundElement = true;
						apicalId = (*itElement)->getCorrecpondingApical(ventralTipIndex); //have the next node
						break;
					}
				}
				if (foundElement){
					Nodes[apicalId]->FixedPos[2] = true;
				}
				else{
					cerr<<"Cannot run zero viscosity simulation, could not find the apical node corresponding to the basal ventral tip, run simulation at your own risk!"<<endl;
				}
*/
        }
    }
}

void Simulation::addCurvatureToColumnar(double h){
	//find the tips:
	double l1 = -1000, l2 = 1000, l3 = -1000;
	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->Position[0] > l1 ){
			l1 = (*itNode)->Position[0]; //displacement from origin to positive x tip
		}
		if ((*itNode)->Position[0] < l2 ){
			l2 = (*itNode)->Position[0]; //displacement from origin to negative x tip
		}
		if ((*itNode)->Position[1] > l3 ){
			l3 = (*itNode)->Position[1];	//displacement from origin to positive y tip
		}
	}
	if (symmetricX){
		l2 = (-1.0)*l1;	//if there is symmetricity in x, then the negative tip will be at zero. This is not correct and the symmetricity should mean the negative displacement will be -l1
	}
	l2 *= -1.0;	//converting the displacement to distance, this is the only negative tip, the rest are already positive
	cout<<"l1: "<<l1 <<" l2: "<<l2<<" l3: "<<l3<<endl;
	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		double x = (*itNode)->Position[0];
		double y = (*itNode)->Position[1];
		double z = (*itNode)->Position[2];
		double a = l1;
		double c = l3;
		if ((*itNode)->Position[0] < 0){
			a = l2;
		}
		double value = (1 - x*x/a/a - y*y/c/c);
		if (value>0){
			double offset = pow((1 - x*x/a/a - y*y/c/c)*h*h,0.5);
			if (h<0){
				offset *= (-1.0);
			}
			if (thereIsPeripodialMembrane && (*itNode)->tissueType != 0){ //node is not columnar
				//there is peripodial membrane, any node above the mid-line of the lumen should be curved in the opposite direction of the columnar layer:
				//But I have a direction choice, if the columnar is curving down, the peripodial shoul curve up
				//If the columnar is curing up, the peripodial should still curve up!

				double heightThreshold = TissueHeight + lumenHeight/2.0;
				if ((*itNode)->Position[2] > heightThreshold) {
					if (offset > 0){
						offset *= (-1.0);
					}
				}
			}
			(*itNode)->Position[2] -= offset;
			if((*itNode)->Id ==8){
				cout<<" a: "<<a<<" c: "<<c<<" x "<<x<<" y "<< y<<" z "<<z<<" value: "<<value<<" offset: "<<offset<<endl;
				cout<<"Node[8] position: "<<(*itNode)->Position[0]<<" "<<(*itNode)->Position[1]<<" "<<(*itNode)->Position[2]<<endl;
			}
		}
	}
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updatePositions(Nodes);
		(*itElement)->updateReferencePositionsToCurentShape();
	}
}

void Simulation::thickecECM(){
	double bbMinX = boundingBox[0][0];
	double bbMaxX = boundingBox[1][0];
	double bbXSize = bbMaxX-bbMinX;
	double thickenDz = 0;//4.17; //make selected regions 3 micron thicker
	double thinnerDz = 2.0;
	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissuePlacement == 0){
			//only basal nodes
			double x = (*itNode)->Position[0];
			double relativeX = (x-bbMinX)/bbXSize;
			if (relativeX < 0.4 || relativeX > 0.65){
				//<0.4 notum, >0.65 pouch
				(*itNode)->Position[2] -= thickenDz;
			}
			//else{
				(*itNode)->Position[2] += thinnerDz;
			//}
		}
	}
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updatePositions(Nodes);
		(*itElement)->updateReferencePositionsToCurentShape();
	}
}

void Simulation::fixNode0InPosition(double x, double y, double z){
	double dx = Nodes[0]->Position[0]-x;
	double dy = Nodes[0]->Position[1]-y;
	double dz = Nodes[0]->Position[2]-z;
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		(*itNode)->Position[0] -=dx;
		(*itNode)->Position[1] -=dy;
		(*itNode)->Position[2] -=dz;
	}
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updatePositions(Nodes);
	}
}

void Simulation::manualPerturbationToInitialSetup(bool deform, bool rotate){
	if(timestep==1){
		 for (vector<Node*>::iterator  itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			 //if ((*itNode)->Id <3 ){
				//fixX((*itNode),0);
				//fixY((*itNode),0);
			 //}
			/*if ((*itNode)->Id != 3 && (*itNode)->Id != 10){
				//fixAllD((*itNode),0);
			}
			else{
				//fixZ((*itNode),0);
				//fixY((*itNode),0);
				(*itNode)->Position[0] -= 5;
			}*/
			//fixAllD((*itNode),0);
		}
		//laserAblateTissueType(1);
		//laserAblate(0.0, 0.0, 0.5);
        double scaleX = 2.0;
        double scaleY = 1.0;
        double scaleZ = 1.0;

        double PI = 3.14159265359;
        double tetX = 0 *PI/180.0;
        double tetY = 45 *PI/180.0;
        double tetZ = 0 *PI/180.0;
    	int axisToRotateOn = 1; //0: rotate around x axis, 1: around y-axis, and 2: around z-axis.
        if(deform){
        	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
        		(*itNode)->Position[0] *=scaleX;
        		(*itNode)->Position[1] *=scaleY;
                (*itNode)->Position[2] *=scaleZ;
                //if ((*itNode)->tissuePlacement ==0 ){
                //	(*itNode)->Position[2]  += 2.0;
                	//double fold = -1.0/(scaleX*scaleY*scaleZ) -1;
                	//(*itNode)->Position[2]  = 12.5/3.0 *scaleZ *fold;
				//}
            }
        }
        if (rotate){
        	cerr<<"Rotating system"<<endl;
            double R[3][3] =  {{1,0,0},{0,1,0},{0,0,1}};
            if (axisToRotateOn == 0){
            	double c = cos(tetX);
            	double s = sin(tetX);
            	//Rx = {{1,0,0},{0,c,-s},{0,s,c}};
            	R[1][1] = c;
				R[1][2] = -s;
				R[2][1] = s;
				R[2][2] = c;
            }
            else if(axisToRotateOn == 1){
            	double c = cos(tetY);
            	double s = sin(tetY);
            	//Ry = {{c,0,s},{0,1,0},{-s,0,c}};
            	R[0][0] = c;
				R[0][2] = s;
				R[2][0] = -s;
				R[2][2] = c;
            }
            else if(axisToRotateOn == 2){
            	double c = cos(tetZ);
            	double s = sin(tetZ);
            	//Rz = {{c,-s,0},{s,c,0},{0,0,1}};
            	R[0][0] = c;
				R[0][1] = -s;
				R[1][0] = s;
				R[1][1] = c;
            }
        	for (vector<Node*>::iterator  itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
                double x = (*itNode)->Position[0]*R[0][0] + (*itNode)->Position[1]*R[0][1] + (*itNode)->Position[2]*R[0][2];
                double y = (*itNode)->Position[0]*R[1][0] + (*itNode)->Position[1]*R[1][1] + (*itNode)->Position[2]*R[1][2];
                double z = (*itNode)->Position[0]*R[2][0] + (*itNode)->Position[1]*R[2][1] + (*itNode)->Position[2]*R[2][2];
                (*itNode)->Position[0]=x;
                (*itNode)->Position[1]=y;
                (*itNode)->Position[2]=z;
            }
        }
    	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    		(*itElement)->updatePositions(Nodes);
        }
    }
}

void Simulation::updateGrowthRotationMatrices(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        if (!(*itElement)->IsAblated){
            (*itElement)->CalculateGrowthRotationByF();
        }
    }
}

void Simulation::calculateDiscretisationLayers(double &hColumnar, int& LumenHeightDiscretisationLayers, double &hLumen, double &peripodialHeight, int& peripodialHeightDiscretisationLayers, double& hPeripodial){
	//The average columnar layer element heigth would be TissueHeight divided by the number of element layers:
	hColumnar = TissueHeight/TissueHeightDiscretisationLayers;
	//Calculating the height of the lumen:
	lumenHeight = TissueHeight*lumenHeightScale;
	//calculating how many layer I need for this lumen:
	LumenHeightDiscretisationLayers = ceil(lumenHeight/hColumnar);
	//calculating the actual element height needed to produce the total height with the calculated layr number:
	hLumen = lumenHeight/LumenHeightDiscretisationLayers;
	//calculating the height of the peripodial membrane:
	peripodialHeight = PeripodialThicnessScale*TissueHeight;
	//calculating how many layers I need for this peripodial membrane:
	peripodialHeightDiscretisationLayers = ceil(peripodialHeight/hColumnar);
	//calculating the actual element height needed to produce the total height with the calculated layr number:
	hPeripodial = peripodialHeight/peripodialHeightDiscretisationLayers;
	//cout<<"Tissue height: "<<TissueHeight<<" columnar element height: "<<hColumnar<<" TissueHeightDiscretisationLayers: "<<TissueHeightDiscretisationLayers<<endl;
	//cout<<"Lumen height:  "<<lumenHeight<<" height of one lumen element: "<<hLumen<<" LumenHeightDiscretisationLayers: "<<LumenHeightDiscretisationLayers<<endl;
	//cout<<"Peripodial height:  "<<peripodialHeight<<" height of one peripodial element: "<<hPeripodial<<" LumenHeightDiscretisationLayers: "<<peripodialHeightDiscretisationLayers<<endl;
}

void Simulation::fillColumnarBasedNodeList(vector< vector<int> > &ColumnarBasedNodeArray, vector <int> &ColumnarCircumferencialNodeList){
	int nCircumference = ColumnarCircumferencialNodeList.size();
	for (int i =0; i<nCircumference; ++i){
    	int currNodeId = ColumnarBasedNodeArray[i][0];
		bool foundElement = true;
		bool finishedTissueThicness  = false;
		while(!finishedTissueThicness && foundElement){ //while the node is not apical, and I could find the next element
			foundElement = false;
			vector<ShapeBase*>::iterator itElement;
			for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
				bool IsBasalOwner = (*itElement)->IsThisNodeMyBasal(currNodeId);
				if (IsBasalOwner){
					foundElement = true;
					break;
				}
			}
			//found the current node as basal node of an element
			//get the corresponding apical node, record in the list:
			currNodeId = (*itElement)->getCorrecpondingApical(currNodeId); //have the next node
			ColumnarBasedNodeArray[i].push_back(currNodeId);
			if (Nodes[currNodeId]->tissuePlacement == 1) {
				//Added the apical node now, I can stop:
				finishedTissueThicness = true;
			}
		}
    }
}

void Simulation::addNodesForPeripodialOnColumnarCircumference (vector< vector<int> > &ColumnarBasedNodeArray, int LumenHeightDiscretisationLayers, double hLumen, int peripodialHeightDiscretisationLayers, double hPeripodial ){
	//add nodes for the lumen and peripodial membrane on top of the columnar nodes:
	int nCircumference = ColumnarBasedNodeArray.size();
	for (int i=0; i<nCircumference; ++i){
		int nodeId0 = ColumnarBasedNodeArray[i][TissueHeightDiscretisationLayers];
		double* pos;
		pos = new double[3];
		for (int j=0; j<Nodes[nodeId0]->nDim; j++){
			pos[j] = Nodes[nodeId0]->Position[j];
		}
		//adding nodes for lumen:
		for (int j=1; j<LumenHeightDiscretisationLayers+1; ++j){
			pos[2] += hLumen;
			int newNodeId = Nodes.size();
			Node* tmp_nd = new Node(newNodeId, 3, pos, 2, 1); //Tissue placement is midlayer (2), tissue type is peripodial membrane (1)
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			ColumnarBasedNodeArray[i].push_back(newNodeId);
		}
		//The last node should also be apical, at it is the bottom of the peripodial membrane, looking into the lumen, change its placement:
		Nodes[Nodes.size()-1]->tissuePlacement = 1; //made tissueplacement apical
		//adding nodes for the peripodial membrane:
		for (int j=1; j<peripodialHeightDiscretisationLayers+1; ++j){
			pos[2] += hPeripodial;
			int newNodeId = Nodes.size();
			Node* tmp_nd = new Node(newNodeId, 3, pos, 2, 1); //Tissue placement is midlayer (2), tissue type is peripodial membrane (1)
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			ColumnarBasedNodeArray[i].push_back(newNodeId);
		}
		//The last node should also be basal, at it is the top of the peripodial membrane, change its placement:
		Nodes[nNodes-1]->tissuePlacement = 0; //made tissueplacement basal
		delete pos;
	}
}
void Simulation::calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, int nodeId2, double* pos, double sideThickness){
	//cout<<"nodeIDs: "<<nodeId0<<" "<<nodeId1<<" "<<nodeId2<<" sideThickness: "<<sideThickness<<endl;
	double* vec1;
	vec1 = new double[3];
	if (symmetricY && nodeId1 == -100){
		vec1[0] = -1.0;
		vec1[1] =  0.0;
		vec1[2] =  0.0;
	}
	else if (symmetricY && nodeId2 == -100){
		vec1[0] =  1.0;
		vec1[1] =  0.0;
		vec1[2] =  0.0;
	}
	else{
		double* vec2;
		vec2 = new double[3];
		for (int j=0; j<Nodes[nodeId0]->nDim; j++){
			vec1[j] = (Nodes[nodeId0]->Position[j] - Nodes[nodeId1]->Position[j]);
			vec2[j] = (Nodes[nodeId0]->Position[j] - Nodes[nodeId2]->Position[j]);
		}
		double dummy = Elements[0]->normaliseVector3D(vec1);
		dummy = Elements[0]->normaliseVector3D(vec2);
		vec1[0] += vec2[0];
		vec1[1] += vec2[1];
		vec1[2] += vec2[2];
		double vec1Mag2 = vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2];
		if (vec1Mag2 < 1E-5){
			//the two nodes are linear, the resulting vector is of zero length.
			//I will rotate one of the vectors 90 degrees and it will be pointing in the correct orientation, the actual direction can be fixed in the below loop:
			//this vector is already normalised, I am skipping the normalisation step
			vec1[0] = -1.0*vec2[1];
			vec1[1] = vec2[0];
		}
		else{
			double dummy = Elements[0]->normaliseVector3D(vec1);
		}
		//now I have the vector to point out from the base node 0. BUT, this will point outwards only if the tissue curvature is convex at all points
		//I need to check if it actually is pointing out, as the experimentally driven tissues can be concave at points.
		//the cross produc of vec[2] and the vector to the cell centre should have the opposite sign with the corss product of my orientation vector and vector 2.
		double* vecCentre;
		vecCentre = new double[3];
		vecCentre[0] =  SystemCentre[0] - Nodes[nodeId0]->Position[0];
		vecCentre[1] =  SystemCentre[1] - Nodes[nodeId0]->Position[1];
		vecCentre[2] =  SystemCentre[2] - Nodes[nodeId0]->Position[2];
		dummy = Elements[0]->normaliseVector3D(vecCentre);
		double* cross1;
		cross1 = new double[3];
		Elements[0]->crossProduct3D(vec2,vecCentre,cross1);
		dummy = Elements[0]->normaliseVector3D(cross1);
		double* cross2;
		cross2 = new double[3];
		Elements[0]->crossProduct3D(vec2,vec1,cross2);
		dummy = Elements[0]->normaliseVector3D(cross2);
		double dotp = Elements[0]->dotProduct3D(cross1,cross2);
		if (dotp >0 ){
			//the vectors are pointing to the same direction! Need to rotate vec1 180 degrees:
			vec1[0] *= -1.0;
			vec1[1] *= -1.0;
			vec1[2] *= -1.0;
		}
		delete[] vecCentre;
		delete[] cross1;
		delete[] cross2;
		delete[] vec2;
	}
	pos[0] = Nodes[nodeId0]->Position[0] + vec1[0]*sideThickness;
	pos[1] = Nodes[nodeId0]->Position[1] + vec1[1]*sideThickness;
	pos[2] = Nodes[nodeId0]->Position[2];
	delete[] vec1;

}


void Simulation::calculateNewNodePosForPeripodialNodeAddition(int nodeId0, int nodeId1, double* pos, double sideThickness){
	double* vec0;
	vec0 = new double[3];
	double midpoint[3] = {0,0,0};
	for (int j=0; j<Nodes[nodeId0]->nDim; j++){
		vec0[j] = (Nodes[nodeId1]->Position[j] - Nodes[nodeId0]->Position[j])/2.0;
		midpoint[j] = Nodes[nodeId0]->Position[j] + vec0[j];
	}
	double dummy = Elements[0]->normaliseVector3D(vec0);
	//The list is sorted counter-cock-wise, to point out, I will rotate normalised vector v0 -90 degrees on z axis:
	// (x,y,z) -> (y,-x,z);
	// then I will add this vector to the calculated mid point to gt the new node's position.
	pos[0] = midpoint[0] + vec0[1]*sideThickness;
	pos[1] = midpoint[1] - vec0[0]*sideThickness;
	pos[2] = midpoint[2];
	delete[] vec0;
}

void Simulation::addNodesForPeripodialOnOuterCircumference (vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray, double hColumnar, int LumenHeightDiscretisationLayers, double hLumen, int peripodialHeightDiscretisationLayers, double hPeripodial ){
	double peripodialSideConnectionThickness =PeripodialLateralThicnessScale*TissueHeight; //in microns
	double avrSide=0.0, dummy =0.0;
	getAverageSideLength(dummy,avrSide);	//first term will get you the average side length of the peripodial membrane elements, second is the columnar elements
	if (avrSide/peripodialSideConnectionThickness > 5 || avrSide/peripodialSideConnectionThickness< 0.2 ){
		cerr<<"WARNING, the lateral connection thickness between the peripodial membrane and the columnar layer is too different than average element size (more than 5 fold diference)"<<endl;
	}
	//Now I need the average side of an element, to add new nodes accordingly:
	int nCircumference = ColumnarBasedNodeArray.size();

	for (int i=0; i<nCircumference; ++i){
		//cout<<"at node in the list: "<<i<<endl;
		//adding 2 point based node:
/*		int nodeId0 = ColumnarBasedNodeArray[i][0];
		int nodeId1;
		if( i == nCircumference - 1){
			nodeId1 = ColumnarBasedNodeArray[0][0];
		}
		else{
			nodeId1 = ColumnarBasedNodeArray[i+1][0];
		}
		double* pos = new double[3];
		calculateNewNodePosForPeripodialNodeAddition(nodeId0, nodeId1, pos, peripodialSideConnectionThickness);
		//cout<<" calculated pos : "<<pos[0] <<" "<<pos[1]<<" "<<pos[2]<<endl;
		//Adding the array of new nodes:
*/
		//adding 3 point based node:
		int nodeId0 = ColumnarBasedNodeArray[i][0];
		int nodeId1;
		int nodeId2;
		int baseIndex0 = i;
		if( i == nCircumference - 1){
			if (symmetricY){
				//the node is the end tip of a tissue with symmetric y. I should not connect it in a loop, it should add
				// the node to +x direction.
				// the node position calculation function will catch this as a flag, and will not calculate
				nodeId1 = -100;
			}
			else{
				nodeId1 = ColumnarBasedNodeArray[0][0];
			}
		}
		else{
			nodeId1 = ColumnarBasedNodeArray[i+1][0];
		}
		if ( i == 0 ){
			if (symmetricY){
				//the node is the beginning tip of a tissue with symmetric y. I should not connect it in a loop, it should add
				// the node to -x direction.
				// the node position calculation function will catch this as a flag, and will not calculate
				nodeId2 = -100;
			}
			else{
				nodeId2 = ColumnarBasedNodeArray[nCircumference-1][0];
			}
		}
		else{
			nodeId2 = ColumnarBasedNodeArray[i-1][0];
		}
		double* pos = new double[3];
		calculateNewNodePosForPeripodialNodeAddition(nodeId0, nodeId1, nodeId2, pos, peripodialSideConnectionThickness);

		//cout<<" calculated pos : "<<pos[0] <<" "<<pos[1]<<" "<<pos[2]<<endl;
		//Adding the array of new nodes:

		//adding the base:
		int newNodeId = Nodes.size();
		Node* tmp_nd = new Node(newNodeId, 3, pos, 0, 2); //Tissue placement basal (0), tissue type is linker zone (2)
		Nodes.push_back(tmp_nd);
		nNodes = Nodes.size();
		OuterNodeArray[i].push_back(newNodeId);
		//adding the nodes for the columnar layer:
		for (int j=1; j<TissueHeightDiscretisationLayers+1; ++j){
			//pos[2] += hColumnar;
			pos[2]  = Nodes[ColumnarBasedNodeArray[baseIndex0][j]]->Position[2];
			//cout<<" pos for columnar aligned new node: "<<pos[0] <<" "<<pos[1]<<" "<<pos[2]<<" hColumnar: "<<hColumnar<<endl;

			int newNodeId = Nodes.size();
			Node* tmp_nd = new Node(newNodeId, 3, pos, 2, 2); //Tissue placement is midlayer (2), tissue type is linker zone (2)
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			OuterNodeArray[i].push_back(newNodeId);
		}
		//adding nodes for lumen:
		for (int j=1; j<LumenHeightDiscretisationLayers+1; ++j){
			pos[2] += hLumen;
			int newNodeId = Nodes.size();
			Node* tmp_nd = new Node(newNodeId, 3, pos, 2, 2); //Tissue placement is midlayer (2), tissue type is linker zone (2)
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			OuterNodeArray[i].push_back(newNodeId);
		}
		//adding nodes for the peripodial membrane:
		for (int j=1; j<peripodialHeightDiscretisationLayers+1; ++j){
			pos[2] += hPeripodial;
			int newNodeId = Nodes.size();
			Node* tmp_nd = new Node(newNodeId, 3, pos, 2, 2); //Tissue placement is midlayer (2), tissue type is linker zone (2)
			Nodes.push_back(tmp_nd);
			nNodes = Nodes.size();
			OuterNodeArray[i].push_back(newNodeId);
		}
		//The last node should also be basal, at it is the top of the peripodial membrane, change its placement:
		Nodes[nNodes-1]->tissuePlacement = 0; //made tissueplacement basal
		delete[] pos;
    }
}

void Simulation::addLateralPeripodialElements(int LumenHeightDiscretisationLayers, int peripodialHeightDiscretisationLayers, vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray){
    // I need to add the elements:
    // Two elements are added for each element on the side:
    // First triangle base will be New node 0, oldNode1, oldnode0, top of the element will be read from the node id stacks
    // Second triangle base will be: New node 0, new node 1, old node 1
	int totalLayers = TissueHeightDiscretisationLayers+LumenHeightDiscretisationLayers+peripodialHeightDiscretisationLayers;
	double peripodialnessFractionStep = 1.0 / (double) (LumenHeightDiscretisationLayers+1.0); //this is the step increment in periopdiallness weight at each element of lumen side. If there is 1 layer, it should be 50% peripodial, 50% columnar etc, if there are two, it should be 33%, 2x33% etc.
	int nCircumference = ColumnarBasedNodeArray.size();
	int nLoop = nCircumference;
	if (symmetricY){
		//in symmetric setup, the last node is not connected to the first node, we simply ignore that step
		nLoop--;
	}
	for (int i=0;i<nLoop; ++i){
    	for (int j=0;j<totalLayers; ++j){
    		//calculating the fraction of peripodialness these two elements should have:
    		//counting which layer in the lumen they are:
    		double peripodialWeight = 1.0; //
    		if (j < TissueHeightDiscretisationLayers) {
    			//the currently added elements are still aligned with columnar, below the lumen.
    			peripodialWeight = 0.0;
    		}
    		else if (j<TissueHeightDiscretisationLayers+LumenHeightDiscretisationLayers){
    			//the currently added elements are at the lumen, below peripodial membrane:
    			peripodialWeight = (j-TissueHeightDiscretisationLayers + 1) * peripodialnessFractionStep;
    		}
    		// if (j>=TissueHeightDiscretisationLayers+LumenHeightDiscretisationLayers), then the added elements are above the lumen, should behave like peripodial, default was 1.0, no change.
    		//preparing the node lists:
    		int indiceTri0Corner0 = i;
			int indiceTri0Corner1 = i+1;
			int indiceTri0Corner2 = i;
			int indiceTri1Corner0 = indiceTri0Corner0;
			int indiceTri1Corner1 = i+1;
			int indiceTri1Corner2 = indiceTri0Corner1;
			if (indiceTri0Corner1 == nCircumference){
				indiceTri0Corner1 = 0;
				indiceTri1Corner2 = 0;
				indiceTri1Corner1 = 0;
			}
			int* NodeIds = new int[6];
			//adding the first element:
			NodeIds[0] = OuterNodeArray[indiceTri0Corner0][j];
			NodeIds[1] = ColumnarBasedNodeArray[indiceTri0Corner1][j];
			NodeIds[2] = ColumnarBasedNodeArray[indiceTri0Corner2][j];
			NodeIds[3] = OuterNodeArray[indiceTri0Corner0][j+1];
			NodeIds[4] = ColumnarBasedNodeArray[indiceTri0Corner1][j+1];
			NodeIds[5] = ColumnarBasedNodeArray[indiceTri0Corner2][j+1];
			Prism* PrismPnt01;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			PrismPnt01->setGrowthWeightsViaTissuePlacement(peripodialWeight);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			if (thereIsExplicitECM){
				PrismPnt01->isECMMimimcingAtCircumference = true;
				PrismPnt01->isECMMimicing = true;
			}
			currElementId++;
			//adding the second element:
			NodeIds[0] = OuterNodeArray[indiceTri1Corner0][j];
			NodeIds[1] = OuterNodeArray[indiceTri1Corner1][j];
			NodeIds[2] = ColumnarBasedNodeArray[indiceTri1Corner2][j];
			NodeIds[3] = OuterNodeArray[indiceTri1Corner0][j+1];
			NodeIds[4] = OuterNodeArray[indiceTri1Corner1][j+1];
			NodeIds[5] = ColumnarBasedNodeArray[indiceTri1Corner2][j+1];
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			PrismPnt01->setGrowthWeightsViaTissuePlacement(peripodialWeight);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;
			if (thereIsExplicitECM){
				PrismPnt01->isECMMimimcingAtCircumference = true;
				PrismPnt01->isECMMimicing = true;
			}
    	}
    }
}

void Simulation::addNodesForPeripodialOnCap(vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &PeripodialCapNodeArray, int TissueHeightDiscretisationLayers, int LumenHeightDiscretisationLayers, int peripodialHeightDiscretisationLayers, double hPeripodial){
    int ncurrNodes = nNodes;
    int nCircumference = ColumnarBasedNodeArray.size();
    for (int i = 0; i<ncurrNodes; ++i){
		if (Nodes[i]->tissuePlacement == 1 && Nodes[i]->tissueType == 0){ //Node is apical and on the columnar layer
			if (Nodes[i]->atCircumference){
				//copy form previous list,
				int lumenCapLayer =  TissueHeightDiscretisationLayers+LumenHeightDiscretisationLayers;
				for (int j=0; j<nCircumference; ++j){
					if(Nodes[i]->Id == ColumnarBasedNodeArray[j][TissueHeightDiscretisationLayers]){
						//found the current node on the existing list
						PeripodialCapNodeArray[j].push_back(ColumnarBasedNodeArray[j][TissueHeightDiscretisationLayers]);
						int n = ColumnarBasedNodeArray[j].size();
						for (int k=lumenCapLayer; k<n; ++k ){
							//copy all the previously added nodes to this list
							//the 0th layer will be the apical surface of the columnar layer, it will be used to read element structure
							//new elements will be added via the nodes on following layers
							PeripodialCapNodeArray[j].push_back(ColumnarBasedNodeArray[j][k]);
						}
						break;
					}
				}
			}
			else{
				//push back a new vector
				int n = PeripodialCapNodeArray.size();
				PeripodialCapNodeArray.push_back(vector<int>(0));
				PeripodialCapNodeArray[n].push_back(Nodes[i]->Id);
				//now I need to add the new nodes:
				//adding nodes for the peripodial membrane:
				double* pos;
				pos = new double[3];
				pos[0] = Nodes[i]->Position[0];
				pos[1] = Nodes[i]->Position[1];
				pos[2] = Nodes[i]->Position[2]+lumenHeight;
				for (int j=0; j<peripodialHeightDiscretisationLayers+1; ++j){
					int newNodeId = Nodes.size();
					int tissuePlacement = 2; //defalut is midlayer
					if (j==0){tissuePlacement = 1;} //The first node should be apical, as it is looking into the lumen, the basal node is corrected outside loop
					Node* tmp_nd = new Node(newNodeId, 3, pos, tissuePlacement, 1); //Tissue placement is midlayer (2), tissue type is peripodial membrane (1)
					Nodes.push_back(tmp_nd);
					nNodes = Nodes.size();
					PeripodialCapNodeArray[n].push_back(newNodeId);
					pos[2] += hPeripodial;
				}
				//The last node should also be basal, at it is the top of the peripodial membrane, change its placement:
				Nodes[nNodes-1]->tissuePlacement = 0; //made tissueplacement basal
				delete[] pos;
			}
		}
	}
 }

void Simulation::constructTriangleCornerListOfApicalSurface( vector< vector<int> > &TriangleList){
	int nTri = 0;
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		vector <int> ApicalTriangles;
		(*itElement)->getApicalTriangles(ApicalTriangles);
		int nList = ApicalTriangles.size();
		for (int k=0; k<nList-2; k+=3){
			if (Nodes[ApicalTriangles[k]]->tissuePlacement == 1 &&
				Nodes[ApicalTriangles[k+1]]->tissuePlacement == 1 &&
				Nodes[ApicalTriangles[k+2]]->tissuePlacement == 1){
				TriangleList.push_back(vector<int>(3));
				TriangleList[nTri][0] = (Nodes[ApicalTriangles[k]]->Id);
				TriangleList[nTri][1] = (Nodes[ApicalTriangles[k+1]]->Id);
				TriangleList[nTri][2] = (Nodes[ApicalTriangles[k+2]]->Id);
				nTri++;
			}
		}
	}
}

void Simulation::addCapPeripodialElements( vector< vector<int> > &TriangleList, vector< vector<int> > &PeripodialCapNodeArray, int peripodialHeightDiscretisationLayers){
	int nTri = TriangleList.size();
	for (int i=0; i<nTri; ++i){
		int indiceTri0Corner0 =-1, indiceTri0Corner1=-1,indiceTri0Corner2=-1;
		int n = PeripodialCapNodeArray.size();
		//cout<<"size of the peripodial cap node array: "<<n<<endl;
		bool found0 = false, found1 = false, found2 = false;
		for (int j =0; j<n; ++j){
		   if (TriangleList[i][0] == PeripodialCapNodeArray[j][0]){
			   indiceTri0Corner0 = j;
			   found0 = true;
		   }
		   else if (TriangleList[i][1] == PeripodialCapNodeArray[j][0]){
			   indiceTri0Corner1 = j;
			   found1 = true;
		   }
		   else if (TriangleList[i][2] == PeripodialCapNodeArray[j][0]){
			   indiceTri0Corner2 = j;
			   found2 = true;
		   }
		   if (found0 && found1 && found2 ){
			   break;
		   }
		}
		if(!(found0 && found1 && found2)){
			cerr<<"could not find all the corners in addCapPeripodialElements, will give error"<<endl;
			cerr.flush();
		}
		//now I have the indices of the corners, adding the elements, as many layers as necessary:
		int* NodeIds = new int[6];
		for (int j =1; j< peripodialHeightDiscretisationLayers+1; ++j){
			NodeIds[0] = PeripodialCapNodeArray[indiceTri0Corner0][j];
			NodeIds[1] = PeripodialCapNodeArray[indiceTri0Corner1][j];
			NodeIds[2] = PeripodialCapNodeArray[indiceTri0Corner2][j];
			NodeIds[3] = PeripodialCapNodeArray[indiceTri0Corner0][j+1];
			NodeIds[4] = PeripodialCapNodeArray[indiceTri0Corner1][j+1];
			NodeIds[5] = PeripodialCapNodeArray[indiceTri0Corner2][j+1];
			Prism* PrismPnt01;
			PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
			Elements.push_back(PrismPnt01);
			nElements = Elements.size();
			currElementId++;
		}
   }
}
void Simulation::correctCircumferentialNodeAssignment(vector< vector<int> > OuterNodeArray){
	//Now I added nodes and elements to the circumference of the tissue.
	//The Nodes that were at the circumference in the columnar layer are embedded in the tissue now,
	// while the new lateral nodes of the peripodial membrane are at the circumference.
	//I will correct this in node booleans:
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		//making all circumferential node flags of the columnar layer false
		if ( (*itNode)->tissueType == 0 && (*itNode)->atCircumference){
			(*itNode)->atCircumference = false;
		}
	}
	//I have a list of the Ids fot the OuterNodes I have added, changing their circumferential node flags to true:
	int n0 = OuterNodeArray.size();
	for (int i =0 ; i<n0; ++i){
		int n1 = OuterNodeArray[i].size();
		for (int j=0; j<n1; ++j){
			int currId = OuterNodeArray[i][j];
			//find the node with the curr id: (I know Ids and index order is the same now, but staying on the safe side for potential future changes
			for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
				//making all circumferential node flags of the columnar layer false
				if ( (*itNode)->Id == currId){
					(*itNode)->atCircumference = true;
					break;
				}
			}
		}
	}
}

bool Simulation::addCurvedPeripodialMembraneToTissue(){
	bool Success = true;
	//here I am calculating the height of the tissue and the discretisation layers used for columnar layer
	double hColumnar; //The average columnar layer element height
	int LumenHeightDiscretisationLayers;  //The number of elements that are used to discretised the lumen height
	double hLumen; //The lumen element height
	double peripodialHeight; 	//Height of the peropodial membrane
	int peripodialHeightDiscretisationLayers; //The number of elements that are used to discretised the peripodial membrane height
	double hPeripodial;  //The peripodial membrane element height
	calculateDiscretisationLayers(hColumnar, LumenHeightDiscretisationLayers, hLumen, peripodialHeight, peripodialHeightDiscretisationLayers, hPeripodial);
	//Now I have all the height information and the number of layers for each subsection of the tissue
	//I can start adding the nodes.
	//Initially, I am starting with the sides.
	//First I want the list of nodes at the basal circumference of the columnar layer:
	vector <int> ColumnarCircumferencialNodeList;
	Success = generateColumnarCircumferenceNodeList(ColumnarCircumferencialNodeList);
	if (!Success){
		cerr<<"Error!! circumferential nodes not extracted properly"<<endl;
	}
	//Now I have the list, I need to sort in to rotate in one direction:
	calculateSystemCentre();
	//I will sort the list counter-clock-wise, while the (+)z is pointing at you
	sortColumnarCircumferenceNodeList(ColumnarCircumferencialNodeList);
	//Now I will create the arrays of vectors that will keep my nodes organised:
	const int nCircumference = ColumnarCircumferencialNodeList.size();
	vector< vector<int> > ColumnarBasedNodeArray( nCircumference , vector<int>(0) );
	for (int i=0;i<nCircumference; ++i){
		ColumnarBasedNodeArray[i].push_back(ColumnarCircumferencialNodeList[i]);
	}
	//Fill the array of node IDs upwards, to cover the whole columnar layer:
	fillColumnarBasedNodeList(ColumnarBasedNodeArray, ColumnarCircumferencialNodeList);
	addNodesForPeripodialOnColumnarCircumference (ColumnarBasedNodeArray, LumenHeightDiscretisationLayers, hLumen, peripodialHeightDiscretisationLayers, hPeripodial );
	//now I will add the positions for outer surface:
	//starting from the base layer of the circumferencial node list, ending at the top of the array
	int divisionSteps = 6;
	vector< vector<int> > OuterNodeArray( nCircumference , vector<int>(0) );
	vector< vector<int> > InnerNodeArray( nCircumference , vector<int>(0) );
	//calculating angle from the centre of the circular edge to each new point:
	//I am thinking of a quarter circle, divided into divisionSteps/2 triangles. The central angle will be divided into segments of:
	double tet = (M_PI /2.0) / (divisionSteps / 2.0);
	double hMidPoint = TissueHeight + 0.5*hLumen;
	for (int i=0;i<nCircumference; ++i){
		int idBase = ColumnarBasedNodeArray[i][0];
		//int idTop = ColumnarBasedNodeArray[i][ColumnarBasedNodeArray[i].size()-1];
		//double posBase[3] = { Nodes[idBase]->Position[0],Nodes[idBase]->Position[1],Nodes[idBase]->Position[2]};
		//double posTop[3] = { Nodes[idTop]->Position[0],Nodes[idTop]->Position[1],Nodes[idTop]->Position[2]};

		//first find the vector pointing out from the base node:
		int nodeId0 = ColumnarBasedNodeArray[i][0];
		int nodeId1;
		int nodeId2;
		if( i == nCircumference - 1){
			if (symmetricY){
				//the node is the end tip of a tissue with symmetric y. I should not connect it in a loop, it should add
				// the node to +x direction.
				// the node position calculation function will catch this as a flag, and will not calculate
				nodeId1 = -100;
			}
			else{
				nodeId1 = ColumnarBasedNodeArray[0][0];
			}
		}
		else{
			nodeId1 = ColumnarBasedNodeArray[i+1][0];
		}
		if ( i == 0 ){
			if (symmetricY){
				//the node is the beginning tip of a tissue with symmetric y. I should not connect it in a loop, it should add
				// the node to -x direction.
				// the node position calculation function will catch this as a flag, and will not calculate
				nodeId2 = -100;
			}
			else{
				nodeId2 = ColumnarBasedNodeArray[nCircumference-1][0];
			}
		}
		else{
			nodeId2 = ColumnarBasedNodeArray[i-1][0];
		}
		double* pos = new double[3];
		double* vec1;
		vec1 = new double[3];
		if (symmetricY && nodeId1 == -100){
			vec1[0] = -1.0;
			vec1[1] =  0.0;
			vec1[2] =  0.0;
		}
		else if (symmetricY && nodeId2 == -100){
			vec1[0] =  1.0;
			vec1[1] =  0.0;
			vec1[2] =  0.0;
		}
		else{
			double* vec2;
			vec2 = new double[3];
			for (int j=0; j<Nodes[nodeId0]->nDim; j++){
				vec1[j] = (Nodes[nodeId0]->Position[j] - Nodes[nodeId1]->Position[j]);
				vec2[j] = (Nodes[nodeId0]->Position[j] - Nodes[nodeId2]->Position[j]);
			}
			double dummy = Elements[0]->normaliseVector3D(vec1);
			dummy = Elements[0]->normaliseVector3D(vec2);
			vec1[0] += vec2[0];
			vec1[1] += vec2[1];
			vec1[2] += vec2[2];
			dummy = Elements[0]->normaliseVector3D(vec1);
			//The list is sorted counter-cock-wise, to point out, I will rotate normalised vector v0 -90 degrees on z axis:
			// (x,y,z) -> (y,-x,z);
			// then I will add this vector to the calculated mid point to gt the new node's position.
			delete[] vec2;
		}
		delete[] pos;

		//first outer point will have two radia, its movement in z direction, and on the direction of the vector pointing out
		// from the base point.
		//similarly, first inner point will have two radia (r3, r4);
		for (int j=1; j<divisionSteps; ++j ){
			double c = cos(j*tet);
			double s = sin(j*tet);
			double r1 = hMidPoint *(1 -c);
			double r2 = hMidPoint * s;
			//outer node:
			drawingPointsX.push_back(Nodes[idBase]->Position[0] + r2 * vec1[0]);
			drawingPointsY.push_back(Nodes[idBase]->Position[1] + r2 * vec1[1]);
			drawingPointsZ.push_back(Nodes[idBase]->Position[2] + r1);
			cout<<"hMidPoint: "<<hMidPoint<<endl;
		}
		delete[] vec1;
	}



	for (int i=0;i<nCircumference; ++i){
		int nColumnarBasedNodeArray = ColumnarBasedNodeArray[i].size();
		for (int j =0; j< nColumnarBasedNodeArray; ++j){
			int id = ColumnarBasedNodeArray[i][j];
			drawingPointsX.push_back( Nodes[id]->Position[0]);
			drawingPointsY.push_back( Nodes[id]->Position[1]);
			drawingPointsZ.push_back( Nodes[id]->Position[2]);
		}
	}

	return Success;
}

bool Simulation::checkIfThereIsPeripodialMembrane(){
	bool Success = true;
	vector<Node*>::iterator itNode;
	for(itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType == 1){
			thereIsPeripodialMembrane = true;
			//I have not added the peripodial membrane as yet,
			//if there is one, it came from the input mesh
			//Then I want to set up the circumferencial input properly.
			setLinkerCircumference();
			break;
		}
	}
	return Success;
}

bool Simulation::addSideECMLayer(){
    bool Success = true;
    //here I am calculating the height of the tissue and the discretisation layers used for columnar layer
    double hColumnar; //The average columnar layer element height
    hColumnar = TissueHeight/TissueHeightDiscretisationLayers;
    //Initially, I am starting with the sides.
    //First I want the list of nodes at the basal circumference of the columnar layer:
    vector <int> ColumnarCircumferencialNodeList;
    Success = generateColumnarCircumferenceNodeList(ColumnarCircumferencialNodeList);
    if (!Success){
        cerr<<"Error!! circumferential nodes not extracted properly"<<endl;
    }
    //Now I have the list, I need to sort in to rotate in one direction:
    calculateSystemCentre();
    //I will sort the list counter-clock-wise, while the (+)z is pointing at you
    sortColumnarCircumferenceNodeList(ColumnarCircumferencialNodeList);
    //Now I will create the arrays of vectors that will keep my nodes organised:
    const int nCircumference = ColumnarCircumferencialNodeList.size();
    vector< vector<int> > ColumnarBasedNodeArray( nCircumference , vector<int>(0) );
    for (int i=0;i<nCircumference; ++i){
        ColumnarBasedNodeArray[i].push_back(ColumnarCircumferencialNodeList[i]);
    }
    //Fill the array of node IDs upwards to cover the whole columnar layer:
    fillColumnarBasedNodeList(ColumnarBasedNodeArray, ColumnarCircumferencialNodeList);
    //Now add nodes for the outer layer:
    vector< vector<int> > OuterNodeArray( nCircumference , vector<int>(0) );
    addNodesForSideECMOnOuterCircumference (ColumnarBasedNodeArray, OuterNodeArray,hColumnar );
    //Now I need to add the elements:
    addSideECMElements(ColumnarBasedNodeArray, OuterNodeArray);
    //now adding the nodes of the central region:
    //Now correct position assignemnts at the circumference, and related node fixing:
    correctCircumferentialNodeAssignment(OuterNodeArray);
    return Success;
}

void Simulation::addNodesForSideECMOnOuterCircumference (vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray , double hColumnar){
    double ECMSideThickness = lateralECMThickness; //in microns
    //Now I need the average side of an element, to add new nodes accordingly:
    int nCircumference = ColumnarBasedNodeArray.size();

    for (int i=0; i<nCircumference; ++i){
        //cout<<"at node in the list: "<<i<<endl;
        //adding 3 point based node:
        int nodeId0 = ColumnarBasedNodeArray[i][0];
        int nodeId1;
        int nodeId2;
        int baseIndex0 = i;
        if( i == nCircumference - 1){
            if (symmetricY){
                //the node is the end tip of a tissue with symmetric y. I should not connect it in a loop, it should add
                // the node to +x direction.
                // the node position calculation function will catch this as a flag, and will not calculate
                nodeId1 = -100;
            }
            else{
                nodeId1 = ColumnarBasedNodeArray[0][0];
            }
        }
        else{
            nodeId1 = ColumnarBasedNodeArray[i+1][0];
        }
        if ( i == 0 ){
            if (symmetricY){
                //the node is the beginning tip of a tissue with symmetric y. I should not connect it in a loop, it should add
                // the node to -x direction.
                // the node position calculation function will catch this as a flag, and will not calculate
                nodeId2 = -100;
            }
            else{
                nodeId2 = ColumnarBasedNodeArray[nCircumference-1][0];
            }
        }
        else{
            nodeId2 = ColumnarBasedNodeArray[i-1][0];
        }
        double* pos = new double[3];
        calculateNewNodePosForPeripodialNodeAddition(nodeId0, nodeId1, nodeId2, pos, ECMSideThickness);

        //cout<<" calculated pos : "<<pos[0] <<" "<<pos[1]<<" "<<pos[2]<<endl;
        //Adding the array of new nodes:

        //adding the base:
		int newNodeId = Nodes.size();
		Node* tmp_nd = new Node(newNodeId, 3, pos, 0, 0); //Tissue placement basal (0), tissue type is columnar
		Nodes.push_back(tmp_nd);
		nNodes = Nodes.size();
		OuterNodeArray[i].push_back(newNodeId);

        //adding the nodes for the columnar layer:
        for (int j=1; j<TissueHeightDiscretisationLayers+1; ++j){
            //pos[2] += hColumnar;
            pos[2] =  Nodes[ColumnarBasedNodeArray[baseIndex0][j]]->Position[2];
            //cout<<" pos for columnar aligned new node: "<<pos[0] <<" "<<pos[1]<<" "<<pos[2]<<" hColumnar: "<<hColumnar<<endl;

            int newNodeId = Nodes.size();
            int tissuePlacement = 2; //midline
            if (j == TissueHeightDiscretisationLayers){
                tissuePlacement=1;//apical
            }

            Node* tmp_nd = new Node(newNodeId, 3, pos, tissuePlacement, 0); //Tissue placement is midlayer (2), tissue type is columnar (0)
            Nodes.push_back(tmp_nd);
            nNodes = Nodes.size();
            OuterNodeArray[i].push_back(newNodeId);
        }
        delete[] pos;
    }
}

void Simulation::addSideECMElements(vector< vector<int> > &ColumnarBasedNodeArray, vector< vector<int> > &OuterNodeArray){
    // I need to add the elements:
    // Two elements are added for each element on the side:
    // First triangle base will be New node 0, oldNode1, oldnode0, top of the element will be read from the node id stacks
    // Second triangle base will be: New node 0, new node 1, old node 1
    int totalLayers = TissueHeightDiscretisationLayers;
    int nCircumference = ColumnarBasedNodeArray.size();
    int nLoop = nCircumference;
    if (symmetricY){
        //in symmetric setup, the last node is not connected to the first node, we simply ignore that step
        nLoop--;
    }

    for (int i=0;i<nLoop; ++i){
        for (int j=0;j<totalLayers; ++j){
            //these are columnar elements, peripodialness is always zero.
            double peripodialWeight = 0.0; //
            //preparing the node lists:
            int indiceTri0Corner0 = i;
            int indiceTri0Corner1 = i+1;
            int indiceTri0Corner2 = i;
            int indiceTri1Corner0 = indiceTri0Corner0;
            int indiceTri1Corner1 = i+1;
            int indiceTri1Corner2 = indiceTri0Corner1;
            if (indiceTri0Corner1 == nCircumference){
                indiceTri0Corner1 = 0;
                indiceTri1Corner2 = 0;
                indiceTri1Corner1 = 0;
            }
            int* NodeIds = new int[6];
            //adding the first element:
            NodeIds[0] = OuterNodeArray[indiceTri0Corner0][j];
            NodeIds[1] = ColumnarBasedNodeArray[indiceTri0Corner1][j];
            NodeIds[2] = ColumnarBasedNodeArray[indiceTri0Corner2][j];
            NodeIds[3] = OuterNodeArray[indiceTri0Corner0][j+1];
            NodeIds[4] = ColumnarBasedNodeArray[indiceTri0Corner1][j+1];
            NodeIds[5] = ColumnarBasedNodeArray[indiceTri0Corner2][j+1];
            Prism* PrismPnt01;
            PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
            PrismPnt01->setGrowthWeightsViaTissuePlacement(peripodialWeight);
            PrismPnt01->setECMMimicing(true);
            PrismPnt01->isECMMimimcingAtCircumference = true;
            Elements.push_back(PrismPnt01);
            nElements = Elements.size();
            currElementId++;
            //adding the second element:
            NodeIds[0] = OuterNodeArray[indiceTri1Corner0][j];
            NodeIds[1] = OuterNodeArray[indiceTri1Corner1][j];
            NodeIds[2] = ColumnarBasedNodeArray[indiceTri1Corner2][j];
            NodeIds[3] = OuterNodeArray[indiceTri1Corner0][j+1];
            NodeIds[4] = OuterNodeArray[indiceTri1Corner1][j+1];
            NodeIds[5] = ColumnarBasedNodeArray[indiceTri1Corner2][j+1];
            PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
            PrismPnt01->setGrowthWeightsViaTissuePlacement(peripodialWeight);
            PrismPnt01->setECMMimicing(true);
            PrismPnt01->isECMMimimcingAtCircumference = true;
            Elements.push_back(PrismPnt01);
            nElements = Elements.size();
            currElementId++;
        }
    }
}

bool Simulation::addStraightPeripodialMembraneToTissue(){
    //cout<<"adding peripodial membrane from scratch" <<endl;
	bool Success = true;
    //here I am calculating the height of the tissue and the discretisation layers used for columnar layer
    double hColumnar; //The average columnar layer element height
    int LumenHeightDiscretisationLayers;  //The number of elements that are used to discretised the lumen height
    double hLumen; //The lumen element height
    double peripodialHeight; 	//Height of the peropodial membrane
    int peripodialHeightDiscretisationLayers; //The number of elements that are used to discretised the peripodial membrane height
    double hPeripodial;  //The peripodial membrane element height
    calculateDiscretisationLayers(hColumnar, LumenHeightDiscretisationLayers, hLumen, peripodialHeight, peripodialHeightDiscretisationLayers, hPeripodial);
    //Now I have all the height information and the number of layers for each subsection of the tissue
    //I can start adding the nodes.
    //Initially, I am starting with the sides.
    //First I want the list of nodes at the basal circumference of the columnar layer:
    vector <int> ColumnarCircumferencialNodeList;
    Success = generateColumnarCircumferenceNodeList(ColumnarCircumferencialNodeList);
    if (!Success){
        cerr<<"Error!! circumferential nodes not extracted properly"<<endl;
    }
    //Now I have the list, I need to sort in to rotate in one direction:
    calculateSystemCentre();
    //I will sort the list counter-clock-wise, while the (+)z is pointing at you
    sortColumnarCircumferenceNodeList(ColumnarCircumferencialNodeList);
    //Now I will create the arrays of vectors that will keep my nodes organised:
    const int nCircumference = ColumnarCircumferencialNodeList.size();
    vector< vector<int> > ColumnarBasedNodeArray( nCircumference , vector<int>(0) );
    for (int i=0;i<nCircumference; ++i){
       	ColumnarBasedNodeArray[i].push_back(ColumnarCircumferencialNodeList[i]);
    }
    //Fill the array of node IDs upwards, to cover the whole columnar layer:
    fillColumnarBasedNodeList(ColumnarBasedNodeArray, ColumnarCircumferencialNodeList);
    addNodesForPeripodialOnColumnarCircumference (ColumnarBasedNodeArray, LumenHeightDiscretisationLayers, hLumen, peripodialHeightDiscretisationLayers, hPeripodial );
    //Now add nodes for the outer layer:
    vector< vector<int> > OuterNodeArray( nCircumference , vector<int>(0) );
    addNodesForPeripodialOnOuterCircumference (ColumnarBasedNodeArray, OuterNodeArray, hColumnar, LumenHeightDiscretisationLayers, hLumen, peripodialHeightDiscretisationLayers, hPeripodial );
    //Now I need to add the elements:
    addLateralPeripodialElements(LumenHeightDiscretisationLayers,peripodialHeightDiscretisationLayers, ColumnarBasedNodeArray, OuterNodeArray);
    //now adding the nodes of the central region:
    // I am initialising the vector as the size of circumference node list, for any node at the circumference, I will copy from previuos lists,
    // for any node that is not at the circumference, I will add a new vector of integers to this list
    vector< vector<int> > PeripodialCapNodeArray( nCircumference , vector<int>(0) );
    addNodesForPeripodialOnCap(ColumnarBasedNodeArray, PeripodialCapNodeArray, TissueHeightDiscretisationLayers, LumenHeightDiscretisationLayers, peripodialHeightDiscretisationLayers, hPeripodial);
    //Now I will generate the triangle list from the apical surface of the columnar layer:
    vector< vector<int> > TriangleList( 0 , vector<int>(0) ); //The IDs of 3 nodes that are supposed to form a triangle
    constructTriangleCornerListOfApicalSurface(TriangleList);
    //Now I have the list of nodes forming triangles
    //For each triangle, find the map of indices in the 0th layer of PeripodialCapNodeArray, and add elements according to the node id stacks:
    addCapPeripodialElements( TriangleList, PeripodialCapNodeArray, peripodialHeightDiscretisationLayers);
    //Peripodial membrane is added, now correct position assignemnts at the circumference, and related node fixing:
    correctCircumferentialNodeAssignment(OuterNodeArray);
    return Success;
}

void Simulation::checkForExperimentalSetupsBeforeIteration(){
	if (stretcherAttached){
		moveClampedNodesForStretcher();
	}
	//moveAFMBead();
}

void Simulation::checkForExperimentalSetupsWithinIteration(){
}

void Simulation::checkForExperimentalSetupsAfterIteration(){
	if (stretcherAttached){
		recordForcesOnClampBorders();
	}
}

void Simulation::calculateStiffnessChangeRatesForActin(int idOfCurrentStiffnessPerturbation){
    startedStiffnessPerturbation[idOfCurrentStiffnessPerturbation] = true;
    const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		bool applyToThisElement = (*itElement)->isActinStiffnessChangeAppliedToElement(ThereIsWholeTissueStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsApicalStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsBasalStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation[idOfCurrentStiffnessPerturbation],ThereIsBasolateralStiffnessPerturbation[idOfCurrentStiffnessPerturbation], stiffnessPerturbationEllipseBandIds[idOfCurrentStiffnessPerturbation], numberOfStiffnessPerturbationAppliesEllipseBands[idOfCurrentStiffnessPerturbation]);
		if (applyToThisElement){
            (*itElement)->calculateStiffnessPerturbationRate(ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation[idOfCurrentStiffnessPerturbation], stiffnessPerturbationBeginTimeInSec[idOfCurrentStiffnessPerturbation],stiffnessPerturbationEndTimeInSec[idOfCurrentStiffnessPerturbation], stiffnessChangedToFractionOfOriginal[idOfCurrentStiffnessPerturbation]);
		}
	}
}

void Simulation::updateStiffnessChangeForActin(int idOfCurrentStiffnessPerturbation){
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		bool applyToThisElement = (*itElement)->isActinStiffnessChangeAppliedToElement(ThereIsWholeTissueStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsApicalStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsBasalStiffnessPerturbation[idOfCurrentStiffnessPerturbation], ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation[idOfCurrentStiffnessPerturbation],  ThereIsBasolateralStiffnessPerturbation[idOfCurrentStiffnessPerturbation], stiffnessPerturbationEllipseBandIds[idOfCurrentStiffnessPerturbation], numberOfStiffnessPerturbationAppliesEllipseBands[idOfCurrentStiffnessPerturbation]);
		if (applyToThisElement){
            (*itElement)->updateStiffnessMultiplier(dt);
            (*itElement)->updateElasticProperties();
            //cout<<" element: "<<(*itElement)->Id<<" stiffness multiplier: "<<(*itElement)->getStiffnessMultiplier()<<endl;
		}
	}
}

void Simulation::checkStiffnessPerturbation(){
	int n = stiffnessPerturbationBeginTimeInSec.size();
	for (int idOfCurrentStiffnessPerturbation=0; idOfCurrentStiffnessPerturbation<n; ++idOfCurrentStiffnessPerturbation){
		if (currSimTimeSec >=stiffnessPerturbationBeginTimeInSec[idOfCurrentStiffnessPerturbation] && currSimTimeSec <stiffnessPerturbationEndTimeInSec[idOfCurrentStiffnessPerturbation]){
			if (startedStiffnessPerturbation[idOfCurrentStiffnessPerturbation] == false){
				calculateStiffnessChangeRatesForActin(idOfCurrentStiffnessPerturbation);
			}
			updateStiffnessChangeForActin(idOfCurrentStiffnessPerturbation);
		}
	}
}


/*
void Simulation::checkStiffnessPerturbation(){
    if (currSimTimeSec >=stiffnessPerturbationBeginTimeInSec && currSimTimeSec <stiffnessPerturbationEndTimeInSec){    
        const int maxThreads = omp_get_max_threads();
        omp_set_num_threads(maxThreads);
        #pragma omp parallel for
        for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
            if((*itElement)->insideEllipseBand){ //this covers the tissue type as well as the position, ellipses are applied to columnar and linker zones only.
                if( ThereIsWholeTissueStiffnessPerturbation  //the whole columnar tissue is perturbed
                   || ((*itElement)->tissuePlacement == 0 && ThereIsBasalStiffnessPerturbation  ) //the basal surface is perturbed and element is basal
                   || ((*itElement)->tissuePlacement == 1 && ThereIsApicalStiffnessPerturbation ) // the apical surface is perturbed and element is apical
                  ){
                    for (int stiffnessPerturbationRangeCounter =0; stiffnessPerturbationRangeCounter<numberOfStiffnessPerturbationAppliesEllipseBands; ++stiffnessPerturbationRangeCounter){
                        if ((*itElement)->coveringEllipseBandId == stiffnessPerturbationEllipseBandIds[stiffnessPerturbationRangeCounter]){
                            if (startedStiffnessPerturbation == false){
                                //this is the first time step I am applying stiffenning.
                                //I need to calculate rates per element first:ß
                                (*itElement)->calculateStiffnessPerturbationRate(stiffnessPerturbationBeginTimeInSec,stiffnessPerturbationEndTimeInSec, stiffnessChangedToFractionOfOriginal);
                            }                        
                            (*itElement)->updateStiffnessMultiplier(dt);
                            (*itElement)->updateElasticProperties();
                        }
                    }
                }
            }
        }
        //I have entered the time loop once, therefore, the stiffness change rates have been calculated.
        startedStiffnessPerturbation = true;
    }
}

*/

void Simulation::updateChangeForExplicitECM(int idOfCurrentECMPerturbation){
    const int maxThreads = omp_get_max_threads();
	#pragma omp parallel for
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		bool applyToThisElement = (*itElement)->isECMChangeAppliedToElement(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
		if (applyToThisElement){
			 (*itElement)->updateStiffnessMultiplier(dt);
			 (*itElement)->updateElasticProperties();
		}
	}
}

void Simulation::updateChangeForViscosityBasedECMDefinition(int idOfCurrentECMPerturbation){
	cout<<"in viscosity update"<<endl;
    //const int maxThreads = omp_get_max_threads();
    #pragma omp parallel for
	for (vector<Node*>::iterator itNode = Nodes.begin(); itNode<Nodes.end(); ++itNode){
		//if (!(*itNode)->allOwnersECMMimicing){
			if(((*itNode)->tissuePlacement == 0 && changeBasalECM[idOfCurrentECMPerturbation] ) || ((*itNode)->tissuePlacement == 1 && changeApicalECM[idOfCurrentECMPerturbation] )){
				if((*itNode)->insideEllipseBand){
					for (int ECMReductionRangeCounter =0; ECMReductionRangeCounter<numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]; ++ECMReductionRangeCounter){
						if ((*itNode)->coveringEllipseBandId == ECMChangeEllipseBandIds[idOfCurrentECMPerturbation][ECMReductionRangeCounter]){
							for (int i =0; i<(*itNode)->nDim; ++i){
								double viscosityChange = (*itNode)->ECMViscosityChangePerHour[i]/3600*dt;
								double newViscosity = (*itNode)->externalViscosity[i] - viscosityChange;
								//avoiding setting negative viscosity!
								if (newViscosity>(*itNode)->maximumExternalViscosity[i]){
									(*itNode)->externalViscosity[i] =(*itNode)->maximumExternalViscosity[i];
								}
								else if(newViscosity<(*itNode)->minimumExternalViscosity[i]){
									(*itNode)->externalViscosity[i] =(*itNode)->minimumExternalViscosity[i];
								}
								else{
									(*itNode)->externalViscosity[i] = newViscosity;
								}
							}
							//cout<<" updated external viscosity: "<<(*itNode)->externalViscosity[0]<<" change per hr: "<<(*itNode)->ECMViscosityChangePerHour[0]<<" min: "<< (*itNode)->minimumExternalViscosity[0]<<" max: "<<(*itNode)->maximumExternalViscosity[0] <<endl;
						}
					}
				}
			}
		//}
	}
	//cout<<"finished viscosity update"<<endl;
}

void Simulation::calculateChangeRatesForECM(int idOfCurrentECMPerturbation){
	changedECM[idOfCurrentECMPerturbation] = true; //this will not be used for emergent ecm perturbations
	//this is the first time step I am changing the ECM stiffness.
	//I need to calculate rates first.
	//If there is explicit ECM, I will calculate the young modulus change via elements.
	//If there is no explicit ECM, the viscosity reflects the ECM stiffness, and I will change the viscosity on a nodal basis.
	if( thereIsExplicitECM){
		cout<<" there is explicit ECM for rate calculation"<<endl;
		#pragma omp parallel for
		for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			bool applyToThisElement = (*itElement)->isECMChangeAppliedToElement(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
			if (applyToThisElement){
				//the first input is used for checking basolateral stiffenning combined with apical relaxation
				//the ECM does not have such options. Will give the boolean as false and continue.
				(*itElement)->calculateStiffnessPerturbationRate(false, ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation],ECMChangeEndTimeInSec[idOfCurrentECMPerturbation], ECMStiffnessChangeFraction[idOfCurrentECMPerturbation]);
			}
		}
		//now I need to check for the viscosity based calculation:
		#pragma omp parallel for
		for (vector<Node*>::iterator itNode = Nodes.begin(); itNode<Nodes.end(); ++itNode){
			bool applyToThisNode = (*itNode)->isECMChangeAppliedToNode(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
			if(applyToThisNode){
				double timeDifferenceInHours = (ECMChangeEndTimeInSec[idOfCurrentECMPerturbation] - ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation])/3600;
				for (int i=0; i<3; ++i){
					(*itNode)->ECMViscosityChangePerHour[i] = (*itNode)->initialExternalViscosity[i]*(1-ECMViscosityChangeFraction[idOfCurrentECMPerturbation])/timeDifferenceInHours;
					if ((*itNode)->ECMViscosityChangePerHour[i]<0){
						(*itNode)->maximumExternalViscosity[i] = (*itNode)->initialExternalViscosity[i]*ECMViscosityChangeFraction[idOfCurrentECMPerturbation];
					}
					else{
						(*itNode)->minimumExternalViscosity[i] = (*itNode)->initialExternalViscosity[i]*ECMViscosityChangeFraction[idOfCurrentECMPerturbation];
					}
				}
			}
		}
	}
	else{
		double timeDifferenceInHours = (ECMChangeEndTimeInSec[idOfCurrentECMPerturbation] - ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation])/3600;
		#pragma omp parallel for
		for (vector<Node*>::iterator itNode = Nodes.begin(); itNode<Nodes.end(); ++itNode){
			bool applyToThisNode = (*itNode)->isECMChangeAppliedToNode(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
			if(applyToThisNode){
				for (int i=0; i<3; ++i){
					(*itNode)->ECMViscosityChangePerHour[i] = (*itNode)->initialExternalViscosity[i]*(1-ECMViscosityChangeFraction[idOfCurrentECMPerturbation])/timeDifferenceInHours;
					if ((*itNode)->ECMViscosityChangePerHour[i]<0){
						(*itNode)->maximumExternalViscosity[i] = (*itNode)->initialExternalViscosity[i]*ECMViscosityChangeFraction[idOfCurrentECMPerturbation];
					}
					else{
						(*itNode)->minimumExternalViscosity[i] = (*itNode)->initialExternalViscosity[i]*ECMViscosityChangeFraction[idOfCurrentECMPerturbation];
					}
				}
			}
		}
	}
}

void Simulation::updateECMRenewalHalflifeMultiplier(int idOfCurrentECMPerturbation){
	if(thereIsExplicitECM){
		if (ECMChangeTypeIsEmergent[idOfCurrentECMPerturbation]){
			double totalTimeChange = ECMChangeEndTimeInSec[idOfCurrentECMPerturbation] - ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation];
			double incrementPerSec = (ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation] - 1.0) / totalTimeChange;
			#pragma omp parallel for
			for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
				if ((*itElement)->isECMMimicing){
					bool applyToThisElement = (*itElement)->isECMChangeAppliedToElement(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
					if (applyToThisElement){
						(*itElement)->plasticDeformationHalfLifeMultiplier +=  incrementPerSec*dt;
						if (incrementPerSec<0 && (*itElement)->plasticDeformationHalfLifeMultiplier < ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation]){
							(*itElement)->plasticDeformationHalfLifeMultiplier = ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation];
						}
						if (incrementPerSec>0 && (*itElement)->plasticDeformationHalfLifeMultiplier > ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation]){
							(*itElement)->plasticDeformationHalfLifeMultiplier = ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation];
						}
					}
				}
			}
		}
		else{
			if (currSimTimeSec>ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation]){
				//perturbation on ECM renewal half life started
				double currECMRenewaHalfLifeMultiplier = 1.0;
				if (currSimTimeSec>=ECMChangeEndTimeInSec[idOfCurrentECMPerturbation]){
					currECMRenewaHalfLifeMultiplier = ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation];
				}
				else{
					double totalTimeChange = ECMChangeEndTimeInSec[idOfCurrentECMPerturbation] - ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation];
					double currTimeChange = currSimTimeSec - ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation];
					currECMRenewaHalfLifeMultiplier = 1 + (ECMRenewalHalfLifeTargetFraction[idOfCurrentECMPerturbation] - 1.0) * currTimeChange / totalTimeChange;
				}
				#pragma omp parallel for
				for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
					if ((*itElement)->isECMMimicing){
						bool applyToThisElement = (*itElement)->isECMChangeAppliedToElement(changeApicalECM[idOfCurrentECMPerturbation], changeBasalECM[idOfCurrentECMPerturbation], ECMChangeEllipseBandIds[idOfCurrentECMPerturbation], numberOfECMChangeEllipseBands[idOfCurrentECMPerturbation]);
						if (applyToThisElement){
							(*itElement)->plasticDeformationHalfLifeMultiplier =  currECMRenewaHalfLifeMultiplier;
						}
					}
				}
			}
		}
	}
}

void Simulation::checkECMChange(){
	//First chak via compartments:
	if( thereIsExplicitECM){
		#pragma omp parallel for
		for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			bool updateStiffness = false;
			if( (*itElement)->isECMMimicing && (*itElement)->tissuePlacement == 0 && (*itElement)->tissueType ==0 ){//columnar basal ecmmimicking element
				//check notum:
				if (notumECMChangeFraction != 1.0 && (*itElement)->compartmentType == 2){//notum:
					if (currSimTimeSec>= notumECMChangeInitTime && currSimTimeSec< notumECMChangeEndTime){
						double currentFraction = 1.0 + (notumECMChangeFraction-1)*(*itElement)->compartmentIdentityFraction;
						(*itElement)->calculateStiffnessPerturbationRate(false, notumECMChangeInitTime,notumECMChangeEndTime, currentFraction);
						updateStiffness=true;
					}
				}
				//check hinge:
				if (hingeECMChangeFraction != 1.0 && (*itElement)->compartmentType == 1){//hinge:
					if (currSimTimeSec>= hingeECMChangeInitTime && currSimTimeSec< hingeECMChangeEndTime){
						double currentFraction = 1.0 + (hingeECMChangeFraction-1)*(*itElement)->compartmentIdentityFraction;
						(*itElement)->calculateStiffnessPerturbationRate(false, hingeECMChangeInitTime,hingeECMChangeEndTime, currentFraction);
						updateStiffness=true;
					}
				}
				//check pouch:
				if (pouchECMChangeFraction != 1.0 && (*itElement)->compartmentType == 0){//pouch:
					if (currSimTimeSec>= pouchECMChangeInitTime && currSimTimeSec< pouchECMChangeEndTime){
						double currentFraction = 1.0 + (pouchECMChangeFraction-1)*(*itElement)->compartmentIdentityFraction;
						(*itElement)->calculateStiffnessPerturbationRate(false, pouchECMChangeInitTime,pouchECMChangeEndTime, currentFraction);
						updateStiffness=true;
					}
				}
			}
			if (updateStiffness){
				(*itElement)->updateStiffnessMultiplier(dt);
				(*itElement)->updateElasticProperties();
			}
		}
	}
	//Now go through ellipses, ellipses overwrite the compartment based definitions
	int n = ECMChangeBeginTimeInSec.size();
	for (int idOfCurrentECMPerturbation=0; idOfCurrentECMPerturbation<n; ++idOfCurrentECMPerturbation){
		if (ECMChangeTypeIsEmergent[idOfCurrentECMPerturbation]){
			if (currSimTimeSec >=ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation]){
				calculateChangeRatesForECM(idOfCurrentECMPerturbation);
				if( thereIsExplicitECM){
					updateECMRenewalHalflifeMultiplier(idOfCurrentECMPerturbation);
					updateChangeForExplicitECM(idOfCurrentECMPerturbation);
				}
				updateChangeForViscosityBasedECMDefinition(idOfCurrentECMPerturbation);
			}
		}
		else{
			if (currSimTimeSec >=ECMChangeBeginTimeInSec[idOfCurrentECMPerturbation] && currSimTimeSec <ECMChangeEndTimeInSec[idOfCurrentECMPerturbation]){
				if (changedECM[idOfCurrentECMPerturbation] == false){
					calculateChangeRatesForECM(idOfCurrentECMPerturbation);
				}
				if( thereIsExplicitECM){
					updateECMRenewalHalflifeMultiplier(idOfCurrentECMPerturbation);
					updateChangeForExplicitECM(idOfCurrentECMPerturbation);
				}
				updateChangeForViscosityBasedECMDefinition(idOfCurrentECMPerturbation);
			}
		}
	}
	cout<<"finished ECM change"<<endl;
}

void Simulation::checkEllipseAllocationWithCurvingNodes(){
	cout<<"checking ellipse allocation with curving nodes, thereIsEmergentEllipseMarking? "<<thereIsEmergentEllipseMarking<<endl;
	int selectedEllipseBandId =100;
	for (vector<Node*>::iterator itNode= Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->onFoldInitiation){
			if((*itNode)->coveringEllipseBandId != 100 && (*itNode)->coveringEllipseBandId != 101){
				if ((*itNode)->tissuePlacement == 1){ //apical collapse: ECM relaxation, cell shortening, volume redistribution to shrink top
					selectedEllipseBandId = 100;
				}
				else{ //basal collapse, volume redistribution to shrink bottom
					selectedEllipseBandId = 101;
				}
				int n = (*itNode)->connectedElementIds.size();
				for (int i=0; i<n; ++i){
					int currElementId = (*itNode)->connectedElementIds[i];
					//cout<<"assigning element "<<currElementId<<" vie owner node "<< (*itNode)->Id<<endl;
					if (Elements[currElementId]->coveringEllipseBandId == 100 || Elements[currElementId]->coveringEllipseBandId == 101 || Elements[currElementId]->isECMMimimcingAtCircumference){
						continue;
					}
					//check if at least two nodes of the element are at curves:
					bool changeEllipseId = Elements[currElementId]->hasEnoughNodesOnCurve(Nodes);
					if (changeEllipseId){
						Elements[currElementId]->coveringEllipseBandId = selectedEllipseBandId;
						cout<<" Id : "<<currElementId<<" assigned "<<selectedEllipseBandId<<" in function checkEllipseAllocationWithCurvingNodes"<<endl;
						Elements[currElementId]->assignEllipseBandIdToWholeTissueColumn(TissueHeightDiscretisationLayers,Nodes,Elements);
					}
				}
			}
		}
	}
}

void Simulation::checkForLeftOutElementsInEllipseAssignment(){
	for (vector<Node*>::iterator itNode = Nodes.begin(); itNode < Nodes.end(); ++itNode){
		if ((*itNode)->checkOwnersforEllipseAsignment){
			int currentEllipseId = (*itNode)->coveringEllipseBandId;
			int n = (*itNode)->connectedElementIds.size();
			for (int i=0; i<n ; ++i){
				int currElementId = (*itNode)->connectedElementIds[i];
				if (Elements[currElementId]->coveringEllipseBandId ==currentEllipseId ){
					continue;
				}
				int nNodes = Elements[currElementId]->getNodeNumber();
				int* nodeIds = Elements[currElementId]->getNodeIds();
				bool hasNodeOutsideEllipse = false;
				for (int j=0; j<nNodes;++j){
					int currNodeId = nodeIds[j];
					if (Nodes[currNodeId]->coveringEllipseBandId != currentEllipseId){
						hasNodeOutsideEllipse = true;
						break;
					}
				}
				if (!hasNodeOutsideEllipse){
					//the element does not have any nodes outside the ellipse
					//should assign to ellipse
					Elements[currElementId]->coveringEllipseBandId = currentEllipseId;
					cout<<" Id : "<<currElementId<<" assigned "<<currentEllipseId<<" in function checkForLeftOutElementsInEllipseAssignment"<<endl;
					Elements[currElementId]->assignEllipseBandIdToWholeTissueColumn(TissueHeightDiscretisationLayers,Nodes,Elements);
				}
			}
		}
	}
}


void Simulation::updateEllipseWithCollapse(){
	for (vector<ShapeBase*>::iterator itEle=Elements.begin(); itEle<Elements.end(); itEle++){
		if ((*itEle)->coveringEllipseBandId == 100 || (*itEle)->coveringEllipseBandId == 101 || (*itEle)->isECMMimimcingAtCircumference){
			continue;
		}
		if((*itEle)->tissuePlacement == 0 || (*itEle)->tissuePlacement == 1 || ((*itEle)->tissuePlacement == 2 && (*itEle)->spansWholeTissue ) ){
			(*itEle)->checkForCollapsedNodes(TissueHeightDiscretisationLayers, Nodes, Elements);
		}
	}
}

void Simulation::checkForEllipseIdUpdateWithECMDegradation(){
	for (vector<ShapeBase*>::iterator itEle=Elements.begin(); itEle<Elements.end(); itEle++){
		if ((*itEle)->isECMMimicing && (*itEle)->coveringEllipseBandId == 100){
			if ((*itEle)->stiffnessMultiplier<shapeChangeECMLimit){
				(*itEle)->coveringEllipseBandId = 102;
				(*itEle)->assignEllipseBandIdToWholeTissueColumn(TissueHeightDiscretisationLayers,Nodes,Elements);
			}
		}
	}
}

void Simulation::updateOnFoldNodesFromCollapse(){
	cout<<"calling updateOnFoldNodesFromCollapse"<<endl;
	if (checkedForCollapsedNodesOnFoldingOnce){
		return;
	}
	//update detection threshold  grid:
	double periAverageSideLength = 0,colAverageSideLength = 0;
	getAverageSideLength(periAverageSideLength, colAverageSideLength);
	if (thereIsPeripodialMembrane){
		colAverageSideLength = (periAverageSideLength+colAverageSideLength)/2.0;
	}
	for (int i=0;i<10;++i){
		for(int j=0;j<5;++j){
			// 0.4*1.2 normal packing detection value
			// 0.25 is stable with clashes
			// 0.333
			packingDetectionThresholdGrid[i][j] = 0.4 * packingDetectionThresholdGrid[i][j];
		}
	}
	//check collapsed nodes, run this code only on initiation!
	for (vector<Node*>::iterator itNode = Nodes.begin(); itNode < Nodes.end(); ++itNode){
		if ( (*itNode)->collapsedWith.size()>0) {
			bool checkForFold = true;
			int collapsedId = (*itNode)->collapsedWith[0];
			if (collapsedId == (*itNode)->Id){
				if ((*itNode)->collapsedWith.size()>1){
					collapsedId = (*itNode)->collapsedWith[1];
				}
				else{
					checkForFold=false;
				}
			}
			if (checkForFold){
				assignFoldRegionAndReleasePeripodial((*itNode), Nodes[collapsedId]);
			}
		}
	}
	checkedForCollapsedNodesOnFoldingOnce = true;
}

void Simulation::checkForEmergentEllipseFormation(){
	updateEllipseWithCollapse();
	updateOnFoldNodesFromCollapse();

	checkEllipseAllocationWithCurvingNodes();
	checkForLeftOutElementsInEllipseAssignment();
	checkForEllipseIdUpdateWithECMDegradation();
}

void Simulation::artificialRelax(){
	for (vector<ShapeBase*>::iterator itEle=Elements.begin(); itEle<Elements.end(); itEle++){
		bool relaxElement = false;
		if ((*itEle)->isECMMimicing){
			if( relaxECMInArtificialRelaxation ){
				relaxElement = true;
			}
		}
		else if((*itEle)->tissueType == 0 ){
			relaxElement = true;
		}
		if (relaxElement){
			(*itEle)->relaxElasticForces();
		}
	}
}


bool Simulation::runOneStep(){

    bool Success = true;
	cout<<"entered run one step, time "<<currSimTimeSec<<endl;
	//ablateSpcific();
    if (currSimTimeSec == -16*3600) {
    	pokeElement(224,0,0,-0.02);
    }
    //cout<<"Position Of Node 187"<<Nodes[187]->Position[2]<<" node 246: "<<Nodes[246]->Position[2]<<endl;
    manualPerturbationToInitialSetup(false,false); //bool deform, bool rotate
    resetForces(true); // reset the packing forces together with all the rest of the forces here
    checkForEmergentEllipseFormation();
    int freq = 10.0/dt ;
    if (freq <1 ){freq =1;}
    if ((timestep - 1)% freq  == 0){
        cout<<"At time -- "<<currSimTimeSec<<" sec ("<<currSimTimeSec/3600<<" hours - "<<timestep<<" timesteps)"<<endl;
        alignTissueDVToXPositive();
        //alignTissueAPToXYPlane();
        calculateBoundingBox();
        calculateDVDistance();
        vector<ShapeBase*>::iterator itElement;
        for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        	(*itElement)->calculateRelativePosInBoundingBox(boundingBox[0][0],boundingBox[0][1],boundingBoxSize[0],boundingBoxSize[1]);
        }
        updateRelativePositionsToApicalPositioning();
    }
    //cout<<"after bounding box"<<endl;
    checkForExperimentalSetupsBeforeIteration();
    if (ThereIsStiffnessPerturbation) {
    	checkStiffnessPerturbation();
    }
    if (thereIsECMChange) {
    	checkECMChange();
    }
    if(nMyosinFunctions > 0){
    	checkForMyosinUpdates();
    }

    //cout<<"after checking for myosin"<<endl;
    bool thereIsRefinement = false;
    if (thereIsRefinement){
    	int oldNodeSize = nNodes;
    	flagElementsThatNeedRefinement(); //this is done in parallel, the actual refinement cannot be parallel
    	refineElements();
    	cout<<"outside refine elements"<<endl;
    	if (nNodes != oldNodeSize){
    		cout<<"re-initiating system forces"<<endl;
    		reInitiateSystemForces(oldNodeSize);
    		NRSolver->reInitiateMatricesAfterRefinement(nNodes);
    	}
    }
    //cout<<"after refinement "<<endl;
    if(nGrowthFunctions>0 || nShapeChangeFunctions >0 || numberOfClones >0 ){
    	checkForPinningPositionsUpdate();
        //outputFile<<"calculating growth"<<endl;
		//if ((timestep - 1)% growthRotationUpdateFrequency  == 0){
			updateGrowthRotationMatrices();
		//}
        if(nGrowthFunctions>0 || numberOfClones > 0){
        	calculateGrowth();
        }
        if(nShapeChangeFunctions>0){
        	calculateShapeChange();
        }
    }
    //cout<<"after growth and shape change"<<endl;
    if (conservingColumnVolumes){
    	conserveColumnVolume();
    }
    //cout<<"after volume conservation"<<endl;
    if(thereIsPlasticDeformation || thereIsExplicitECM){
    	updatePlasticDeformation();
    }
    //cout<<"after plastic deformation"<<endl;
    if(thereIsCellMigration){
    	cellMigrationTool->updateMigratingElements(Elements);
    	cellMigrationTool->updateMigrationLists(Elements, dt);
    	cellMigrationTool->updateVolumesWithMigration(Elements);
    }
    //cout<<"after migration"<<endl;
    checkForVolumeRedistributionInTissue();
   // cout<<" after volume redistribution"<<endl;
    for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        if (!(*itElement)->IsAblated){
        	//(*itElement)->updateInternalViscosityTest();
        	(*itElement)->growShapeByFg();
        	//(*itElement)->defineFgByGrowthTemplate();
        	//(*itElement)->changeShapeByFsc(dt);
        }
    }
    //cout<<"after grow by Fg"<<endl;
    updateNodeMasses();
    updateNodeViscositySurfaces();
    //cout<<"updating connected node lists"<<endl;
    updateElementToConnectedNodes(Nodes);
    //cout<<"calculate Myosin Forces"<<endl;
    calculateMyosinForces();
    //cout<<"before update position of collapsing nodes"<<endl;
    updatePositionsOfNodesCollapsingInStages();
    //cout<<"after update position of collapsing nodes"<<endl;
    detectPacingNodes();
	//detectPackingToAFMBead();
    //cout<<"after detect packing"<<endl;
    if (encloseTissueBetweenSurfaces){
    	detectPacingToEnclosingSurfacesNodes();
    }
    if (implicitPacking == false){
    	calculatePackingForcesExplicit3D();
    }
    if (addingRandomForces){
    	calculateRandomForces();
    }
    //cout<<"collapse check "<<endl;
    bool thereIsBinding = false;
    if (thereNodeCollapsing){
    	thereIsBinding = checkEdgeLenghtsForBindingPotentiallyUnstableElements();
    }
    //cout<<"binding check "<<endl;
    if (thereIsBinding){
		NRSolver->boundNodesWithSlaveMasterDefinition = true;
	}
   // cout<<"adhesion check "<<endl;
    if (thereIsAdhesion){
    	thereIsBinding = adhereNodes();
    }
    //cout<<"bind nodes on matrix "<<endl;
    if (thereIsBinding){
    	NRSolver->boundNodesWithSlaveMasterDefinition = true;
    }
    //cout<<" beginning NR "<<endl;
    updateStepNR();
    //cout<<" after NR "<<endl;

    /* double xx, yy, zz;
    xx = gsl_matrix_get(Elements[0]->Strain,0,0);
    yy = gsl_matrix_get(Elements[0]->Strain,1,0);
    zz = gsl_matrix_get(Elements[0]->Strain,2,0);
	cout<<"Element[0] strains: " <<xx<<" "<<yy<<" "<<zz<<endl;
    */

    calculateBoundingBox();

    //calculateColumnarLayerBoundingBox();
    //cout<<" step: "<<timestep<<" Pressure: "<<SuctionPressure[2]<<" Pa, maximum z-height: "<<boundingBox[1][2]<<" L/a: "<<(boundingBox[1][2]-50)/(2.0*pipetteRadius)<<endl;
    Success = checkFlip();
    if (thereIsArtificaialRelaxation && artificialRelaxationTime == currSimTimeSec){
		artificialRelax();
	}
    timestep++;
    currSimTimeSec += dt;
    if (Success){
    	processDisplayDataAndSave();
    }
    return Success;
}

void Simulation::assignIfElementsAreInsideEllipseBands(){
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->IsAblated && ((*itElement)->tissueType == 0 || (*itElement)->tissueType == 2)){
			//The element is not ablated, it is either a columnar or a linker element.
			//Peripodial elements are not counted in the ellipses, as the perturbations of
			// interest are not applied to them			
			(*itElement)->checkIfInsideEllipseBands(nMarkerEllipseRanges, markerEllipseBandXCentres,markerEllipseBandR1Ranges, markerEllipseBandR2Ranges, Nodes);
		}
	}
}

void Simulation::updateRelativePositionsToApicalPositioning(){
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ( (*itElement)->tissuePlacement == 1 && (*itElement)->tissueType == 0){//apical element of columnar layer
			//I have all the elements on this column stored in "elementsIdsOnSameColumn", first id of the stack
			//being the basal element and the last being the apical element. This is an apical element.
			//I will equate tissue positioning of this element to all other elements on the list:
			double* ReletivePos = new double[2];
			(*itElement)->getRelativePosInBoundingBox(ReletivePos);

			for (int i=0; i < TissueHeightDiscretisationLayers-2;++i){
				int idOfElementBelowThisApicalElement = (*itElement)->elementsIdsOnSameColumn[i];
				Elements[idOfElementBelowThisApicalElement]->setRelativePosInBoundingBox(ReletivePos[0],ReletivePos[1]);
			}
			delete[] ReletivePos;
		}
	}
}

void Simulation::checkForPinningPositionsUpdate(){
	//cout<<"checking for pinning update, currTime: "<<currSimTimeSec<<" nGrowthPinning: "<<nGrowthPinning<<" GridGrowthsPinnedOnInitialMesh? "<<GridGrowthsPinnedOnInitialMesh<<endl;
	if (GridGrowthsPinnedOnInitialMesh){
		//cout<<" grid is pinned "<<endl;
		for (int i=0; i<nGrowthPinning; ++i){
			//cout<<"growthPinUpdateTime["<<i<<"]: "<<growthPinUpdateTime[i]<<" growthPinUpdateBools["<<i<<"]: "<<growthPinUpdateBools[i]<<endl;
			if (currSimTimeSec >= growthPinUpdateTime[i] &&  !growthPinUpdateBools[i]){
				updatePinningPositions();
				growthPinUpdateBools[i] = true;
			}
		}
	}
}

void Simulation::updatePinningPositions(){
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->setInitialRelativePosInBoundingBox();
	}
}

void Simulation::flagElementsThatNeedRefinement(){
	//cout<<"inside flag for refinement"<<endl;
	double thresholdSide = 5;
	double thresholdArea = thresholdSide/2.0 * thresholdSide/2.0 *pow(3,0.5);
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for //private(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces, outputFile, dt)
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->IsAblated && (*itElement)->tissueType ==0){//Refine only the columnar layer elements.
			if ((*itElement)->tissuePlacement == 1){ // element is apical, I should check with apical area
				(*itElement)->calculateApicalArea();
				(*itElement)->doesElementNeedRefinement(thresholdArea, 1); //checking refinement with apical surface;
			}
			else if ((*itElement)->tissuePlacement == 0){ // element is basal, I should check with basal area
				(*itElement)->calculateBasalArea();
				(*itElement)->doesElementNeedRefinement(thresholdArea, 0); //checking refinement with basal surface;
			}
		}
	}

}

void Simulation::refineElements(){
	int n = Elements.size();
	for (int i=0; i<n; ++i){
		if (Elements[i]->willBeRefined){
			//I will refine all the elements in this column, they are
			//stored on elementsIdsOnSameColumn attribute of each element
			//This list starts from basal elements and moves to apical;
			//Clearing pu the flags that I may have on the other end of the tissue first:
			for (int j=0; j<TissueHeightDiscretisationLayers; ++j){
				Elements[Elements[i]->elementsIdsOnSameColumn[j]]->willBeRefined = false;
			}
			//I will add discretisation layers + 1 number of nodes, I want to keep the ids recorded for now:
			int* newNodeIdList = new int [TissueHeightDiscretisationLayers+1];
			addNodesForRefinement(Elements[i],newNodeIdList);
			addElementsForRefinement(Elements[i]->elementsIdsOnSameColumn,newNodeIdList);
			delete[] newNodeIdList;
		}
	}


}

void Simulation::addNodesForRefinement(ShapeBase* currElement, int* newNodeIdList){
	//cout<<" inside add nodes for refinement"<<endl;
	Node* tmp_nd;
	//now I will add nodes:
	//Adding the basal node first:
	int dividedElementId = currElement->elementsIdsOnSameColumn[0];
	double* pos = new double[3];
	Elements[dividedElementId]->getBasalCentre(pos);
	double* externalViscosity = Elements[dividedElementId]->getBasalMinViscosity(Nodes);
	//cout<<" basal centre: "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
	int tissuePos = 0; //at basal position
	int tissueType = 0; //columnar layer tissue, I only refine columnar layer
	int atCircumference = 0; //newly added nodes will never be at circumferenece or at symmetry borders
	tmp_nd = new Node(nNodes, 3, pos,tissuePos, tissueType);
	tmp_nd->atCircumference = atCircumference;
	for (int j=0; j<3; ++j){
		tmp_nd->externalViscosity[j] = externalViscosity[j];
	}
	Nodes.push_back(tmp_nd);
	nNodes = Nodes.size();
	newNodeIdList[0] = tmp_nd->Id;
	//Now adding the apical nodes:
	tissuePos = 2;//rest of the nodes are at mid-line, except for the last node to be added
	for (int i=0; i< TissueHeightDiscretisationLayers ;++i){
		if (i == TissueHeightDiscretisationLayers-1){
			tissuePos = 1; // the last node to be added is apical
		}
		dividedElementId = currElement->elementsIdsOnSameColumn[i];
		Elements[dividedElementId]->getApicalCentre(pos);
		//cout<<" apical centre: "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
		double* externalViscosity = Elements[dividedElementId]->getApicalMinViscosity(Nodes);
		tmp_nd = new Node(nNodes, 3, pos,tissuePos, tissueType);
		tmp_nd->atCircumference = atCircumference;
		for (int j=0; j<3; ++j){
			tmp_nd->externalViscosity[j] = externalViscosity[j];
		}
		Nodes.push_back(tmp_nd);
		nNodes = Nodes.size();
		newNodeIdList[i+1] = tmp_nd->Id;
		/*cout<<" Element: "<<(*itElement)->Id<<" column : ";
		for (int j=0; j<TissueHeightDiscretisationLayers; ++j){
			cout<<" "<<(*itElement)->elementsIdsOnSameColumn[j]<<" ";
		}
		cout<<endl;*/
	}
	delete[] pos;
}

void Simulation::addElementsForRefinement(int* elementsIdsOnColumn, int* newNodeIdList){
	//I will convert the existing element for each layer to the element using
	//corner 0-1-new basal node, 3 - 4 - new apical node. Then I will create two new
	//elements, that will use corners  1-2 and 2-0 (basally).

	for (int i=0; i< TissueHeightDiscretisationLayers ;++i){
		//cout<<"Started the layer: "<<i<<endl;
		ShapeBase* currElement = Elements[elementsIdsOnColumn[i]];
		const int dim = currElement->getDim();
		const int n = currElement->getNodeNumber();
		double** referenceOfShapeToBeRefined = currElement->getReferencePos();
		double * basalReferenceCentre = new double [dim];
		double * apicalReferenceCentre = new double [dim];
		currElement->getReferenceBasalCentre(basalReferenceCentre);
		currElement->getReferenceApicalCentre(apicalReferenceCentre);
		//cout<<" basal reference centre: "<<basalReferenceCentre[0]<<" "<<basalReferenceCentre[1]<<" "<<basalReferenceCentre[2]<<endl;
		//cout<<" apical reference centre: "<<apicalReferenceCentre[0]<<" "<<apicalReferenceCentre[1]<<" "<<apicalReferenceCentre[2]<<endl;

		//Removing the current element from the nodes it is attached to:
		for (int j=0; j<n; ++j){
			Nodes[currElement->NodeIds[j]]->removeFromConnectedElements(currElement->Id, currElement->VolumePerNode);
		}
		//adding the immediate neighbours before I change the node list of elemetn to be divided.
		//I am adding the new basal node to the lists of basal nodes of the element (0-2);
		Nodes[currElement->NodeIds[0]]->addToImmediateNeigs(newNodeIdList[i]);
		Nodes[currElement->NodeIds[1]]->addToImmediateNeigs(newNodeIdList[i]);
		Nodes[currElement->NodeIds[2]]->addToImmediateNeigs(newNodeIdList[i]);
		//I am adding the new apical node to the lists of apical nodes of the element (3-5);
		Nodes[currElement->NodeIds[3]]->addToImmediateNeigs(newNodeIdList[i+1]);
		Nodes[currElement->NodeIds[4]]->addToImmediateNeigs(newNodeIdList[i+1]);
		Nodes[currElement->NodeIds[5]]->addToImmediateNeigs(newNodeIdList[i+1]);
		//Now moving on to manipulating elements:
		int* NodeIds;
		NodeIds = new int[n];
		double** referencePos;
		referencePos = new double*[n];
		for (int j=0; j<n; ++j){
			referencePos[j] = new double[dim];
			for(int k=0; k<dim; ++k){
				referencePos[j][k] = 0.0;
			}
		}
		for (int k =0; k<3; ++k){
			int basalIdIndice0;
			int basalIdIndice1;
			int apicalIdIndice0;
			int apicalIdIndice1;
			if ( k == 0 ){
				basalIdIndice0 = 1;
				basalIdIndice1 = 2;
				apicalIdIndice0 = 4;
				apicalIdIndice1 = 5;
			}
			else if ( k == 1 ){
				basalIdIndice0 = 2;
				basalIdIndice1 = 0;
				apicalIdIndice0 = 5;
				apicalIdIndice1 = 3;
			}
			else if ( k == 2 ){
				//This is converting the main element into a smaller one!
				basalIdIndice0 = 0;
				basalIdIndice1 = 1;
				apicalIdIndice0 = 3;
				apicalIdIndice1 = 4;
			}
			NodeIds[0] = currElement->NodeIds[basalIdIndice0];
			NodeIds[1] = currElement->NodeIds[basalIdIndice1];
			NodeIds[2] = newNodeIdList[i];
			NodeIds[3] = currElement->NodeIds[apicalIdIndice0];
			NodeIds[4] = currElement->NodeIds[apicalIdIndice1];
			NodeIds[5] = newNodeIdList[i+1];
			for (int j=0; j<dim; ++j){
				referencePos[0][j] = referenceOfShapeToBeRefined[basalIdIndice0][j];
				referencePos[1][j] = referenceOfShapeToBeRefined[basalIdIndice1][j];
				referencePos[2][j] = basalReferenceCentre[j];
				referencePos[3][j] = referenceOfShapeToBeRefined[apicalIdIndice0][j];
				referencePos[4][j] = referenceOfShapeToBeRefined[apicalIdIndice1][j];
				referencePos[5][j] = apicalReferenceCentre[j];
			}
			if (k == 0 || k == 1){
				//cout<<" k is "<<k<<endl;
				Prism* PrismPnt01;
				PrismPnt01 = new Prism(NodeIds, Nodes, currElementId,thereIsPlasticDeformation);
				PrismPnt01->updateReferencePositionMatrixFromInput(referencePos);
				PrismPnt01->checkRotationConsistency3D();
				PrismPnt01->copyElementInformationAfterRefinement(currElement,TissueHeightDiscretisationLayers,thereIsPlasticDeformation);
				PrismPnt01->calculateElementShapeFunctionDerivatives();
				//Add the weight of element to its nodes:
				for (int j = 0; j<n;++j){
					Nodes[NodeIds[j]]->addToConnectedElements(PrismPnt01->Id,PrismPnt01->VolumePerNode);
				}
				Elements.push_back(PrismPnt01);
				nElements = Elements.size();
				currElementId++;
			}
			else if (k == 2){
				//cout<<" k is two, manipulating old element"<<endl;
				//converting the current prism to a smaller version.
				currElement->updateNodeIdsForRefinement(NodeIds);
				currElement->updateReferencePositionMatrixFromInput(referencePos);
				currElement->checkRotationConsistency3D();
				currElement->updatePositions(Nodes);
				currElement->GrownVolume /=3.0;
				currElement->VolumePerNode /= 3.0;
				currElement->calculateElementShapeFunctionDerivatives();
				//Add the weight of element to its nodes, I have removed this element from its nodes above:
				for (int j = 0; j<n;++j){
					Nodes[currElement->NodeIds[j]]->addToConnectedElements(currElement->Id,currElement->VolumePerNode);
				}
				//cout<<"finalised manipulating old element: "<<currElement->Id<<endl;
			}
		}
		//clean up:
		delete[] NodeIds;
		for (int j=0; j<n; ++j){
			delete[] referencePos[j];
		}
		delete[] referencePos;
		delete[] apicalReferenceCentre;
		delete[] basalReferenceCentre;
	}
}

bool Simulation::checkFlip(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->isFlipped){
			//there is a flipped element:
			outputFile<<"There is A flipped element: "<<(*itElement)->Id<<endl;
			return false;
		}
	}
	return true;
}
void Simulation::wrapUpAtTheEndOfSimulation(){
    alignTissueDVToXPositive();
    //alignTissueAPToXYPlane();
}


void Simulation::clearProjectedAreas(){
    for (int i=0;i<nNodes; ++i){
        Nodes[i]->zProjectedArea = 0.0;
    }
}

void Simulation::correctzProjectedAreaForMidNodes(){
    for (int i=0;i<nNodes; ++i){
        if (Nodes[i]->tissuePlacement == 2 ){ // the node is on midlayer
            //For the midline nodes, the area is added from apical and basal surfaces of elemetns on both sides.
            Nodes[i]->zProjectedArea /= 2.0;
        }
    }
}

void Simulation::calculateZProjectedAreas(){
    clearProjectedAreas();
    vector<ShapeBase*>::iterator itElement;
    for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	(*itElement)->calculateZProjectedAreas();
        (*itElement)->assignZProjectedAreas(Nodes);
    }
    correctzProjectedAreaForMidNodes();
}

void Simulation::updatePlasticDeformation(){
	//double rate = plasticDeformationRate/3600.0*dt; //convert from per hour to per de(in sec)
	#pragma omp parallel for //private(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces, outputFile, dt)
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->IsAblated ){
			if (thereIsExplicitECM &&  (*itElement)->isECMMimicing){
				//The simulation is defining an explicit ECM layer. All basal elements
				//should be subject to non-volum econserving plastic deformation.
				//The volume conservation should be false (first inout)
				//the half time used should be the ecm remodelling half time
				//there should not be thresholds for z remodelling, there is no
				//volume conservation, therefore these are not relevant anyway.
				(*itElement)->calculatePlasticDeformation3D(false,dt,ECMRenawalHalfLife, 0.1, 10.0);
			}
			else if (thereIsPlasticDeformation){
				//The ECM will always have plastic deformation, hence this function will be called.
				//Then for all other elements, I should first check if there is plastic deformation in the first place
				//Then I will check the tissue types
				//the lateral elements will have plastic deformation if either columnar or peripodial elemetns have the deformation.
				if(( ((*itElement)->tissueType == 0 || (*itElement)->tissueType == 2) && plasticDeformationAppliedToColumnar) || ( ((*itElement)->tissueType == 1 || (*itElement)->tissueType == 2) && plasticDeformationAppliedToPeripodial))
				{
					(*itElement)->calculatePlasticDeformation3D(volumeConservedInPlasticDeformation,dt,plasticDeformationHalfLife, zRemodellingLowerThreshold, zRemodellingUpperThreshold);
				}
			}
		}
		else{
			(*itElement)->setPlasticDeformationIncrement(1.0,1.0,1.0);
		}
	}
}

void Simulation::calculateNumericalJacobian(bool displayMatricesDuringNumericalCalculation){
	int dim = 3;
	gsl_matrix_set_zero(NRSolver->Knumerical);
	NRSolver->calculateDisplacementMatrix(dt);
	//PACKING SHOULD BE ADDED HERE If using this nuerical calculation!!!

	//Trying to see the manual values:
	resetForces(true); // reset the packing forces together with all the rest of the forces here
	//No perturbation:
	gsl_matrix* ge_noPerturb = gsl_matrix_calloc(dim*nNodes,1);
	gsl_matrix* gvInternal_noPerturb = gsl_matrix_calloc(dim*nNodes,1);
	NRSolver->calculateForcesAndJacobianMatrixNR(Nodes, Elements, dt, recordForcesOnFixedNodes, FixedNodeForces);
	NRSolver->writeForcesTogeAndgvInternal(Nodes, Elements, SystemForces);
	gsl_matrix_memcpy(ge_noPerturb, NRSolver->ge);
	gsl_matrix_memcpy(gvInternal_noPerturb, NRSolver->gvInternal);
	//perturbation loop:
	gsl_matrix* uk_original = gsl_matrix_calloc(dim*nNodes,1);
	gsl_matrix_memcpy(uk_original,NRSolver->uk);
	for (int i=0; i<nNodes; ++i){
		for (int j=0; j<3; ++j){
	        if (implicitPacking){
	            resetForces(true);	// reset packing forces
	        }
	        else{
	            resetForces(false);	// do not reset packing forces
	        }
	        gsl_matrix_set(NRSolver->uk,i*3+j,0,gsl_matrix_get(NRSolver->uk,i*3+j,0)+1E-6);
			NRSolver->calculateDisplacementMatrix(dt);
			Nodes[i]->Position[j] += 1E-6;
			updateElementPositions();
			NRSolver->calculateForcesAndJacobianMatrixNR(Nodes, Elements, dt, recordForcesOnFixedNodes, FixedNodeForces);
			gsl_matrix* ge_withPerturb = gsl_matrix_calloc(dim*nNodes,1);
			gsl_matrix* gvInternal_withPerturb = gsl_matrix_calloc(dim*nNodes,1);
			NRSolver->writeForcesTogeAndgvInternal(Nodes, Elements, SystemForces);
			gsl_matrix_memcpy(ge_withPerturb, NRSolver->ge);
			gsl_matrix_memcpy(gvInternal_withPerturb, NRSolver->gvInternal);
			//Calculate dg/dx:
			gsl_matrix_sub(ge_withPerturb,ge_noPerturb);
			gsl_matrix_sub(gvInternal_withPerturb,gvInternal_noPerturb);
			gsl_matrix_scale(ge_withPerturb,1.0/1E-6);
			gsl_matrix_scale(gvInternal_withPerturb,1.0/1E-6);
			for (int k=0; k<nNodes*3; ++k){
				double valueElastic =   0;//gsl_matrix_get(ge_withPerturb,k,0);
				double valueViscous = 	gsl_matrix_get(gvInternal_withPerturb,k,0);
				double value = valueElastic + valueViscous;
				value *= -1.0;
				gsl_matrix_set(NRSolver->K,i*3+j,k,value);
			}
			gsl_matrix_memcpy(NRSolver->uk,uk_original);
			//gsl_matrix_set(uk,i*3+j,0,gsl_matrix_get(uk,i*3+j,0)-1E-6);
			Nodes[i]->Position[j] -= 1E-6;
			updateElementPositions();
			gsl_matrix_free(ge_withPerturb);
			gsl_matrix_free(gvInternal_withPerturb);
		}
	}
	NRSolver->calcutateFixedK(Nodes);
	if (displayMatricesDuringNumericalCalculation){
		Elements[0]->displayMatrix(NRSolver->K,"numericalK");
	}
	gsl_matrix_memcpy(NRSolver->Knumerical,NRSolver->K);
	NRSolver->setMatricesToZeroInsideIteration();
	gsl_matrix_free(ge_noPerturb);
	gsl_matrix_free(gvInternal_noPerturb);
	gsl_matrix_free(uk_original);
}

void Simulation::conserveColumnVolume(){
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	//#pragma omp parallel for //private(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces, outputFile, dt)
    for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	vector <int> elementIdsForRedistribution;
    	if ((*itElement)->tissueType == 0 && (*itElement)->tissuePlacement == 1){
    		//columnar apical element, take the elements on same column
    		for (int i=0; i < TissueHeightDiscretisationLayers;++i){
    			int idOfElementOnSameColumn = (*itElement)->elementsIdsOnSameColumn[i];
    			if (Elements[idOfElementOnSameColumn]->tissueType != 0){
    				//must be columnar
					continue;
				}
    			if (Elements[idOfElementOnSameColumn]->isECMMimicing){
    				//exclude ECM
    				continue;
    			}
    			if (Elements[idOfElementOnSameColumn]->isActinMimicing){
    				//exclude top actin layer if there is one
					continue;
				}
    			elementIdsForRedistribution.push_back(idOfElementOnSameColumn);
    		}
    		//now I have a list of elements that I would like to redistribute the volumes of:
    		//get sum of ideal volumes:
    		const int n = elementIdsForRedistribution.size();
    		double ratioOfCurrentVolumeToIdeal[n];
    		double sumIdealVolumes = 0;
    		double sumCurrentVolumes = 0;
    		for (int i=0; i < n;++i){
    			double idealVolume = Elements[elementIdsForRedistribution[i]]->GrownVolume;
    			sumIdealVolumes += idealVolume;
    			//get current volume from average detF:
    			double currentVolume = Elements[elementIdsForRedistribution[i]]->getCurrentVolume();
    			if (currentVolume < 10E-10){
    				currentVolume = idealVolume;
    			}
    			ratioOfCurrentVolumeToIdeal[i]= currentVolume/idealVolume;
    			sumCurrentVolumes += currentVolume;
    			if (elementIdsForRedistribution[i] == 7009 || elementIdsForRedistribution[i] == 5207 || elementIdsForRedistribution[i] == 4622 || elementIdsForRedistribution[i] == 2820){
    				cout<<" redistributing element: "<<elementIdsForRedistribution[i]<<" idealVolume: "<<idealVolume<<" currentVolume: "<<currentVolume<<" ratio: "<<ratioOfCurrentVolumeToIdeal[i]<<endl;
    			}
    		}
    		double redistributionFactors[n];
    		for (int i=0; i < n;++i){
    			//calculate redistribution factors:
    			redistributionFactors[i] = sumIdealVolumes/sumCurrentVolumes *ratioOfCurrentVolumeToIdeal[i];
    			if (elementIdsForRedistribution[i] == 7009 || elementIdsForRedistribution[i] == 5207 || elementIdsForRedistribution[i] == 4622 || elementIdsForRedistribution[i] == 2820){
    				cout<<" redistributing element: "<<elementIdsForRedistribution[i]<<" redist factor: "<<redistributionFactors[i]<<endl;
    			}
    			//calculateGrowth:
    			double uniformGrowthIncrement = pow(redistributionFactors[i],1.0/3.0);
    			gsl_matrix* columnarFgIncrement = gsl_matrix_calloc(3,3);
    			gsl_matrix* peripodialFgIncrement = gsl_matrix_calloc(3,3);
    			gsl_matrix_set_identity(columnarFgIncrement);
    			gsl_matrix_set_identity(peripodialFgIncrement);
    			gsl_matrix_set(columnarFgIncrement,0,0,uniformGrowthIncrement);
    			gsl_matrix_set(columnarFgIncrement,1,1,uniformGrowthIncrement);
    			gsl_matrix_set(columnarFgIncrement,2,2,uniformGrowthIncrement);
    			Elements[elementIdsForRedistribution[i]]->updateGrowthIncrement(columnarFgIncrement,peripodialFgIncrement);
    			gsl_matrix_free(columnarFgIncrement);
    			gsl_matrix_free(peripodialFgIncrement);
    		}
    	}
    }
}

void Simulation::updateStepNR(){
    int iteratorK = 0;
    int maxIteration =20;
    bool converged = false;

    bool numericalCalculation = false;
    bool displayMatricesDuringNumericalCalculation = false;
    bool useNumericalKIncalculation = false;

    //double clock0 = ( std::clock() - simulationStartClock ) / (double) CLOCKS_PER_SEC;
    NRSolver->setMatricesToZeroAtTheBeginningOfIteration(numericalCalculation);
    NRSolver->constructUnMatrix(Nodes);
    //NRSolver->constructLumpedMassExternalViscosityMatrix(Nodes);
    NRSolver->initialteUkMatrix();
    while (!converged){
        cout<<"iteration: "<<iteratorK<<endl;
        if (implicitPacking){
            resetForces(true);	// reset packing forces
        }
        else{
            resetForces(false);	// do not reset packing forces
        }
        NRSolver->setMatricesToZeroInsideIteration();
        if (numericalCalculation){
        	calculateNumericalJacobian(displayMatricesDuringNumericalCalculation);
        }
        NRSolver->calculateDisplacementMatrix(dt);
        NRSolver->calculateForcesAndJacobianMatrixNR(Nodes, Elements, dt, recordForcesOnFixedNodes, FixedNodeForces);
	    //Writing elastic Forces and elastic Ke:
		NRSolver->writeForcesTogeAndgvInternal(Nodes, Elements, SystemForces);
	    NRSolver->writeImplicitElementalKToJacobian(Elements);
	    if (numericalCalculation){
			NRSolver->calculateDifferenceBetweenNumericalAndAnalyticalJacobian(Nodes, displayMatricesDuringNumericalCalculation);
			if(useNumericalKIncalculation){
				NRSolver->useNumericalJacobianInIteration();
			}
		}
		NRSolver->calculateExternalViscousForcesForNR(Nodes);
	    NRSolver->addImplicitKViscousExternalToJacobian(Nodes,dt);
	    //These are the calculation of packing forces that would work if I wanted implicit packing
	    if (implicitPacking){
	    	//if (!thereIsAdhesion){
	    		calculatePackingForcesImplicit3D();
	    		calculatePackingJacobian3D(NRSolver->K);
	    		//calculatePackingToAFMBeadJacobian3D(NRSolver->K);
	    	//}
	        if (encloseTissueBetweenSurfaces){
	        	calculatePackingForcesToEnclosingSurfacesImplicit3D();
	        	calculatePackingToEnclosingSurfacesJacobian3D(NRSolver->K);
	        }
	    }
	    //End of packing forces.
        //Now I will check if there are any nodes with zero mass, then I will be able to fill in the zero K matrix with identity if necessary.
        NRSolver->checkJacobianForAblatedNodes(AblatedNodes);
        NRSolver->calculateSumOfInternalForces();
        if (PipetteSuction && timestep >= PipetteInitialStep){
			packToPipetteWall();
			calculateZProjectedAreas();
			addPipetteForces(NRSolver->gExt);
		}
		addMyosinForces(NRSolver->gExt);
		//packing can come from both encapsulation and tissue-tissue packing. I ad the forces irrespective of adhesion.
		addPackingForces(NRSolver->gExt);

		//------------------
/*
		int nFibres = 1;
		vector<double> orientationAngles;
		vector <double> fibreForceMagnitudes;
		vector <double> initTimesInHr;
		vector <double > matureTimesinHr;
		vector<double> initAngle;
		vector<double> endAngle;
		orientationAngles.push_back(30.0);
		fibreForceMagnitudes.push_back(1000.0);
		initTimesInHr.push_back(2.0);
		matureTimesinHr.push_back(4.0);
		initAngle.push_back(0.0);
		endAngle.push_back(45.0);
		MuscleFibres* fibres;
		fibres = new MuscleFibres(nFibres, orientationAngles, fibreForceMagnitudes,initTimesInHr, matureTimesinHr, initAngle, endAngle);
		fibres->assignNodesForFibreAttachement(Nodes, symmetricX, symmetricY, SystemCentre[0], SystemCentre[1]);
		fibres->addFibreForces(Nodes, NRSolver->gExt, currSimTimeSec/3600.0);
		delete fibres;
*/
		//------------------


		if (addingRandomForces){
			addRandomForces(NRSolver->gExt);
		}
        NRSolver->addExernalForces();
        checkForExperimentalSetupsWithinIteration();
    	NRSolver->calcutateFixedK(Nodes);
        //Elements[0]->displayMatrix(NRSolver->NSlaveMasterBinding,"NSlaveMasterBinding before fixing");
        //Elements[0]->displayMatrix(NRSolver->ISlaveMasterBinding,"ISlaveMasterBinding before fixing");
        //Elements[0]->displayMatrix(NRSolver->K,"K before fixing");
        //Elements[0]->displayMatrix(NRSolver->gSum,"gSum before fixing");
    	NRSolver->calculateBoundKWithSlavesMasterDoF();
        //Elements[0]->displayMatrix(NRSolver->K,"K after fixing");
        //Elements[0]->displayMatrix(NRSolver->gSum,"gSum after fixing");
        //cout<<"displaying the jacobian after all additions"<<endl;
        //Elements[0]->displayMatrix(NRSolver->K,"theJacobian");
        //cout<<"checking convergence with forces"<<endl;
        //converged = NRSolver->checkConvergenceViaForce();
        if (converged){
            break;
        }
        //cout<<"solving for deltaU"<<endl;
        //Elements[0]->displayMatrix(NRSolver->K,"K");
        NRSolver->solveForDeltaU();
        //cout<<"checking convergence"<<endl;
        converged = NRSolver->checkConvergenceViaDeltaU();
        //Elements[0]->displayMatrix(NRSolver->deltaU,"deltaU after fixing");

        NRSolver->updateUkInIteration();
        updateElementPositionsinNR(NRSolver->uk);
        updateNodePositionsNR(NRSolver->uk);
        iteratorK ++;
        if (!converged && iteratorK > maxIteration){
            cerr<<"Error: did not converge!!!"<<endl;
            converged = true;
        }
    }
    checkForExperimentalSetupsAfterIteration();
    //Now the calculation is converged, I update the node positions with the latest positions uk:
    updateNodePositionsNR(NRSolver->uk);
     //Element positions are already up to date.
    cout<<"finished run one step"<<endl;
    if (PipetteSuction){
    	//find the max z:
    	double zMax = -10000;
    	int idMax = -10;
    	for (vector<Node*>::iterator itNode = Nodes.begin(); itNode < Nodes.end(); ++itNode){
    		if ((*itNode)->Position[2] > zMax){
    			zMax = (*itNode)->Position[2];
    			idMax = (*itNode)->Id;
    		}
    	}
    	cout<<"Pipette suction: "<<SuctionPressure[2]<<" max suction: "<<zMax<<" from node "<<idMax<<endl;
    }
}

void Simulation::calculateRandomForces(){
	randomForces.clear();
	randomForces=RandomGenerator::Obj().getNormRV( randomForceMean,randomForceVar, 3*nNodes );
	//making the sum of forces zero:
	double sumRandomForceX = 0.0;
	double sumRandomForceY = 0.0;
	double sumRandomForceZ = 0.0;

	vector<double>::iterator itDouble;
	for (itDouble =randomForces.begin(); itDouble < randomForces.end(); itDouble = itDouble+3){
		sumRandomForceX+= (*itDouble);
		sumRandomForceY+= (*(itDouble+1));
		sumRandomForceZ+= (*(itDouble+2));
	}
	sumRandomForceX /= nNodes;
	sumRandomForceY /= nNodes;
	sumRandomForceZ /= nNodes;

	for (itDouble =randomForces.begin(); itDouble < randomForces.end(); itDouble = itDouble+3){
		(*itDouble) -= sumRandomForceX;
		(*(itDouble+1)) -= sumRandomForceY;
		(*(itDouble+2)) -= sumRandomForceZ;
	}
}

void Simulation::addRandomForces(gsl_matrix* gExt){
		for (int j=0; j<3*nNodes; ++j){
			double F = randomForces[j];
			F += gsl_matrix_get(gExt,j,0);
			gsl_matrix_set(gExt,j,0,F);
		}
}

void Simulation::updateNodePositionsNR(gsl_matrix* uk){
    int dim = 3;
    for (int i = 0; i<nNodes; ++i){
    	for (int j=0; j<dim; ++j){
            Nodes[i]->Position[j]=gsl_matrix_get(uk,dim*i+j,0);
        }
    }
    //cout<<"finised node pos update"<<endl;
}

void Simulation::updateElementPositionsinNR(gsl_matrix* uk){
    int dim = 3;
    for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        int* nodeIds = (*itElement)->getNodeIds();
        int nNodes= (*itElement)->getNodeNumber();
        for (int j=0; j<nNodes; ++j){
            double x = gsl_matrix_get(uk,dim*nodeIds[j],0);
            double y = gsl_matrix_get(uk,dim*nodeIds[j]+1,0);
            double z = gsl_matrix_get(uk,dim*nodeIds[j]+2,0);
            (*itElement)->Positions[j][0] = x;
            (*itElement)->Positions[j][1] = y;
            (*itElement)->Positions[j][2] = z;
        }
    }
}

void Simulation::calculatePackingForcesExplicit3D(){
	//cout<<"inside calculatePackingForcesExplicit3D, size of packing node couples: "<<pacingNodeCouples0.size()<<endl;
	int n = pacingNodeCouples0.size();
	for(int i = 0 ; i<n; ++i){
		if (pacingNodeCouplesHaveAdhered[i]){
			continue;
		}
		int id0 = pacingNodeCouples0[i];
		int id1 = pacingNodeCouples1[i];

		double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];

		double d = pow((dx*dx + dy*dy + dz*dz),0.5);

		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		double sigmoidSaturation = 5;

		double F = 5000.0 * averageMass / (1 + exp(sigmoidSaturation* d / packingThreshold));

		double Fx = F * initialWeightPointx[i];
		double Fy = F * initialWeightPointy[i];
		double Fz = F * initialWeightPointz[i];
		PackingForces[id0][0] += Fx;
		PackingForces[id0][1] += Fy;
		PackingForces[id0][2] += Fz;
		PackingForces[id1][0] -= Fx;
		PackingForces[id1][1] -= Fy;
		PackingForces[id1][2] -= Fz;
		//if (id0 == 3444 || id0 == 3444 || id1 == 3444 || id1 == 3444){
			//cout<<"id0: "<<id0<<" id1: "<<id1<<endl;
			//cout<<" dx: "<<dx<<" dy: "<<dy<<" dz: "<<dz<<" d: "<<d<<endl;
			//cout<<" Fz: "<<Fx<<" Fy: "<<Fy<<" Fz: "<<Fz<<" F: "<<F<<endl;
			//cout<<" initialWeightPoins: "<<initialWeightPointx[i]<<" "<<initialWeightPointy[i]<<" "<<initialWeightPointz[i]<<endl;
			//cout<<" PackingForces["<<id0<<"]: "<<PackingForces[id0][0]<<" "<<PackingForces[id0][1]<<" "<<PackingForces[id0][2]<<endl;
			//cout<<" PackingForces["<<id1<<"]: "<<PackingForces[id1][0]<<" "<<PackingForces[id1][1]<<" "<<PackingForces[id1][2]<<endl;
		//}
	}
}

void Simulation::calculatePackingForcesToEnclosingSurfacesImplicit3D(){
	int nPositive=nodesPackingToPositiveSurface.size();
	int nNegative=nodesPackingToNegativeSurface.size();
	double multiplier = packingMultiplier;
	double sigmoidSaturation = sigmoidSaturationForPacking;
	for(int i = 0 ; i<nNegative; ++i){
		int id0 = nodesPackingToNegativeSurface[i];
		double dz = Nodes[id0]->Position[2] - zEnclosementBoundaries[0];
		if (initialWeightPackingToNegativeSurface[i]>0){
			dz *= -1.0;
		}
		double mass = Nodes[id0]->mass;
		double Fz = multiplier * mass / (1 + exp(sigmoidSaturation / packingToEnclosingSurfacesThreshold * (-1.0*dz)));
		Fz *= initialWeightPackingToNegativeSurface[i];
		PackingForces[id0][2] += Fz;
	}
	for(int i = 0 ; i<nPositive; ++i){
		int id0 = nodesPackingToPositiveSurface[i];
		double dz = Nodes[id0]->Position[2] - zEnclosementBoundaries[1];
		if (initialWeightPackingToPositiveSurface[i]>0){
			dz *= -1.0;
		}
		double mass = Nodes[id0]->mass;
		double Fz = multiplier * mass / (1 + exp(sigmoidSaturation / packingToEnclosingSurfacesThreshold * (-1.0*dz)));
		Fz *= initialWeightPackingToPositiveSurface[i];
		PackingForces[id0][2] += Fz;
	}
}

void Simulation::calculatePackingToEnclosingSurfacesJacobian3D(gsl_matrix* K){
	int nPositive=nodesPackingToPositiveSurface.size();
	int nNegative=nodesPackingToNegativeSurface.size();
	double sigmoidSaturation = sigmoidSaturationForPacking;
	double multiplier = packingMultiplier;
	#pragma omp parallel for
	for(int i = 0 ; i<nPositive; ++i){
		int id0 =  nodesPackingToPositiveSurface[i];
		//sigmoid test:
		double dz = Nodes[id0]->Position[2] - zEnclosementBoundaries[1];
		double mass = Nodes[id0]->mass;
		if (initialWeightPackingToPositiveSurface[i]>0){
			dz *= -1.0;
		}
		double sigmoidz =  1 / (1 + exp(sigmoidSaturation/ packingToEnclosingSurfacesThreshold * (-1.0 * dz) ));
		double dFzdz0 = sigmoidz * (1 - sigmoidz) * multiplier * mass * initialWeightPackingToPositiveSurface[i] * (sigmoidSaturation/packingToEnclosingSurfacesThreshold);
		if (initialWeightPackingToPositiveSurface[i]>0){
			dFzdz0 *= -1.0;
		}
		//z values:
		double value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value -= dFzdz0;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);
	}
	#pragma omp parallel for
	for(int i = 0 ; i<nNegative; ++i){
		int id0 =  nodesPackingToNegativeSurface[i];
		//sigmoid test:
		double dz = Nodes[id0]->Position[2] - zEnclosementBoundaries[0];
		double mass = Nodes[id0]->mass;
		if (initialWeightPackingToNegativeSurface[i]>0){
			dz *= -1.0;
		}
		double sigmoidz =  1 / (1 + exp(sigmoidSaturation/ packingToEnclosingSurfacesThreshold * (-1.0 * dz) ));
		double dFzdz0 = sigmoidz * (1 - sigmoidz) * multiplier * mass * initialWeightPackingToNegativeSurface[i] * (sigmoidSaturation/packingToEnclosingSurfacesThreshold);
		if (initialWeightPackingToNegativeSurface[i]>0){
			dFzdz0 *= -1.0;
		}
		//z values:
		double value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value -= dFzdz0;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);
	}
}

void Simulation::calculatePackingToAFMBeadJacobian3D(gsl_matrix* K){
	int n=nodesPackingToBead.size();
	double sigmoidSaturation = sigmoidSaturationForPacking;
	double multiplier = packingMultiplier;
	#pragma omp parallel for
	for(int i = 0 ; i<n; ++i){
		int id0 =  nodesPackingToBead[i];
		double dGap = distanceToBead[i];
		double dx = initialWeightPackingToBeadx[i]*dGap;
		double dy = initialWeightPackingToBeady[i]*dGap;
		double dz = initialWeightPackingToBeadz[i]*dGap;
		double mass = Nodes[id0]->mass;
		if (initialWeightPackingToBeadx[i]>0){
			dx *= -1.0;
		}
		if (initialWeightPackingToBeady[i]>0){
			dy *= -1.0;
		}
		if (initialWeightPackingToBeadz[i]>0){
			dz *= -1.0;
		}
		double sigmoidx =  1 / (1 + exp(sigmoidSaturation/ packingToBeadThreshold * (-1.0 * dx) ));
		double dFxdx0 = sigmoidx * (1 - sigmoidx) * multiplier * mass * initialWeightPackingToBeadx[i] * (sigmoidSaturation/packingToBeadThreshold);
		if (initialWeightPackingToBeadx[i]>0){
			dFxdx0 *= -1.0;
		}
		double sigmoidy =  1 / (1 + exp(sigmoidSaturation/ packingToBeadThreshold * (-1.0 * dy) ));
		double dFydy0 = sigmoidy * (1 - sigmoidy) * multiplier * mass * initialWeightPackingToBeady[i] * (sigmoidSaturation/packingToBeadThreshold);
		if (initialWeightPackingToBeady[i]>0){
			dFydy0 *= -1.0;
		}
		double sigmoidz =  1 / (1 + exp(sigmoidSaturation/ packingToBeadThreshold * (-1.0 * dz) ));
		double dFzdz0 = sigmoidz * (1 - sigmoidz) * multiplier * mass * initialWeightPackingToBeadz[i] * (sigmoidSaturation/packingToBeadThreshold);
		if (initialWeightPackingToBeadz[i]>0){
			dFzdz0 *= -1.0;
		}
		//x values:
		double value = gsl_matrix_get(K,3*id0,3*id0);
		value -= dFxdx0;
		gsl_matrix_set(K,3*id0,3*id0,value);

		//y values:
		value = gsl_matrix_get(K,3*id0+1,3*id0+1);
		value -= dFydy0;
		gsl_matrix_set(K,3*id0+1,3*id0+1,value);

		//z values:
		value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value -= dFzdz0;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);
	}
}

bool Simulation::areNodesOnNeighbouingElements(int masterNoeId, int slaveNodeId){
	bool neigElements = false;
	int nOwnersMaster = Nodes[masterNoeId]->connectedElementIds.size();
	int nOwnersSlave = Nodes[slaveNodeId]->connectedElementIds.size();
	for (int i =0; i<nOwnersMaster; ++i){
		int idOwnerMaster = Nodes[masterNoeId]->connectedElementIds[i];
		int nNodesMasterOwner = Elements[idOwnerMaster]->getNodeNumber();
		int*  nodeIdsMasterOwner = Elements[idOwnerMaster]->getNodeIds();
		for (int j =0; j<nOwnersSlave; ++j){
			int idOwnerSlave = Nodes[slaveNodeId]->connectedElementIds[j];
			int counter = 0;
			int nNodesSlaveOwner = Elements[idOwnerSlave]->getNodeNumber();
			int*  nodeIdsSlaveOwner = Elements[idOwnerSlave]->getNodeIds();
			for (int iterator_nodesMaster=0; iterator_nodesMaster<nNodesMasterOwner; ++iterator_nodesMaster){
				for (int iterator_nodesSlave=0; iterator_nodesSlave<nNodesSlaveOwner; ++iterator_nodesSlave){
					if (nodeIdsMasterOwner[iterator_nodesMaster] == nodeIdsSlaveOwner[iterator_nodesSlave]){
						counter++;
					}
					if (counter>2){
						neigElements = true;
						return neigElements;
					}
				}
			}
		}
	}
	return neigElements;
}

void Simulation::manualAdhesion(int masterNodeId,int slaveNodeId){
	pacingNodeCouples0.push_back(masterNodeId);
	pacingNodeCouples1.push_back(slaveNodeId);
	pacingNodeCouplesHaveAdhered.push_back(false);
}

bool Simulation::isAdhesionAllowed(int masterNodeId, int slaveNodeId){
	//if node has an ongoing movement due to staged collapse, do not adhere further:
	if(Nodes[slaveNodeId]->positionUpdateOngoing || Nodes[masterNodeId]->positionUpdateOngoing){
		//do not adhere columnar to peripodial sections
		//cout<<"position update ongoing"<<endl;
		return false;
	}
	//Do not adhere element at circumference, they have further limitations on them in most cases
	if(Nodes[slaveNodeId]->tissueType != Nodes[masterNodeId]->tissueType){
		//do not adhere columnar to peripodial sections
		//cout<<"tissue types different!"<<endl;
		return false;
	}
	//the ecm mimicking nodes at the circumference should only adhere to other circumrefernce nodes
	//otherwise, the thin nodes will adhere to nodes further away. This is to allow for adhesion of two outer surfaces on the lateral side.
	if(  (Nodes[slaveNodeId]->atCircumference && !Nodes[masterNodeId]->atCircumference)
		||(!Nodes[slaveNodeId]->atCircumference && Nodes[masterNodeId]->atCircumference) )
	{
		//cout<<"one is circumferential other is not"<<endl;
		return false;
	}
	//do not adhere lateral elements, they are too thin in most cases.
	//if(Nodes[slaveNodeId]->hasLateralElementOwner || Nodes[masterNodeId]->hasLateralElementOwner){
	//	continue;
	//}
	//one to one adhesion is linked to collapse now, if I dont have collapse, I dont limit to one-to-one adhesion
	if(collapseNodesOnAdhesion){
		if (Nodes[slaveNodeId]->adheredTo > -1 || Nodes[masterNodeId]->adheredTo > -1) {
			//cout<<"one is already adhered"<<endl;
			return false;
		}
	}
	//Need to check these again as the lists are updated since the last step
	//check if the master node is in the collapsed list?
	if (binary_search(Nodes[masterNodeId]->collapsedWith.begin(), Nodes[masterNodeId]->collapsedWith.end(),slaveNodeId)){
		//the couple is already collapsed:
		//cout<<"Master/slave couple: "<<masterNodeId<<"/"<<slaveNodeId<<" already collapsed"<<endl;
		//cout<<"master node is in the collapsed list"<<endl;
		return false;
	}
	//does this node belong to a neighbour element of mine?
	bool nodeIsOnNeigElement = areNodesOnNeighbouingElements(masterNodeId,slaveNodeId);
	if (nodeIsOnNeigElement){
		//cout<<"node are neighbours"<<endl;
		return false;
	}
	bool collapedNodeIsNeig = false;
	//check if the collapsed nodes of master are neigs of slave
	collapedNodeIsNeig = Nodes[masterNodeId]->isNeigWithMyCollapsedNodes(slaveNodeId,Nodes);
	if (collapedNodeIsNeig){
		//cout<<"Master/slave couple: "<<masterNodeId<<"/"<<slaveNodeId<<" already neig (masters collapsed is neig of slave)"<<endl;
		//cout<<"node are collapsed with neighbours (master)"<<endl;
		return false;
	}
	collapedNodeIsNeig = Nodes[slaveNodeId]->isNeigWithMyCollapsedNodes(masterNodeId,Nodes);
	if (collapedNodeIsNeig){
		//cout<<"Master/slave couple: "<<masterNodeId<<"/"<<slaveNodeId<<" already neig (slaves collapsed is neig of master)"<<endl;
		//cout<<"node are collapsed with neighbours (slave)"<<endl;
		return false;
	}
	return true;
}

double Simulation::distanceSqBetweenNodes(int id0, int id1){
	double dx = Nodes[id0] ->Position[0] - Nodes[id1] ->Position[0];
	double dy = Nodes[id0] ->Position[1] - Nodes[id1] ->Position[1];
	double dz = Nodes[id0] ->Position[2] - Nodes[id1] ->Position[2];
	double dSq = dx*dx + dy*dy + dz*dz;
	return dSq;
}

bool Simulation::checkForElementFlippingUponNodeCollapse(vector<int> &newCollapseList, double* avrPos){
	//check for element flipping after adhesion collapse:
	//loop over elements:
	int nCollapsedList = newCollapseList.size();
	bool elementWillCollapse = false;
	for (int nodeIterator=0; nodeIterator<nCollapsedList;++nodeIterator){
		int currNodeId = newCollapseList[nodeIterator];
		int nElements=  Nodes[currNodeId]->connectedElementIds.size();
		for (int i=0; i<nElements;++i){
			//cout<<"checking element : "<< Nodes[currNodeId]->connectedElementIds[i]<<" from node : "<< currNodeId<<endl;
			elementWillCollapse = Elements[Nodes[currNodeId]->connectedElementIds[i]]->isElementFlippedInPotentialNewShape(currNodeId, avrPos[0], avrPos[1], avrPos[2]);
			if (elementWillCollapse){
				cout<<" element "<<Nodes[currNodeId]->connectedElementIds[i]<<" will flip if adhered"<<endl;
				return false;
			}
		}
	}
	return true;
}


bool Simulation::updatePositionsOfNodesCollapsingInStages(){
	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->positionUpdateOngoing){
			//cout<<"updating node "<<(*itNode)->Id<<endl;
			double* avrPos = new double[3];
			bool* fix = new bool[3];
			avrPos[0] = 0; avrPos[1]=0; avrPos[2]=0;
			fix[0]=false; fix[1]=false; fix[2]=false;
			int n = (*itNode)->collapsedWith.size();
			for (int i = 0; i<n; ++i){
				for (int j=0; j<3; ++j){
					avrPos[j] += Nodes[(*itNode)->collapsedWith[i]]->Position[j]/n;
				}
			}
			//cout<<" average position: "<<avrPos[0]<<" "<<avrPos[1]<<" "<<avrPos[2]<<endl;
			for (int i = 0; i<n; ++i){
				for (int j=0; j<3; ++j){
					if(Nodes[(*itNode)->collapsedWith[i]]->FixedPos[j]){
						avrPos[j] = Nodes[(*itNode)->collapsedWith[i]]->Position[j];
						fix[j] = true;
					}
				}
			}
			//cout<<" fix pos : "<<fix[0]<<" "<<fix[1]<<" "<<fix[2]<<endl;
			(*itNode)->updatePositionTowardsPoint(avrPos,fix);
			//cout<<"updated towards average"<<endl;
			delete[] avrPos;
			delete[] fix;
		}

	}
}

bool Simulation::adhereNodes(){
	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		 (*itNode)->clearDuplicatesFromCollapseList();
	}
	int n = pacingNodeCouples0.size();
	int dim = 3;
	bool thereIsBinding = false;
	vector<bool> adhesionAllowedList;
	//I am constructing a list to track if the adhesion between the recorded couples would be allowed
	for(int nodeCoupleIterator = 0 ; nodeCoupleIterator<n; ++nodeCoupleIterator){
		bool adhesionAllowed = isAdhesionAllowed(pacingNodeCouples0[nodeCoupleIterator], pacingNodeCouples1[nodeCoupleIterator]);
		adhesionAllowedList.push_back(adhesionAllowed);
		//cout<<" pair : "<<pacingNodeCouples0[nodeCoupleIterator]<<" "<<pacingNodeCouples1[nodeCoupleIterator]<<" adhesion feasible? "<<adhesionAllowed<<endl;
	}
	for(int nodeCoupleIterator = 0 ; nodeCoupleIterator<n; ++nodeCoupleIterator){
		if (!adhesionAllowedList[nodeCoupleIterator]){
			cout<<" pair : "<<pacingNodeCouples0[nodeCoupleIterator]<<" "<<pacingNodeCouples1[nodeCoupleIterator]<<" adhesion not feasible"<<endl;
			continue;
		}
		//cout<<" checking pair for adhesion : "<<pacingNodeCouples0[nodeCoupleIterator]<<" "<<pacingNodeCouples1[nodeCoupleIterator]<<endl;
		int masterNodeId = pacingNodeCouples0[nodeCoupleIterator];
		int slaveNodeId = pacingNodeCouples1[nodeCoupleIterator];
		//Is there a closer pair for any of the selected master and slave?
		//I will not allow multi-pair adhesion if there is collapse, therefore I need to check the closest if possible
		if (collapseNodesOnAdhesion){
			for (int nodeCoupleIteratorFromHereOn = nodeCoupleIterator+1; nodeCoupleIteratorFromHereOn<n; ++nodeCoupleIteratorFromHereOn){
				if (masterNodeId == pacingNodeCouples0[nodeCoupleIteratorFromHereOn] ||
					masterNodeId == pacingNodeCouples1[nodeCoupleIteratorFromHereOn] ||
					slaveNodeId  == pacingNodeCouples0[nodeCoupleIteratorFromHereOn] ||
					slaveNodeId  == pacingNodeCouples1[nodeCoupleIteratorFromHereOn]
					){
					//one of the node couples occurs somewhere else, is this occurance viable?
					if (adhesionAllowedList[nodeCoupleIteratorFromHereOn]){
						//yes this adhesion is possible, compare distances:
						double dSqOriginal = distanceSqBetweenNodes(masterNodeId, slaveNodeId);
						double dSqNext = distanceSqBetweenNodes(pacingNodeCouples0[nodeCoupleIteratorFromHereOn],pacingNodeCouples1[nodeCoupleIteratorFromHereOn]);
						cout<<"two occurrences, distances are: "<<dSqOriginal<<" "<<dSqNext<<endl;
						if (dSqNext < dSqOriginal){
							//the next couple I will reach is closer, adhere them (or another one that may be further down the list and even closer)
							adhesionAllowedList[nodeCoupleIterator] = false;
							cout<<"I will adhere the next pair";
							break;
						}
						else{
							//this couple is further away than my original couple, adhere my original, this is not viable any more
							adhesionAllowedList[nodeCoupleIteratorFromHereOn] = false;
							cout<<"I will adhere this pair";
						}
					}
				}
			}
			//I may have decided to check the next adhesion couple first, this one is eliminated this time
			if (!adhesionAllowedList[nodeCoupleIterator]){
				continue;
			}
		}

		//check if collapse is feasible from elements perspective, i.e. will any element flip if I move the node(s)
		bool isAdhesionCollapseFeasible = true;
		if (collapseNodesOnAdhesion){
			vector<int> newCollapseList;
			double* avrPos = new double[3];
			bool* fix = new bool[3];
			Nodes[slaveNodeId]->getNewCollapseListAndAveragePos(newCollapseList, avrPos, fix,Nodes, masterNodeId);
			isAdhesionCollapseFeasible = checkForElementFlippingUponNodeCollapse(newCollapseList, avrPos);
			/*if (isAdhesionCollapseFeasible){
				//Adhering the nodes
				Nodes[slaveNodeId]->adheredTo = masterNodeId;
				Nodes[masterNodeId]->adheredTo = slaveNodeId;
				//collapse nodes if needed
				if (collapseNodesOnAdhesion){
					Nodes[slaveNodeId]->collapseOnNode(newCollapseList, avrPos, fix,Nodes, masterNodeId);
				}
			}
			delete[] avrPos;
			delete[] fix;
			if (!isAdhesionCollapseFeasible){
				cout<<"adhesion would cause element flipping, no adhesion"<<endl;
				adhesionAllowedList[nodeCoupleIterator] = false;
				continue;
			}*/
			Nodes[slaveNodeId]->adheredTo = masterNodeId;
			Nodes[masterNodeId]->adheredTo = slaveNodeId;
			if (collapseNodesOnAdhesion){
				//if (isAdhesionCollapseFeasible){
				//	Nodes[slaveNodeId]->collapseOnNode(newCollapseList, avrPos, fix,Nodes, masterNodeId);
				//}
				//else{
					Nodes[slaveNodeId]->collapseOnNodeInStages(newCollapseList, avrPos, fix,Nodes, masterNodeId);
				//}
			}
			delete[] avrPos;
			delete[] fix;
		}
		else{
			//Adhering the nodes
			Nodes[slaveNodeId]->adheredTo = masterNodeId;
			Nodes[masterNodeId]->adheredTo = slaveNodeId;
		}
		pacingNodeCouplesHaveAdhered[nodeCoupleIterator] = true;
		//cout<<"after main checks"<<endl;
		for (int i =0 ; i<3; ++i){
			//if the dimension is fixed in space, fix the other and move on.
			if (Nodes[masterNodeId]->FixedPos[i]){
				Nodes[slaveNodeId]->FixedPos[i]=true;
			}
			else if (Nodes[slaveNodeId]->FixedPos[i]){
				Nodes[masterNodeId]->FixedPos[i]=true;
			}
			//not using an else, as the slave could be fixed in given dimension independent of the master
			if (!Nodes[slaveNodeId]->FixedPos[i]){
				int dofmaster = masterNodeId*dim+i;
				int dofslave  = slaveNodeId*dim+i;
				//if slave is already slave of another node
				//   make the master of the slave the new slave
				//   algorithm will take care of the rest.
				if (Nodes[slaveNodeId]->slaveTo[i] > -1){
					//Current slave has a master. If this master is the current master, or the master of current master, than dont do anything, all is fine:
					if(Nodes[slaveNodeId]->slaveTo[i] == masterNodeId || Nodes[slaveNodeId]->slaveTo[i] == Nodes[masterNodeId]->slaveTo[i] ){
						pacingNodeCouplesHaveAdhered[nodeCoupleIterator] = true;
						//cout<<" eliminated, same numbers"<<endl;
						continue;
					}
					dofslave = Nodes[slaveNodeId]->slaveTo[i]*dim+i;
					slaveNodeId = Nodes[slaveNodeId]->slaveTo[i];
				}
				//check if the master dof is already bound to something:
				//cout<<"check master update"<<endl;
				NRSolver->checkMasterUpdate(dofmaster,masterNodeId);
				//It may have been that the slave was the master of the master node.
				//Now I have updated the master to the slave, and they are equal.
				//I do not need to add anything, as the master-slave relation is already implemented.
				//cout<<"DOF not fixed, after master update   : "<<dofmaster<<" "<<dofslave<<endl;
				if (dofmaster != dofslave){
					bool continueAddition =  NRSolver->checkIfCombinationExists(dofslave,dofmaster);
					//cout<<"DOF not fixed, continueAddition? : "<<continueAddition<<endl;
					if (continueAddition){
						bool madeChange = NRSolver->checkIfSlaveIsAlreadyMasterOfOthers(dofslave,dofmaster);
						if (madeChange){
							for (int nodeIt = 0 ; nodeIt<Nodes.size(); ++nodeIt){
								if(Nodes[nodeIt]->slaveTo[i]==slaveNodeId){
									Nodes[nodeIt]->slaveTo[i]=masterNodeId;
								}
							}
						}
						vector <int> fixDOF;
						fixDOF.push_back(dofslave);
						fixDOF.push_back(dofmaster);
						NRSolver->slaveMasterList.push_back(fixDOF);
						Nodes[slaveNodeId]->slaveTo[i] = masterNodeId;
						Nodes[masterNodeId]->isMaster[i] = true;
						//cout<<" adhereing nodes: "<<masterNodeId<<" "<<slaveNodeId<<" in dof "<<i<<endl;
						thereIsBinding = true;
						pacingNodeCouplesHaveAdhered[nodeCoupleIterator] = true;
						if (adherePeripodialToColumnar && Nodes[masterNodeId]->attachedToPeripodial){
							NRSolver->cleanPeripodialBindingFromMaster(dofmaster, Nodes);
							if (i==2){
								//freed z, not bound to peripodial anymore
								Nodes[masterNodeId]->attachedToPeripodial=false;
							}
						}
					}
				}
			}
		}
	}
	return thereIsBinding;
}


void Simulation::calculatePackingForcesImplicit3D(){
	//cout<<"inside calculatePackingForcesImplicit3D, size of packing node couples: "<<pacingNodeCouples0.size()<<endl;
	int n = pacingNodeCouples0.size();
	for(int i = 0 ; i<n; ++i){
		if (pacingNodeCouplesHaveAdhered[i]){
			continue;
		}
		int id0 = pacingNodeCouples0[i];
		int id1 = pacingNodeCouples1[i];
		double multiplier = packingMultiplier;

		//sigmoid test:
		double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];

		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		double sigmoidSaturation = sigmoidSaturationForPacking;

		if (initialWeightPointx[i]>0){
			dx *= -1.0;
		}
		if (initialWeightPointy[i]>0){
			dy *= -1.0;
		}
		if (initialWeightPointz[i]>0){
			dz *= -1.0;
		}

		double Fx = multiplier * averageMass / (1 + exp(sigmoidSaturation / packingThreshold * (-1.0*dx)));
		Fx *= initialWeightPointx[i];

		double Fy = multiplier * averageMass / (1 + exp(sigmoidSaturation / packingThreshold * (-1.0*dy)));
		Fy *= initialWeightPointy[i];

		double Fz = multiplier * averageMass / (1 + exp(sigmoidSaturation / packingThreshold * (-1.0*dz)));
		Fz *= initialWeightPointz[i];

		PackingForces[id0][0] += Fx;
		PackingForces[id0][1] += Fy;
		PackingForces[id0][2] += Fz;
		PackingForces[id1][0] -= Fx;
		PackingForces[id1][1] -= Fy;
		PackingForces[id1][2] -= Fz;
		if (isnan(Fx)){
		      cout<<" packing force Fx is nan for nodes "<<pacingNodeCouples0[i]<<" - "<<pacingNodeCouples1[i]<<endl;
		}
		if (isnan(Fy)){
		      cout<<" packing force Fy is nan for nodes "<<pacingNodeCouples0[i]<<" - "<<pacingNodeCouples1[i]<<endl;
		}
		if (isnan(Fz)){
		      cout<<" packing force Fz is nan for nodes "<<pacingNodeCouples0[i]<<" - "<<pacingNodeCouples1[i]<<endl;
		}
	}
}

void Simulation::calculatePackingJacobian3D(gsl_matrix* K){
	int n = pacingNodeCouples0.size();
	for(int i = 0 ; i<n; ++i){
		if (pacingNodeCouplesHaveAdhered[i]){
			continue;
		}
		int id0 = pacingNodeCouples0[i];
		int id1 = pacingNodeCouples1[i];
		double multiplier = packingMultiplier;

		//sigmoid test:
		double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];
		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		double sigmoidSaturation = sigmoidSaturationForPacking;
		if (initialWeightPointx[i]>0){
			dx *= -1.0;
		}
		if (initialWeightPointy[i]>0){
			dy *= -1.0;
		}
		if (initialWeightPointz[i]>0){
			dz *= -1.0;
		}
		double sigmoidx =  1 / (1 + exp(sigmoidSaturation/ packingThreshold * (-1.0 * dx) ));
		double dFxdx0 = sigmoidx * (1 - sigmoidx) * multiplier * averageMass * initialWeightPointx[i] * (sigmoidSaturation/packingThreshold);
		double dFxdx1 = -1.0*dFxdx0;
		if (initialWeightPointx[i]>0){
			dFxdx0 *= -1.0;
			dFxdx1 *= -1.0;
		}

		double sigmoidy =  1 / (1 + exp(sigmoidSaturation/ packingThreshold * (-1.0 * dy) ));
		double dFydy0 = sigmoidy * (1 - sigmoidy) * multiplier * averageMass * initialWeightPointy[i] * (sigmoidSaturation/packingThreshold);
		double dFydy1 = -1.0*dFydy0;
		if (initialWeightPointy[i]>0){
			dFydy0 *= -1.0;
			dFydy1 *= -1.0;
		}

		double sigmoidz =  1 / (1 + exp(sigmoidSaturation/ packingThreshold * (-1.0 * dz) ));
		double dFzdz0 = sigmoidz * (1 - sigmoidz) * multiplier * averageMass * initialWeightPointz[i] * (sigmoidSaturation/packingThreshold);
		double dFzdz1 = -1.0*dFzdz0;
		if (initialWeightPointz[i]>0){
			dFzdz0 *= -1.0;
			dFzdz1 *= -1.0;
		}
		//x values:
		double value = gsl_matrix_get(K,3*id0,3*id0);
		value -= dFxdx0;
		gsl_matrix_set(K,3*id0,3*id0,value);
		value = gsl_matrix_get(K,3*id0,3*id1);
		value -= dFxdx1;
		gsl_matrix_set(K,3*id0,3*id1,value);
		value = gsl_matrix_get(K,3*id1,3*id1);
		value -= dFxdx0;
		gsl_matrix_set(K,3*id1,3*id1,value);
		value = gsl_matrix_get(K,3*id1,3*id0);
		value -= dFxdx1;
		gsl_matrix_set(K,3*id1,3*id0,value);

		//y values:
		value = gsl_matrix_get(K,3*id0+1,3*id0+1);
		value -= dFydy0;
		gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		value = gsl_matrix_get(K,3*id0+1,3*id1+1);
		value -= dFydy1;
		gsl_matrix_set(K,3*id0+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id1+1,3*id1+1);
		value -= dFydy0;
		gsl_matrix_set(K,3*id1+1,3*id1+1,value);
		value = gsl_matrix_get(K,3*id1+1,3*id0+1);
		value -= dFydy1;
		gsl_matrix_set(K,3*id1+1,3*id0+1,value);

		//z values:
		value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value -= dFzdz0;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);
		value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		value -= dFzdz1;
		gsl_matrix_set(K,3*id0+2,3*id1+2,value);
		value = gsl_matrix_get(K,3*id1+2,3*id1+2);
		value -= dFzdz0;
		gsl_matrix_set(K,3*id1+2,3*id1+2,value);
		value = gsl_matrix_get(K,3*id1+2,3*id0+2);
		value -= dFzdz1;
		gsl_matrix_set(K,3*id1+2,3*id0+2,value);
	}
}

void Simulation::addValueToMatrix(gsl_matrix* K, int i, int j , double value){
	value += gsl_matrix_get(K,i,j);
	gsl_matrix_set(K,i,j,value);
}

void Simulation::calculatePackingNumerical(gsl_matrix* K){
	int n = pacingNodeCouples0.size();
	for(int i = 0 ; i<n; ++i){
		if (pacingNodeCouplesHaveAdhered[i]){
			continue;
		}
		int id0 = pacingNodeCouples0[i];
		int id1 = pacingNodeCouples1[i];
		double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		double multiplier = 1;
		double c = multiplier * averageMass;
		double p = 2; //the power of the division (d/t);
		double threshold = 6.0; //packing threshold;

		double perturbation = 0.0000001;
		double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];
		double Fmag = c * (1 - pow(dz/threshold,p));
		double F0[3] = {0, 0, Fmag};
		
		dz += perturbation;
		Fmag = c * (1 - pow(dz/threshold,p));
		double F3[3] = {0 , 0 , Fmag};

		double dgczz = (F3[2]-F0[2])/ perturbation;

		//cout<<"dgczz: "<<dgczz<<endl;//" dgcxy: "<<dgcxy<<" dgcyx: "<<dgcyx<<" c: "<<c<<endl;
		//cout<<"K matrix "<<3*id0+2<<" , "<<3*id0+2<<" : "<<gsl_matrix_get(K,3*id0+2,3*id0+2)<<endl;
		double value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value += dgczz;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);
		//cout<<"K matrix "<<3*id0+2<<" , "<<3*id0+2<<" : "<<gsl_matrix_get(K,3*id0+2,3*id0+2)<<endl;
		//cout<<"K matrix "<<3*id0+2<<" , "<<3*id1+2<<" : "<<gsl_matrix_get(K,3*id0+2,3*id1+2)<<endl;
		value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id0+2,3*id1+2,value);
		//cout<<"K matrix "<<3*id0+2<<" , "<<3*id1+2<<" : "<<gsl_matrix_get(K,3*id0+2,3*id1+2)<<endl;

		value = gsl_matrix_get(K,3*id1+2,3*id1+2);
		value += dgczz;
		gsl_matrix_set(K,3*id1+2,3*id1+2,value);

		value = gsl_matrix_get(K,3*id1+2,3*id0+2);
		value -= dgczz;
		gsl_matrix_set(K,3*id1+2,3*id0+2,value);

		/*
		double value = gsl_matrix_get(K,3*id0,3*id0);
		value +=dgcxx;
		gsl_matrix_set(K,3*id0,3*id0,value);
		//y_i y_i
		value = gsl_matrix_get(K,3*id0+1,3*id0+1);
		value +=dgcyy;
		gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		//z_i z_i
		value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		value +=dgczz;
		cout<<"K matrix 92,92: "<<gsl_matrix_get(K,92,92)<<endl;
		gsl_matrix_set(K,3*id0+2,3*id0+2,value);
		cout<<"K matrix 92,92: "<<gsl_matrix_get(K,92,92)<<endl;
		//x_i y_i  
		value = gsl_matrix_get(K,3*id0,3*id0+1);
		value +=dgcxy;
		gsl_matrix_set(K,3*id0,3*id0+1,value);
		//y_i x_i
		value = gsl_matrix_get(K,3*id0+1,3*id0);
		value +=dgcyx;
		gsl_matrix_set(K,3*id0+1,3*id0,value);
		//x_i z_i  
		value = gsl_matrix_get(K,3*id0,3*id0+2);
		value +=dgcxz;
		gsl_matrix_set(K,3*id0,3*id0+2,value);
		//z_i x_i
		value = gsl_matrix_get(K,3*id0+2,3*id0);
		value +=dgczx;
		gsl_matrix_set(K,3*id0+2,3*id0,value);
		//y_i z_i  
		value = gsl_matrix_get(K,3*id0+1,3*id0+2);
		value +=dgcyz;
		gsl_matrix_set(K,3*id0+1,3*id0+2,value);
		//z_i y_i
		value = gsl_matrix_get(K,3*id0+2,3*id0+1);
		value +=dgczy;
		gsl_matrix_set(K,3*id0+2,3*id0+1,value);

		cout<<"K matrix 92,92: before adding to slave "<<gsl_matrix_get(K,92,92)<<endl;
		//add the negative values to the slave:
		//x_i x_slave
		value = gsl_matrix_get(K,3*id0,3*id1);
		value -=dgcxx;
		gsl_matrix_set(K,3*id0,3*id1,value);
		//y_i y_slave
		value = gsl_matrix_get(K,3*id0+1,3*id1+1);
		value -=dgcyy;
		gsl_matrix_set(K,3*id0+1,3*id1+1,value);
		//z_i z_slave
		value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		value -=dgczz;
		gsl_matrix_set(K,3*id0+2,3*id1+2,value);
		//x_i y_slave  
		value = gsl_matrix_get(K,3*id0,3*id1+1);
		value -=dgcxy;
		gsl_matrix_set(K,3*id0,3*id1+1,value);
		//y_i x_slave
		value = gsl_matrix_get(K,3*id0+1,3*id1);
		value -=dgcyx;
		gsl_matrix_set(K,3*id0+1,3*id1,value);
		//x_i z_slave  
		value = gsl_matrix_get(K,3*id0,3*id1+2);
		value -=dgcxz;
		gsl_matrix_set(K,3*id0,3*id1+2,value);
		//z_i x_slave
		value = gsl_matrix_get(K,3*id0+2,3*id1);
		value -=dgczx;
		gsl_matrix_set(K,3*id0+2,3*id1,value);
		//y_i z_slave  
		value = gsl_matrix_get(K,3*id0+1,3*id1+2);
		value -=dgcyz;
		gsl_matrix_set(K,3*id0+1,3*id1+2,value);
		//z_i y_slave
		value = gsl_matrix_get(K,3*id0+2,3*id1+1);
		value -=dgczy;
		gsl_matrix_set(K,3*id0+2,3*id1+1,value);
*/
	}
}

void Simulation::calculatePackingK(gsl_matrix* K){
	int n = pacingNodeCouples0.size();
	for(int i = 0 ; i<n; ++i){
		 if (pacingNodeCouplesHaveAdhered[i]){
			continue;
		 }
		 int id0 = pacingNodeCouples0[i];
		 int id1 = pacingNodeCouples1[i];
		 double multiplier = 1.0;
		 double p = 0.3; //the power of the division (d/t);
		 double threshold = 10.0; //packing threshold;
		 double tp = pow(threshold,p);

		 double dx = Nodes[id0]->Position[0] - Nodes[id1]->Position[0];
		 double dy = Nodes[id0]->Position[1] - Nodes[id1]->Position[1];
		 double dz = Nodes[id0]->Position[2] - Nodes[id1]->Position[2];
		 double d = dx*dx + dy*dy + dz*dz;  //distance between the two nodes at current itertion
		 d = pow(d,0.5);
		 double averageMass = 0.5 *( Nodes[id0]->mass + Nodes[id1]->mass );
		 double c = multiplier * averageMass;
		 //cout<<"C: "<<c<<endl;
		 double dp = pow(d,p);

		 double ddx = dx / d;
		 double ddy = dy / d;
		 double ddz = dz / d;

		 //double dgcxx = ( -1.0 * c /tp *p *dp /d * ddx) * dx/d + (c * (1 - dp/tp)) * ( 1.0 / d - dx /d/d * ddx );
		 //double dgcxy = ( -1.0 * c /tp *p *dp /d * ddy) * dx/d - (c * (1 - dp/tp)) * dx /d/d * ddy;
		 //double dgcxz = ( -1.0 * c /tp *p *dp /d * ddz) * dx/d - (c * (1 - dp/tp)) * dx /d/d * ddz;
		 //double dgcyy = ( -1.0 * c /tp *p *dp /d * ddy) * dy/d + (c * (1 - dp/tp)) * (1.0 / d - dy /d/d * ddy);
		 //double dgcyz = ( -1.0 * c /tp *p *dp /d * ddz) * dy/d - (c * (1 - dp/tp)) * dy /d/d * ddz;
		 //double dgczz = ( -1.0 * c /tp *p *dp /d * ddz) * dz/d + (c * (1 - dp/tp)) * (1.0 / d - dz /d/d * ddz);

		 double dgcxx = (-1.0 * c/d/d * ddx - (p-1)/tp*dp/d/d*ddx) * dx + c/d*(1.0 + dp/tp);
		 double dgcxy = (-1.0 * c/d/d * ddy - (p-1)/tp*dp/d/d*ddy) * dx;
		 double dgcxz = (-1.0 * c/d/d * ddz - (p-1)/tp*dp/d/d*ddz) * dx;
		 double dgcyy = (-1.0 * c/d/d * ddy - (p-1)/tp*dp/d/d*ddy) * dy + c/d*(1.0 + dp/tp);
		 double dgcyz = (-1.0 * c/d/d * ddz - (p-1)/tp*dp/d/d*ddz) * dy;
		 double dgczz = (-1.0 * c/d/d * ddz - (p-1)/tp*dp/d/d*ddz) * dz + c/d*(1.0 + dp/tp);
		 // cout<<"dgcxx: "<<dgcxx<<" dgcyy: "<<dgcyy<<" dgcxy: "<<dgcxy<<endl;
		 //x_i x_i
		 double value = gsl_matrix_get(K,3*id0,3*id0);
		 value +=dgcxx;
		 gsl_matrix_set(K,3*id0,3*id0,value);
		 //y_i y_i
		 value = gsl_matrix_get(K,3*id0+1,3*id0+1);
		 value +=dgcyy;
		 gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		 //z_i z_i
		 value = gsl_matrix_get(K,3*id0+2,3*id0+2);
		 value +=dgczz;
		 gsl_matrix_set(K,3*id0+2,3*id0+2,value);
		 //x_i y_i  && y_i x_i
		 value = gsl_matrix_get(K,3*id0,3*id0+1);
		 value +=dgcxy;
		 gsl_matrix_set(K,3*id0,3*id0+1,value);
		 value = gsl_matrix_get(K,3*id0+1,3*id0);
		 value +=dgcxy;
		 gsl_matrix_set(K,3*id0+1,3*id0,value);
		 //x_i z_i  && z_i x_i
		 value = gsl_matrix_get(K,3*id0,3*id0+2);
		 value +=dgcxz;
		 gsl_matrix_set(K,3*id0,3*id0+2,value);
		 value = gsl_matrix_get(K,3*id0+2,3*id0);
		 value +=dgcxz;
		 gsl_matrix_set(K,3*id0+2,3*id0,value);
		 //y_i z_i  && z_i y_i
		 value = gsl_matrix_get(K,3*id0+1,3*id0+2);
		 value +=dgcyz;
		 gsl_matrix_set(K,3*id0+1,3*id0+2,value);
		 value = gsl_matrix_get(K,3*id0+2,3*id0+1);
		 value +=dgcyz;
		 gsl_matrix_set(K,3*id0+2,3*id0+1,value);

		 //add the negative values to the slave:
		 //x_i x_slave
		 value = gsl_matrix_get(K,3*id0,3*id1);
		 value -=dgcxx;
		 gsl_matrix_set(K,3*id0,3*id1,value);
		 //y_i y_slave
		 value = gsl_matrix_get(K,3*id0+1,3*id1+1);
		 value -=dgcyy;
		 gsl_matrix_set(K,3*id0+1,3*id0+1,value);
		 //z_i z_slave
		 value = gsl_matrix_get(K,3*id0+2,3*id1+2);
		 value -=dgczz;
		 gsl_matrix_set(K,3*id0+2,3*id0+2,value);
		 //x_i y_slave  && y_i x_slave
		 value = gsl_matrix_get(K,3*id0,3*id1+1);
		 value -=dgcxy;
		 gsl_matrix_set(K,3*id0,3*id1+1,value);
		 value = gsl_matrix_get(K,3*id0+1,3*id1);
		 value -=dgcxy;
		 gsl_matrix_set(K,3*id0+1,3*id1,value);
		 //x_i z_slave  && z_i x_slave
		 value = gsl_matrix_get(K,3*id0,3*id1+2);
		 value -=dgcxz;
		 gsl_matrix_set(K,3*id0,3*id1+2,value);
		 value = gsl_matrix_get(K,3*id0+2,3*id1);
		 value -=dgcxz;
		 gsl_matrix_set(K,3*id0+2,3*id1,value);
		 //y_i z_slave  && z_i y_slave
		 value = gsl_matrix_get(K,3*id0+1,3*id1+2);
		 value -=dgcyz;
		 gsl_matrix_set(K,3*id0+1,3*id1+2,value);
		 value = gsl_matrix_get(K,3*id0+2,3*id1+1);
		 value -=dgcyz;
		 gsl_matrix_set(K,3*id0+2,3*id1+1,value);
	}
}

void Simulation::processDisplayDataAndSave(){
    //if (displayIsOn && !DisplaySave){
        //The simulation is not displaying a saved setup, it is running and displaying
        //I need to correct the values to be displayed, and store averages, otherwise
        //the displayed values will be from artificial setups of different RK steps. (RK1 or RK4 depending on parameter)
        //updateDisplaySaveValuesFromRK();
    //}
	if (saveData && ( (int) (currSimTimeSec/dt) )% dataSaveInterval == 0){ //timestep % dataSaveInterval == 0){
		saveStep();
    }
}

void Simulation::updateNodeMasses(){
    for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
    	(*itNode)->mass = 0;
    }
    for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        if (!(*itElement)->IsAblated){
        	(*itElement)->assignVolumesToNodes(Nodes);
        }
    }
}

void Simulation::updateNodeViscositySurfaces(){
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for
    for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->IsAblated){
			(*itElement)->calculateViscositySurfaces();
		}
	}
    //calculated the areas, now assigning them
	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		(*itNode)->viscositySurface = 0;
	}
    for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	if (!(*itElement)->IsAblated){
    		(*itElement)->assignViscositySurfaceAreaToNodes(Nodes);
    	}
    }
}

void 	Simulation:: updateElementToConnectedNodes(vector <Node*>& Nodes){
	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
	    if ((*itNode)->mass > 0){//an ablated node will have this as zero
			int n = (*itNode)->connectedElementIds.size();
			for (int i=0; i<n; ++i){
				//if(!Elements[Nodes[j]->connectedElementIds[i]] -> IsAblated){
				(*itNode)->connectedElementWeights[i] = Elements[(*itNode)->connectedElementIds[i]]->VolumePerNode/(*itNode)->mass;
				//}
			}
	    }
    }
}

void Simulation::fillInNodeNeighbourhood(){
	vector<ShapeBase*>::iterator itEle;
	for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
		(*itEle)->fillNodeNeighbourhood(Nodes);
	}
}

void Simulation::setBasalNeighboursForApicalElements(){
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for //private(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces, outputFile, dt)
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->setBasalNeigElementId(Elements);
	}
}

void Simulation::fillInElementColumnLists(){
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for //private(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces, outputFile, dt)
	for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->IsAblated && (*itElement)->tissueType ==0 ){//Check only the columnar layer elements.
			if ((*itElement)->tissuePlacement == 0){
				//start from the basal element and move up
				//then you can copy the element numbers to the other elements on the list
				(*itElement)->constructElementStackList(TissueHeightDiscretisationLayers, Elements);
			}
		}
	}
}

void Simulation::getApicalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner){
	if (!ElementPointer->ApicalNormalForPackingUpToDate){
		ElementPointer->calculateNormalForPacking(1);
	}
	normalForPacking[0] = ElementPointer->ApicalNormalForPacking[0];
	normalForPacking[1] = ElementPointer->ApicalNormalForPacking[1];
	normalForPacking[2] = ElementPointer->ApicalNormalForPacking[2];
	ElementPointer->getApicalNodePos(posCorner);
}

void Simulation::getBasalNormalAndCornerPosForPacking(ShapeBase* ElementPointer, double* normalForPacking,double* posCorner){
	if (!ElementPointer->BasalNormalForPackingUpToDate){
		ElementPointer->calculateNormalForPacking(0);
	}
	normalForPacking[0] = ElementPointer->BasalNormalForPacking[0];
	normalForPacking[1] = ElementPointer->BasalNormalForPacking[1];
	normalForPacking[2] = ElementPointer->BasalNormalForPacking[2];
	ElementPointer->getBasalNodePos(posCorner);
}

void Simulation::getNormalAndCornerPosForPacking(Node* NodePointer, ShapeBase* ElementPointer, double* normalForPacking,double* posCorner){
	if (NodePointer->tissuePlacement == 1){
		//node is apical, should pack to apical nodes only
		getApicalNormalAndCornerPosForPacking(ElementPointer, normalForPacking, posCorner);
	}
	else if (NodePointer->tissuePlacement == 0){
		//node is basal, should pack to basal nodes only
		getBasalNormalAndCornerPosForPacking(ElementPointer, normalForPacking, posCorner);
	}
}

inline void Simulation::CapPackingForce(double& Fmag){
	double Fcap = 1e2;
	if (Fmag > Fcap){
		Fmag = Fcap;
	}
}

void Simulation::checkPackingToPipette(bool& packsToPip, double* pos, double* pipF, double mass){
	double pipThickness = 2; //microns
	double multiplier = 0.05; //just a term to scale the forces down
	double threshold = 1.0;	 //pipette packing forces start at 2 microns
	double t2 = threshold*threshold;
	double pipRange2[2] = {pipetteInnerRadius-threshold, pipetteInnerRadius+pipThickness+threshold};
	pipRange2[0] *= pipRange2[0];
	pipRange2[1] *= pipRange2[1];
	double dist[3] = {pos[0]-pipetteCentre[0], pos[1]-pipetteCentre[1], pos[2]-pipetteCentre[2]};
	double dist2[3] = {dist[0]*dist[0], dist[1]*dist[1],dist[2]*dist[2]};
	double xydist2 = dist2[0]+dist2[1];
	if (xydist2 > pipRange2[0] && xydist2<pipRange2[1]){
		//the node is within the packing ring range
		//is it close enough to bottom?
		if (dist2[2]<t2){
			//yes the node should be pushe (-ve)z for apical suciton and (+)ve z for basal suction:
			if (dist2[2]<1e-2){dist2[2] = 1e-2;}
			double Fmag = mass * multiplier * (1.0/dist2[2] - 1.0/t2);
			pipF[0] = 0.0;
			pipF[1] = 0.0;
			pipF[2] = Fmag;
			if (ApicalSuction){
				pipF[2] *= -1.0;
			}
			packsToPip = true;
			//cout<<"Node is pushed from bottom -  packing force: "<<pipF[0]<<" "<<pipF[1]<<" "<<pipF[2]<<endl;
			return;
		}
		//the node is not close to the pipette tip in x&y but not in range to be pushed by tip bottom, it may be too far away
		//now I need to check the distance in z and then the walls:
		if( (ApicalSuction && dist[2]>0) || (!ApicalSuction && dist[2]<0) ) {
			double midWallDist2 =  pipetteInnerRadius + 0.5*pipThickness;
			midWallDist2 *= midWallDist2;
			double d2 = 0.0;
			if (midWallDist2>xydist2){
				//the node should be pushed inside:
				multiplier *= -1.0;
				if(xydist2>pipetteInnerRadius*pipetteInnerRadius){
					//the node is inside the pipette: maximum force, d2 is zero:
					d2 =1e-2;
				}
				else{
					double dx = pos[0]-pipetteInnerRadius;
					double dy = pos[1]-pipetteInnerRadius;
					d2 = dx*dx+dy*dy;
				}
			}
			else{
				if(xydist2<(pipetteInnerRadius+pipThickness)*(pipetteInnerRadius+pipThickness)){
					//the node is inside the pipette: maximum force, d2 is zero:
					d2 =1e-2;
				}
				else{
					double dx = pos[0]-pipetteInnerRadius-pipThickness;
					double dy = pos[1]-pipetteInnerRadius-pipThickness;
					d2 = dx*dx+dy*dy;
				}
			}
			if (d2<1e-2){d2 = 1e-2;}
			double Fmag = mass * multiplier * (1.0/d2 - 1.0/t2);
			double xydist = pow(xydist2,0.5);
			Fmag /= xydist;
			pipF[0] = Fmag*dist[0];
			pipF[1] = Fmag*dist[1];
			pipF[2] = 0.0;
			//cout<<"Node "<<id<<" is pushed from side wall - packing force: "<<pipF[0]<<" "<<pipF[1]<<" "<<pipF[2]<<endl;
			packsToPip = true;
			return;
		}
	}
}
void Simulation::cleanUpPacingCombinations(){
	//each edge should be recorded once for a base node:
	//int m = pacingNodeEdgeList0.size();
	//for (int a = 0; a<m ; ++a){
	//	cout<<"edgelist: "<<pacingNodeEdgeList0[a]<<" "<<pacingNodeEdgeList1[a]<<" "<<pacingNodeEdgeList2[a]<<endl;
	//}
	//m = pacingNodePointList0.size();
	//for (int a = 0; a<m ; ++a){
	//	cout<<"nodelist: "<<pacingNodePointList0[a]<<" "<<pacingNodePointList1[a]<<endl;
	//}
	int n = pacingNodeEdgeList0.size();
	int i =0;
	while ( i<n){
		int currBase = pacingNodeEdgeList0[i];
		int edge1 = pacingNodeEdgeList1[i];
		int edge2 = pacingNodeEdgeList2[i];
		int j = i+1;
		while(j<n){
			if (pacingNodeEdgeList0[j] == currBase){
				if ( (edge1 == pacingNodeEdgeList1[j] && edge2 == pacingNodeEdgeList2[j]) || (edge1 == pacingNodeEdgeList2[j] && edge2 == pacingNodeEdgeList1[j])){
					//the edge is recorded again, remove these points:
					pacingNodeEdgeList0.erase(pacingNodeEdgeList0.begin()+j);
					pacingNodeEdgeList1.erase(pacingNodeEdgeList1.begin()+j);
					pacingNodeEdgeList2.erase(pacingNodeEdgeList2.begin()+j);
					initialSignsEdgex.erase(initialSignsEdgex.begin()+j);
					initialSignsEdgey.erase(initialSignsEdgey.begin()+j);
					initialSignsEdgez.erase(initialSignsEdgez.begin()+j);
					initialWeightEdgex.erase(initialWeightEdgex.begin()+j);
					initialWeightEdgey.erase(initialWeightEdgey.begin()+j);
					initialWeightEdgez.erase(initialWeightEdgez.begin()+j);
					j--;
					n--;
				}
			}
			j++;
		}
		i++;
	}
	//cleaning up nodes, each node should back to another node once:
	n = pacingNodePointList0.size();
	i =0;
	while (i<n){
		int node0 = pacingNodePointList0[i];
		int node1 = pacingNodePointList1[i];
		int j = i+1;
		while(j<n){
			if ((pacingNodePointList0[j] == node0 && pacingNodePointList1[j] == node1) || (pacingNodePointList1[j] == node0 && pacingNodePointList0[j] == node1)) {
				pacingNodePointList0.erase(pacingNodePointList0.begin()+j);
				pacingNodePointList1.erase(pacingNodePointList1.begin()+j);
				initialSignsPointx.erase(initialSignsPointx.begin()+j);
				initialSignsPointy.erase(initialSignsPointy.begin()+j);
				initialSignsPointz.erase(initialSignsPointz.begin()+j);
				initialWeightPointx.erase(initialWeightPointx.begin()+j);
				initialWeightPointy.erase(initialWeightPointy.begin()+j);
				initialWeightPointz.erase(initialWeightPointz.begin()+j);
				j--;
				n--;
			}
			j++;
		}
		i++;
	}
	//m = pacingNodeEdgeList0.size();
	//for (int a = 0; a<m ; ++a){
	//	cout<<"edgelist: "<<pacingNodeEdgeList0[a]<<" "<<pacingNodeEdgeList1[a]<<" "<<pacingNodeEdgeList2[a]<<endl;
	//}
	//m = pacingNodePointList0.size();
	//for (int a = 0; a<m ; ++a){
	//	cout<<"nodelist: "<<pacingNodePointList0[a]<<" "<<pacingNodePointList1[a]<<endl;
	//}
}

void Simulation::detectPacingCombinations(){
	double threshold = 7;	 //packing forces start at 4 microns - keep this lower than the force threshold!!
	double t2 = threshold*threshold;	//threshold square for rapid calculation
	//TO DO: make this the function  emptyPackingVectors();
	pacingNodeSurfaceList0.clear();
	pacingNodeSurfaceList1.clear();
	pacingNodeSurfaceList2.clear();
	pacingNodeSurfaceList3.clear();
	initialSignsSurfacex.clear();
	initialSignsSurfacey.clear();
	initialSignsSurfacez.clear();
	pacingNodeEdgeList0.clear();
	pacingNodeEdgeList1.clear();
	pacingNodeEdgeList2.clear();
	initialSignsEdgex.clear();
	initialSignsEdgey.clear();
	initialSignsEdgez.clear();
	pacingNodePointList0.clear();
	pacingNodePointList1.clear();
	initialSignsPointx.clear();
	initialSignsPointy.clear();
	initialSignsPointz.clear();
	initialWeightSurfacex.clear();
	initialWeightSurfacey.clear();
	initialWeightSurfacez.clear();
	initialWeightEdgex.clear();
	initialWeightEdgey.clear();
	initialWeightEdgez.clear();
	initialWeightPointx.clear();
	initialWeightPointy.clear();
	initialWeightPointz.clear();

	//
	vector<Node*>::iterator itNode;
	vector<ShapeBase*>::iterator itEle;
	// end of function
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		bool NodeHasPacking = (*itNode)->checkIfNodeHasPacking();
		if (NodeHasPacking){
			double* pos;
			pos = new double[3];
			(*itNode)->getCurrentPosition(pos);
			bool packedToSurface = false;
			bool packedToEdge = false;
			for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
				//excluding elements that own this node
				bool PackingToThisElement =  (*itEle)->checkPackingToThisNodeViaState(TissueHeightDiscretisationLayers, (*itNode));
				if (PackingToThisElement){
					int id1 = -1, id2 = -1, id3 = -1;
					(*itEle)->getRelevantNodesForPacking((*itNode)->tissuePlacement, id1, id2, id3);
					//make a separate function to detect packing to surface, filling the vecotrs if necessary and returning packedToSurface bool;
					//get mid point:
					if (!packedToSurface && !packedToEdge){
						//did not pack to surface or edge, packing to nodes:
						double dx = pos[0] - Nodes[id1]->Position[0];
						double dy = pos[1] - Nodes[id1]->Position[1];
						double dz = pos[2] - Nodes[id1]->Position[2];
						double d2 = dx*dx + dy*dy + dz*dz;
						if (d2<t2){
							pacingNodePointList0.push_back((*itNode)->Id);
							pacingNodePointList1.push_back(id1);
							if (dx >0) {initialSignsPointx.push_back(1);}
							else {initialSignsPointx.push_back(-1);}
							if (dy >0) {initialSignsPointy.push_back(1);}
							else {initialSignsPointy.push_back(-1);}
							if (dz >0) {initialSignsPointz.push_back(1);}
							else {initialSignsPointz.push_back(-1);}
							double d = pow(d2,0.5);
							initialWeightPointx.push_back(dx/d);
							initialWeightPointy.push_back(dy/d);
							initialWeightPointz.push_back(dz/d);
							//cout<<"Node: "<<(*itNode)->Id<<" "<<" slave edges: "<<id1<<" from element: "<<(*itEle)->Id<<endl;
						}
						dx = pos[0] - Nodes[id2]->Position[0];
						dy = pos[1] - Nodes[id2]->Position[1];
						dz = pos[2] - Nodes[id2]->Position[2];
						d2 = dx*dx + dy*dy + dz*dz;
						if (d2<t2){
							pacingNodePointList0.push_back((*itNode)->Id);
							pacingNodePointList1.push_back(id2);
							if (dx >0) {initialSignsPointx.push_back(1);}
							else {initialSignsPointx.push_back(-1);}
							if (dy >0) {initialSignsPointy.push_back(1);}
							else {initialSignsPointy.push_back(-1);}
							if (dz >0) {initialSignsPointz.push_back(1);}
							else {initialSignsPointz.push_back(-1);}
							double d = pow(d2,0.5);
							initialWeightPointx.push_back(dx/d);
							initialWeightPointy.push_back(dy/d);
							initialWeightPointz.push_back(dz/d);
							//cout<<"Node: "<<(*itNode)->Id<<" "<<" slave edges: "<<id2<<" from element: "<<(*itEle)->Id<<endl;
						}
						dx = pos[0] - Nodes[id3]->Position[0];
						dy = pos[1] - Nodes[id3]->Position[1];
						dz = pos[2] - Nodes[id3]->Position[2];
						d2 = dx*dx + dy*dy + dz*dz;
						if (d2<t2){
							pacingNodePointList0.push_back((*itNode)->Id);
							pacingNodePointList1.push_back(id3);
							if (dx >0) {initialSignsPointx.push_back(1);}
							else {initialSignsPointx.push_back(-1);}
							if (dy >0) {initialSignsPointy.push_back(1);}
							else {initialSignsPointy.push_back(-1);}
							if (dz >0) {initialSignsPointz.push_back(1);}
							else {initialSignsPointz.push_back(-1);}
							double d = pow(d2,0.5);
							initialWeightPointx.push_back(dx/d);
							initialWeightPointy.push_back(dy/d);
							initialWeightPointz.push_back(dz/d);
							//cout<<"Node: "<<(*itNode)->Id<<" "<<" slave edges: "<<id3<<" from element: "<<(*itEle)->Id<<endl;
						}
					}

				}
			}
			delete[] pos;
		}
	}
	cleanUpPacingCombinations();
}

void Simulation::setUpAFM(){
	packingToBeadThreshold = 5.0;
	beadR = 7.5;
	beadPos[0] = 0.0;
	beadPos[1] = 0.0;
	beadPos[2] = -7.75;
}

void Simulation::detectPackingToAFMBead(){
	double packingDetectionToBeadThreshold = 1.2 * packingToBeadThreshold;
	nodesPackingToBead.clear();
	initialWeightPackingToBeadx.clear();
	initialWeightPackingToBeady.clear();
	initialWeightPackingToBeadz.clear();
	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->mass >0){ //node is not ablated
			if ((*itNode)->tissuePlacement == 0){
				double* pos;
				pos = new double[3];
				(*itNode)->getCurrentPosition(pos);
				double dx = pos[0] - beadPos[0];
				double dy = pos[1] - beadPos[1];
				double dz = pos[2] - beadPos[2];
				double d2 = dx*dx + dy*dy +  dz*dz;
				double d = pow (d2,0.5);
				double dGap = d - beadR;
				dx = dx/d*dGap;
				dy = dy/d*dGap;
				dz = dz/d*dGap;
				if (dGap<packingDetectionToBeadThreshold){
					//node is close enough to pack to surface:
					nodesPackingToBead.push_back((*itNode)->Id);
					cout<<"packing to bead: "<<(*itNode)->Id<<" dGap: "<<dGap<<" threshold: "<< packingDetectionToBeadThreshold<<" d: "<<dx<<" "<<dy<<" "<<dz<<endl;
					initialWeightPackingToBeadx.push_back(dx/dGap);
					initialWeightPackingToBeady.push_back(dy/dGap);
					initialWeightPackingToBeadz.push_back(dz/dGap);
					distanceToBead.push_back(dGap);
				}
				delete[] pos;
			}
		}
	}
}

void Simulation::detectPacingToEnclosingSurfacesNodes(){
	packingToEnclosingSurfacesThreshold = 3;  //pack to the boundary at 3 microns distance
	packingDetectionToEnclosingSurfacesThreshold = 1.2 * packingToEnclosingSurfacesThreshold;
	double t2 = packingDetectionToEnclosingSurfacesThreshold*packingDetectionToEnclosingSurfacesThreshold;	//threshold square for rapid calculation
	nodesPackingToPositiveSurface.clear();
	nodesPackingToNegativeSurface.clear();
	initialWeightPackingToPositiveSurface.clear();
	initialWeightPackingToNegativeSurface.clear();
	//calculate the current boundaries:
	if (currSimTimeSec < initialTimeToEncloseTissueBetweenSurfacesSec){
		return;
	}
	else if (currSimTimeSec>=finalTimeToEncloseTissueBetweenSurfacesSec){
		zEnclosementBoundaries[0] = finalZEnclosementBoundaries[0];
		zEnclosementBoundaries[1] = finalZEnclosementBoundaries[1];
	}
	else{
		double totalTime = finalTimeToEncloseTissueBetweenSurfacesSec - initialTimeToEncloseTissueBetweenSurfacesSec;
		double currTimeDiff = currSimTimeSec - initialTimeToEncloseTissueBetweenSurfacesSec;
		double zNegDifference = finalZEnclosementBoundaries[0] - initialZEnclosementBoundaries[0];
		double zPosDifference = finalZEnclosementBoundaries[1] - initialZEnclosementBoundaries[1];
		zEnclosementBoundaries[0] = initialZEnclosementBoundaries[0] + currTimeDiff*zNegDifference / totalTime;
		zEnclosementBoundaries[1] = initialZEnclosementBoundaries[1] + currTimeDiff*zPosDifference / totalTime;

	}
	cout<<"curr time in sec: "<<currSimTimeSec<<" z enclosement boundaries: "<<zEnclosementBoundaries[0]<<" "<<zEnclosementBoundaries[1]<<" initial boundaries: "<<initialZEnclosementBoundaries[0]<<" "<<initialZEnclosementBoundaries[1]<<" final boundaries: "<<finalZEnclosementBoundaries[0]<<" "<<finalZEnclosementBoundaries[1]<<endl;

	//go over nodes to detect packing:
	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->mass >0){ //node is not ablated
			bool checkAgainstPositiveSurface = false;
			bool checkAgainstNegativeSurface = false;
			if ((*itNode)->atCircumference){
				//Node is at the ciscumference, with sufficient rotation, it can pack anywhere
				//chack against both surfaces.
				checkAgainstPositiveSurface = true;
				checkAgainstNegativeSurface = true;
			}
			else {
				//basal node of the columnar layer is checked against the negative surface:
				if ( (*itNode)->tissuePlacement == 0){
					//node is basal
					//if it is columnar layer, check against negative surface:
					if( (*itNode)->tissueType ==0 ){
						checkAgainstNegativeSurface = true;
					}
					else if((*itNode)->tissueType ==1 ){
						//if it is peripodial layer, check against positive surface:
						checkAgainstPositiveSurface = true;
					}
				}
				//if tissue placement is apical, and the tissue tyoe is columnar,
				//check against the positive surface:
				if( (*itNode)->tissuePlacement == 1 && (*itNode)->tissueType ==0 ){	//Node is apical, check against positive border
					checkAgainstPositiveSurface = true;
				}
			}
			if( checkAgainstNegativeSurface ){
				double* pos;
				pos = new double[3];
				(*itNode)->getCurrentPosition(pos);
				double dz = pos[2] - zEnclosementBoundaries[0];
				double d2 = dz*dz;
				if (d2<t2){
					//node is close enough to pack to surface:
					nodesPackingToNegativeSurface.push_back((*itNode)->Id);
					double d = pow (d2,0.5);
					initialWeightPackingToNegativeSurface.push_back(dz/d);
				}
				delete[] pos;
			}
			if( checkAgainstPositiveSurface){
				double* pos;
				pos = new double[3];
				(*itNode)->getCurrentPosition(pos);
				double dz = pos[2] - zEnclosementBoundaries[1];
				double d2 = dz*dz;
				if (d2<t2){
					//node is close enough to pack to surface:
					nodesPackingToPositiveSurface.push_back((*itNode)->Id);
					double d = pow (d2,0.5);
					initialWeightPackingToPositiveSurface.push_back(dz/d);
				}
				delete[] pos;
			}
		}
	}
}


void 	Simulation::assignFoldRegionAndReleasePeripodial(Node* NodeMaster, Node* NodeSlave ){
	int nDim  = 3;
	int masterXGridIndex = floor( (NodeMaster->Position[0] - boundingBox[0][0])/boundingBoxSize[0] );
	int masterYGridIndex = floor( (NodeMaster->Position[1] - boundingBox[0][1])/boundingBoxSize[1] );

	NodeMaster->onFoldInitiation = true;
	NodeSlave->onFoldInitiation = true;
	//cout<<"assigning to fold initiation, nodes "<<(*itNode)->Id <<" and "<<(*itNodeSlave)->Id<<endl;
	//check for other nodes in the vicinity:
	//if I am on the apical side, I will check the borders from basal sides
	//If I am on the basal side, I will chack from apical sides:
	int masterCorrespondingX = NodeMaster->Position[0];
	int slaveCorrespoindingX = NodeSlave->Position[0];
	int masterCorrespondingY = NodeMaster->Position[1];
	int slaveCorrespoindingY = NodeSlave->Position[1];

	for (vector<int>::iterator itInt = NodeMaster->immediateNeigs.begin(); itInt < NodeMaster->immediateNeigs.end(); ++itInt){
		if(Nodes[(*itInt)]->tissuePlacement != NodeMaster->tissuePlacement){
			int commonOwnerId = NodeMaster->getCommonOwnerId(Nodes[(*itInt)]);
			//cout<<"master Node : "<<NodeMaster->Id<<" checked node : "<<Nodes[(*itInt)]->Id<<" common owner: "<<commonOwnerId<<endl;
			if (commonOwnerId > -1){
				//there is common owner
				//are nodes on top of each other on this element?
				bool nodesConnected = Elements[commonOwnerId]->areNodesDirectlyConnected(NodeMaster->Id, Nodes[(*itInt)]->Id);
				if (nodesConnected){
					masterCorrespondingX = Nodes[(*itInt)]->Position[0];
					masterCorrespondingY = Nodes[(*itInt)]->Position[1];
					//cout<<"master Node : "<<NodeMaster->Id<<" corresponding; "<<Nodes[(*itInt)]->Id<<endl;
					break;
				}
			}
		}
	}
	for (vector<int>::iterator itInt = NodeSlave->immediateNeigs.begin(); itInt < NodeSlave->immediateNeigs.end(); ++itInt){
		if(Nodes[(*itInt)]->tissuePlacement != NodeSlave->tissuePlacement){
			int commonOwnerId = NodeSlave->getCommonOwnerId(Nodes[(*itInt)]);
			if (commonOwnerId > -1){
				//there is common owner
				//are nodes on top of each other on this element?
				bool nodesConnected = Elements[commonOwnerId]->areNodesDirectlyConnected(NodeSlave->Id, Nodes[(*itInt)]->Id);
				if (nodesConnected){
					slaveCorrespoindingX = Nodes[(*itInt)]->Position[0];
					slaveCorrespoindingY = Nodes[(*itInt)]->Position[1];
					//cout<<"slave Node : "<<NodeSlave->Id<<" corresponding; "<<Nodes[(*itInt)]->Id<<endl;
					break;
				}
			}
		}
	}
	double xmin = min(masterCorrespondingX, slaveCorrespoindingX);
	double xmax = max(masterCorrespondingX, slaveCorrespoindingX);
	double ymin = min(masterCorrespondingY, slaveCorrespoindingY);
	double ymax = max(masterCorrespondingY, slaveCorrespoindingY);
	//double xmin = min((*itNode)->Position[0], (*itNodeSlave)->Position[0]);
	//double xmax = max((*itNode)->Position[0], (*itNodeSlave)->Position[0]);
	//double ymin = min((*itNode)->Position[1], (*itNodeSlave)->Position[1]);
	//double ymax = max((*itNode)->Position[1], (*itNodeSlave)->Position[1]);
	double zmax = max(NodeMaster->Position[2], NodeSlave->Position[2]);
	double zmin = min(NodeMaster->Position[2], NodeSlave->Position[2]);
	ymin -= 4.0*packingDetectionThresholdGrid[masterXGridIndex][masterYGridIndex];
	ymax += 4.0*packingDetectionThresholdGrid[masterXGridIndex][masterYGridIndex];
	for (vector<Node*>::iterator itNodeInBetween = Nodes.begin(); itNodeInBetween < Nodes.end(); itNodeInBetween++){
		if((*itNodeInBetween)->tissuePlacement == NodeMaster->tissuePlacement ){
		if((*itNodeInBetween)->attachedToPeripodial || !(*itNodeInBetween)->onFoldInitiation ){
			//for apical curves, chack for zmax, for basal , check for z min
			if( ((*itNodeInBetween)->tissuePlacement == 1 && (*itNodeInBetween)->Position[2] <= zmax ) ||((*itNodeInBetween)->tissuePlacement == 0 && (*itNodeInBetween)->Position[2] >= zmin )   ){
				if((*itNodeInBetween)->Position[0] >= xmin && (*itNodeInBetween)->Position[0] <= xmax){
					if((*itNodeInBetween)->Position[1] >= ymin && (*itNodeInBetween)->Position[1] <= ymax){
						(*itNodeInBetween)->onFoldInitiation = true;
						//cout<<"adding node in between "<<(*itNodeInBetween)->Id<<endl;
						if((*itNodeInBetween)->attachedToPeripodial){
							cout<<" Releasing peripodial from node "<<(*itNodeInBetween)->Id<<" as in between"<<endl;
							for (int j=0;j<nDim;++j){
								double dof =0;
								dof = (*itNodeInBetween)->Id * nDim+j;
								NRSolver->cleanPeripodialBindingFromMaster(dof, Nodes);
							}
							(*itNodeInBetween)->attachedToPeripodial = false;
						}
					}
				}
			}
		}}
	}

	for (int j=0;j<nDim;++j){
		double dof =0;
		if (NodeMaster->attachedToPeripodial){
			cout<<"Releasing peripodial from node "<< NodeMaster->Id<<" via "<<NodeSlave->Id <<endl;
			dof = NodeMaster->Id * nDim+j;
			NRSolver->cleanPeripodialBindingFromMaster(dof, Nodes);
		}
		if (NodeSlave->attachedToPeripodial){
			cout<<"Releasing peripodial from node "<< NodeSlave->Id<<" via "<<NodeMaster->Id <<endl;
			dof = NodeSlave->Id * nDim+j;
			NRSolver->cleanPeripodialBindingFromMaster(dof, Nodes);
		}
	}
	NodeMaster->attachedToPeripodial = false;
	NodeSlave->attachedToPeripodial = false;
}



void Simulation::detectPacingNodes(){
	cout<<"in detect packing"<<endl;
	double periAverageSideLength = 0,colAverageSideLength = 0;
	getAverageSideLength(periAverageSideLength, colAverageSideLength);

	if (thereIsPeripodialMembrane){
		colAverageSideLength = (periAverageSideLength+colAverageSideLength)/2.0;
	}
	double packingDetectionThresholdGridSq[10][5];
	for (int i=0;i<10;++i){
		for(int j=0;j<5;++j){
			// 0.4*1.2 normal packing detection value
			// 0.25 is stable with clashes
			// 0.333
			packingDetectionThresholdGrid[i][j] = 0.4 * packingDetectionThresholdGrid[i][j];
			packingDetectionThresholdGridSq[i][j] = packingDetectionThresholdGrid[i][j]*packingDetectionThresholdGrid[i][j];
		}
	}
	packingThreshold = 0.4*colAverageSideLength;  //0.6 was old value
	packingDetectionThreshold = 1.2 * packingThreshold; //0.9 * packingThreshold;
	packingMultiplier = 5000; //1000 in original, increasing 5 fold to test run54 on 16 april 2018
	sigmoidSaturationForPacking = 5;
	double t2 = packingDetectionThreshold*packingDetectionThreshold;	//threshold square for rapid calculation
	pacingNodeCouples0.clear();
	pacingNodeCouples1.clear();
	pacingNodeCouplesHaveAdhered.clear();
	initialWeightPointx.clear();
	initialWeightPointy.clear();
	initialWeightPointz.clear();
	#pragma omp parallel for
	for (vector<ShapeBase*>::iterator itEle = Elements.begin(); itEle<Elements.end(); ++itEle){
		if ((*itEle)->tissueType == 0 && ((*itEle)->tissuePlacement == 1 || (*itEle)->spansWholeTissue)){
			//columnar element at the apical surface or spans whole tissue
			(*itEle)->calculateApicalNormalCurrentShape();
		}
	}
	//Added for parallelisation:
	//node based only:
	const int nArray = 16; //the size of array that I will divide the elements into.
	int parallellisationSegmentSize = ceil(nNodes/nArray);
	if (parallellisationSegmentSize<1){
		parallellisationSegmentSize =1;
	}
	const int nSegments = ceil(nNodes/parallellisationSegmentSize); //this does not need to be the same as desired size
	//if there are 5 nodes and I want 16 segments, the actual nSegments will be 5.
	const int segmentBoundariesArraySize = nSegments+1;
	int parallellisationSegmentBoundaries[segmentBoundariesArraySize];
	parallellisationSegmentBoundaries[0] = 0;
	parallellisationSegmentBoundaries[segmentBoundariesArraySize-1] = nNodes;
	for (int i=1; i<segmentBoundariesArraySize-1; ++i){
		parallellisationSegmentBoundaries[i] = i*parallellisationSegmentSize;
	}
	vector <vector<int> > arrayForParallelisationPacingNodeCouples0(nSegments, vector <int>(0));
	vector <vector<int> > arrayForParallelisationPacingNodeCouples1(nSegments, vector <int>(0));
	vector <vector<double> > arrayForParallelisationInitialWeightPointx(nSegments, vector <double>(0));
	vector <vector<double> > arrayForParallelisationInitialWeightPointy(nSegments, vector <double>(0));
	vector <vector<double> > arrayForParallelisationInitialWeightPointz(nSegments, vector <double>(0));
	//parallelise this loop:
	//cout<<" at parallelisation loop for packing, within detectPacingNodes"<<endl;
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	/*cout<<"parallellisationSegmentBoundaries: "<<endl;
	for (int i = 0; i<segmentBoundariesArraySize; ++i){
		cout<<" "<<parallellisationSegmentBoundaries[i]<<" ";
	}
	cout<<endl;*/
	//#pragma omp parallel for //private(arrayForParallelisationPacingNodeCouples0, arrayForParallelisationPacingNodeCouples1, arrayForParallelisationInitialWeightPointx, arrayForParallelisationInitialWeightPointy, arrayForParallelisationInitialWeightPointz)
	for (int a =0; a<nSegments; ++a){
		int initialpoint = parallellisationSegmentBoundaries[a];
		int breakpoint =  parallellisationSegmentBoundaries[a+1];
		//cout<<" a = "<<a<<" initialpoint: "<<initialpoint<<" breakpoint "<<breakpoint<<" nNodes: "<<nNodes<<endl;
		for (vector<Node*>::iterator itNode=Nodes.begin()+ initialpoint; itNode<Nodes.begin()+breakpoint; ++itNode){
			bool NodeHasPacking = (*itNode)->checkIfNodeHasPacking();
			if (NodeHasPacking){//702 - 596
				double* pos;
				pos = new double[3];
				(*itNode)->getCurrentPosition(pos);
				int masterXGridIndex = floor( (pos[0] - boundingBox[0][0])/boundingBoxSize[0] );
				int masterYGridIndex = floor( (pos[1] - boundingBox[0][1])/boundingBoxSize[1] );
				t2 = packingDetectionThresholdGridSq[masterXGridIndex][masterYGridIndex];
				for (vector<Node*>::iterator itNodeSlave=itNode+1; itNodeSlave<Nodes.end(); ++itNodeSlave){
					bool SlaveNodeHasPacking = (*itNodeSlave)->checkIfNodeHasPacking();
					//bool nodesOnSeperateSurfaces = (*itNodeSlave)->tissueType != (*itNode)->tissueType;
					//if (SlaveNodeHasPacking && nodesOnSeperateSurfaces){
					if (SlaveNodeHasPacking){
						if ((*itNode)->tissuePlacement == (*itNodeSlave)->tissuePlacement){
							//nodes can pack, are they connected?
							bool neigbours = (*itNode)->isMyNeig((*itNodeSlave)->Id);
							if (neigbours){
								continue;
							}
							//check if the master node is in the collapsed list?
							//the vector is sorted, as I sort it to remove duplicates
							if (binary_search((*itNode)->collapsedWith.begin(), (*itNode)->collapsedWith.end(),(*itNodeSlave)->Id)){
								//the couple is already collapsed:
								continue;
							}
							//check if the collapsed nodes of master are neigs of slave
							neigbours = (*itNode)->isNeigWithMyCollapsedNodes((*itNodeSlave)->Id,Nodes);
							if (neigbours){
								continue;
							}
							neigbours = (*itNodeSlave)->isNeigWithMyCollapsedNodes((*itNode)->Id,Nodes);
							if (neigbours){
								continue;
							}
							//the nodes can potentially pack, are they close enough?
							double* posSlave;
							posSlave = new double[3];
							(*itNodeSlave)->getCurrentPosition(posSlave);
							double dx = pos[0] - posSlave[0];
							double dy = pos[1] - posSlave[1];
							double dz = pos[2] - posSlave[2];
							double d2 = dx*dx + dy*dy + dz*dz;
							if (d2<t2){
								//close enough for packing , add to list:
								arrayForParallelisationPacingNodeCouples0[a].push_back((*itNode)->Id);
								arrayForParallelisationPacingNodeCouples1[a].push_back((*itNodeSlave)->Id);
								double d = pow (d2,0.5);
								arrayForParallelisationInitialWeightPointx[a].push_back(dx/d);
								arrayForParallelisationInitialWeightPointy[a].push_back(dy/d);
								arrayForParallelisationInitialWeightPointz[a].push_back(dz/d);
							}
							bool checkForCurvedRegions = false;

							double curveIdentificationThresholdSq = packingDetectionThresholdGrid[masterXGridIndex][masterYGridIndex];
							curveIdentificationThresholdSq *= curveIdentificationThresholdSq;
							curveIdentificationThresholdSq *=9;//0.625;
							//if ((*itNode)->Id == 5207 && (*itNodeSlave)->Id == 5293 ){
							//	cout<<" checking couple :"<<(*itNode)->Id <<" "<<(*itNodeSlave)->Id<<" detech threshold: "<< packingDetectionThresholdGrid[masterXGridIndex][masterYGridIndex]<<" d2 "<<d2<<" thres sq: "<<curveIdentificationThresholdSq<<endl;
							//}
							if (d2<curveIdentificationThresholdSq){
								if (thereIsEmergentEllipseMarking && (!(*itNode)->onFoldInitiation || !(*itNodeSlave)->onFoldInitiation) ){
									//only check if the nodes are not covered by emergent ellipse bands
									checkForCurvedRegions=true;
								}
								else if (adherePeripodialToColumnar && ((*itNode)->attachedToPeripodial || (*itNodeSlave)->attachedToPeripodial)){
									//only check if one of the nodes was attached to peripodial
									checkForCurvedRegions=true;
								}
							}
							if (checkForCurvedRegions){
								//cout<<" check for curved regions active"<<endl;
								int elementIdMaster = (*itNode)->connectedElementIds[0];
								int elementIdSlave = (*itNodeSlave)->connectedElementIds[0];

								double dotP = Elements[elementIdSlave]->dotProduct3D(Elements[elementIdMaster]->apicalNormalCurrentShape,Elements[elementIdSlave]->apicalNormalCurrentShape);
								//cout<<"dotp: "<<dotP<<endl;
								if (dotP <0){
									//normals point opposite directions
									//do they face each other?
									bool releaseNodes = false;
									double *connectingVec = new double[3];
									for (int dimIter=0; dimIter<3; ++dimIter){
										connectingVec[dimIter] = (*itNode)->Position[dimIter] - (*itNodeSlave)->Position[dimIter];
									}
									double dotPmaster = Elements[elementIdMaster]->dotProduct3D(connectingVec,Elements[elementIdMaster]->apicalNormalCurrentShape);
									if (dotPmaster >0 ){
										releaseNodes = true;
									}
									else{
										double dotPslave = Elements[elementIdSlave]->dotProduct3D(connectingVec,Elements[elementIdSlave]->apicalNormalCurrentShape);
										if (dotPslave>0){
											releaseNodes = true;
										}
									}
									if (!releaseNodes){
										continue;
									}
									assignFoldRegionAndReleasePeripodial((*itNode),(*itNodeSlave));
/*									(*itNode)->onFoldInitiation = true;
									(*itNodeSlave)->onFoldInitiation = true;
									//cout<<"assigning to fold initiation, nodes "<<(*itNode)->Id <<" and "<<(*itNodeSlave)->Id<<endl;
									//check for other nodes in the vicinity:
									//if I am on the apical side, I will check the borders from basal sides
									//If I am on the basal side, I will chack from apical sides:
									int masterCorrespondingX = (*itNode)->Position[0];
									int slaveCorrespoindingX = (*itNodeSlave)->Position[0];
									int masterCorrespondingY = (*itNode)->Position[1];
									int slaveCorrespoindingY = (*itNodeSlave)->Position[1];

									for (vector<int>::iterator itInt = (*itNode)->immediateNeigs.begin(); itInt < (*itNode)->immediateNeigs.end(); ++itInt){
										if(Nodes[(*itInt)]->tissuePlacement != (*itNode)->tissuePlacement){
											masterCorrespondingX = Nodes[(*itInt)]->Position[0];
											masterCorrespondingY = Nodes[(*itInt)]->Position[1];
											cout<<"Node : "<<(*itNode)->Id<<" corresponding; "<<(*itInt)<<endl;
											break;
										}
									}

									for (vector<int>::iterator itInt = (*itNodeSlave)->immediateNeigs.begin(); itInt < (*itNodeSlave)->immediateNeigs.end(); ++itInt){
										if(Nodes[(*itInt)]->tissuePlacement != (*itNodeSlave)->tissuePlacement){
											slaveCorrespoindingX = Nodes[(*itInt)]->Position[0];
											slaveCorrespoindingY = Nodes[(*itInt)]->Position[1];
											cout<<"Node : "<<(*itNode)->Id<<" corresponding; "<<(*itInt)<<endl;
											break;
										}
									}
									double xmin = min(masterCorrespondingX, slaveCorrespoindingX);
									double xmax = max(masterCorrespondingX, slaveCorrespoindingX);
									double ymin = min(masterCorrespondingY, slaveCorrespoindingY);
									double ymax = max(masterCorrespondingY, slaveCorrespoindingY);
									//double xmin = min((*itNode)->Position[0], (*itNodeSlave)->Position[0]);
									//double xmax = max((*itNode)->Position[0], (*itNodeSlave)->Position[0]);
									//double ymin = min((*itNode)->Position[1], (*itNodeSlave)->Position[1]);
									//double ymax = max((*itNode)->Position[1], (*itNodeSlave)->Position[1]);
									double zmax = max((*itNode)->Position[2], (*itNodeSlave)->Position[2]);
									double zmin = min((*itNode)->Position[2], (*itNodeSlave)->Position[2]);
									ymin -= 4.0*packingDetectionThresholdGrid[masterXGridIndex][masterYGridIndex];
									ymax += 4.0*packingDetectionThresholdGrid[masterXGridIndex][masterYGridIndex];
									for (vector<Node*>::iterator itNodeInBetween = Nodes.begin(); itNodeInBetween < Nodes.end(); itNodeInBetween++){
										if((*itNodeInBetween)->tissuePlacement == (*itNode)->tissuePlacement ){
										if((*itNodeInBetween)->attachedToPeripodial || !(*itNodeInBetween)->onFoldInitiation ){
											//for apical curves, chack for zmax, for basal , check for z min
											if( ((*itNodeInBetween)->tissuePlacement == 1 && (*itNodeInBetween)->Position[2] < zmax ) ||((*itNodeInBetween)->tissuePlacement == 0 && (*itNodeInBetween)->Position[2] > zmin )   ){
												if((*itNodeInBetween)->Position[0] > xmin && (*itNodeInBetween)->Position[0] < xmax){
													if((*itNodeInBetween)->Position[1] > ymin && (*itNodeInBetween)->Position[1] < ymax){
														(*itNodeInBetween)->onFoldInitiation = true;
														//cout<<"adding node in between "<<(*itNodeInBetween)->Id<<endl;
														if((*itNodeInBetween)->attachedToPeripodial){
															cout<<" Releasing peripodial from node "<<(*itNodeInBetween)->Id<<" as in between"<<endl;
															for (int j=0;j<nDim;++j){
																double dof =0;
																dof = (*itNodeInBetween)->Id * nDim+j;
																NRSolver->cleanPeripodialBindingFromMaster(dof, Nodes);
															}
															(*itNodeInBetween)->attachedToPeripodial = false;
														}
													}
												}
											}
										}}
									}

									for (int j=0;j<nDim;++j){
										double dof =0;
										if ((*itNode)->attachedToPeripodial){
											cout<<"Releasing peripodial from node "<< (*itNode)->Id<<" via "<<(*itNodeSlave)->Id <<endl;
											dof = (*itNode)->Id * nDim+j;
											NRSolver->cleanPeripodialBindingFromMaster(dof, Nodes);
										}
										if ((*itNodeSlave)->attachedToPeripodial){
											cout<<"Releasing peripodial from node "<< (*itNodeSlave)->Id<<" via "<<(*itNode)->Id <<endl;
											dof = (*itNodeSlave)->Id * nDim+j;
											NRSolver->cleanPeripodialBindingFromMaster(dof, Nodes);
										}
									}
									(*itNode)->attachedToPeripodial = false;
									(*itNodeSlave)->attachedToPeripodial = false;*/
								}
							}
							delete[] posSlave;
						}
					}
				}
				delete[] pos;
			}
		}
		//calculate packing here
	}
	//outside parallellisation, need to combine the nodes:
	for (int a =0; a<nSegments-1; ++a){
		int n = arrayForParallelisationPacingNodeCouples0[a].size();
		for (int i=0; i<n; ++i){
			pacingNodeCouples0.push_back(arrayForParallelisationPacingNodeCouples0[a][i]);
			pacingNodeCouples1.push_back(arrayForParallelisationPacingNodeCouples1[a][i]);
			pacingNodeCouplesHaveAdhered.push_back(false);
			initialWeightPointx.push_back(arrayForParallelisationInitialWeightPointx[a][i]);
			initialWeightPointy.push_back(arrayForParallelisationInitialWeightPointy[a][i]);
			initialWeightPointz.push_back(arrayForParallelisationInitialWeightPointz[a][i]);
			//cout<<"packing nodes: "<<pacingNodeCouples0.back()<<" "<<pacingNodeCouples1.back()<<endl;
		}
	}
}

void Simulation::calculatePacking(){
	//cout<<"in calculate packing: "<<endl;
	double multiplier = 0.3; //just a term to scale the forces down
	double thresholdForNodeDetection = 20; //check every node within 5 microns of distance;
	double threshold = 10;	 //packing forces start at 2 microns
	double thresholdNeg = (-1.0)*threshold;
	double t2 = threshold*threshold;	//threshold square for rapid calculation
	vector<Node*>::iterator itNode;
	vector<ShapeBase*>::iterator itEle;
	//resetting all the normal update flags
	for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
		(*itEle)->ApicalNormalForPackingUpToDate = false;
		(*itEle)->BasalNormalForPackingUpToDate = false;
	}
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		/*double sumPack[3] = {0,0,0};
		for (int pp=0; pp<nNodes; ++pp){
			sumPack[0] += PackingForces[pp][0];
			sumPack[1] += PackingForces[pp][1];
			sumPack[2] += PackingForces[pp][2];
		}
		cout<<"checking node: "<<(*itNode)->Id<<" sum of Packing forces: "<<sumPack[0]<<" "<<sumPack[1]<<" "<<sumPack[2]<<endl;
		*/
		//if ((*itNode)->Id == 60) {cout<<" checking node : "<<(*itNode)->Id<<endl;}
		bool NodeHasPacking = (*itNode)->checkIfNodeHasPacking();
		if (NodeHasPacking){
			//cout<<"	node has packing"<<endl;
			double* pos;
			pos = new double[3];
			(*itNode)->getCurrentPosition(pos);
			vector <int> edgeNodeData0, edgeNodeData1;
			vector <double> edgeNormalDataX, edgeNormalDataY, edgeNormalDataZ;
			bool pushedBySurface = false;
			bool pushedByEdge = false;
			for (itEle=Elements.begin(); itEle<Elements.end(); ++itEle){
				//excluding elements that own this element
				//if ((*itNode)->Id == 60 && ( (*itEle)->Id == 184 ||  (*itEle)->Id == 195 ) ) {cout<<" checking 1st to element  : "<<(*itEle)->Id<<endl;}
				bool PackingToThisElement =  (*itEle)->checkPackingToThisNodeViaState(TissueHeightDiscretisationLayers, (*itNode));
				if (PackingToThisElement){
					//cout<<"	node may pack to element: "<<(*itEle)->Id<<endl;
					//the node does not belong to the element, and placed correctly. Lets have a preliminary distance check:
					//if ((*itNode)->Id == 60 && ( (*itEle)->Id == 184 ||  (*itEle)->Id == 195 ) ) {cout<<" checking 2nd to element  : "<<(*itEle)->Id<<endl;}
					PackingToThisElement = (*itEle)->IsPointCloseEnoughForPacking((*itNode)->Position, thresholdForNodeDetection, (*itNode)->tissuePlacement);
				}
				if (PackingToThisElement){
					//if ((*itNode)->Id == 60 && ( (*itEle)->Id == 184 ||  (*itEle)->Id == 195 ) ) {cout<<" continuing packing to element  : "<<(*itEle)->Id<<endl;}
					//cout<<"	node can pack to element: "<<(*itEle)->Id<<endl;
					//all positions and stated are correct, this node can pack to this element
					double* normalForPacking;
					normalForPacking = new double[3];
					normalForPacking[0]=0.0;normalForPacking[1]=0.0;normalForPacking[2]=0.0;
					double* posCorner;
					posCorner = new double[3];
					//getting the normal to the element surface in the direction of the pushing force, and one corner of the surface
					getNormalAndCornerPosForPacking((*itNode),(*itEle),normalForPacking,posCorner);
					double dInNormalDir = 0.0;
					//if ((*itNode)->Id == 78 ) {cout<<" checking distance to surface for  : "<<(*itEle)->Id<<endl;}
					for (int i=0; i<3; ++i){
						dInNormalDir += (pos[i] - posCorner[i])*normalForPacking[i];
					}
					//cout<<"	distance to surface in direction of normal: "<<dInNormalDir<<endl;
					if ((dInNormalDir > 0 && dInNormalDir < threshold) || (dInNormalDir < 0 && dInNormalDir > thresholdNeg)){
						//if the edge is not composed of any of my neigs, then add to list:
						addToEdgeList((*itNode), (*itEle), edgeNodeData0, edgeNodeData1);
						edgeNormalDataX.push_back(normalForPacking[0]);
						edgeNormalDataY.push_back(normalForPacking[1]);
						edgeNormalDataZ.push_back(normalForPacking[2]);
						//if ((*itNode)->Id == 60 && ( (*itEle)->Id == 184 ||  (*itEle)->Id == 195 ) ) {cout<<" checking projection to surface for  : "<<(*itEle)->Id<<endl;}
						//the node is close enough in the direction of the normal, now I should check if the projection falls within hte triangle:
						double projectedPoint[3];
						for (int i=0; i<3; ++i){
							projectedPoint[i] = pos[i] - dInNormalDir*normalForPacking[i];
						}
						//cout<<"	pos: " <<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
						//cout<<"	normalForPacking: "<<normalForPacking[0]<<" "<<normalForPacking[1]<<" "<<normalForPacking[2]<<endl;
						//cout<<"	projected point: "<<projectedPoint[0]<<" "<<projectedPoint[1]<<" "<<projectedPoint[2]<<endl;
						bool pointInsideTriangle = (*itEle)->IspointInsideTriangle((*itNode)->tissuePlacement,projectedPoint[0],projectedPoint[1],projectedPoint[2]);
						//cout<<"	pointInsideTriangle: "<<pointInsideTriangle<<endl;
						if (pointInsideTriangle){
							//if ((*itNode)->Id == 60 && ( (*itEle)->Id == 184 ||  (*itEle)->Id == 195 ) ) {cout<<" calculating force for  : "<<(*itEle)->Id<<endl;}
							//the point is within packing distance, calculate packing force now:
							//if dInNormalDir is negative, then the node penetrated the surface
							// if this is the case, I do not calculate force via distance, I set distance to zero:
							float d2;
							if (dInNormalDir<0){
								d2=0;
							}
							else{
								d2 = dInNormalDir*dInNormalDir;
							}
							double averageMass = 0.5*((*itEle)->VolumePerNode + (*itNode)->mass);
							//double Fmag = averageMass* multiplier * (1.0/d2 - 1.0/t2);
							//double Fmag = averageMass* multiplier * (1.0 - (d2/t2));
							double Fmag = averageMass* multiplier * (1.0 - pow((d2/t2),0.5));
							//cout<<"	Fmag before capping: "<<Fmag<<endl;
							CapPackingForce(Fmag);
							//cout<<"	Fmag after capping: "<<Fmag<<endl;
							double F[3] = {Fmag*normalForPacking[0],  Fmag*normalForPacking[1], Fmag*normalForPacking[2]};
							//cout<<" packing force: "<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
							//if (PackingForces[126][2]> 0 ){
							//	cout<<" force on 126 became non-zero via node, node Id: "<<(*itNode)->Id<<" prisim: "<<(*itEle)->Id<<endl;
							//}
							for(int i=0; i<3; ++i){
								if ((*itNode)->FixedPos[i]){
									F[i] = 0.0;
									//if the node is fixed on  a certain axis, then do not reflect the force on the surface
								}
							}

							bool allCornersFixed[3] = {false,false,false};
							(*itEle)->AddPackingToSurface((*itNode)->tissuePlacement, F[0],F[1],F[2], PackingForces, Nodes, allCornersFixed[0], allCornersFixed[1], allCornersFixed[2]);
							for(int i=0; i<3; ++i){
								if (!allCornersFixed[i] && !(*itNode)->FixedPos[i]){
									//only adding the force if I could reflect it on a surface and the node is not fixed (already made it zero above, keeping here as a reminder
									PackingForces[(*itNode)->Id][i] += F[i];
								}
							}

							//if (PackingForces[126][2]> 0 ){
							///	cout<<" force on 126 is non-zero after surface addition, node Id: "<<(*itNode)->Id<<" node tissue placement: "<<(*itNode)->tissuePlacement<<" has lateral owner? "<<(*itNode)->hasLateralElementOwner<<" prisim: "<<(*itEle)->Id<<" element tissue type: "<<(*itEle)->tissueType<<endl;
							//}
							pushedBySurface = true;
						}
					}
					delete[] normalForPacking;
					delete[] posCorner;
				}
			}
			int closestNode  = -1000;
			double currMinDist2 = 1000;
			double closestNodeVec[3] = {0.0,0.0,0.0};
			if (!pushedBySurface){
				//find closest node, if any edge is closer than the node, it will pack:
				int n = edgeNodeData0.size();
				bool isNodeChecked[(const int) n];
				bool isEdgeChecked[(const int) n];
				for( int i = 0; i<n; ++i){
					isNodeChecked[i] = false;
					isEdgeChecked[i] = false;
				}
				for( int i = 0; i<n; ++i){
					bool isNeig = (*itNode)->checkIfNeighbour(edgeNodeData0[i]);
					if (isNeig){
						isNodeChecked[i] = true;
						isEdgeChecked[i] = true;
					}
					else{
						//check the other end of the edge:
						isNeig = (*itNode)->checkIfNeighbour(edgeNodeData1[i]);
						if (isNeig){
								isEdgeChecked[i] = true;
						}
					}
				}
				//clear all the nodes and corresponding edges that are my immeditate neigs:

				for (int i=0; i<n; ++i){
					//looping over the nodes
					if (!isNodeChecked[i]){ //do not double check the node if it has been checked previously
						int node0 = edgeNodeData0[i];
						//go over all nodes in the list (after this one) and mark all same id as checked:
						for (int j=i; j<n; ++j){
							if (edgeNodeData0[i] == edgeNodeData0[j]){
								isNodeChecked[j] = true;
							}
						}
						//calculate distace:
						double dx = (*itNode)->Position[0] - Nodes[node0]->Position[0];
						double dy = (*itNode)->Position[1] - Nodes[node0]->Position[1];
						double dz = (*itNode)->Position[2] - Nodes[node0]->Position[2];
						double d2 = dx*dx +  dy*dy  + dz*dz ;
						if (d2 < currMinDist2){
							currMinDist2 = d2;
							closestNode = node0;
							closestNodeVec[0] = dx;
							closestNodeVec[1] = dy;
							closestNodeVec[2] = dz;
						}
					}
				}
				//found closest node
				//if ((*itNode)->Id == 78 ) {cout<<" checking edge packing for node : "<<(*itNode)->Id<<endl;}
				//the node have not been pushed by nay of the surfaces, I need to check if it is being pushed by edges:
				n = edgeNormalDataX.size();
				for (int i=0; i<n; ++i){
					//looping over the normals vectors
					for (int j =0; j<3; ++j){
						//three node couples recorded for each node data point
						if (!isEdgeChecked[3*i+j]){ //do not double check the edge if it has been checked previously
							double normal1[3] = {0.0,0.0,1.0};
							int node0 = edgeNodeData0[3*i+j];
							int node1 = edgeNodeData1[3*i+j];
							double normal0[3] = {edgeNormalDataX[i],edgeNormalDataY[i],edgeNormalDataZ[i]};
							bool secondNormalFound = false;
							isEdgeChecked[3*i+j] = true;
							//I need to start looing at the nodes that are after this triplet, to find the second occurace of this edge:
							for (int k=i+1; k<n; ++k){
								if (!isEdgeChecked[3*k] && (edgeNodeData0[3*k] == node0 || edgeNodeData0[3*k] == node1)){
									if (edgeNodeData1[3*k] == node0 || edgeNodeData1[3*k] == node1){
										normal1[0] = edgeNormalDataX[k];
										normal1[1] = edgeNormalDataY[k];
										normal1[2] = edgeNormalDataZ[k];
										isEdgeChecked[3*k] = true;
										secondNormalFound = true;
										break;
									}
								}
								if (!isEdgeChecked[3*k+1] && (edgeNodeData0[3*k+1] == node0 || edgeNodeData0[3*k+1] == node1)){
									if (edgeNodeData1[3*k+1] == node0 || edgeNodeData1[3*k+1] == node1){
										normal1[0] = edgeNormalDataX[k];
										normal1[1] = edgeNormalDataY[k];
										normal1[2] = edgeNormalDataZ[k];
										isEdgeChecked[3*k+1] = true;
										secondNormalFound = true;
										break;
									}
								}
								if (!isEdgeChecked[3*k+2] && (edgeNodeData0[3*k+2] == node0 || edgeNodeData0[3*k+2] == node1)){
									if (edgeNodeData1[3*k+2] == node0 || edgeNodeData1[3*k+2] == node1){
										normal1[0] = edgeNormalDataX[k];
										normal1[1] = edgeNormalDataY[k];
										normal1[2] = edgeNormalDataZ[k];
										isEdgeChecked[3*k+2] = true;
										secondNormalFound = true;
										break;
									}
								}
							}
							if (!secondNormalFound){
								//the edge is recorded only once, the second time is not recorded.
								//the edge is either on the edge of the tissue, or it is at the symmetricity border
								normal1[0] = normal0[0];
								normal1[2] = normal0[2];
								if (Nodes[node0]->atSymmetricityBorder && Nodes[node1]->atSymmetricityBorder){
									normal1[1] = (-1.0)*normal0[1];
								}
								else{
									normal1[1] = normal0[1];
								}
							}
							//if ((*itNode)->Id == 78 ) {
							//	cout<<" checking edge couple : "<<node0<<" "<<node1<<endl;
							//	cout<<"normals: "<<normal0[0]<<" "<<normal0[1]<<" "<<normal0[2]<<" / "<<normal1[0]<<" "<<normal1[1]<<" "<<normal1[2]<<endl;
							//}
							//now I have normals from both sides for this edge
							double v[3]; //vector from the edge to the point;
							bool isPointOnEdge = checkIfPointIsOnEdge(node0, node1, pos[0], pos[1], pos[2], v[0], v[1], v[2]);
							//if ((*itNode)->Id == 78 ) {
							//	cout<<" is point on the edge? : "<<isPointOnEdge<<endl;
							//}
							if (isPointOnEdge){
								double d2 = v[0]*v[0] +  v[1]*v[1] +  v[2]*v[2];
								//if ((*itNode)->Id == 78 ) {
								//	cout<<" point on the edge, distance squared is : "<<d2<<endl;
								//}
								if (d2 < t2 && d2 <currMinDist2){ // the edge is close enough to pack, and is also closer than the closest node:
									//if ((*itNode)->Id == 78 ) {
									//	cout<<"distance was close enough, nodes were: "<<node0<<" "<<node1<<endl;
									//}
									double d = pow(d2,0.5);
									v[0] /= d;
									v[1] /= d;
									v[2] /= d;
									//check if the direction is correct:
									double dot = (normal0[0] + normal1[0])*v[0] + (normal0[1] + normal1[1])*v[1] + (normal0[2] + normal1[2])*v[2];
									if (dot < 0 ){
										//point penetrated into the elemetns, the pushing force should be inverted:
										v[0] *= (-1.0);
										v[1] *= (-1.0);
										v[2] *= (-1.0);
									}
									//point close enough for packing:
									double averageMass = 1.0/3.0 *( (*itNode)->mass + Nodes[node0]->mass + Nodes[node1]->mass );
									//double Fmag = averageMass* multiplier * (1.0 - (d2/t2));
									double Fmag = averageMass* multiplier * (1.0 - pow((d2/t2),0.5));
									CapPackingForce(Fmag);
									double F[3] = {Fmag*v[0],  Fmag*v[1], Fmag*v[2]};
									//cout<<" packing force: "<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
									for(int i=0; i<3; ++i){
										/*if (!(*itNode)->FixedPos[i]){
											PackingForces[(*itNode)->Id][i] += F[i];
										}
										//add to the edges
										if (!Nodes[node0]->FixedPos[i]){
											PackingForces[node0][i] -= 0.5*F[i];
										}
										if (!Nodes[node1]->FixedPos[i]){
											PackingForces[node1][i] -= 0.5*F[i];
										}*/
										if (!(*itNode)->FixedPos[i]){
											//the node being checked for packing is not fixed on this axis, may add forces:
											bool edgesFixed[2] = {Nodes[node0]->FixedPos[i], Nodes[node1]->FixedPos[i]};
											if (!edgesFixed[0] && !edgesFixed[1]){
												//both edges  are not fixed, add forces normally:
												PackingForces[(*itNode)->Id][i] += F[i];
												PackingForces[node0][i] -= 0.5*F[i];
												PackingForces[node1][i] -= 0.5*F[i];
											}
											else if(!edgesFixed[0] || !edgesFixed[1]){
												//at least one of the edges are not fixed, find which one, and add all the force on it:
												PackingForces[(*itNode)->Id][i] += F[i];
												if (!edgesFixed[0]){
													PackingForces[node0][i] -= F[i];
												}else {
													PackingForces[node1][i] -= F[i];
												}
											}
											//the other condition is that both edges are fixed, and the node should not be applying force in that direction
										}
									}
									pushedByEdge = true;
								}
							}
						}
					}
				}
			}
			//if ((*itNode)->Id == 78 ) {
			//	cout<<" pushedByEdge? : "<<pushedByEdge<<endl;
			//}
			if(!pushedBySurface && !pushedByEdge){
				//if ((*itNode)->Id == 78 ) {cout<<" checking node based packing for node : "<<(*itNode)->Id<<endl;}
				//checking point to point:
				if (currMinDist2 < t2 ){
					//if ((*itNode)->Id == 211 ) {
					//	cout<<" node 211 is packing to node : "<<closestNode<<endl;
					//}
					double d = pow(currMinDist2,0.5);
					closestNodeVec[0] /= d;
					closestNodeVec[1] /= d;
					closestNodeVec[2] /= d;
					//point close enough for packing:
					double averageMass = 1.0/2.0 *( (*itNode)->mass + Nodes[closestNode]->mass );
					//double Fmag = averageMass* multiplier * (1.0 - (currMinDist2/t2));
					double Fmag = averageMass* multiplier * (1.0 - pow((currMinDist2/t2),0.5));
					CapPackingForce(Fmag);
					double F[3] = {Fmag*closestNodeVec[0],  Fmag*closestNodeVec[1], Fmag*closestNodeVec[2]};
							//cout<<" packing force: "<<F[0]<<" "<<F[1]<<" "<<F[2]<<endl;
					for(int i=0; i<3; ++i){
						//if (!(*itNode)->FixedPos[i]){
						//	PackingForces[(*itNode)->Id][i] += F[i];
						//}
						//add to the edges
						//if (!Nodes[closestNode]->FixedPos[i]){
						//	PackingForces[closestNode][i] -= F[i];
						//}
						//if both the packing node and the corresponding node are not fixed:
						if (!(*itNode)->FixedPos[i] && !Nodes[closestNode]->FixedPos[i]){
							PackingForces[closestNode][i] -= F[i];
						}
					}
				}
			}
			//if ((*itNode)->Id == 211 ) {cout<<"node 211: packing to surf: "<<pushedBySurface<<" packing to edge: "<<pushedByEdge<<endl;}
			delete[] pos;
		}
	}
}

bool Simulation::checkIfPointIsOnEdge(int node0, int node1, double x, double y, double z, double& vx, double& vy, double& vz){
	double vec0[3] = {Nodes[node1]->Position[0] - Nodes[node0]->Position[0], Nodes[node1]->Position[1] - Nodes[node0]->Position[1], Nodes[node1]->Position[2]- Nodes[node0]->Position[2]};
	double vec1[3] = {x - Nodes[node0]->Position[0], y - Nodes[node0]->Position[1], z- Nodes[node0]->Position[2]};
	//if ((node0 == 202 || node0 == 220) && (node1==202 || node1 == 220) ){
	//	cout<<"vec of edge: "<<vec0[0]<<" "<<vec0[1]<<" "<<vec0[2]<<" vec from "<<node0<<" to curr node: " <<vec1[0]<<" "<<vec1[1]<<" "<<vec1[2]<<endl;
	//}
	//get length of edge and normalise edge vector
	double v0Mag = pow(vec0[0]*vec0[0] + vec0[1]*vec0[1]+vec0[2]*vec0[2],0.5);
	vec0[0] /= v0Mag;
	vec0[1] /= v0Mag;
	vec0[2] /= v0Mag;
	//if ((node0 == 202 || node0 == 220) && (node1==202 || node1 == 220) ){
	//	cout<<"vec of edge ("<<node0<<" - "<<node1<<") after normalisaiton: "<<vec0[0]<<" "<<vec0[1]<<" "<<vec0[1]<<" mag was:  "<<v0Mag<<endl;
	//}
	//dot product of normalised vec0 and vec1 should be positive, and should be smaller than the length of vec0:
	double dot = vec0[0]*vec1[0] + vec0[1]*vec1[1] + vec0[2]*vec1[2];
	//if ((node0 == 202 || node0 == 220) && (node1==202 || node1 == 220) ){
	//	cout<<"dotp: "<<dot<<endl;
	//}
	if (dot > 0 && dot < v0Mag){
		//the point falls within the edge,
		//now record the normal vector from edge to point:
		vx = vec1[0] - vec0[0]*dot;
		vy = vec1[1] - vec0[1]*dot;
		vz = vec1[2] - vec0[2]*dot;
		return true;
	}
	vx = 0.0;
	vy = 0.0;
	vz = 0.0;
	return false;
}

void Simulation::addToEdgeList(Node* nodePointer, ShapeBase* elementPointer, vector <int> & edgeNodeData0, vector <int> & edgeNodeData1 ){
	int  E0Index = -1, E1Index = -1, E2Index = -1;
	if (nodePointer->tissuePlacement == 0 ){ //checking basal surface,
		if (elementPointer->tissueType == 0){ //element is columnar
			E0Index = 0;
			E1Index = 1;
			E2Index = 2;
		}
		else{//element is periodial:
			E0Index = 3;
			E1Index = 4;
			E2Index = 5;
		}
	}
	else if (nodePointer->tissuePlacement  == 1 ){ //checking apical surface
		if (elementPointer->tissueType == 0){ //element is columnar
			E0Index = 3;
			E1Index = 4;
			E2Index = 5;
		}
		else{//element is periodial:
			E0Index = 0;
			E1Index = 1;
			E2Index = 2;
		}
	}

	edgeNodeData0.push_back(elementPointer->getNodeId(E0Index));
	edgeNodeData1.push_back(elementPointer->getNodeId(E1Index));
	edgeNodeData0.push_back(elementPointer->getNodeId(E0Index));
	edgeNodeData1.push_back(elementPointer->getNodeId(E2Index));
	edgeNodeData0.push_back(elementPointer->getNodeId(E1Index));
	edgeNodeData1.push_back(elementPointer->getNodeId(E2Index));
}

void Simulation::addPackingForces(gsl_matrix* gExt){
	double sumPack[3] = {0.0,0.0,0.0};
	//double sumPackPre[3] = {0.0,0.0,0.0};
	//double sumAv[3] = {0.0,0.0,0.0};
	//double sumgExt[3] = {0.0,0.0,0.0};
	//cout<<"in add packing forces"<<endl;
	for (int j=0; j<nNodes; ++j){
		//cout<<"calculating node: "<<j<<" of "<<nNodes<<endl;
		//double Fx = 0.333* ( PackingForces[j][0] + PackingForcesPreviousStep[j][0] + PackingForcesTwoStepsAgoStep[j][0] );
		//double Fy = 0.333* ( PackingForces[j][1] + PackingForcesPreviousStep[j][1] + PackingForcesTwoStepsAgoStep[j][1] );
		//double Fz = 0.333* ( PackingForces[j][2] + PackingForcesPreviousStep[j][2] + PackingForcesTwoStepsAgoStep[j][2] );
		double Fx = PackingForces[j][0];
		double Fy = PackingForces[j][1];
		double Fz = PackingForces[j][2];
		sumPack[0] += PackingForces[j][0];
		//sumPackPre[0] += PackingForcesPreviousStep[j][0];
		//sumAv[0] += Fx;
		sumPack[1] += PackingForces[j][1];
		//sumPackPre[1] += PackingForcesPreviousStep[j][1];
		//sumAv[1] += Fy;
		sumPack[2] += PackingForces[j][2];
		//sumPackPre[2] += PackingForcesPreviousStep[j][2];
		//sumAv[2] += Fz;
		//cout<<"node: "<<j<<" after summation"<<endl;
		int indice = j*3;
		Fx += gsl_matrix_get(gExt,indice,0);
		gsl_matrix_set(gExt,indice,0,Fx);
		Fy += gsl_matrix_get(gExt,indice+1,0);
		gsl_matrix_set(gExt,indice+1,0,Fy);
		Fz += gsl_matrix_get(gExt,indice+2,0);
		gsl_matrix_set(gExt,indice+2,0,Fz);
		//sumgExt[0] += gsl_matrix_get(gExt,indice,0);
		//sumgExt[1] += gsl_matrix_get(gExt,indice+1,0);
		//sumgExt[2] += gsl_matrix_get(gExt,indice+2,0);
	}
	//cout<<" sum pack: "<<sumPack[0]<<" "<<sumPack[1]<<" "<<sumPack[2]<<endl;
	//cout<<" sum gExt: "<<sumgExt[0]<<" "<<sumgExt[1]<<" "<<sumgExt[2]<<endl;
}

void Simulation::updateElementPositions(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updatePositions(Nodes);
	}
}


void Simulation::updateElementPositionsSingle(int i ){
	Elements[i]->updatePositions(Nodes);
}

void Simulation::assignTips(){
	double xTips[2] ={ 10000, -10000};
	double yTips[2] ={ 10000, -10000};
	for (int i =0; i<nNodes; ++i){
		if (Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0 ) { //taking the basal nodes of the columnar layer only
			if (Nodes[i]->Position[0] < xTips[0] ){
				dorsalTipIndex = i;
				xTips[0] = Nodes[i]->Position[0];
			}
			if (Nodes[i]->Position[1] < yTips[0] ){
				anteriorTipIndex = i;
				yTips[0] = Nodes[i]->Position[1];
			}
			if (Nodes[i]->Position[0] > xTips[1] ){
				ventralTipIndex = i;
				xTips[1] = Nodes[i]->Position[0];
			}
			if (Nodes[i]->Position[1] > yTips[1] ){
				posteriorTipIndex = i;
				yTips[1] = Nodes[i]->Position[1];
			}
		}
		//cout<<"DV Tip node indexes: "<<dorsalTipIndex<<" "<<ventralTipIndex<<endl;
		//cout<<"AP Tip node indexes: "<<anteriorTipIndex<<" "<<posteriorTipIndex<<endl;
	}
	//I have identified the tips assuming the tissue is represented in full.
	//If there is symmetry in the input, then I will need to correct this.
	//In the case of x symmetry, the minimum x will be correct in DV, but not the maximum x.
	//In the case of y symmetry, the maximum y will be correct in AP. but not the minimum.
	if (symmetricX){
		double delta = 0.02;
		//the ventralTipIndex is correct, it is at the minimum of x axis. The dorsalTipIndex is at x=0. but there are a series of nodes at x=0,
		//and the selected one should be at y=0 too.
		for (int i =0; i<nNodes; ++i){
				if (Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0 ) { //taking the basal nodes of the columnar layer only
					if (Nodes[i]->Position[0]< delta &&  Nodes[i]->Position[0]> -delta){
						//this node is at x = 0, is it at y = 0 too?
						if (Nodes[i]->Position[1]< delta &&  Nodes[i]->Position[1]> -delta){
							dorsalTipIndex = i;
							break;
						}
					}
				}
		}
	}
	if (symmetricY){
		double delta = 0.02;
		//the posteriorTipIndex is correct, it is at the maximum of y axis. The anteriorTipIndex is at y=0. but there are a series of nodes at y=0,
		//and the selected one should be at x=0 too.
		for (int i =0; i<nNodes; ++i){
				if (Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0 ) { //taking the basal nodes of the columnar layer only
					if (Nodes[i]->Position[1]< delta &&  Nodes[i]->Position[1]> -delta){
						//this node is at x = 0, is it at y = 0 too?
						if (Nodes[i]->Position[0]< delta &&  Nodes[i]->Position[0]> -delta){
							anteriorTipIndex = i;
							break;
						}
					}
				}
		}
	}
	cout<<"Tip node indexes - dorsalTipIndex: "<<dorsalTipIndex<<" ventralTipIndex: "<<ventralTipIndex<<endl;
	cout<<"Tip node indexes - anteriorTipIndex: "<<anteriorTipIndex<<" posteriorTipIndex: "<<posteriorTipIndex<<endl;
}

void Simulation::alignTissueDVToXPositive(){
	if (!symmetricX){
		double* u = new double[3];
		double* v = new double[3];
		//For simulations with no external viscosity, the position of Dorsal tip is fixed, Ventral tip is fixed in y and z, another node is fixed in z
		for (int i=0;i<3;++i){
			u[i] = Nodes[ventralTipIndex]->Position[i] - Nodes[dorsalTipIndex]->Position[i];
		}
		//cout<<" ventralTipIndex: "<<ventralTipIndex<<" dorsalTipIndex "<<dorsalTipIndex<<endl;
		double dummy = Elements[0]->normaliseVector3D(u);
		v[0]=1;v[1]=0;v[2]=0;
		double c, s;
		Elements[0]->calculateRotationAngleSinCos(u,v,c,s);
		double *rotAx;
		rotAx = new double[3];
		double *rotMat;
		rotMat = new double[9]; //matrix is written in one row
		Elements[0]->calculateRotationAxis(u,v,rotAx,c);	//calculating the rotation axis that is perpendicular to both u and v
		Elements[0]->constructRotationMatrix(c,s,rotAx,rotMat);
		for(int i=1;i<nNodes;++i){
			for(int j = 0; j< Nodes[i]->nDim; ++j){
				u[j] = Nodes[i]->Position[j] - Nodes[0]->Position[j];
			}
			Elements[0]->rotateVectorByRotationMatrix(u,rotMat);

			for(int j = 0; j< Nodes[i]->nDim; ++j){
				Nodes[i]->Position[j] = u[j] + Nodes[0]->Position[j];
			}
		}
		vector<ShapeBase*>::iterator itElement;
		for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			(*itElement)->updatePositions(Nodes);
		}
		delete[] rotAx;
		delete[] rotMat;
		delete[] u;
		delete[] v;
	}
}



void Simulation::alignTissueAPToXYPlane(){
	double* u = new double[3];
	double* v = new double[3];
	for (int i=0;i<3;++i){
		u[i] = Nodes[anteriorTipIndex]->Position[i] - Nodes[posteriorTipIndex]->Position[i];
	}
	double dummy = Elements[0]->normaliseVector3D(u);
	v[0]=u[0];v[1]=u[1];v[2]=0;
	dummy = Elements[0]->normaliseVector3D(v);
	double c, s;
	Elements[0]->calculateRotationAngleSinCos(u,v,c,s);
	double *rotAx;
	rotAx = new double[3];
	double *rotMat;
	rotMat = new double[9]; //matrix is written in one row
	Elements[0]->calculateRotationAxis(u,v,rotAx,c);	//calculating the rotation axis that is perpendicular to both u and v
	Elements[0]->constructRotationMatrix(c,s,rotAx,rotMat);
	for(int i=0;i<nNodes;++i){
		for(int j = 0; j< Nodes[i]->nDim; ++j){
			u[j] = Nodes[i]->Position[j];
		}
		Elements[0]->rotateVectorByRotationMatrix(u,rotMat);
		for(int j = 0; j< Nodes[i]->nDim; ++j){
			Nodes[i]->Position[j] = u[j];
		}
	}
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updatePositions(Nodes);
	}
	delete[] rotAx;
	delete[] rotMat;
	delete[] u;
	delete[] v;
}

void Simulation::calculateDVDistance(){
	double d[3];
	if (symmetricX){
		for (int i=0;i<3;++i){
			d[i] = Nodes[ventralTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	else{
		for (int i=0;i<3;++i){
			d[i] = Nodes[dorsalTipIndex]->Position[i] - Nodes[ventralTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	double dmag = d[0]+d[1]+d[2];
	dmag = pow(dmag,0.5);
	outputFile<<"time: "<<currSimTimeSec<<" DV distance is: "<<dmag<<" ";
	cout<<"time: "<<currSimTimeSec<<" DV distance is: "<<dmag<<" ";
	if (symmetricY){
		for (int i=0;i<3;++i){
			d[i] = Nodes[posteriorTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	else{
		for (int i=0;i<3;++i){
			d[i] = Nodes[anteriorTipIndex]->Position[i] - Nodes[posteriorTipIndex]->Position[i];
			d[i] *=d[i];
		}
	}
	dmag = d[0]+d[1]+d[2];
	dmag = pow(dmag,0.5);
	outputFile<<" AP distance is: "<<dmag<<endl;
	cout<<" AP distance is: "<<dmag<<endl;
	//cout<<"time: "<<currSimTimeSec<<" Node 0 position: "<<Nodes[0]->Position[0]<<" "<<Nodes[0]->Position[1]<<" "<<Nodes[0]->Position[2]<<endl;
	//cout<<"TIPS: "<<dorsalTipIndex<<" "<<ventralTipIndex<<" "<<anteriorTipIndex<<" "<<posteriorTipIndex<<endl;
}

void Simulation::resetForces(bool resetPacking){
	int dim = 3;
	//n nodes
	for (int j=0;j<nNodes;++j){
		//3 dimensions
		for (int k=0;k<dim;++k){
			SystemForces[j][k]=0.0;
			//packing forces should be reset at the beginning of the step, but not during NR iterations
			if (resetPacking){
				//PackingForcesTwoStepsAgoStep[j][k] = PackingForcesPreviousStep[j][k];
				//PackingForcesPreviousStep[j][k] = PackingForces[j][k];
				PackingForces[j][k]=0.0;
			}
			if (recordForcesOnFixedNodes){
				FixedNodeForces[j][k]=0.0;
			}
		}
	}
}


void Simulation::calculateApicalSize(){
	double sizeLim[2][2] = {{0.0,0.0},{0.0,0.0}};
	bool found[2][2] = {{false,false},{false,false}};
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType == 0 && (*itNode)->tissuePlacement == 1){ //checking only columnar apical layer nodes
			for (int i=0; i<2; ++i){
				if ( (*itNode)->Position[i] < sizeLim[0][i] ){
					sizeLim[0][i] = (*itNode)->Position[i];
					found[0][i] = true;
				}
				else if((*itNode)->Position[i]>sizeLim[1][i]){
					sizeLim[1][i] = (*itNode)->Position[i];
					found[1][i] = true;
				}
			}
		}
	}
	if (!found[0][0] && !found[0][1] && !found[1][0] && !found[1][1]){
		cerr<<" error in apical bounding  box calculation! Found? :"<<found[0][0]<<" "<<found[0][1]<<" "<<found[1][0]<<" "<<found[1][1]<<endl;
	}
	double DV = sizeLim[1][0] - sizeLim[0][0];
	double AP = sizeLim[1][1] - sizeLim[0][1];
	outputFile<<"At time: "<<currSimTimeSec<<" apical bounding box size: "<<DV<<" "<<AP<<endl;
}

void Simulation::calculateBoundingBox(){
	boundingBox[0][0] =  100000.0;	//lower left x
	boundingBox[0][1] =  100000.0;	//lower left y
	boundingBox[0][2] =  100000.0;	//lower z
	boundingBox[1][0] = -100000.0;	//upper right x
	boundingBox[1][1] = -100000.0;	//upper right y
	boundingBox[1][2] = -100000.0;	//upper z
	bool found[2][3] = {{false,false,false},{false,false,false}};
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		//Do not count node if it is part of an explicit ECM:
		if(!(*itNode)->allOwnersECMMimicing){
		//if (!(*itNode)->allOwnersAblated){
			//There is at least one element owning this node that is not ablated
		//if ((*itNode)->tissueType == 0 || (*itNode)->tissueType==1) {	//only consider peripodial or columnar nodes, not the linkers
			for (int i=0; i<(*itNode)->nDim; ++i){
				if ( (*itNode)->Position[i] < boundingBox[0][i] ){
					boundingBox[0][i] = (*itNode)->Position[i];
					found[0][i] = true;
				}
				else if((*itNode)->Position[i]>boundingBox[1][i]){
					boundingBox[1][i] = (*itNode)->Position[i];
					found[1][i] = true;
				}
			}
		}
		//}
		//}
	}
	//cout<<"calculating bounding box, symmetricity x & y: "<<symmetricX<<" "<<symmetricY<<endl;
	//cout<<"bounding box before update: "<<boundingBox[0][0]<<" "<<boundingBox[0][1]<<" "<<boundingBox[1][0]<<" "<<boundingBox[1][1]<<endl;
	if (symmetricY){
		boundingBox[0][1] = (-1.0)*boundingBox[1][1]; //if there is Y symmetricity, then the bounding box is extended to double the size in y
	}
	cout<<"bounding box after update: "<<boundingBox[0][0]<<" "<<boundingBox[0][1]<<" "<<boundingBox[1][0]<<" "<<boundingBox[1][1]<<endl;
	for (int i=0; i<3; ++i){
		boundingBoxSize[i] = boundingBox[1][i] - boundingBox[0][i];
	}
	if (!found[0][0] && !found[0][1] && !found[0][2] && !found[1][0] && !found[1][1] && !found[1][2]){
		cerr<<" error in bounding box calculation! Found? :"<<found[0][0]<<" "<<found[0][1]<<" "<<found[0][2]<<" "<<found[1][0]<<" "<<found[1][1]<<" "<<found[1][2]<<endl;
	}
}


void Simulation::bringMyosinStimuliUpToDate(){
	//need to go throug all myosin functions and apply them to here:
	//for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
	//	(*itElement)->bringMyosinStimuliUpToDate();
	//}
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->adjustCMyosinFromSave();
	}
}

/*
void Simulation::calculateColumnarLayerBoundingBox(){
	columnarBoundingBox[0][0] =  100000.0;	//lower left x
	columnarBoundingBox[0][1] =  100000.0;	//lower left y
	columnarBoundingBox[0][2] =  100000.0;	//lower z
	columnarBoundingBox[1][0] = -100000.0;	//upper right x
	columnarBoundingBox[1][1] = -100000.0;	//upper right y
	columnarBoundingBox[1][2] = -100000.0;	//upper z
	bool found[2][3] = {{false,false,false},{false,false,false}};
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType == 0){ //checking only columnar layer nodes
			for (int i=0; i<(*itNode)->nDim; ++i){
				if ( (*itNode)->Position[i] < columnarBoundingBox[0][i] ){
					columnarBoundingBox[0][i] = (*itNode)->Position[i];
					found[0][i] = true;
				}
				else if((*itNode)->Position[i]>columnarBoundingBox[1][i]){
					columnarBoundingBox[1][i] = (*itNode)->Position[i];
					found[1][i] = true;
				}
			}
		}
	}
	for (int i=0; i<3; ++i){
		columnarBoundingBoxSize[i] = columnarBoundingBox[1][i] - columnarBoundingBox[0][i];
	}
	if (!found[0][0] && !found[0][1] && !found[0][2] && !found[1][0] && !found[1][1] && !found[1][2]){
		cerr<<" error in bounding box calculation! Found? :"<<found[0][0]<<" "<<found[0][1]<<" "<<found[0][2]<<" "<<found[1][0]<<" "<<found[1][1]<<" "<<found[1][2]<<endl;
	}
}

void Simulation::calculatePeripodialBoundingBox(){
	peripodialBoundingBox[0][0] =  100000.0;	//lower left x
	peripodialBoundingBox[0][1] =  100000.0;	//lower left y
	peripodialBoundingBox[0][2] =  100000.0;	//lower z
	peripodialBoundingBox[1][0] = -100000.0;	//upper right x
	peripodialBoundingBox[1][1] = -100000.0;	//upper right y
	peripodialBoundingBox[1][2] = -100000.0;	//upper z
	bool found[2][3] = {{false,false,false},{false,false,false}};
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		if ((*itNode)->tissueType != 0){ //checking any node that is not columnar
			for (int i=0; i<(*itNode)->nDim; ++i){
				if ( (*itNode)->Position[i] < peripodialBoundingBox[0][i] ){
					peripodialBoundingBox[0][i] = (*itNode)->Position[i];
					found[0][i] = true;
				}
				else if((*itNode)->Position[i]>peripodialBoundingBox[1][i]){
					peripodialBoundingBox[1][i] = (*itNode)->Position[i];
					found[1][i] = true;
				}
			}
		}
	}
	for (int i=0; i<3; ++i){
		peripodialBoundingBoxSize[i] = peripodialBoundingBox[1][i] - peripodialBoundingBox[0][i];
	}
	if (!found[0][0] && !found[0][1] && !found[0][2] && !found[1][0] && !found[1][1] && !found[1][2]){
		cerr<<" error in bounding box calculation! Found? :"<<found[0][0]<<" "<<found[0][1]<<" "<<found[0][2]<<" "<<found[1][0]<<" "<<found[1][1]<<" "<<found[1][2]<<endl;
	}
}
*/
void Simulation::saveStep(){
	outputFile<<"Saving step: "<< timestep<<" this is :"<<currSimTimeSec<<" sec ("<<currSimTimeSec/3600<<" hr )"<<endl;
	writeSaveFileStepHeader();
	writeNodes();
	writeElements();
	writeSaveFileStepFooter();
	writeTensionCompression();
    writeGrowth();
	writeForces();
	writePacking();
	writeProteins();
	writePhysicalProp();
	writeGrowthRedistribution();
	writeNodeBinding();
	writeCollapseAndAdhesion();
}

void Simulation::writeSaveFileStepHeader(){
    saveFileMesh<<"=============== TIME: ";
	saveFileMesh.precision(6);
	saveFileMesh.width(10);
	saveFileMesh<<currSimTimeSec;
	saveFileMesh<<"==================================================="<<endl;
}

void Simulation::writeSaveFileStepFooter(){
	saveFileMesh<<"=============== END OF TIME: ";
	saveFileMesh.precision(6);
	saveFileMesh.width(10);
	saveFileMesh<<currSimTimeSec;
	saveFileMesh<<"============================================"<<endl;
}

void Simulation::writeNodes(){
	saveFileMesh<<nNodes<<endl;
	for (int i = 0; i<nNodes; ++i){
		if(Nodes[i]->nDim==2){
			saveFileMesh.precision(10);saveFileMesh.width(20);
			saveFileMesh<<Nodes[i]->Position[0];
			saveFileMesh.precision(10);saveFileMesh.width(20);
			saveFileMesh<<Nodes[i]->Position[1];
			saveFileMesh.precision(10);saveFileMesh.width(20);
			saveFileMesh<<0.0;

		}
		else{
			int ndim = Nodes[i]->nDim;
			for (int j=0;j<ndim;++j){
				saveFileMesh.precision(10);saveFileMesh.width(20);
				saveFileMesh<<Nodes[i]->Position[j];
			}

		}
		saveFileMesh.width(4);
		saveFileMesh<<Nodes[i]->tissuePlacement;
		saveFileMesh.width(4);
		saveFileMesh<<Nodes[i]->tissueType;
		saveFileMesh.width(4);
		saveFileMesh<<Nodes[i]->atCircumference;
		saveFileMesh<<endl;
	}
}

void Simulation::writeElements(){
	saveFileMesh<<nElements<<endl;
	for (int i = 0; i<nElements; ++i){
		int shapetype = Elements[i]->getShapeType();
		saveFileMesh.width(4);
		saveFileMesh<<shapetype;
		if(shapetype == 4 ){
			double height = Elements[i]->getElementHeight();
			saveFileMesh.precision(5);saveFileMesh.width(12);
			saveFileMesh<<height;
		}
		saveFileMesh.width(4);
		saveFileMesh<<Elements[i]->IsAblated;
		int nodeNumber = Elements[i]->getNodeNumber();
		int*  NodeIds = Elements[i]->getNodeIds();
		for (int j = 0; j<nodeNumber; ++j ){
			saveFileMesh.width(9);saveFileMesh<<NodeIds[j];
		}
		int dim  = Elements[i]->getDim();
		double** refPos = Elements[i]->getReferencePos();
		for (int j = 0; j<nodeNumber; ++j ){
			for (int k = 0; k<dim; ++k ){
				saveFileMesh.precision(4);saveFileMesh.width(15);
				saveFileMesh<<refPos[j][k];
			}
		}
		saveFileMesh<<endl;
	}
}

void Simulation::writeTensionCompression(){
	//for (int i=0;i<6;++i){
	//	cout<<" at timestep :"<< timestep<<" the plastic strains of element 0:	"<<Elements[0]->PlasticStrain(i)<<"	normal strain: 	"<<Elements[i]->Strain(0)<<endl;
	//}
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        for (int j=0; j<6; ++j){
            double S = gsl_matrix_get((*itElement)->Strain,j,0);
            saveFileTensionCompression.write((char*) &S, sizeof S);
        }
    }
	saveFileTensionCompression.flush();
}



void Simulation::writeGrowthRedistribution(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		bool thereIsDistribution = (*itElement)->thereIsGrowthRedistribution;
		bool shrinksElement = (*itElement)->growthRedistributionShrinksElement;
		saveFileGrowthRedistribution.write((char*) &thereIsDistribution, sizeof thereIsDistribution);
		saveFileGrowthRedistribution.write((char*) &shrinksElement, sizeof shrinksElement);
	}
	saveFileGrowthRedistribution.flush();
}

void Simulation::writeNodeBinding(){
	//count slave dof:
	//cout<<"writing node binding"<<endl;
	vector<int> slaveDOFs, masterDOFs;
	int counter = 0;
	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		for (int i=0; i<3; ++i){
			if ((*itNode)->slaveTo[i]>0){
				int slavedof = 3*(*itNode)->Id + i;
				int masterdof = 3*(*itNode)->slaveTo[i] + i;
				slaveDOFs.push_back(slavedof);
				masterDOFs.push_back(masterdof);
				counter++;
			}
		}
	}
	/*saveFileNodeBinding.write((char*) &counter, sizeof counter);
	for(int i=0;i<counter;++i){
		int slave  =slaveDOFs[i];
		int master = masterDOFs[i];
		saveFileNodeBinding.write((char*) &slave, sizeof slave);
		saveFileNodeBinding.write((char*) &master, sizeof master);
		cout<<"slaveDOFs :"<<slaveDOFs[i]<<" masterDOFs "<<masterDOFs[i]<<endl;
	}*/
	saveFileNodeBinding<<counter<<" "<<endl;
	for(int i=0;i<counter;++i){
		int slave  =slaveDOFs[i];
		int master = masterDOFs[i];
		saveFileNodeBinding<<slave<<" ";
		saveFileNodeBinding<<master<<" "<<endl;
		//cout<<"slaveDOFs :"<<slaveDOFs[i]<<" masterDOFs "<<masterDOFs[i]<<endl;
	}
}

void Simulation::writeCollapseAndAdhesion(){
	//[adhered to] [# of collapsed partners] [ids of collapsed partners]
	for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		saveFileCollapseAndAdhesion.write((char*) &(*itNode)->adheredTo, sizeof (*itNode)->adheredTo);
		saveFileCollapseAndAdhesion.write((char*) &(*itNode)->positionUpdateOngoing, sizeof (*itNode)->positionUpdateOngoing);
		saveFileCollapseAndAdhesion.write((char*) &(*itNode)->positionUpdateCounter, sizeof (*itNode)->positionUpdateCounter);
		int n = (*itNode)->collapsedWith.size();
		saveFileCollapseAndAdhesion.write((char*) &n, sizeof n);
		for(int i=0;i<n;++i){
			saveFileCollapseAndAdhesion.write((char*) &(*itNode)->collapsedWith[i], sizeof (*itNode)->collapsedWith[i]);
		}
	}
}

void Simulation::writeGrowth(){
    //for (int i=0;i<6;++i){
    //	cout<<" at timestep :"<< timestep<<" the plastic strains of element 0:	"<<Elements[0]->PlasticStrain(i)<<"	normal strain: 	"<<Elements[i]->Strain(0)<<endl;
    //}
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        gsl_matrix* currFg = (*itElement)->getFg();
        double* growthRate = (*itElement)->getGrowthRate();
        for (int j=0; j<3; ++j){
            for (int k=0; k<3; ++k){
                double Fgjk = gsl_matrix_get(currFg,j,k);
                saveFileGrowth.write((char*) &Fgjk, sizeof Fgjk);
            }
            saveFileGrowthRate.write((char*) &growthRate[j], sizeof growthRate[j]);
        }
        gsl_matrix_free(currFg);
    }
    saveFileGrowth.flush();
    saveFileGrowthRate.flush();
}



void Simulation::writeForces(){
	//Write system forces first, on a nodal basis
	for (int i=0;i<nNodes;++i){
		saveFileForces.write((char*) &SystemForces[i][0], sizeof SystemForces[i][0]);
		saveFileForces.write((char*) &SystemForces[i][1], sizeof SystemForces[i][1]);
		saveFileForces.write((char*) &SystemForces[i][2], sizeof SystemForces[i][2]);
	}
	//Then write myosin forces
	for (int i=0;i<nElements;++i){
		int n=Elements[i]->getNodeNumber();
		for (int j=0; j<n; ++j){
			saveFileForces.write((char*) &Elements[i]->MyoForce[j][0], sizeof Elements[i]->MyoForce[j][0]);
			saveFileForces.write((char*) &Elements[i]->MyoForce[j][1], sizeof Elements[i]->MyoForce[j][1]);
			saveFileForces.write((char*) &Elements[i]->MyoForce[j][2], sizeof Elements[i]->MyoForce[j][2]);
		}
	}
	saveFileForces.flush();
}

void Simulation::writePacking(){
	int n = pacingNodeCouples0.size();
	vector <int>::iterator itId0;
	vector <int>::iterator itId1;
	itId1=pacingNodeCouples1.begin();
	saveFilePacking.write((char*) &n, sizeof n);
	for (itId0 = pacingNodeCouples0.begin(); itId0 < pacingNodeCouples0.end(); ++itId0){
		saveFilePacking.write((char*) &(*itId0), sizeof (*itId0));
		saveFilePacking.write((char*) &(*itId1), sizeof (*itId1));
		itId1++;
	}
	for (int i=0;i<n;++i){
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples0[i]][0], sizeof PackingForces[pacingNodeCouples0[i]][0]);
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples0[i]][1], sizeof PackingForces[pacingNodeCouples0[i]][1]);
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples0[i]][2], sizeof PackingForces[pacingNodeCouples0[i]][2]);
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples1[i]][0], sizeof PackingForces[pacingNodeCouples1[i]][0]);
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples1[i]][1], sizeof PackingForces[pacingNodeCouples1[i]][1]);
		saveFilePacking.write((char*) &PackingForces[pacingNodeCouples1[i]][2], sizeof PackingForces[pacingNodeCouples1[i]][2]);
	}
	saveFilePacking.flush();
}

void Simulation::writeProteins(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double* cMyo = new double[4];
		double* cMyoEq = new double[4];
		(*itElement)->getMyosinLevels(cMyo);
		(*itElement)->getEquilibriumMyosinLevels(cMyoEq);
		saveFileProteins.write((char*) &cMyo[0], sizeof cMyo[0]);
		saveFileProteins.write((char*) &cMyo[1], sizeof cMyo[1]);
		saveFileProteins.write((char*) &cMyo[2], sizeof cMyo[2]);
		saveFileProteins.write((char*) &cMyo[3], sizeof cMyo[3]);
		saveFileProteins.write((char*) &cMyoEq[0], sizeof cMyoEq[0]);
		saveFileProteins.write((char*) &cMyoEq[1], sizeof cMyoEq[1]);
		saveFileProteins.write((char*) &cMyoEq[2], sizeof cMyoEq[2]);
		saveFileProteins.write((char*) &cMyoEq[3], sizeof cMyoEq[3]);
		delete[] cMyo;
		delete[] cMyoEq;
	}
	saveFileProteins.flush();
}

void Simulation::writePhysicalProp(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double E = (*itElement)->getYoungModulus();
		double interalVisc	= (*itElement)->getInternalViscosity();
		double zRemodellingSoFar = (*itElement)->getZRemodellingSoFar();
		saveFilePhysicalProp.write((char*) &E, sizeof E);
		saveFilePhysicalProp.write((char*) &interalVisc, sizeof interalVisc);
		saveFilePhysicalProp.write((char*) &zRemodellingSoFar, sizeof zRemodellingSoFar);
	}
	vector<Node*>::iterator itNode;
	for (itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		saveFilePhysicalProp.write((char*) &(*itNode)->externalViscosity[0], sizeof (*itNode)->externalViscosity[0]);
		saveFilePhysicalProp.write((char*) &(*itNode)->externalViscosity[1], sizeof (*itNode)->externalViscosity[1]);
		saveFilePhysicalProp.write((char*) &(*itNode)->externalViscosity[2], sizeof (*itNode)->externalViscosity[2]);
	}
	saveFilePhysicalProp.flush();
}

void Simulation::calculateMyosinForces(){
	//cout<<"Entered calculateMyosinForces"<<endl;
	cleanUpMyosinForces();
	bool basedOnArea = false;
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updateMyosinConcentration(dt, kMyo, thereIsMyosinFeedback, MyosinFeedbackCap);
		if (basedOnArea){
			//The forces are based on the total apical/basal area of the element, it may cause oscillations with rapid shape changes driven by myosin itself.
			(*itElement)->calculateMyosinForcesAreaBased(forcePerMyoMolecule);
		}
		else{
			//The forces are based on the total size of the element, it does not fluctuate with rapid shape changes, based on the grown volume
			(*itElement)->calculateMyosinForcesTotalSizeBased(forcePerMyoMolecule);
		}
	}
}

void Simulation::cleanUpMyosinForces(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->cleanMyosinForce();
	}
}

void Simulation::checkForMyosinUpdates(){
	for (int i=0; i<nMyosinFunctions; ++i){
		if(currSimTimeSec == myosinFunctions[i]->initTime  ){ //the application time of the signal is given in seconds.
			updateEquilibriumMyosinsFromInputSignal(myosinFunctions[i]);
		}
	}
}


void Simulation::updateEquilibriumMyosinsFromInputSignal(MyosinFunction* currMF){
	//cout<<"inside updateEquilibriumMyosinsFromInputSignal "<<endl;
	int nGridX = currMF->getGridX();
	int nGridY = currMF->getGridY();
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((currMF->applyToColumnarLayer && (*itElement)->tissueType == 0) || (currMF->applyToPeripodialMembrane && (*itElement)->tissueType == 1) ){//|| (*itElement)->tissueType == 2){
			if ( (*itElement)->spansWholeTissue
			     || (currMF->isApical && (*itElement)->tissuePlacement == 1)
			     || (!currMF->isApical && !(*itElement)->isECMMimicing && ( (*itElement)->tissuePlacement == 0 || (*itElement)->atBasalBorderOfECM) )
			     || currMF->isLateral
				){
				// 1) The element spans the whole tissue therefore both apical and basal responses should be applied
				// 2) The myosin response is applicable to apical surface, and the tissue placement of the element is apical,
				// 3) The myosin response is applicable to basal surface, and the tissue placement of the lement is basal.
				// 4) The myosin response is applied laterally, all elements should be checked
				if (currMF->manualStripes){
					//Here, I need to check if all the nodes of the element fall into the stripe of myosin up-regulations.
					//double stripeSize1, stripeSize2;
					//double initialPoint, endPoint;
					bool inActiveZone = (*itElement)->calculateIfInsideActiveStripe(currMF->initialPoint,currMF->endPoint,currMF->stripeSize1,currMF->stripeSize2);
					if (inActiveZone){
						if (currMF->isPolarised){
							(*itElement)->updateUnipolarEquilibriumMyosinConcentration(currMF->isApical,currMF->manualCMyoEq, currMF->manualOrientation[0], currMF->manualOrientation[1]);
						}
						else{
							(*itElement)->updateUniformEquilibriumMyosinConcentration(currMF->isApical, currMF->manualCMyoEq);
						}
					}
				}
				else if(currMF->useEllipses){
					bool applyToThisElement = (*itElement)->isMyosinViaEllipsesAppliedToElement(currMF->isApical, currMF->isLateral, myosinEllipseBandIds, numberOfMyosinAppliedEllipseBands);
					if (applyToThisElement){
						(*itElement)->updateUniformEquilibriumMyosinConcentration(currMF->isApical, currMF->manualCMyoEq);
					}
				}
				else{
					//calculating the grid indices:
					double* ReletivePos = new double[2];
					//normalising the element centre position with bounding box
					(*itElement)->getRelativePosInBoundingBox(ReletivePos);
					/*if ((*itElement)->tissueType == 0){
						(*itElement)->getRelativePosInColumnarBoundingBox(ReletivePos);
					}
					else{
						(*itElement)->getRelativePosInPeripodialBoundingBox(ReletivePos);
					}*/
					int indexX, indexY;
					double fracX, fracY;
					(*itElement)->convertRelativePosToGridIndex(ReletivePos, indexX, indexY, fracX, fracY, nGridX, nGridY);
					//reading the equilibrium myosin value
					double cEqYmid[2]= {0.0,0.0};
					double cEq;
					cEqYmid[0] = currMF->getEquilibriumMyoMatrixElement(indexX,indexY)*(1.0-fracX) + currMF->getEquilibriumMyoMatrixElement(indexX+1,indexY)*fracX;
					cEqYmid[1] = currMF->getEquilibriumMyoMatrixElement(indexX,indexY+1)*(1.0-fracX) + currMF->getEquilibriumMyoMatrixElement(indexX+1,indexY+1)*fracX;
					cEq = cEqYmid[0]*(1.0-fracY) + cEqYmid[1]*fracY;
					if (currMF->isPolarised){
						//If the function is polarised, reading the orientation
						double orientation[2];
						for (int axis=0; axis<2; ++axis){
							cEqYmid[0] = currMF->getOrientationMatrixElement(indexX,indexY,axis)*(1.0-fracX) + currMF->getOrientationMatrixElement(indexX+1,indexY,axis)*fracX;
							cEqYmid[1] = currMF->getOrientationMatrixElement(indexX,indexY+1,axis)*(1.0-fracX) + currMF->getOrientationMatrixElement(indexX+1,indexY+1,axis)*fracX;
							orientation[axis] = cEqYmid[0]*(1.0-fracY) + cEqYmid[1]*fracY;
					}
					//updating the values of the shape for polarised myosin:
					(*itElement)->updateUnipolarEquilibriumMyosinConcentration(currMF->isApical, cEq, orientation[0], orientation[1]);
					}
					else{
						//updating the values of the shape for uniform contractile:
						(*itElement)->updateUniformEquilibriumMyosinConcentration(currMF->isApical, cEq);
					}
					delete[] ReletivePos;
				}
			}
		}
	}
	//cout<<"finalised updateEquilibriumMyosinsFromInputSignal "<<endl;
}


void Simulation::calculateGrowth(){
	//cout<<"Calculating Growth"<<endl;
	cleanUpGrowthRates();
	for (int i=0; i<nGrowthFunctions; ++i){
		if (GrowthFunctions[i]->Type == 1){
			//cout<<"Calculating Uniform Growth, function: "<<i<<endl;
			calculateGrowthUniform(GrowthFunctions[i]);
		}
		else if(GrowthFunctions[i]->Type == 2){
			calculateGrowthRing(GrowthFunctions[i]);
		}
		else if(GrowthFunctions[i]->Type == 3){
			calculateGrowthGridBased(GrowthFunctions[i]);
		}
	}
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->isMutated){
			(*itElement)->updateGrowthByMutation(dt);
		}
	}
}

void Simulation::calculateShapeChange(){
	cleanUpShapeChangeRates();
	//cout<<"outside cleanUpShapeChangeRates"<<endl;
	for (int i=0; i<nShapeChangeFunctions; ++i){
		if (ShapeChangeFunctions[i]->Type == 1){
			//cout<<"calling calculate Shape change"<<endl;
			calculateShapeChangeUniform(ShapeChangeFunctions[i]);
		}
		if (ShapeChangeFunctions[i]->Type == 2){
			//cout<<"calling calculate Shape change"<<endl;
			calculateShapeChangeMarkerEllipseBased(ShapeChangeFunctions[i]);
		}
	}
}


void Simulation::cleanUpGrowthRates(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->setGrowthRate(dt,0.0,0.0,0.0);
		(*itElement)->updateGrowthIncrementFromRate();
	}
}

void Simulation::cleanUpShapeChangeRates(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->setShapeChangeRate(0.0,0.0,0.0,0.0,0.0,0.0);
		(*itElement)->setShapeChangeInrementToIdentity();
	}
}

void Simulation::calculateShapeChangeMarkerEllipseBased (GrowthFunctionBase* currSCF){
	if(currSimTimeSec >= currSCF->initTime && currSimTimeSec < currSCF->endTime ){
			gsl_matrix* columnarShapeChangeIncrement = gsl_matrix_calloc(3,3);
			double *growthRates = new double[3];
			currSCF->getShapeChangeRateRate(growthRates);
	    	vector<ShapeBase*>::iterator itElement;
	    	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
				bool appliedToElement = (*itElement)->isShapeChangeAppliedToElement(currSCF->appliedEllipseBandIds, currSCF->applyToBasalECM, currSCF->applyToLateralECM, currSCF->applyTissueApical, currSCF->applyTissueBasal, currSCF->applyTissueMidLine );
				if (appliedToElement){
					gsl_matrix_set_identity(columnarShapeChangeIncrement);
					(*itElement)->calculateShapeChangeIncrementFromRates(dt, growthRates[0],growthRates[1],growthRates[2],columnarShapeChangeIncrement);
					(*itElement)->updateShapeChangeIncrement(columnarShapeChangeIncrement);
				}
	    	}
			delete[] growthRates;
			gsl_matrix_free(columnarShapeChangeIncrement);
		}
}

void Simulation::calculateShapeChangeUniform (GrowthFunctionBase* currSCF){
	if(currSimTimeSec >= currSCF->initTime && currSimTimeSec < currSCF->endTime ){
		//cout<<"calculating shape change uniform"<<endl;
		double *maxValues;
        maxValues = new double[3];
        currSCF->getGrowthRate(maxValues);
    	vector<ShapeBase*>::iterator itElement;
    	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			//tissue type == 0 is columnar layer, ==1 is peripodial membrane, ==2 id linker zone
			if ( currSCF->applyToColumnarLayer){
				if ((*itElement)->tissueType == 0){ //columnar layer, grow directly
					(*itElement)->updateShapeChangeRate(maxValues[0],maxValues[1],maxValues[2],0,0,0);
				}
				else if ((*itElement)->tissueType == 2){ //Linker zone, need to weight the growth
					double weight = (*itElement)->getColumnarness();
					(*itElement)->updateShapeChangeRate(weight*maxValues[0],weight*maxValues[1],weight*maxValues[2],0,0,0);
				}
			}
			if ( currSCF->applyToPeripodialMembrane){
				if ((*itElement)->tissueType == 1){ //peripodial membrane, grow directly
					(*itElement)->updateShapeChangeRate(maxValues[0],maxValues[1],maxValues[2],0,0,0);
				}
				else if ((*itElement)->tissueType == 2){ //Linker zone, need to weight the growth
					double weight = (*itElement)->getPeripodialness();
					(*itElement)->updateShapeChangeRate(weight*maxValues[0],weight*maxValues[1],weight*maxValues[2],0,0,0);
				}
			}
		}
    	//cout<<"finalised shape change uniform"<<endl;
		delete[] maxValues;
	}
}

void Simulation::calculateGrowthUniform(GrowthFunctionBase* currGF){
	//cout<<"inside uniform growth function, initTime: "<<currGF->initTime <<" endtime: "<<currGF->endTime<<" simTime"<<simTime<<endl;
	if(currSimTimeSec >= currGF->initTime && currSimTimeSec < currGF->endTime ){
		gsl_matrix* columnarFgIncrement = gsl_matrix_calloc(3,3);
		gsl_matrix* peripodialFgIncrement = gsl_matrix_calloc(3,3);
		double *growthRates = new double[3];
		currGF->getGrowthRate(growthRates);
    	vector<ShapeBase*>::iterator itElement;
    	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			gsl_matrix_set_identity(columnarFgIncrement);
			gsl_matrix_set_identity(peripodialFgIncrement);
    		//if (!(*itElement)->isMutated){
				if (currGF->applyToBasalECM || currGF->applyToLateralECM){
					if (currGF->applyToBasalECM){
						if ((*itElement)->isECMMimicing && (*itElement)->tissuePlacement == 0){
							(*itElement)->calculateFgFromRates(dt, growthRates[0],growthRates[1],growthRates[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
						}
					}
					if (currGF->applyToLateralECM){
						if ((*itElement)->isECMMimimcingAtCircumference && !(*itElement)->tissuePlacement == 0){ //do not grow the basal element twice
							(*itElement)->calculateFgFromRates(dt, growthRates[0],growthRates[1],growthRates[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
						}
					}
				}
				if (!thereIsExplicitECM || !(*itElement)->isECMMimicing ){
					//There is either no explicit ECM definition, or the element is not ECM mimicing.
					//If there is explicit ECM, the basal elements should not grow, all others should proceed as usual
					//If there is no explicit ecm, then all should proceed as usual.
					if (currGF->applyToColumnarLayer){
						(*itElement)->calculateFgFromRates(dt, growthRates[0],growthRates[1],growthRates[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
					}
					if (currGF->applyToPeripodialMembrane){
						(*itElement)->calculateFgFromRates(dt, growthRates[0],growthRates[1],growthRates[2], currGF->getShearAngleRotationMatrix(), peripodialFgIncrement, 1, currGF->zMin, currGF->zMax);
					}
				}
				(*itElement)->updateGrowthIncrement(columnarFgIncrement,peripodialFgIncrement);
			//}
    	}
		delete[] growthRates;
		gsl_matrix_free(columnarFgIncrement);
		gsl_matrix_free(peripodialFgIncrement);
	}
}

void Simulation::calculateGrowthRing(GrowthFunctionBase* currGF){
	if(currSimTimeSec >= currGF->initTime && currSimTimeSec < currGF->endTime ){
		//The growth function is active at current time, now I will grow the elements.
		//First get the remaining data from the growth function parameters
		float centre[2];
		currGF->getCentre(centre[0], centre[1]);
		float innerRadius = currGF->getInnerRadius();
		float outerRadius = currGF->getOuterRadius();
		double* maxValues;
        maxValues = new double[3];
		currGF->getGrowthRate(maxValues);
		float innerRadius2 = innerRadius*innerRadius;
		float outerRadius2 = outerRadius*outerRadius;
		gsl_matrix* columnarFgIncrement = gsl_matrix_calloc(3,3);
		gsl_matrix* peripodialFgIncrement = gsl_matrix_calloc(3,3);
    	vector<ShapeBase*>::iterator itElement;
    	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			double* Elementcentre = new double[3];
			Elementcentre = (*itElement)->getCentre();
			//the distance is calculated in the x-y projection
			double d[2] = {centre[0] - Elementcentre[0], centre[1] - Elementcentre[1]};
			double dmag2 = d[0]*d[0] + d[1]*d[1];
			if (dmag2 > innerRadius2 && dmag2 < outerRadius2){
				//the element is within the growth zone.
				float distance = pow(dmag2,0.5);
				//calculating the growth rate: as a fraction increase within this time point
				double sf = (1.0 - (distance - innerRadius) / (outerRadius - innerRadius) );

				/*sf = 0.0;
				int arrayOfInterest[30];
				arrayOfInterest[0] = 3509;
				arrayOfInterest[1] = 3511;
				arrayOfInterest[2] = 3838;
				arrayOfInterest[3] = 3747;
				arrayOfInterest[4] = 3751;
				arrayOfInterest[5] = 3427;
				for (int aa=0;aa<4;++aa){
					for (int bb=0;bb<6; ++bb){
						arrayOfInterest[(aa+1)*6+bb] = arrayOfInterest[aa+bb] - 848;
					}
				}
				for (int aa=0;aa<30;++aa){
					cout<<"arrayOfInterest["<<aa<<"]: "<<arrayOfInterest[aa]<<endl;
				}
				for (int aa=0;aa<30;++aa){
					if((*itElement)->Id == arrayOfInterest[aa]){
						sf = 1.0;
						break;
					}
				}
				*/
				double growthscale[3] = {maxValues[0]*sf,maxValues[1]*sf,maxValues[2]*sf};
				gsl_matrix_set_identity(columnarFgIncrement);
				gsl_matrix_set_identity(peripodialFgIncrement);
				//if (!(*itElement)->isMutated){
					if (currGF->applyToBasalECM || currGF->applyToLateralECM){
						if (currGF->applyToBasalECM){
							if ((*itElement)->isECMMimicing && (*itElement)->tissuePlacement == 0){
								(*itElement)->calculateFgFromRates(dt, growthscale[0],growthscale[1],growthscale[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
							}
						}
						if (currGF->applyToLateralECM){
							if ((*itElement)->isECMMimimcingAtCircumference && !(*itElement)->tissuePlacement == 0){ //do not grow the basal element twice
								(*itElement)->calculateFgFromRates(dt, growthscale[0],growthscale[1],growthscale[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
							}
						}
					}
					if (!thereIsExplicitECM || !(*itElement)->isECMMimicing ){
						//There is either no explicit ECM definition, or the element is not ECM mimicing.
						//If there is explicit ECM, the basal elements should not grow, all others should proceed as usual
						//If there is no explicit ecm, then all should proceed as usual.
						if (currGF->applyToColumnarLayer){
							(*itElement)->calculateFgFromRates(dt, growthscale[0],growthscale[1],growthscale[2], currGF->getShearAngleRotationMatrix(), columnarFgIncrement, 0, currGF->zMin, currGF->zMax);
						}
						if (currGF->applyToPeripodialMembrane){
							(*itElement)->calculateFgFromRates(dt, growthscale[0],growthscale[1],growthscale[2], currGF->getShearAngleRotationMatrix(), peripodialFgIncrement, 1, currGF->zMin, currGF->zMax);
						}
					}
					(*itElement)->updateGrowthIncrement(columnarFgIncrement,peripodialFgIncrement);
				//}
			}
			delete[] Elementcentre;
		}
		gsl_matrix_free(columnarFgIncrement);
		gsl_matrix_free(peripodialFgIncrement);
		delete[] maxValues;
	}
}

void Simulation::setUpECMMimicingElements(){
	for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
		(*itNode)->allOwnersECMMimicing = true;
		//will be cleared in the loop for assigning the ECM mimicking elements
	}

	for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ( (*itElement)->tissuePlacement == 0 ){
			if((*itElement)->tissueType == 0){
				//basal elements are mimicing the ECM:
				(*itElement)->setECMMimicing(true);
			}
		}
		else if ( (*itElement)->tissuePlacement == 1 ){
			if((*itElement)->tissueType == 1){
				//apical elements are mimicing the ECM on peripodial:
				(*itElement)->setECMMimicing(true);
			}
		}
		else if ( (*itElement)->tissuePlacement == 2 && (*itElement)->spansWholeTissue ){
			//elements that  are called mid -line, as they
			//span the whole tissue, are treated as ECM mimicing
			(*itElement)->setECMMimicing(true);
		}
		if (AddPeripodialMembrane){
			/*double peripodialSideConnectionThickness =PeripodialLateralThicnessScale*TissueHeight; //in microns
			double ECMThicknessSquare = peripodialSideConnectionThickness*peripodialSideConnectionThickness;

			//element centre

			//If there is no peripodial membrane, I would like to have a circumferential region as ECM as well:
			//If I am adding peripodial membrane in the simulation (not the sophisticated round shape but simple,
			//rectangular addition, then I also want to assign the curcumference. Here, I need to ensure the
			//actin mimicking element attributes are cleared once I assign an element to be ECM,
			//bool isElementAtCircumference = (*itElement)->areanyOfMyNodesAtCircumference(Nodes);
			//if (isElementAtCircumference){
			//	(*itElement)->setECMMimicing(true);
			//	(*itElement)->isECMMimimcingAtCircumference = true;
			//	if ((*itElement)->isActinMimicing){
			//		(*itElement)->setActinMimicing(false);
			//	}
			//}
			double *c = (*itElement)->getCentre();
			for(vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
				if ( (*itNode)->atCircumference ){
					//The node is at circumference. Calculate distance of the element to this node.
					//If The distance is lower than the ECM thickness, then assign the element as
					//ECM mimicking.
					double dx = c[0]-(*itNode)->Position[0];
					double dy = c[1]-(*itNode)->Position[1];
					double d = dx*dx + dy*dy;
					if (d<ECMThicknessSquare){
						(*itElement)->setECMMimicing(true);
						(*itElement)->isECMMimimcingAtCircumference = true;
						if ((*itElement)->isActinMimicing){
							(*itElement)->setActinMimicing(false);
						}
						break;
					}
				}
			}
			delete[] c;*/
			if ((*itElement)->isECMMimimcingAtCircumference){
				//I have already assigned this in setting up peripodial membrane
				(*itElement)->setECMMimicing(true);
				(*itElement)->isECMMimimcingAtCircumference = true;
				if ((*itElement)->isActinMimicing){
					(*itElement)->setActinMimicing(false);
				}
			}
		}
		if (!(*itElement)->isECMMimicing){
			//the element is not ecm mimicking, its nodes will not be ecm mimicking:
			int n = (*itElement)->getNodeNumber();
			int* nodeIds;
			nodeIds = (*itElement)->getNodeIds();
			for (int i =0; i<n; ++i){
				Nodes[nodeIds[i]]->allOwnersECMMimicing = false;
				//cout<<"Node ["<<nodeIds[i]<<" is not ECM mimicking"<<endl;
			}
		}
	}
	assigneElementsAtTheBorderOfECM();
}


void Simulation::assigneElementsAtTheBorderOfECM(){
	vector<int> nodeListToCheck;
	for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ( (*itElement)->isECMMimicing  && !(*itElement)->isECMMimimcingAtCircumference){ //I do not want to track the elements at the side borders
			//The element is ECM mimicking. All elements that share nodes with this element, but are
			//not ECMMimicking themselves, should be treated as basal (or whole tissue spanning) elements:
			int n = (*itElement)->getNodeNumber();
			for (int i=0; i<n; ++i){
				bool alreadyRecorded = false;
				int nVector = nodeListToCheck.size();
				for (int j=0 ;j< nVector; ++j){
					if (nodeListToCheck[j] ==(*itElement)->NodeIds[i] ){
						alreadyRecorded = true;
						break;
					}
				}
				if (!alreadyRecorded){
					nodeListToCheck.push_back((*itElement)->NodeIds[i]);
				}
			}
		}
	}
	int nVector = nodeListToCheck.size();
	for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->isECMMimicing && !(*itElement)->atBasalBorderOfECM){
			for (int j=0 ;j< nVector; ++j){
				if ((*itElement)->DoesPointBelogToMe(nodeListToCheck[j])){
					(*itElement)->atBasalBorderOfECM = true;
					break;
				}
			}
		}
	}
}

void Simulation::assigneElementsAtTheBorderOfActin(){
	vector<int> nodeListToCheck;
	for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ( (*itElement)->isActinMimicing  ){
			//The element is actin mimicking. All elements that share nodes with this element, but are
			//not ActinMimicing themselves, and are not ECM mimicing, should be treated as basal (or whole tissue spanning) elements:
			int n = (*itElement)->getNodeNumber();
			for (int i=0; i<n; ++i){
				bool alreadyRecorded = false;
				int nVector = nodeListToCheck.size();
				for (int j=0 ;j< nVector; ++j){
					if (nodeListToCheck[j] ==(*itElement)->NodeIds[i] ){
						alreadyRecorded = true;
						break;
					}
				}
				if (!alreadyRecorded){
					nodeListToCheck.push_back((*itElement)->NodeIds[i]);
				}
			}
		}
	}
	int nVector = nodeListToCheck.size();
	for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if (!(*itElement)->isActinMimicing && !(*itElement)->isECMMimicing && !(*itElement)->atBasalBorderOfECM){
			for (int j=0 ;j< nVector; ++j){
				if ((*itElement)->DoesPointBelogToMe(nodeListToCheck[j])){
					(*itElement)->atApicalBorderOfActin = true;
					break;
				}
			}
		}
	}
}

void Simulation::setUpActinMimicingElements(){
	for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ( (*itElement)->tissuePlacement == 1 ){
			//apical elements are mimicing actin:
			if (!(*itElement)->isECMMimimcingAtCircumference){
				(*itElement)->setActinMimicing(true);
			}
		}
		if ( (*itElement)->tissuePlacement == 2 && (*itElement)->spansWholeTissue ){
			//elemetns that  are called mid -line, as they
			//span the whole tissue, are treated as Actin mimicing
			if (!(*itElement)->isECMMimimcingAtCircumference){
				(*itElement)->setActinMimicing(true);
			}
		}
	}
	assigneElementsAtTheBorderOfActin();
}

void Simulation::calculateGrowthGridBased(GrowthFunctionBase* currGF){
	int nGridX = currGF->getGridX();
	int nGridY = currGF->getGridY();
	if(currSimTimeSec >= currGF->initTime && currSimTimeSec < currGF->endTime ){
		gsl_matrix* columnarFgIncrement = gsl_matrix_calloc(3,3);
		gsl_matrix* peripodialFgIncrement = gsl_matrix_calloc(3,3);
		for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			gsl_matrix_set_identity(columnarFgIncrement);
			gsl_matrix_set_identity(peripodialFgIncrement);

			int IndexX = 0.0, IndexY = 0.0;
			double FracX = 1.0,  FracY = 1.0;
			//(*itElement)->getTissuePositionWeigths(columnarnessWeight, peripodialnessWeight);
			if (GridGrowthsPinnedOnInitialMesh){
				(*itElement)->getInitialRelativePositionInTissueInGridIndex(nGridX, nGridY, IndexX, IndexY, FracX, FracY);
			}
			else{
				(*itElement)->getRelativePositionInTissueInGridIndex(nGridX, nGridY, IndexX, IndexY, FracX, FracY);
			}
			//if (!(*itElement)->isMutated){
				if (currGF->applyToBasalECM || currGF->applyToLateralECM){
					if (currGF->applyToBasalECM){
						if ((*itElement)->isECMMimicing && (*itElement)->tissuePlacement == 0){
							(*itElement)->calculateFgFromGridCorners(gridGrowthsInterpolationType, dt, currGF, columnarFgIncrement, 0, IndexX,  IndexY, FracX, FracY); 	//sourceTissue is 0 for columnar Layer
						}
					}
					if (currGF->applyToLateralECM){
						if ((*itElement)->isECMMimimcingAtCircumference && !(*itElement)->tissuePlacement == 0){ //do not grow the basal element twice
							(*itElement)->calculateFgFromGridCorners(gridGrowthsInterpolationType, dt, currGF, columnarFgIncrement, 0, IndexX,  IndexY, FracX, FracY); 	//sourceTissue is 0 for columnar Layer
						}
					}
				}
				if (!thereIsExplicitECM || !(*itElement)->isECMMimicing ){
					//There is either no explicit ECM definition, or the element is not ECM mimicing.
					//If there is explicit ECM, the basal elements should not grow, all others should proceed as usual
					//If there is no explicit ecm, then all should proceed as usual.
					if (currGF->applyToColumnarLayer){
						(*itElement)->calculateFgFromGridCorners(gridGrowthsInterpolationType, dt, currGF, columnarFgIncrement, 0, IndexX,  IndexY, FracX, FracY); 	//sourceTissue is 0 for columnar Layer
					}
					if (currGF->applyToPeripodialMembrane){
						(*itElement)->calculateFgFromGridCorners(gridGrowthsInterpolationType, dt, currGF, peripodialFgIncrement, 1, IndexX,  IndexY, FracX, FracY); 	//sourceTissue is 1 for peripodial membrane
					}
				}
				(*itElement)->updateGrowthIncrement(columnarFgIncrement,peripodialFgIncrement);
			//}
		}
		gsl_matrix_free(columnarFgIncrement);
		gsl_matrix_free(peripodialFgIncrement);
	}
}

void Simulation::TissueAxisPositionDisplay(){
	cerr<<"DV border: "<<endl;
	for (int i=0;i<nNodes;++i){
		double x= Nodes[i]->Position[0];
		if (x < 0.2 && x > -0.2 ){
			cout<<Nodes[i]->Position[0]<<" "<<Nodes[i]->Position[1]<<" "<<Nodes[i]->Position[2]<<endl;
		}
	}
	cerr<<"AP border: "<<endl;
	for (int i=0;i<nNodes;++i){
		double y= Nodes[i]->Position[1];
		if (y < 0.2 && y > -0.2 ){
			cout<<Nodes[i]->Position[0]<<" "<<Nodes[i]->Position[1]<<" "<<Nodes[i]->Position[2]<<endl;
		}
	}
}

void Simulation::coordinateDisplay(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		int type =(*itElement)-> getShapeType();
		if(type ==1){
			for(int j=0;j<3;++j){
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[0]<<" ";
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[1]<<" ";
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[2]<<" ";
			}
		}
	}
	cout<<endl;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		int type =(*itElement)-> getShapeType();
		if(type ==1){
			for(int j=3;j<6;++j){
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[0]<<" ";
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[1]<<" ";
				cout<<Nodes[(*itElement)->NodeIds[j]]->Position[2]<<" ";
			}
		}
	}
	cout<<endl;
}

void Simulation::setStretch(){
	cout<<"setting the stretcher"<<endl;
	recordForcesOnFixedNodes = true;
	for (int i=0;i<3;i++){
		leftClampForces[i] = 0.0;
		rightClampForces[i] = 0.0;
	}
	vector <int> clampedNodeIds;
	double distance = 0;


	if (DVClamp){
		distance = fabs(Nodes[ventralTipIndex]->Position[0] - Nodes[dorsalTipIndex]->Position[0]);
		cerr<<"Total DV distance: "<<distance<<" ";
		distanceIndex = 0; //if the clamp is on DV axis, then the direction of interest is x, index is 0;
	}
	else{
		distance = fabs(Nodes[anteriorTipIndex]->Position[1] - Nodes[posteriorTipIndex]->Position[1]);
		cerr<<"Total AP distance: "<<distance<<" ";
		distanceIndex = 1; //if the clamp is on AP axis, then the direction of interest is y, index is 1.
	}
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->Position[distanceIndex]> StretchMax || Nodes[i]->Position[distanceIndex] < StretchMin){
			Nodes[i]->FixedPos[0]=1;
			Nodes[i]->FixedPos[1]=1;
			Nodes[i]->FixedPos[2]=1;
			clampedNodeIds.push_back(Nodes[i]->Id);
		}
	}
	setUpClampBorders(clampedNodeIds);
	//the distance that is to be moved:
	distance *= StretchStrain;
	cerr<<"the distance that is to be moved: "<<distance<<" ";
	//the time steps that the stretch operation should take place in:
	double StretchTimeSteps = (StretchEndTime - StretchInitialTime)/dt;
	cerr<<"stretchTimeSteps: "<<StretchTimeSteps<<" ";
	StretchDistanceStep = 0.5* (distance / StretchTimeSteps);
	cerr<<"StretchDistanceStep: "<<StretchDistanceStep<<endl;
}

void Simulation::setUpClampBorders(vector<int>& clampedNodeIds){
	int* nodeIds;
	int n = clampedNodeIds.size();
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		nodeIds = (*itElement)->getNodeIds();
		bool hasClampedNode = false;
		bool hasNonClampedNode = false;
		bool leftHandSide = false;
		vector <int> clampedBorderNodes;
		int nNodes = (*itElement)->getNodeNumber();
		for (int j=0; j< nNodes; j++){
			bool lastNodeWasClamped = false;
			for (int k=0; k<n; k++){
				if (nodeIds[j] == clampedNodeIds[k]){
					hasClampedNode = true;
					lastNodeWasClamped = true;
					//check if node is recorded:
					bool alreadyRecorded = false;
					int nLeftClamp = leftClampBorder.size();
					for (int m=0; m<nLeftClamp; m++){
						if (leftClampBorder[m] == nodeIds[j]){
							alreadyRecorded = true;
							break;
						}
					}
					if (!alreadyRecorded){
						int nRightClamp = rightClampBorder.size();
						for (int m=0; m<nRightClamp; m++){
							if (rightClampBorder[m] == nodeIds[j]){
								alreadyRecorded = true;
								break;
							}
						}
					}
					if (!alreadyRecorded){
						clampedBorderNodes.push_back(nodeIds[j]);
					}
					if (Nodes[nodeIds[j]]->Position[distanceIndex] < StretchMin){
						leftHandSide = true;
					}
					break;
				}
			}
			if (lastNodeWasClamped == false){
				hasNonClampedNode = true;
			}
		}
		if (hasClampedNode && hasNonClampedNode){
			cout<<" Element "<<(*itElement)->Id<<" is at the border"<<endl;
			int nClamp = clampedBorderNodes.size();
			if(leftHandSide){
				cout<<"Element is on left hand side"<<endl;
				for (int k=0; k<nClamp; k++){
					leftClampBorder.push_back(clampedBorderNodes[k]);
				}
			}
			else {
				cout<<"Element is on right hand side"<<endl;
				for (int k=0; k<nClamp; k++){
					rightClampBorder.push_back(clampedBorderNodes[k]);
				}
			}
		}
	}int nLeftClamp = leftClampBorder.size();
	for (int k=0; k<nLeftClamp; k++){
		cout<<"left clamp border nodes: "<<leftClampBorder[k]<<endl;
	}
	int nRightClamp = rightClampBorder.size();
	for (int k=0; k<nRightClamp; k++){
		cout<<"right clamp border nodes: "<<rightClampBorder[k]<<endl;
	}
}

void Simulation::recordForcesOnClampBorders(){
	if (recordForcesOnFixedNodes){
		for (int i=0;i<3;i++){
			leftClampForces[i] = 0.0;
			rightClampForces[i] = 0.0;
		}
		int nLeftClamp = leftClampBorder.size();
		for (int k=0; k<nLeftClamp; k++){
			for (int j=0; j<3; ++j){
				leftClampForces[j] += FixedNodeForces[leftClampBorder[k]][j];
			}
		}
		int nRightClamp = rightClampBorder.size();
		for (int k=0; k<nRightClamp; k++){
			for (int j=0; j<3; ++j){
				rightClampForces[j] += FixedNodeForces[rightClampBorder[k]][j];
			}
		}
	}
	outputFile<<"Forces on clamps lhs: "<<leftClampForces[0]<<" "<<leftClampForces[1]<<" "<<leftClampForces[2]<<" rhs: "<<rightClampForces[0]<<" "<<rightClampForces[1]<<" "<<rightClampForces[2]<<endl;
}

void Simulation::moveAFMBead(){
	if (beadPos[2] != -1000){
		beadPos[2] += 0.05;
		cout<<"new bead pos:" <<beadPos[2]<<endl;
	}
}

void Simulation::moveClampedNodesForStretcher(){
	if (currSimTimeSec>=StretchInitialTime && currSimTimeSec<StretchEndTime){
		for (int i=0; i<nNodes; ++i){
			if (Nodes[i]->Position[distanceIndex]> StretchMax){
				Nodes[i]->Position[distanceIndex] += StretchDistanceStep;
			}
			else if( Nodes[i]->Position[distanceIndex] < StretchMin ){
				Nodes[i]->Position[distanceIndex] -= StretchDistanceStep;
			}
		}
		vector<ShapeBase*>::iterator itElement;
		for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			(*itElement)->updatePositions(Nodes);
        }
	}
}

void Simulation::setupPipetteExperiment(){
	pipetteInnerRadiusSq = pipetteInnerRadius*pipetteInnerRadius;
	effectLimitsInZ[0] = pipetteCentre[2] - pipetteDepth;
	effectLimitsInZ[1] = pipetteCentre[2] + pipetteDepth;
	//Now I set up the upper and lower boundaries to cover both apical and basal suction, but I need to extend the boundary
	//inside the tube, such that the nodes will not stop being suced in after they go inside a certain length into the tube
	//This means increasing the upper boundary for apical suction, and reducing the lower boundary for basal suction:
	if(ApicalSuction){
		effectLimitsInZ[1] += 1000;
	}
	else{
		effectLimitsInZ[0] -= 1000;
	}
	//cout<<"set the system, fixing nodes:"<<endl;
    //Now I am sticking the other side of the tissue to a surface
	if (TissueStuckOnGlassDuringPipetteAspiration){
		if (ApicalSuction){
			for (int i=0; i<nNodes; ++i){
				//fix basal nodes of columnar layer:
				if(Nodes[i]->tissuePlacement == 0 && Nodes[i]->tissueType == 0) {
					fixAllD(i, false); //this is fixing with adhesives, should be a hard fix at all times
				}
			}
		}
		else{
			for (int i=0; i<nNodes; ++i){
				//fix apical nodes and all peripodial membrane nodes:
				if(Nodes[i]->tissuePlacement == 1 || Nodes[i]->tissueType == 1) {
					fixAllD(i, false); //this is fixing with adhesives, should be a hard fix at all times
				}
			}
		}
	}
}

void Simulation::addMyosinForces(gsl_matrix* gExt){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
	    int* nodeIds = (*itElement)->getNodeIds();
	    int nNodes= (*itElement)->getNodeNumber();
	    for (int j=0; j<nNodes; ++j){
			double Fx = (*itElement)->MyoForce[j][0];
			double Fy = (*itElement)->MyoForce[j][1];
			double Fz = (*itElement)->MyoForce[j][2];
			int indice = nodeIds[j]*3;
			Fx += gsl_matrix_get(gExt,indice,0);
			gsl_matrix_set(gExt,indice,0,Fx);
			Fy += gsl_matrix_get(gExt,indice+1,0);
			gsl_matrix_set(gExt,indice+1,0,Fy);
			Fz += gsl_matrix_get(gExt,indice+2,0);
			gsl_matrix_set(gExt,indice+2,0,Fz);
	    }
	}
}

void Simulation::addPipetteForces(gsl_matrix* gExt){
    //cout<<"in add pipette forces, pipette pos: "<<pipetteCentre[0]<<" "<<pipetteCentre[1]<<endl;
	for (int i=0; i<nPipetteSuctionSteps; ++i){
		if (currSimTimeSec >= pipetteSuctionTimes[i]){
			SuctionPressure[2] = pipetteSuctionPressures[i];
		}
		else{
			break;
		}
	}
	int dim = 3;
	for (int i=0; i<nNodes; ++i){
        //cout<<"Node "<<i<<" z pos: "<<Nodes[i]->Position[2]<<" effectLimitsInZ: "<<effectLimitsInZ[0]<<" "<<effectLimitsInZ[1]<<endl;
		//The node is apical and the suction is apical OR the suction is basal and the node is basal!
		//Should not fix midline nodes just because they are coming close to the pipette border
		if( ( ApicalSuction && Nodes[i]->tissuePlacement == 1) || ( !ApicalSuction && Nodes[i]->tissuePlacement == 0)){
			if (Nodes[i]->Position[2]> effectLimitsInZ[0] &&  Nodes[i]->Position[2]< effectLimitsInZ[1]){
				//cout<<"Node "<<i<<" is within z range"<<endl;
				double dx = pipetteCentre[0] - Nodes[i]->Position[0];
				double dy = pipetteCentre[1] - Nodes[i]->Position[1];
				//cout<<"dx: "<<dx<<" dy: "<<dy<<" dx*dx+dy*dy: "<<dx*dx+dy*dy<<endl;
				double d2 = dx*dx+dy*dy ;
				if (d2 < pipetteInnerRadiusSq){
					//cout<<"Node "<<i<<" is within planar range"<<endl;
					double multiplier = 1.0;
					bool scalePressure = false;
					if(scalePressure){ multiplier = (1 - d2/pipetteInnerRadiusSq);}
					double LocalPressure[3] = { multiplier*SuctionPressure[0], multiplier*SuctionPressure[1],multiplier*SuctionPressure[2]};
					//cout<<"Node: "<<i<<" being sucked into pipette, force: "<<SuctionPressure[0]*Nodes[i]->surface<<" "<<SuctionPressure[1]*Nodes[i]->surface<<" "<<SuctionPressure[2]*Nodes[i]->surface<<endl;
					//node is within suciton range
					if(!Nodes[i]->FixedPos[0]){
						double value = gsl_matrix_get(gExt,i*dim,0);
						value +=LocalPressure[0]*Nodes[i]->zProjectedArea;
						gsl_matrix_set(gExt,i*dim,0,value);
					}
					if(!Nodes[i]->FixedPos[1]){
						double value = gsl_matrix_get(gExt,i*dim+1,0);
						value +=LocalPressure[1]*Nodes[i]->zProjectedArea;
						gsl_matrix_set(gExt,i*dim+1,0,value);
					}
					if(!Nodes[i]->FixedPos[2]){
						double value = gsl_matrix_get(gExt,i*dim+2,0);
						value +=LocalPressure[2]*Nodes[i]->zProjectedArea;
						gsl_matrix_set(gExt,i*dim+2,0,value);
					}
					//Elements[0]->displayMatrix(gExt,"gExtInsidePipetteFunc");
				}
			}
		}
	}
}

void Simulation::packToPipetteWall(){
	//Fixing the z-height of nodes near the pipette wall:
	double pipLim[2] = {0.95*pipetteInnerRadius, 1.2*(pipetteInnerRadius+pipetteThickness)};
	double pipLim2[2] = {pipLim[0] *pipLim[0], pipLim[1]* pipLim[1]};
	bool croppedTipPipette = false;
	int nZfix = TransientZFixListForPipette.size();
	for (int zfixiter =0; zfixiter <nZfix; zfixiter++){
		Nodes[TransientZFixListForPipette[zfixiter]]->FixedPos[2] = 0;
	}
	TransientZFixListForPipette.clear();
	for (int currId=0; currId<nNodes;currId++){
		if (Nodes[currId]->FixedPos[2] == 0){
			//The node is apical and the suction is apical OR the suction is basal and the node is basal!
			//Should not fix midline nodes just because they are coming close to the pipette border
			if( ( ApicalSuction && Nodes[currId]->tissuePlacement == 1) || ( !ApicalSuction && Nodes[currId]->tissuePlacement == 0)){
				bool freezeNodeZ = false;
				//the node is not already fixed in z
				//I will check if it is in the range of the tip of pipette:
				double zLimits[2] = {-2.0,2.0}; //the range in z to freeze nodes
				double v0[3] = {Nodes[currId]->Position[0]-pipetteCentre[0], Nodes[currId]->Position[1]-pipetteCentre[1], Nodes[currId]->Position[2]-pipetteCentre[2]};
				if (v0[2] < zLimits[1] && v0[2]> zLimits[0]){
					//node is in Z range, is it in x-y range?
					double xydist2 = v0[0]*v0[0]+v0[1]*v0[1];
					if (xydist2>pipLim2[0] && xydist2<pipLim2[1]){
						//node is within x-y range,
						//Now I want to know if I am using cropped pipette tip, and if so,
						//I will check z limit again:
						if (croppedTipPipette){
							double R = pow(xydist2,0.5);
							double CroppedZAtCurrentR =  zLimits[0] + ((zLimits[1] - zLimits[0])/(pipLim[1] - pipLim[0]) * (R - pipLim[0] ));
							if (v0[2] < CroppedZAtCurrentR){
								freezeNodeZ = true;
							}
						}
						else{
							freezeNodeZ = true;
						}
						if (freezeNodeZ){
							Nodes[currId]->FixedPos[2] = 1;
							TransientZFixListForPipette.push_back(currId);
						}
					}
				}
			}
		}
	}
}

void Simulation::laserAblateTissueType(int ablationType){
	vector <int> AblatedElements;
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->tissueType == ablationType && !Nodes[i]->hasLateralElementOwner){
			//I want to ablate a whole tisse type, such as ablating the disc proper at the beginning of simulation.
			//BUT, I want the nodes that are owned by the linker nodes to stay intact. Otherwise, I will loose part of the unintended tissue.
			AblatedNodes.push_back(i);
		}
		else if (Nodes[i]->Position[2] > 14.25){//I do not want any node above 14.25;
			AblatedNodes.push_back(i);
		}
	}

	int nAN = AblatedNodes.size();
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if(!(*itElement)->IsAblated){
			for (int j =0; j<nAN; ++j){
				bool IsAblatedNow = (*itElement)->DoesPointBelogToMe(AblatedNodes[j]);
				if (IsAblatedNow){
					(*itElement)->removeMassFromNodes(Nodes);
					(*itElement)->IsAblated = true;
					//cerr<<"Ablating element:" <<Elements[i]->Id<<endl;
					break;
				}
			}
		}
	}
	//some nodes are left with zero mass, which will cause problems in later calculations:
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->mass <=0){
			Nodes[i]->mass = 0.1;
		}
	}
}

void Simulation::laserAblate(double OriginX, double OriginY, double Radius){
	vector <int> AblatedNodes;
	vector <int> AblatedElements;
	double thres2 = Radius*Radius;
	for (int i=0; i<nNodes; ++i){
		double dx = Nodes[i]->Position[0]- OriginX;
		double dy = Nodes[i]->Position[1]- OriginY;
		double d2 = dx *dx + dy*dy;
		if (d2 < thres2){
			AblatedNodes.push_back(i);
		}
		//if (i != 0 && i != 1 && i != 22 && i != 88 && i != 66 && i != 67  ){
			//AblatedNodes.push_back(i);
		//}
	}

	int nAN = AblatedNodes.size();
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if(!(*itElement)->IsAblated){
			for (int j =0; j<nAN; ++j){
				bool IsAblatedNow = (*itElement)->DoesPointBelogToMe(AblatedNodes[j]);
				if (IsAblatedNow){
					(*itElement)->removeMassFromNodes(Nodes);
					(*itElement)->IsAblated = true;
					cerr<<"Ablating element:" <<(*itElement)->Id<<endl;
					break;
				}
			}
		}
	}
	//some nodes are ledt with zero mass, which will cause problems in later calculations:
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->mass <=0){
			Nodes[i]->mass = 0.1;
		}
	}
}

void Simulation::updateElementVolumesAndTissuePlacements(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		cout<<"updating element: "<<(*itElement)->Id<<endl;
		(*itElement)->updateElementVolumesAndTissuePlacementsForSave(Nodes);
	}
}

void Simulation::updateElasticPropertiesForAllNodes(){
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updateElasticProperties();
	}
}

void Simulation::clearNodeMassLists(){
	for (int i=0 ;i<nNodes;++i){
		Nodes[i]->connectedElementIds.size();
		Nodes[i]->connectedElementIds.clear();
		Nodes[i]->connectedElementWeights.clear();
		Nodes[i]->mass=0.0;
		//Nodes[i]->surface=0.0;
	}
}

void Simulation::clearLaserAblatedSites(){
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if ((*itElement)->IsAblated){
			(*itElement)->removeMassFromNodes(Nodes);
		}
	}
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->mass <=0){
			Nodes[i]->mass = 0.1;
		}
	}
}

void Simulation::setupYsymmetricity(){
	double yLimPos = 0.1;
	double yLimNeg = (-1.0)*yLimPos;
	vector <int> AblatedNodes;
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->Position[1]< yLimPos){
			if (Nodes[i]->Position[1] > yLimNeg){
				//symmetricYBoundaryNodes.push_back(Nodes[i]);
				Nodes[i]->atSymmetricityBorder = true;
				fixY(Nodes[i],false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
			}
			else{
				AblatedNodes.push_back(i);
				fixAllD(Nodes[i], false); //this is fixing for ablated nodes, no need for calculations
				//setSymmetricNode(Nodes[i],yLimPos);
			}
		}
	}
	int nAN = AblatedNodes.size();
	//fix the position of all ablated nodes for effective Newton Raphson calculation:
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if(!(*itElement)->IsAblated){
			for (int j =0; j<nAN; ++j){
				bool IsAblatedNow = (*itElement)->DoesPointBelogToMe(AblatedNodes[j]);
				if (IsAblatedNow){
					(*itElement)->removeMassFromNodes(Nodes);
					(*itElement)->IsAblated = true;
					break;
				}
			}
		}
	}
}

void Simulation::setupXsymmetricity(){
	double xLimPos = 0.1;
	double xLimNeg = (-1.0)*xLimPos;
	vector <int> AblatedNodes;
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->Position[0]< xLimPos){
			if (Nodes[i]->Position[0] > xLimNeg){
				//symmetricXBoundaryNodes.push_back(Nodes[i]);
				Nodes[i]->atSymmetricityBorder = true;
				fixX(Nodes[i],false); //this is for symmetricity, the fixing has to be hard fixing, not with external viscosity under any condition
			}
			else{
				AblatedNodes.push_back(i);
				fixAllD(Nodes[i], false); //this is fixing for ablated nodes, no need for calculations
				//setSymmetricNode(Nodes[i],yLimPos);
			}
		}
	}
	int nAN = AblatedNodes.size();
	//fix the position of all ablated nodes for effective Newton Raphson calculation:
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		if(!(*itElement)->IsAblated){
			for (int j =0; j<nAN; ++j){
				bool IsAblatedNow = (*itElement)->DoesPointBelogToMe(AblatedNodes[j]);
				if (IsAblatedNow){
					(*itElement)->removeMassFromNodes(Nodes);
					(*itElement)->IsAblated = true;
					break;
				}
			}
		}
	}
}

void Simulation::ablateSpcific(){
	vector <int> AblatedNodes;
	for (int i=0; i<nNodes; ++i){
		fixAllD(Nodes[i],  false); //this is fixing for ablated nodes, no need for calculations);
	}
	//int nAN = AblatedNodes.size();
	//fix the position of all ablated nodes for effective Newton Raphson calculation:
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		double* c = new double[3];
		c = (*itElement)->getCentre();
		//if ((*itElement)->Id != 38 && (*itElement)->Id != 39 /*&& (*itElement)->Id != 30 && (*itElement)->Id != 41*/){
		if (c[0]<20 || c[1] >10 || c[1]<-10){
			(*itElement)->removeMassFromNodes(Nodes);
			(*itElement)->IsAblated = true;
		}
		else{
			for (int i=1; i<nNodes; ++i){
				if ((*itElement)->DoesPointBelogToMe(i)){
					Nodes[i]->FixedPos[0]= false;
					Nodes[i]->FixedPos[1]= false;
					Nodes[i]->FixedPos[2]= false;
				}
			}
		}
		delete[] c;
	}
	Nodes[1]->FixedPos[2]= true;
	Nodes[18]->FixedPos[1]= true;
	Nodes[18]->FixedPos[2]= true;
}

void Simulation::pokeElement(int elementId, double dx, double dy, double dz){
	cout<<" poking element: "<<elementId<<endl;
	int* nodeIds = Elements[elementId]->getNodeIds();
	int nNodes= Elements[elementId]->getNodeNumber();
	for (int j=0; j<nNodes; ++j){
		Nodes[nodeIds[j]]->Position[0] += dx;
		Nodes[nodeIds[j]]->Position[1] += dy;
		Nodes[nodeIds[j]]->Position[2] += dz;
	}
	vector<ShapeBase*>::iterator itElement;
	for(itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->updatePositions(Nodes);
    }
}

void Simulation::writeMeshRemovingAblatedRegions(){
	cout<<"writing non-ablated mesh"<<endl;
	int nonAblatedNodeMap[( const int ) nNodes];
	int nNonAblatedNode = 0;
	for (int i=0; i<nNodes; ++i){
		nonAblatedNodeMap[i] = -10;
		if (Nodes[i]->mass > 0){
			nonAblatedNodeMap[i] = nNonAblatedNode;
			nNonAblatedNode++;
		}
	}
	cout<<"got non-ablated node number and map "<<endl;
	int nNonAblatedElements = 0;
	for (int i=0; i<nElements; ++i){
		if (!Elements[i]->IsAblated){
			nNonAblatedElements++;
		}
	}
	cout<<"opening file"<<endl;
	string meshSaveString = saveDirectory +"/MeshFromNonAblated.mesh";
	const char* name_meshSaveString = meshSaveString.c_str();;
	ofstream file;
	file.open(name_meshSaveString, ofstream::out);
	file<<nNonAblatedNode;
	file<<endl;
	for (int i=0; i<nNodes; ++i){
		if (Nodes[i]->mass > 0){
			file << Nodes[i]->Position[0];
			file<<" \t";
			file << Nodes[i]->Position[1];
			file<<" \t";
			file << Nodes[i]->Position[2];
			file<<" \t";
			file << Nodes[i]->tissuePlacement;
			file<<" \t";
			file << Nodes[i]->tissueType;
			file<<" \t";
			file << Nodes[i]->atCircumference;
			file << endl;
		}
	}
	file<<nNonAblatedElements;
	file<<endl;
	for (int i=0; i<nElements; ++i){
		if (!Elements[i]->IsAblated){
			file<< Elements[i]->getShapeType();
			file<<" \t";
			const int n = Elements[i]->getNodeNumber();
			int* NodeIds;
			NodeIds = new int[n];
			NodeIds = Elements[i]->getNodeIds();
			for (int j=0; j<n; ++j){
				file<< nonAblatedNodeMap[NodeIds[j]];
				file<<" \t";
			}
			int dim  = Elements[i]->getDim();
			double** refPos = Elements[i]->getReferencePos();
			for (int j = 0; j<6; ++j ){
				for (int k = 0; k<dim; ++k ){
					file.precision(5);file.width(12);
					file<<refPos[j][k];
				}
			}
			file << endl;
		}
	}
	//recording tissue weights
	file <<1<<endl;
	for (int i=0; i<nElements; ++i){
			if (!Elements[i]->IsAblated){
				file <<Elements[i]->getPeripodialness() << endl;
			}
	}
	file.close();
}

void Simulation::checkForVolumeRedistributionInTissue(){
	clearScaleingDueToApikobasalRedistribution();
	for (int i=0; i< nApikobasalVolumeRedistributionFunctions; ++i){
		if (currSimTimeSec >= apikobasalVolumeRedistributionBeginTimeInSec[i] && currSimTimeSec < apikobasalVolumeRedistributionEndTimeInSec[i]){
			bool thisFunctionShrinksApical = apikobasalVolumeRedistributionFunctionShrinksApical[i];
			double scale =apikobasalVolumeRedistributionScales[i];
			const int maxThreads = omp_get_max_threads();
			omp_set_num_threads(maxThreads);
			#pragma omp parallel for
			for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
				(*itElement)->updateGrowthWillBeScaledDueToApikobasalRedistribution(thisFunctionShrinksApical, scale, apikobasalVolumeRedistributionFunctionEllipseBandIds[i]);
			}
		}
	}
}

void Simulation::clearScaleingDueToApikobasalRedistribution(){
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for
	for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
		(*itElement)->thereIsGrowthRedistribution = false;
		(*itElement)->growthRedistributionScale = 0.0;
	}
}

void Simulation::assignCompartment(){
	double bufferLength = 6.0; //microns;
	double notumPrisms[2] = {2183, 3505};
	double pouchPrisms[3] = {1494,305,368};
	if (Elements.size()<2184){
		for (vector<ShapeBase*>::iterator itEle = Elements.begin(); itEle<Elements.end(); ++itEle){
			(*itEle)->compartmentIdentityFraction = 1.0;
			(*itEle)->compartmentType = -1;
		}
		return;
	}
	//notum slope and constant: m_notum x + b_notum = y
	double* notum_p0 = Elements[notumPrisms[0]]->getCentre();
	double* notum_p1 = Elements[notumPrisms[1]]->getCentre();

	double* pouch_p0 = Elements[pouchPrisms[0]]->getCentre();
	double* pouch_p1 = Elements[pouchPrisms[1]]->getCentre();
	double* pouch_p2 = Elements[pouchPrisms[2]]->getCentre();

	double m_notum = (notum_p0[1]-notum_p1[1])/(notum_p0[0] - notum_p1[0]);
	double b_notum = notum_p0[1] - m_notum * notum_p0[0];
	double m_pouch0 = (pouch_p0[1]-pouch_p1[1])/(pouch_p0[0] - pouch_p1[0]);
	double b_pouch0 = pouch_p0[1] - m_pouch0 * pouch_p0[0];
	double m_pouch1 = (pouch_p1[1]-pouch_p2[1])/(pouch_p1[0] - pouch_p2[0]);
	double b_pouch1 = pouch_p1[1] - m_pouch1 * pouch_p1[0];
	//cout<<" notum0: "<<notum_p0[0]<<" "<<notum_p0[1]<<endl;
	//cout<<" notum1: "<<notum_p1[0]<<" "<<notum_p1[1]<<endl;
	//cout<<" pouch0: "<<pouch_p0[0]<<" "<<pouch_p0[1]<<endl;
	//cout<<" pouch1: "<<pouch_p1[0]<<" "<<pouch_p1[1]<<endl;
	//cout<<" pouch2: "<<pouch_p2[0]<<" "<<pouch_p2[1]<<endl;
	//cout<<" m_notum, b_notum, m_pouch0, b_pouch0, m_pouch1, b_pouch1"<<endl;
	//cout<< m_notum<<" "<< b_notum<<" "<< m_pouch0<<" "<< b_pouch0<<" "<< m_pouch1<<" "<< b_pouch1<<endl;

	for (vector<ShapeBase*>::iterator itEle = Elements.begin(); itEle<Elements.end(); ++itEle){
		if ((*itEle) -> tissueType == 0){ //columnar layer
			(*itEle)->compartmentType = 1; //default is set to hinge
			double* p = (*itEle)->getCentre();
			//check if notum:
			double xOnLine = (p[1]-b_notum)/m_notum;
			double d = p[0]-xOnLine; //correct for hinge
			//cout<<"(*itEle)->Id: "<< (*itEle)->Id<<" notum_p0 "<<notum_p0[0]<<" "<<notum_p0[1]<<" notum_p1: "<<notum_p1[0]<<" "<<notum_p1[1]<<" m & b: "<<m_notum<<" "<<b_notum<<" p "<<p[0]<<" "<<p[1]<<" xOnLine "<<xOnLine<<endl;
			if (xOnLine > p[0]){
				//element is on the notum side!
				(*itEle)->compartmentType = 2; //notum
				d *= -1.0; //corrected distance o reflect the notum side.
				//threshold at 3 microns, this will make a total of 6 microns on each side
			}
			else{
				//check pouch:
				if (p[1] < pouch_p1[1]){
					//check first segment
					xOnLine = (p[1]-b_pouch0)/m_pouch0;
					if (xOnLine < p[0]){
						//on pouch, correct d:
						d = p[0]- xOnLine;
						//element is on the pouch side!
						(*itEle)->compartmentType = 0; //notum
					}
					else{
						//on hinge, is it closer to this boundary?
						double dpouch = xOnLine- p[0];
						if( dpouch < d){
							d = dpouch;
						}
					}

				}
				else{
					//check second segment
					xOnLine = (p[1]-b_pouch1)/m_pouch1;
					if (xOnLine < p[0]){
						//on pouch, correct d:
						d = p[0]- xOnLine;
						//element is on the pouch side!
						(*itEle)->compartmentType = 0; //notum
					}
					else{
						//on hinge, is it closer to this boundary?
						double dpouch = xOnLine- p[0];
						if( dpouch < d){
							d = dpouch;
						}
					}
				}
			}
			if (d>bufferLength){
				(*itEle)->compartmentIdentityFraction = 1.0;
			}
			else{
				(*itEle)->compartmentIdentityFraction = d/bufferLength;
			}
			delete[] p;
		}
		else{
			(*itEle)->compartmentType = -1;
		}
	}
	delete[] notum_p0;
	delete[] notum_p1;
	delete[] pouch_p0;
	delete[] pouch_p1;
	delete[] pouch_p2;
}
