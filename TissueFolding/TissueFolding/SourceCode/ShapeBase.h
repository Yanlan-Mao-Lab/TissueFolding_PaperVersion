#ifndef ShapeBase_H
#define ShapeBase_H


#include <iostream>
#include <ostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>
//this was the working version in linux. It should be working with correct addition of the path to INCLUDEPATH in .pro file
//#include </usr/include/gsl/gsl_matrix.h>
//#include </usr/include/gsl/gsl_linalg.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "Node.h"
#include "ReferenceShapeBase.h"
#include "GrowthFunctionBase.h"
#include "GrowthFunctionTypes.h"


/* TO DO
 *  The variables are defined for prisms only, convert them to be generic:
 *  	exposedApicalSurfaceNodeIds[3];				///< The int array of size 3, listing the node IDs of element that form the exposed apical surface. The IDs are the node IDs on the element (0-5 for prism), not the actual Node::Id.
 *  	exposedBasalSurfaceNodeIds[3];				///< The int array of size 3, listing the node IDs of element that form the exposed basal surface. The IDs are the node IDs on the element (0-5 for prism), not the actual Node::Id.
 *  	exposedLateralAreaApicalSideNodeIds[4];		///< The int array of size 4, listing the node IDs of element that form the lateral surface exposed apically. The IDs are the node IDs on the element (0-5 for prism), not the actual Node::Id.
 *  	exposedLateralAreaBasalSideNodeIds[4];		///< The int array of size 4, listing the node IDs of element that form the lateral surface exposed basally. The IDs are the node IDs on the element (0-5 for prism), not the actual Node::Id.
 * 		The actual numbers 3 & 4 are already stored in variables and used as such (nLateralSurfaceAreaNodeNumber, nSurfaceAreaNodeNumber);
 */
using namespace std;

class ShapeBase{
private:
	void ParentErrorMessage(string functionName);					 	///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string inputs.
	bool ParentErrorMessage(string functionName, bool returnValue); 	///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string and bool inputs, returning bool.
	double ParentErrorMessage(string functionName, double returnValue);	///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string and double inputs, returning double.
	int ParentErrorMessage(string functionName, int returnValue);		///<Error message displayed when a virtual function is called through the ShapeBase(parent), while it should have been called through a child (eg. Prism). For functions taking in string and int inputs, returning int.
protected:
	int 		ShapeType;							///< The integer defining the type of the shape, Prisms shape type = 1;
	int 		nNodes;								///< The number of nodes of the element, it is based on ShapeBase#ShapeType
	int 		nDim;								///< The number of dimensions for the positions of each of the nodes of the element
	int* 		IdentifierColour;					///< The unique identifier colour of the element, this is used for "picking" in the visual interface.
	double* 	GrowthRate;							///< Growth rate recording for display purposes only. The recorded growth rate in x, y, and z  coordinates, does not record shear deformation induced in growth. Recorded in exponential form through time step, converted to rate per hour for display within the visual interface
	gsl_matrix* growthIncrement;					///< The matrix (3,3) representing the incremental growth in current time step. Reset to identity at the beginning of each time step, updated in growth functions, and utilised to update Fg.
	gsl_matrix* plasticDeformationIncrement;		///< The matrix (3,3) representing the incremental plastic deformation (treated as growth) in current time step. Set in plastic deformation calculation at each step, and utilised to update Fg.
	gsl_matrix* shapeChangeIncrement;				///< The matrix (3,3) representing the incremental shape change in current time step. Reset to identity at the beginning of each time step, updated in shape change functions, and utilised to update Fg.
	double 		zRemodellingSoFar;
	double  	columnarGrowthWeight;				///< The fraction defining how close to the columnar layer the element is. 1.0 for columnar layer, 0.0 for peripodial membrane elements, and scaled according to position in the elements surrounding the lumen.
	double  	peripodialGrowthWeight;				///< The fraction defining how close to the peripodial membrane the element is. 0.0 for columnar layer, 1.0 for peripodial membrane elements, and scaled according to position in the elements surrounding the lumen.
	double* 	ShapeChangeRate;					///< Shape change rate of the elements, only orthagonal shape changes are allowed (x, y, z). Shape changes will be scaled to conserve volume, thus three values will not be independent.
    bool    	rotatedGrowth;						///< The boolean stating if the element has rotated from the growth axis, hence the calculated growth requires further rotation to follow tissue axes.
    double* 	relativePosInBoundingBox;  			///< The relative position on x-y plane, within the bounding box of the tissue(x,y).
    double* 	initialRelativePosInBoundingBox; 	///< The relative position on x-y plane, within the bounding box of the tissue(x,y) at the beginning of simulation. This is used when growth rates are pinned to the initial structure of the tissue.
    double 		initialRelativePositionInZ;			///<< The relative position on z-height of tissue, taken not in z direction but in tissue layers, 0 being on the apical surface and 1 being on the basal surface.
    int			numberOfGaussPoints;				///<< The number of Gauss points used in numerical deforamtion calculation.
    double**	gaussPoints;						///<< The array contianing all the Gauss points for element.
    double* 	gaussWeights;					///<< The array for storing the weights of each Gauss point for element.
    gsl_matrix 	**ShapeFuncDerivatives;				///< The array of matrices for shape function derivatives. The array stores a ShapeBase#nDim by ShapeBase#nNodes matrix for each gauss point (there are 3 Gauss points for prisms).
    gsl_matrix 	**ShapeFuncDerStacks;				///< The array of matrices of shape function derivatives in stacked format for ease of matrix operations. The array stores a (ShapeBase#nDim * ShapeBase#nDim) by (ShapeBase#nDim * ShapeBase#nNodes) matrix for each gauss point (there are 3 Gauss points for prisms).
    gsl_matrix 	**InvdXdes;							///< The array stores inverse of the matrix for derivatives of world coordinates with respect to barycentric coordinates (dX / de). The array stores an ShapeBase#nDim by ShapeBase#nDim  matrix for each gauss point (there are 3 Gauss points for prisms).
    double* detdXdes;								///< The array stores the determinants of the matrices for derivatives of world coordinates with respect to barycentric coordinates (dX / de). The array stores a double value for each gauss point (there are 3 Gauss points for prisms).
    gsl_matrix 	**Bmatrices;						///< The array stores the B matrix for the calculation of stiffness matrix, see for ShapeBase#calculateBTforNodalForces calculation. The array stores an ShapeBase#nNodes by (ShapeBase#nDim*ShapeBase#nNodes)  matrix for each Gauss point (there are 3 Gauss points for prisms).
    gsl_matrix 	**FeMatrices;						///< The array stores the elastic part of the deformation matrix. The array stores an ShapeBase#nDim by ShapeBase#nDim  matrix for each Gauss point (there are 6 Gauss points for prisms).
    gsl_matrix 	**invJShapeFuncDerStack;			///< The array stores the shape function derivatives multiplied by the inverse Jacobian stack, for each Gauss point. See ShapeBase#calculateBTforNodalForces for calculation.
    gsl_matrix 	**invJShapeFuncDerStackwithFe;		///< See ShapeBase#calculateInvJShFuncDerSWithFe for calculation.
    gsl_matrix 	**elasticStress;					///< The array of matrices for elastic stress of the element. The array stores a 6 by 6 matrix for each Gauss point (there are 6 Gauss points for prisms).
    gsl_matrix 	**viscousStress;					///< The array of matrices for internal viscous stress of the element. The array stores a 6 by 6 matrix for each Gauss point (there are 6 Gauss points for prisms).
    gsl_matrix*	TriPointF;							///< The deformation matrix of the element resulting from iteration over all Gauss points. The dimensions of the matrix is ShapeBase#nDim by ShapeBase#nDim.
	gsl_matrix*	ElementalElasticSystemForces;   	///< The matrix stores the elemental elastic forces. The dimensions of the matrix is ShapeBase#nNodes by ShapeBase#nDim.
	gsl_matrix*	ElementalInternalViscousSystemForces; ///< The matrix stores the elemental internal viscous forces. The dimensions of the matrix is ShapeBase#nNodes by ShapeBase#nDim.
    double* 	detFs;								///< The array stores the determinant of the deformation matrix for each Gauss point.
    double 		ZProjectedBasalArea;				///< The z-projected area of the basal surface of the element.
    double 		ZProjectedApicalArea;				///< The z-projected area of the apical surface of the element.
    double 		BasalArea;							///< The area of the basal surface of the element.
    double 		ApicalArea;							///< The area of the apical surface of the element.
    double		exposedLateralAreaApicalSide;		///< The area of the element on a linker position, and has lateral sides exposed to outside of the tissue, on the apical side, therefore should feel external viscosity.
    double		exposedLateralAreaBasalSide;		///< The area of the element on a linker position, and has lateral sides exposed to outside of the tissue, on the basal side, therefore should feel external viscosity.
    double 		cMyoUniform[2]; 					///< Myosin concentration in the uniformly distributed pool, array of 2: [apical][basal]
    double 		cMyoUnipolar[2]; 					///< Myosin concentration in the polarised pool, array of 2: [apical][basal]
    double 		cMyoUniformEq[2]; 					///< The equilibrium level of myosin concentration in uniformly distributed pool, array of 2: [apical][basal]
    double 		cMyoUnipolarEq[2]; 					///< The equilibrium level of myosin concentration in polarised pool, array of 2: [apical][basal]
    gsl_matrix*	myoPolarityDir;						///< The orientation of myosin polarity, unit vector in world coordinates.
    bool		cellsMigrating;						///< The boolean stating if the cells inside the element are migrating


    bool elementHasExposedApicalSurface;			///< The boolean stating if the element has any apical surface exposed to the environment
    bool elementHasExposedBasalSurface;				///< The boolean stating if the element has any basal surface exposed to the environment
    bool elementHasExposedLateralApicalSurface;		///< The boolean stating if the element has any lateral surfaces exposed to the environment on the apical side of the tissue
    bool elementHasExposedLateralBasalSurface;		///< The boolean stating if the element has any lateral surfaces exposed to the environment on the basal side of the tissue
    int exposedApicalSurfaceNodeIds[3];				///< The int array of size 3, listing the node IDs of element that form the exposed apical surface. The IDs are the node IDs on the element (0-5 for prism), not the actual Node#Id.
    int exposedBasalSurfaceNodeIds[3];				///< The int array of size 3, listing the node IDs of element that form the exposed basal surface. The IDs are the node IDs on the element (0-5 for prism), not the actual Node#Id.
    int exposedLateralAreaApicalSideNodeIds[4];		///< The int array of size 4, listing the node IDs of element that form the lateral surface exposed apically. The IDs are the node IDs on the element (0-5 for prism), not the actual Node::Id.
    int exposedLateralAreaBasalSideNodeIds[4];		///< The int array of size 4, listing the node IDs of element that form the lateral surface exposed basally. The IDs are the node IDs on the element (0-5 for prism), not the actual Node::Id.
    int nLateralSurfaceAreaNodeNumber;				///< Number of nodes that form the lateral surfaces for the element.
    int nSurfaceAreaNodeNumber;						///< Number of nodes that form the apical/basal surfaces for the element.

    double	stiffnessPerturbationRateInSec;	///< The rate at which the stiffness ofthe element will be perturbed, used with the model inputs from "Stiffness_Perturbation:" header in model input file
    double minimumValueOfStiffnessMultiplier;
    double maximumValueOfStiffnessMultiplier;

    double mutationGrowthRatePerSec;
    double mutationGrowthFold;
    void 	setShapeType(string TypeName);				///< The function sets the type of the shape.
    void 	readNodeIds(int* tmpNodeIds);				///< The function sets the Node#Id array that constructs the shape.
    void 	setPositionMatrix(vector<Node*>& Nodes);	///< The function sets the ShapeBase#Positions matrix to define the locations of each constructing node.
    void 	setTissuePlacement(vector<Node*>& Nodes);	///< The function sets the placement of the element within the tissue
    void 	setTissueType(vector<Node*>& Nodes);		///< The function sets the tissue type of the element
    void 	setReferencePositionMatrix();				///< The function sets the RefereneceShapeBase#Positions matrix to define the reference positions of the element.
    void 	setIdentificationColour();					///< The function sets the unique ShapeBase#IdentifierColour colour for the element, which is used in element picking from the user interface.
    void 	rotateReferenceElementByRotationMatrix(double* rotMat); ///< The function rotates the reference of the element (ShapeBase#ReferenceShape) by input rotation matrix, provided as a double pointer of 9 doubles.
    bool 	InvertMatrix(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse); ///< The function takes the first input matrix, and writes the inverse on the second input. False is returned if the matrix is not inverted. Input format is ublas matrices (slow).
    bool 	InvertMatrix(gsl_matrix* input, gsl_matrix* inverse); ///< The function takes the first input matrix, and writes the inverse on the second input. False is returned if the matrix is not inverted. Input format is gsl matrices (fast).

	void 		updateNodeIdsFromSave(ifstream& file);	///< The function reads the ShapeBase#NodeIds of the current shape from save file provided as input.
	void 		updateReferencePositionMatrixFromSave(ifstream& file); ///< The function reads and updates the ShapeBase#ReferenceShape positions (ReferenceShapeBase#Positions) of the current shape from save file provided as input.
	virtual void calculateReferenceVolume(){ParentErrorMessage("calculateReferenceVolume");};  ///<Virtual function of the ShapeBase class to calculate volume of the ShapeBase#ReferenceShape

	bool 		calculateGrowthStrainsRotMat(double* v);	///< The function calculates the rotation matrix to apply on growth strains to align growth with the current x axis of the tissue.
	void		calculateForces3D(vector <Node*>& Nodes,  gsl_matrix* displacementPerDt, bool recordForcesOnFixedNodes, double **FixedNodeForces); ///< The function calculates the viscous and elastic forces generated by the element.
	gsl_matrix* calculateEForNodalForcesKirshoff(gsl_matrix* C);
	gsl_matrix* calculateCauchyGreenDeformationTensor(gsl_matrix* Fe);
	gsl_matrix* calculateSForNodalForcesKirshoff(gsl_matrix* E);
	gsl_matrix* calculateSForNodalForcesNeoHookean(gsl_matrix* invC, double lnJ);
	void 		updateLagrangianElasticityTensorNeoHookean(gsl_matrix* invC,double lnJ, int pointNo);
	gsl_matrix* calculateCompactStressForNodalForces(double detFe,gsl_matrix* Fe, gsl_matrix* S, gsl_matrix *Stress);
    gsl_matrix* calculateInverseJacobianStackForNodalForces(gsl_matrix* Jacobian);
    gsl_matrix* calculateBTforNodalForces(gsl_matrix* InvJacobianStack, gsl_matrix* ShapeFuncDerStack, gsl_matrix *B, gsl_matrix* invJShFuncDerS);
    void		calculateInvJShFuncDerSWithFe(gsl_matrix * currFe, gsl_matrix * InvDXde, gsl_matrix* ShapeFuncDerStack, gsl_matrix *invJShFuncDerSWithF);
    gsl_matrix* calculateVelocityGradientTensor(gsl_matrix* B, gsl_matrix* displacementPerDt);
    gsl_matrix* constructElementalDisplacementMatrix(gsl_matrix* displacement); 		///< The function will assemble elemental node displacement matrix from the input displacement matrix for the whole system.
    gsl_matrix* calculateRateOfDeformationTensor(gsl_matrix* l);						///< The function will calculate rate of deformation tensor from velocity gradient tensor
    void 		calculateViscousStress(gsl_matrix* d, gsl_matrix* viscousStress);		///< The function will calculate internal viscous stress of the element from rate of deformation matrix.
    void 		calculateViscousForces(gsl_matrix*  gv, gsl_matrix*  BTdetFdetdXde, gsl_matrix* viscousStress); ///< The function  will calculate the elemental viscous forces from viscous stress.

    void    	consturctBaTBb(gsl_matrix* B, gsl_matrix* BaT, gsl_matrix* Bb, int a, int b);
    void   		calculateElasticKIntegral1(gsl_matrix* currElementalK,int pointNo);
    void		calculateElasticKIntegral2(gsl_matrix* currElementalK,int pointNo);
    void		calculateViscousKIntegral1(gsl_matrix* currElementalK, gsl_matrix* paranthesisTermForKv1, int pointNo);
    void		calculateViscousKIntegral2(gsl_matrix* currElementalK,int pointNo);
    void		calculateVelocityGradient( gsl_matrix* velocityGradient, gsl_matrix* displacementPerDt, int pointNo);
    void		calculateOuterProduct(gsl_matrix* a, gsl_matrix* b, gsl_matrix* outerProduct);
    gsl_matrix* calculateSymmetricisedTensorProduct(gsl_matrix* a, gsl_matrix* b);

    bool 	disassembleRotationMatrixForZ(gsl_matrix* rotMat);
    bool 	calculate3DRotMatFromF(gsl_matrix* rotMat);

    gsl_matrix* D;
    gsl_matrix* CoeffMat;
    //double D81[3][4][4][4][4];
    double***** D81;
    //boost::numeric::ublas::vector<double> Forces;

	double E, v;
	double internalViscosity;
	double originalInternalViscosity;
    double lambda, mu;

    gsl_matrix* InvFg;		///< Inverse of growth matrix
    gsl_matrix* Fsc;		///< Shape change matrix
    gsl_matrix* InvFsc;		///< Inverse of shape change matrix

    gsl_matrix* TriPointKe;
    gsl_matrix* TriPointKv;

public:
    double 	stiffnessMultiplier;						///< The double for the multiplier that will define Young's modulus stress stiffening.

    gsl_matrix*	remodellingPlaneRotationMatrix;			///< The rotation matrix converting the xyz coordinate system to the plane of remodelling for the lateral elements.
    gsl_matrix* Fg;			///< Growth matrix

    int 	Id;
    int		ShapeDim;
    int* 	NodeIds;
    virtual ~ShapeBase(){
        //while deleting a ShapeBase* that happens to point a child, this destructor will be called after the child destructor
    };
    double** Positions;
    ReferenceShapeBase* ReferenceShape;
    gsl_matrix* Strain;

    //bool 	IsGrowing;
    bool 	isFlipped;
    bool 	IsChangingShape;
    bool	willBeRefined;
    bool	ApicalNormalForPackingUpToDate;
    bool	BasalNormalForPackingUpToDate;
    int 	tissuePlacement; //1 -> apical, 0 -> basal, 2->middle, 3 -> lateral
    int 	tissueType;	///< The tissue type is 0 for columnar layer, 1 for peripodial membrane, and 2 for linker zone
    bool	spansWholeTissue; ///< Boolean staing is the element spans the whole tissue. This is used to identify mid-layer tagged tissues (tissuePlacement = 2), that should still have apical abd basal responses (such as myosin).
    int 	compartmentType; ///< integer identifying the compartment of the tissue in DV axis, 0 pouch, 1 hinge, 2 notum
    double  compartmentIdentityFraction;
    bool	isECMMimicing;
    bool	isECMMimimcingAtCircumference;
    bool	atBasalBorderOfECM;
    bool	isActinMimicing;
    bool	atApicalBorderOfActin;
	bool	IsAblated;
	bool	atSymetricityBoundary;
	bool	IsClippedInDisplay;
	bool 	IsXSymmetricClippedInDisplay;
	bool	IsYSymmetricClippedInDisplay;
	double 	CurrShapeChangeToAdd[3];
	double* ApicalNormalForPacking;
	double* BasalNormalForPacking;
    double  GrownVolume;
	double  VolumePerNode;
	bool 	capElement;
	double** MyoForce;
	int* 	elementsIdsOnSameColumn;
	int 	basalNeigElementId; 	///<This is recorded only for apical nodes of the columnar layer. If not recorded, id is -1.
	bool 	insideEllipseBand;
	int 	coveringEllipseBandId;
    gsl_matrix* ECMThicknessPlaneRotationalMatrix;

	double emergentShapeLongAxis[2];
	double emergentShapeShortAxis[2];

	double plasticDeformationHalfLifeMultiplier;
    bool isMutated;

    bool thereIsGrowthRedistribution;
    bool growthRedistributionShrinksElement;
    double growthRedistributionScale;

    double* apicalNormalCurrentShape;
	int 	getId();
	string 	getName();
	int 	getShapeType();
	int 	getNodeNumber();
	int* 	getNodeIds();
	int		getNodeId(int i);
	int 	getDim();
	int* 	getIdentifierColour();
	double* getCentre();
	double	getPeripodialness();
	double	getColumnarness();
	void	getRelativePositionInTissueInGridIndex(int nGridX, int nGridY, int& IndexX, int& IndexY, double& FracX, double& FracY);
	void	getInitialRelativePositionInTissueInGridIndex(int nGridX, int nGridY, int& IndexX, int& IndexY, double& FracX, double& FracY);
	double	getStiffnessMultiplier();
	double 	getCurrentVolume();
	gsl_matrix* getCurrentFe();
	void 	relaxElasticForces();
	bool 	isGrowthRateApplicable(int sourceTissue, double& weight, double zmin, double zmax);
	void 	updateGrowthWillBeScaledDueToApikobasalRedistribution(bool thisFunctionShrinksApical, double scale, vector<int>& ellipseBandIdsForGrowthRedistribution);
	void 	scaleGrowthForZRedistribution( double& x, double& y, double& z,int sourceTissue);
	void 	calculateFgFromRates(double dt, double x, double y, double z, gsl_matrix* rotMat, gsl_matrix* increment, int sourceTissue, double zMin, double zMax);
	void 	calculateFgFromGridCorners(int gridGrowthsInterpolationType, double dt, GrowthFunctionBase* currGF, gsl_matrix* increment, int sourceTissue, int IndexX, int IndexY, double FracX, double dFracY);
	gsl_matrix* getGrowthIncrement();
	void 	updateGrowthIncrement(gsl_matrix* columnar, gsl_matrix* peripodial);
	void 	updateGrowthByMutation(double dt);
	void	scaleGrowthIncrement(double multiuplier);
	double 	getGrowthMutationMultiplier();
	void 	calculateShapeChangeIncrementFromRates(double dt, double rx, double ry, double rz, gsl_matrix* increment);
	void 	updateShapeChangeIncrement(gsl_matrix* columnarShapeChangeIncrement);
	void	calculateRelativePosInBoundingBox(double boundingBoxXMin, double boundingBoxYMin, double boundingBoxLength, double boundingBoxWidth);
	void	mutateElement(double growthFold, double growthRatePerHour);
	//void	calculateRelativePosInBoundingBox(double columnarBoundingBoxXMin, double columnarBoundingBoxYMin, double columnarBoundingBoxLength, double columnarBoundingBoxWidth, double peipodialBoundingBoxXMin, double peipodialBoundingBoxYMin, double peipodialBoundingBoxLength, double peipodialBoundingBoxWidth);
	void	updateReferencePositionMatrixFromInput(double** input);
	void	displayRelativePosInBoundingBox();
	void	getRelativePosInBoundingBox(double* relativePos);
	void 	setRelativePosInBoundingBox(double x, double y);
	void	setInitialRelativePosInBoundingBox();
	void 	setInitialZPosition(double zMax, double TissueHeight);
	void	getInitialRelativePosInBoundingBox(double* relativePos);
	//void	getRelativePosInColumnarBoundingBox(double* relativePos);
	//void	getRelativePosInPeripodialBoundingBox(double* relativePos);
	void 	convertRelativePosToGridIndex(double* relpos, int& indexX, int &indexY, double &fracX, double &fracY, int nGridX, int nGridY);
	void 	getStrain(int type, float &StrainMag);
	void 	getNodeBasedPysProp(int type, int NodeNo, vector<Node*>& Nodes, float& PysPropMag);
    void 	getPysProp(int type, float &PysPropMag, double dt);
    double	getInternalViscosity();
    double	getOriginalInternalViscosity();
    void updateInternalViscosityTest();
    double 	getYoungModulus();
	double 	getPoissonRatio();
	double* getGrowthRate();
	double* getShapeChangeRate();
	double** getReferencePos();
    void    getPos(gsl_matrix* Pos);
    gsl_matrix* getFg();
    gsl_matrix* getInvFg();
    gsl_matrix* getFsc();
    gsl_matrix* getInvFsc();
    gsl_matrix* getFe();
    double 	getZRemodellingSoFar();
    void 	setZRemodellingSoFar(double zRemodellingSoFar);
	void 	displayName();
	void	displayNodeIds();
	void 	displayPositions();
	void 	displayReferencePositions();
	void 	displayIdentifierColour();
    void    setFg(gsl_matrix* currFg);
	void 	setGrowthWeightsViaTissuePlacement (double periWeight);
	void 	setYoungsModulus(double E);
    virtual void setElasticProperties(double /*EApical*/,double /*EBasal*/, double /*EMid*/, double /*EECM*/, double /*v*/){ParentErrorMessage("setElasticProperties");};
    virtual void checkEdgeLenghtsForBinding(vector<int>& /*masterIds*/,vector<int>& /*slaveIds*/){ParentErrorMessage("checkEdgeLenghtsForBinding");}
    void 	setViscosity(double viscosityApical,double viscosityBasal, double viscosityMid);
    void 	setViscosity(double viscosityApical,double viscosityBasal);
    void 	setViscosity(double viscosity);
    void	setCellMigration(bool migratingBool);
	double	calculateEmergentShapeOrientation();

    bool isMyosinViaEllipsesAppliedToElement(bool isApical, bool isLateral, vector <int> & myosinEllipseBandIds, int numberOfMyosinAppliedEllipseBands);
    bool isActinStiffnessChangeAppliedToElement(bool ThereIsWholeTissueStiffnessPerturbation, bool ThereIsApicalStiffnessPerturbation, bool ThereIsBasalStiffnessPerturbation, bool ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation, bool ThereIsBasolateralStiffnessPerturbation, vector <int> &stiffnessPerturbationEllipseBandIds, int numberOfStiffnessPerturbationAppliesEllipseBands );
    bool isECMChangeAppliedToElement(bool changeApicalECM, bool changeBasalECM, vector<int> &ECMChangeEllipseBandIds, int numberOfECMChangeEllipseBands);
    bool isShapeChangeAppliedToElement(vector<int> &ellipseBandIds, bool applyBasalECM, bool applyToLateralECM, bool applyApically, bool applyBasally, bool applyMidLayer );
    void 	calculateStiffnessPerturbationRate(bool ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation, double stiffnessPerturbationBeginTimeInSec, double stiffnessPerturbationEndTimeInSec, double stiffnessChangedToFractionOfOriginal);
    void 	updateStiffnessMultiplier(double dt); ///< The funciton will update the actin multiplier as a result of stiffness perturbations.
//    virtual void	fillLateralNeighbours();
    bool	getCellMigration();
    virtual void calculateBasalNormal(double * /*normal*/){ParentErrorMessage("calculateBasalNormal");};
    virtual void calculateApicalNormalCurrentShape(){ParentErrorMessage("calculateApicalNormal");};
    virtual void AlignReferenceBaseNormalToZ(){ParentErrorMessage("AlignReferenceBaseNormalToZ");};
    void 	calculateCurrentGrowthIncrement(gsl_matrix* resultingGrowthIncrement, double dt, double growthx, double growthy, double growthz, gsl_matrix* ShearAngleRotationMatrix);
    void 	updateShapeChangeRate(double x, double y, double z, double xy, double yz, double xz);
    virtual void calculateReferenceStiffnessMatrix(){ParentErrorMessage("calculateReferenceStiffnessMatrix");};
    virtual void calculateElementShapeFunctionDerivatives(){ParentErrorMessage("calculateElementShapeFunctionDerivatives");};
    virtual void calculateCurrNodalForces(gsl_matrix */*gslcurrge*/, gsl_matrix */*gslcurrgv*/, gsl_matrix */*gslcurrF*/, gsl_matrix* /*displacementPerDt*/, int /*pointNo*/){ParentErrorMessage("calculateCurrNodalForces");};
    virtual void calculateCurrTriPointFForRotation(gsl_matrix */*currF*/,int /*pointNo*/){ParentErrorMessage("calculateCurrTriPointFForRotation");};
    virtual void copyElementInformationAfterRefinement(ShapeBase* /*baseElement*/){ParentErrorMessage("copyElementInformationAfterRefinement");};
    virtual void calculateApicalArea(){ParentErrorMessage("calculateApicalArea");};
    virtual void calculateBasalArea(){ParentErrorMessage("calculateBasalArea");};
    double 		calculateCurrentGrownAndEmergentVolumes();
    void 	updateNodeIdsForRefinement(int* tmpNodeIds);
    virtual void updateElasticProperties(){ParentErrorMessage("updateElasticProperties");};
    virtual void  fillLateralNeighbours(vector<Node*>& /*Nodes*/, vector<int>& /*lateralNeigbours*/ ){ParentErrorMessage("fillInLateralNeigbours");};
    void	calculateStiffnessFeedback(double dt);
    void 	writeInternalForcesTogeAndgv(gsl_matrix* ge, gsl_matrix* gvInternal, double** SystemForces, vector <Node*>& Nodes);


    void 	calculateForces(vector <Node*>& Nodes, gsl_matrix* displacementPerDt, bool recordForcesOnFixedNodes, double** FixedNodeForces);
    void 	updatePositions(vector<Node*>& Nodes);
    void	updateReferencePositionsToCurentShape();
    void 	setGrowthRate(double dt, double rx, double ry, double rz);
    void 	setGrowthRateExpFromInput(double x, double y, double z);
    void 	updateGrowthIncrementFromRate();
    void 	cleanMyosinForce();
    void	updateUniformEquilibriumMyosinConcentration(bool isApical, double cEqUniform);
	void	updateUnipolarEquilibriumMyosinConcentration(bool isApical, double cEqUnipolar, double orientationX, double orientationY);
	void	adjustCMyosinFromSave();
	void	updateMyosinConcentration(double dt, double kMyo, bool thereIsMyosinFeedback, double MyosinFeedbackCap);
	double  calculateVolumeForInputShapeStructure(double** shapePositions, int nTriangularFaces, int** triangularFaces, double* midPoint );
	void 	calculatePrincipalStrains3D(double& e1, double &e2,  double &e3, gsl_matrix* eigenVec);
	void 	calculatePrincipalStrains2D(double& e1, double &e2,  double &e3, gsl_matrix* eigenVec);
	void 	calculatePrincipalStrainAxesOnXYPlane(double& e1, double &e2, double& tet);
	void	updateEquilibriumMyoWithFeedbackFromZero(double MyosinFeedbackCap);
	void	updateEquilibriumMyoWithFeedbackFromFixedTotal(double totalMyosinLevel);
	bool	checkIfXYPlaneStrainAboveThreshold(double thres);
	bool 	calculateIfInsideActiveStripe(double initialPoint,double endPoint, double stripeSize1, double stripeSize2);
	double	getCmyosinUniformForNode (int TissuePlacement);
	double	getCmyosinUnipolarForNode (int TissuePlacement);
	void 	getMyosinLevels (double *cMyo);
	void 	getEquilibriumMyosinLevels (double *cMyoEq);
	gsl_matrix* getMyosinDirection();
	void 	setMyosinLevels (double cUni0, double cUni1, double cPol0, double cPol1);
	void 	setEquilibriumMyosinLevels (double cUni0, double cUni1, double cPol0, double cPol1);
	virtual void calculateMyosinForcesAreaBased(double /*forcePerMyoMolecule*/){ParentErrorMessage("calculateMyosinForces");};
	virtual void calculateMyosinForcesTotalSizeBased(double /*forcePerMyoMolecule*/){ParentErrorMessage("calculateMyosinForces");};
    virtual void distributeMyosinForcesAreaBased(bool /*isIsotropic*/, bool /*apical*/, double /*forcePerMyoMolecule*/){ParentErrorMessage("distributeMyosinAreaBased");};
    virtual void distributeMyosinForcesTotalSizeBased(bool /*isIsotropic*/, bool /*apical*/, double /*forcePerMyoMolecule*/){ParentErrorMessage("distributeMyosinToitalSizeBased");};
    void 	setShapeChangeRate(double x, double y, double z, double xy, double yz, double xz);
    void  	setShapeChangeInrementToIdentity();
    void 	updateElementVolumesAndTissuePlacementsForSave(vector<Node*>& Nodes);
    bool 	readNodeIdData(ifstream& file);
    bool	readReferencePositionData(ifstream& file);
    void 	convertPlasticStrainToGrowthStrain();
    void 	setLateralElementsRemodellingPlaneRotationMatrix(double systemCentreX, double systemCentreY);
    void	setECMMimicingElementThicknessGrowthAxis();

    bool areanyOfMyNodesAtCircumference(vector<Node*>& Nodes);

    virtual void  checkHealth(){ParentErrorMessage("checkHealth");};
    void 	resetCurrStepShapeChangeData();
    void    writeKelasticToMainKatrix(gsl_matrix* K);
    void    writeKviscousToMainKatrix(gsl_matrix* K);
    void    calculateImplicitKElastic();
    void 	calculateImplicitKViscous(gsl_matrix* displacementPerDt, double dt);
    void	calculateForceFromStress(int nodeId, gsl_matrix* Externalstress, gsl_matrix* ExternalNodalForces);


	void 	updateShapeFromSave(ifstream& file);
	void 	displayMatrix(boost::numeric::ublas::matrix<double>& mat, string matname);
	void 	displayMatrix(boost::numeric::ublas::matrix<int>& mat, string matname);
	void 	displayMatrix(boost::numeric::ublas::vector<double>& vec, string matname);
    void 	displayMatrix(gsl_matrix* mat, string matname);
    void 	displayMatrix(gsl_vector* mat, string matname);
    void 	createMatrixCopy(gsl_matrix *dest, gsl_matrix* src);
	double	calculateMagnitudeVector3D(double* v);
	double	normaliseVector3D(double* v);
	void	normaliseVector3D(gsl_vector* v);
	double	getNormVector3D(gsl_vector* v);
	double 	determinant3by3Matrix(double* rotMat);
	double 	determinant3by3Matrix(boost::numeric::ublas::matrix<double>& Mat);
    double 	determinant3by3Matrix(gsl_matrix* Mat);
	double 	determinant2by2Matrix(boost::numeric::ublas::matrix<double>& Mat);
	void	calculateRotationAngleSinCos(double* u, double* v, double& c, double& s);
	void	calculateRotationAxis(double* u, double* v,double* rotAx, double c);
	void	constructRotationMatrix(double c, double s, double* rotAx, double* rotMat);
	void	rotateVectorByRotationMatrix(double* u,double* rotMat);
	void	rotateVectorByRotationMatrix(double* u,gsl_matrix* rotMat);

	void    CalculateGrowthRotationByF();
	void 	calculateTriPointFForRatation();
	void 	setPlasticDeformationIncrement(double xx, double yy, double zz);
    void 	growShapeByFg();
    void 	changeShapeByFsc(double dt);
    void	checkIfInsideEllipseBands(int nMarkerEllipseRanges, vector<double> markerEllipseBandXCentres, vector<double> markerEllipseBandR1Ranges, vector<double> markerEllipseBandR2Ranges, vector<Node*>& Nodes);
    bool	checkZCappingInRemodelling(bool volumeConserved, double zRemodellingLowerThreshold, double zRemodellingUpperThreshold, gsl_matrix* increment, gsl_matrix* eigenVec);

    bool	assignSoftHinge(double lowHingeLimit, double highHingeLimit,double softnessLevel);

    void	calculatePlasticDeformation3D(bool volumeConserved, double dt, double plasticDeformationHalfLife, double zRemodellingLowerThreshold, double zRemodellingUpperThreshold);
    void 	addMigrationIncrementToGrowthIncrement(gsl_matrix* migrationIncrement);
    void 	displayDebuggingMatrices();
	virtual double getApicalSideLengthAverage(){return ParentErrorMessage("getApicalSideLengthAverage",0.0);};
	virtual double getBasalSideLengthAverage(){return ParentErrorMessage("getBasalSideLengthAverage",0.0);};
	virtual void getApicalTriangles(vector <int> &/*ApicalTriangles*/){ParentErrorMessage("getApicalTriangles");};
	virtual int getCorrecpondingApical(int /*currNodeId*/){return ParentErrorMessage("getCorrecpondingApical", -100);};
	virtual bool IsThisNodeMyBasal(int /*currNodeId*/){return ParentErrorMessage("IsThisNodeMyBasal", false);};
	virtual bool IsThisNodeMyApical(int /*currNodeId*/){return ParentErrorMessage("IsThisNodeMyApical", false);};
	virtual double getElementHeight(){return ParentErrorMessage("getElementHeight", 0.0);};
	virtual void getRelevantNodesForPacking(int /*TissuePlacementOfPackingNode*/, int& /*id1*/, int& /*id2*/, int& /*id3*/){return ParentErrorMessage("getRelevantNodesForPacking");}
	virtual bool IsPointCloseEnoughForPacking(double* /*Pos*/,  float /*threshold*/, int /*TissuePlacementOfPackingNode*/){return ParentErrorMessage("IsPointCloseEnoughForPacking", false);};
	virtual void calculateNormalForPacking(int /*tissuePlacementOfNormal*/){ParentErrorMessage("calculateNormalForPacking");};
	virtual void AddPackingToSurface(int /*tissueplacement*/, double /*Fx*/, double /*Fy*/,double /*Fz*/, double **/*PackingForces*/, vector<Node*> &/*Nodes*/, bool& /*allCornersFixedX*/, bool& /*allCornersFixedY*/, bool& /*allCornersFixedZ*/){ParentErrorMessage("AddPackingToApicalSurface");};
	virtual void getApicalNodePos(double* /*posCorner*/){ParentErrorMessage("getApicalNodePos");};
	virtual void getBasalNodePos(double* /*posCorner*/){ParentErrorMessage("getBasalNodePos");};
	virtual bool IspointInsideTriangle(int /*tissueplacement*/, double /*x*/, double /*y*/,double /*z*/){return ParentErrorMessage("IspointInsideTriangle",false );};
	virtual void constructElementStackList(const int /*discretisationLayers*/, vector<ShapeBase*>& /*elementsList*/){ParentErrorMessage("constructElementStackList");};
	virtual void getApicalCentre(double* /*centre*/){ParentErrorMessage("getApicalCentre");};
	virtual void getBasalCentre(double* /*centre*/){ParentErrorMessage("getBasalCentre");};
	virtual void getReferenceApicalCentre(double* /*centre*/){ParentErrorMessage("getReferenceApicalCentre");};
	virtual void getReferenceBasalCentre(double* /*centre*/){ParentErrorMessage("getReferenceBasalCentre");};
	virtual double* getApicalMinViscosity(vector<Node*> /*Nodes*/){ParentErrorMessage("getApicalMinViscosity");double* dummy; return dummy;};
	virtual double* getBasalMinViscosity(vector<Node*> /*Nodes*/){ParentErrorMessage("getBasalMinViscosity");double* dummy; return dummy;};
	virtual void copyElementInformationAfterRefinement(ShapeBase* /*baseElement*/, int /*layers*/, bool /*thereIsPlasticDeformation*/){ParentErrorMessage("copyElementInformationAfterRefinement");};
	virtual void checkRotationConsistency3D(){ParentErrorMessage("checkRotationConsistency3D");};
	virtual void getLumenFacingNodeIds(int* /*nodeIds*/, int& /*numberOfTriangles*/){ParentErrorMessage("getLumenFacingNodeIds");};
	virtual bool areNodesDirectlyConnected(int /*node0*/, int /*node1*/){ParentErrorMessage("areNodesDirectlyConnected");};
	bool checkPackingToThisNodeViaState(int ColumnarLayerDiscretisationLAyers, Node* NodePointer);
	bool DoesPointBelogToMe(int IdNode);
	void growShape();
	void assignVolumesToNodes(vector <Node*>& Nodes);
	//void assignSurfaceAreaToNodes(vector <Node*>& Nodes);
    void calculateZProjectedAreas();
    void assignZProjectedAreas(vector <Node*> Nodes);
	void assignElementToConnectedNodes(vector <Node*>& Nodes);
	void removeMassFromNodes(vector <Node*>& Nodes);
	void setECMMimicing(bool IsECMMimicing);
	void setActinMimicing(bool isActinMimicing);

	void 	convertLocalStrainToTissueStrain(double* strainsToAdd);
	virtual void assignExposedSurfaceAreaIndices(vector <Node*>& /*Nodes*/){ParentErrorMessage("assignExposedSurfaceAreaIndices");};
	void	calculateExposedLateralAreaApicalSide();
	void 	calculateExposedLateralAreaBasalSide();
	void 	calculateViscositySurfaces();
	void 	assignViscositySurfaceAreaToNodes(vector <Node*>& Nodes);

	bool RotatedElement;
    gsl_matrix* GrowthStrainsRotMat;

    void calculateEmergentRotationAngles();//this is for test and display purposes, calculates the rotation of the element in plane
    void doesElementNeedRefinement(double areaThreshold, int surfaceidentifier);

	bool 	calculateAlignmentScore(double** RefNormalised);
	void 	bringShapePositionsToOrigin(double** RefNormalised, double* refCentre);
	void 	updateElementsNodePositions(int RKId, double ***SystemForces, vector <Node*>& Nodes, double dt);
	void 	updateReferencePositionMatrixFromMeshInput(ifstream& file);
	void	fillNodeNeighbourhood(vector<Node*>& Nodes);
	void 	checkDisplayClipping(double xClip, double yClip, double zClip);
	double  dotProduct3D(double* u, double* v);
	void	crossProduct3D(double* u, double* v, double* cross);
	void	crossProduct3D(gsl_vector* u, gsl_vector* v, gsl_vector* cross);
	void	alignGrowthCalculationOnReference();
	void	readNewGrowthRate(double* NewGrowth, double& ex, double&ey, double& ez, double& exy, double& exz, double& eyz);
	void	updateUniformOrRingGrowthRate(double* NewGrowth, int GrowthId);
	void	updateGridBasedGrowthRate(double* NewGrowth, int GrowthId, int i, int j);
	virtual void setBasalNeigElementId(vector<ShapeBase*>& /*elementsList*/){ParentErrorMessage("setBasalNeigElementId");};
	bool 	isElementFlippedInPotentialNewShape(int nodeId, double newX, double newY, double newZ);
	void 	checkForCollapsedNodes(int TissueHeightDiscretisationLayers, vector<Node*>& Nodes, vector<ShapeBase*>& Elements);
	bool 	hasEnoughNodesOnCurve(vector<Node*>&Nodes);
	void 	assignEllipseBandIdToWholeTissueColumn(int TissueHeightDiscretisationLayers, vector<Node*>& Nodes, vector<ShapeBase*>& Elements);
	void 	assignEllipseBandId(vector<Node*>& Nodes, int selectedEllipseBandId);
	void 	assignEllipseBandIdToNodes(vector<Node*>& Nodes);
};

#endif
