#include "ShapeBase.h"
#include "Node.h"
#include <sstream>

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


void 	ShapeBase::ParentErrorMessage(string functionName){
	cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<endl;
}

bool 	ShapeBase::ParentErrorMessage(string functionName, bool returnValue){
	cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<endl;
	return returnValue;
}

double 	ShapeBase::ParentErrorMessage(string functionName, double returnValue){
	cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<endl;
	return returnValue;
}

int 	ShapeBase::ParentErrorMessage(string functionName, int returnValue){
	cerr<<"You are calling the function: "<<functionName<<" from a parent here, check declaration is via pointers"<<endl;
	return returnValue;
}

void	ShapeBase::setShapeType(string TypeName){
	/**
	 *  The function will set the shape type of the element, stored in variable ShapeBase#ShapeType.
	 *  The mapping is as follows:
	 *  	- 1 	: Prism
	 *  	- 2 	: PrismLateral
	 *  	- 3 	: Tetrahedron
	 *  	- 4 	: Triangle
	 *  	- (-100): Default number if input name is not recognised.
	 *
	 */
	if (TypeName == "Prism"){
		this->ShapeType = 1;
	}
	else if (TypeName == "PrismLateral"){
		this->ShapeType = 2;
	}
	else if (TypeName == "Tetrahedron"){
		this->ShapeType = 3;
	}
	else if (TypeName == "Triangle"){
		//cout<<"set shape type to triangle"<<endl;
		this->ShapeType = 4;
	}
	else{
		this->ShapeType= -100;
	};
	//cout<<"finalised set shape type"<<endl;
}

void	ShapeBase::setIdentificationColour(){
	/**
	 *  The function sets the unique ShapeBase#IdentifierColour, which on click, will allow the
	 *  user interface to identify this element as the selected element. ShapeBase#IdentifierColour
	 *  is an integer array of size 3, and corresponds to rgb channels, in 0-255 range. Therefore, the
	 *  selection tool works for 255^3 = 16581375 elements.
	 *  The colour storing is started from b channel moving towards r, and done with element Id
	 *
	 */
	IdentifierColour[2] = Id % 255;
	int a = (Id - IdentifierColour[2]) / 255;
	IdentifierColour[1] = ( a ) % 255;
	if (a>255){
		IdentifierColour[0] = (a - IdentifierColour[1]) / 255;
	}
	else{
		IdentifierColour[0] = 0;
	}
	//cout<<"IdentifierColour: "<<IdentifierColour[0]<<" "<<IdentifierColour[1]<<" "<<IdentifierColour[2]<<endl;
}

int		ShapeBase::getShapeType(){
	return ShapeType;
}

int 	ShapeBase::getId(){
	return Id;
}

int 	ShapeBase::getNodeNumber(){
	return nNodes;
}

int* 	ShapeBase::getNodeIds(){
	return NodeIds;
}

int		ShapeBase::getNodeId(int i){
	return NodeIds[i];
}

int 	ShapeBase::getDim(){
	return nDim;
}

string 	ShapeBase::getName(){
	string name;
	if (ShapeType == 1){
		name = "Prism";
	}
	else if (ShapeType == 2){
		name = "PrismLateral";
	}
	else if (ShapeType == 3){
		name = "Tetrahedron";
	}
	else if (ShapeType == 4){
		name = "Triangle";
	}
	else{
		name = "Unknown";
	}
	stringstream inter;
	inter.fill('0');
	inter.width(4);
	inter<<Id;
	name = name + inter.str();
	return name;
}

double** ShapeBase::getReferencePos(){
	return ReferenceShape->Positions;
}

void ShapeBase::getPos(gsl_matrix* Pos){
    for (int i=0; i<nNodes; ++i){
        for (int j =0; j<nDim; ++j){
            gsl_matrix_set (Pos, i, j, Positions[i][j]);
        }
    }
}

double 	ShapeBase::getInternalViscosity(){
	return internalViscosity;
}

double 	ShapeBase::getOriginalInternalViscosity(){
	return originalInternalViscosity;
}

double 	ShapeBase::getZRemodellingSoFar(){
	return zRemodellingSoFar;
}

double 	ShapeBase::getStiffnessMultiplier(){
	return stiffnessMultiplier;
}

void 	ShapeBase::setZRemodellingSoFar(double zRemodellingSoFar){
	this -> zRemodellingSoFar = zRemodellingSoFar;
}

void ShapeBase::updateInternalViscosityTest(){
	double d[2] = {0.0,0.0};
	for (int i = 0; i<nNodes; ++i ){
		d[0] += Positions[i][0];
		d[1] += Positions[i][1];
	}
	d[0] /= nNodes; d[1] /= nNodes;
	double dmag = d[0]*d[0] + d[1]*d[1];
	dmag = pow(dmag,0.5);
	internalViscosity = dmag*originalInternalViscosity;
}


double 	ShapeBase::getYoungModulus(){
	return (stiffnessMultiplier*E);
}

double 	ShapeBase::getPoissonRatio(){
	return v;
}

double* ShapeBase::getGrowthRate(){
	//cout<<"Element "<<Id<<" Growth rate: "<<GrowthRate[0]<<" "<<GrowthRate[1]<<" "<<GrowthRate[2]<<endl;
	return GrowthRate;
}

gsl_matrix* ShapeBase::getFg(){
    gsl_matrix* tmpFg =gsl_matrix_calloc(nDim, nDim);
    createMatrixCopy(tmpFg,Fg);
    return tmpFg;
}

gsl_matrix* ShapeBase::getFe(){
    gsl_matrix* tmpFe =gsl_matrix_calloc(nDim, nDim);
    for (int iter =0; iter<3;++iter){
		gsl_matrix_add(tmpFe, FeMatrices[iter]);
	}
    gsl_matrix_scale(tmpFe,1.0/3.0);
    return tmpFe;
}

gsl_matrix* ShapeBase::getInvFg(){
    gsl_matrix* tmpInvFg =gsl_matrix_calloc(nDim, nDim);
    createMatrixCopy(tmpInvFg,InvFg);
    return tmpInvFg;
}


gsl_matrix* ShapeBase::getFsc(){
    gsl_matrix* tmpFsc =gsl_matrix_calloc(nDim, nDim);
    createMatrixCopy(tmpFsc,Fsc);
    return tmpFsc;
}

gsl_matrix* ShapeBase::getInvFsc(){
    gsl_matrix* tmpInvFsc =gsl_matrix_calloc(nDim, nDim);
    createMatrixCopy(tmpInvFsc,InvFsc);
    return tmpInvFsc;
}

void ShapeBase::createMatrixCopy(gsl_matrix* dest, gsl_matrix* src){
    int m = src->size1;
    int n = src->size2;
    gsl_matrix_set_zero(dest);
    double tmp= 0.0;
    for (int i=0; i<m; ++i){
        for (int j=0 ; j<n; ++j){
            tmp = gsl_matrix_get(src,i,j);
            gsl_matrix_set(dest,i,j,tmp);
        }
    }
}

double* ShapeBase::getShapeChangeRate(){
	return ShapeChangeRate;
}

double* ShapeBase::getCentre(){
	double* d = new double[3];
	d[0]= 0.0; d[1]= 0.0; d[2]=0.0;
	for (int i = 0; i<nNodes; ++i ){
		for (int j = 0; j< nDim; ++j){
			d[j] += Positions[i][j];
		}
	}
	d[0] /= nNodes; d[1] /= nNodes; d[2] /= nNodes;
	return d;
}

double ShapeBase::getPeripodialness(){
	return peripodialGrowthWeight;
}

double ShapeBase::getColumnarness(){
	return columnarGrowthWeight;
}


void ShapeBase::getRelativePositionInTissueInGridIndex(int nGridX, int nGridY , int& IndexX, int& IndexY, double& FracX, double& FracY){
	//cout<<"inside getRelativePositionInTissueInGridIndex"<<endl;
	double* reletivePos = new double[2];
	getRelativePosInBoundingBox(reletivePos);
	convertRelativePosToGridIndex(reletivePos, IndexX, IndexY, FracX, FracY, nGridX, nGridY);
	delete[] reletivePos;
}

void ShapeBase::getInitialRelativePositionInTissueInGridIndex(int nGridX, int nGridY, int& IndexX, int& IndexY, double& FracX, double& FracY){
	double* reletivePos = new double[2];
	getInitialRelativePosInBoundingBox(reletivePos);
	convertRelativePosToGridIndex(reletivePos, IndexX, IndexY, FracX, FracY, nGridX, nGridY);
	delete[] reletivePos;

}

bool ShapeBase::isGrowthRateApplicable( int sourceTissue, double& weight, double zmin, double zmax){
	//weight is the weight of the current tissue in linker sites
	if (initialRelativePositionInZ < zmin ||  initialRelativePositionInZ > zmax){
		return false;
	}
	if (sourceTissue == 0){//columnar layer growth
		if (tissueType == 0){ //columnar
			weight = 1.0;
			return true;
		}
		else if(tissueType == 2){ //linker
			weight = columnarGrowthWeight;
			return true;
		}
	}
	else if (sourceTissue == 1){//peripodial membrane growth
		if (tissueType == 1){ //peripodial
			weight = 1.0;
			return  true;
		}
		else if ( tissueType == 2) { //linker
			weight = peripodialGrowthWeight;
			return true;
		}
	}
	return false;
}

void ShapeBase::updateGrowthWillBeScaledDueToApikobasalRedistribution(bool thisFunctionShrinksApical, double scale, vector<int>& ellipseBandIdsForGrowthRedistribution){
	if(tissueType == 0){ //columnar tissue
		bool insideEllipseBandWithRedistribution = false;
		if (!isECMMimicing && insideEllipseBand){
			for (int i=0; i<ellipseBandIdsForGrowthRedistribution.size(); ++i){
				if (coveringEllipseBandId == ellipseBandIdsForGrowthRedistribution[i]){
					insideEllipseBandWithRedistribution =  true;
					break;
				}
			}
		}
		if (insideEllipseBandWithRedistribution){
			thereIsGrowthRedistribution = true;
			if ((tissuePlacement == 1 ) || atApicalBorderOfActin) { //apical
				//this element is apically positioned, it should shrink if the function is shrinking
				//the apical layer, or expand it if not.
				growthRedistributionShrinksElement = thisFunctionShrinksApical;
			}
			else{
				//the element is basally positioned, if the function is shrinking apical, this element should
				//expand, and vice versa.
				growthRedistributionShrinksElement = !thisFunctionShrinksApical;
			}
		}
	}
}

void ShapeBase::scaleGrowthForZRedistribution( double& x, double& y, double& z, int sourceTissue){
	if(thereIsGrowthRedistribution){ //columnar tissue
		double growthRedistributionTime = 24*3600; //apply in terms of the effect in 24 hours
		double shrunkRates[3] = {0,0,0}, increasedRates[3] = {0,0,0};
		double xyGrowth = exp((x+y)*growthRedistributionTime);
		if (growthRedistributionShrinksElement){
			double scaleFactorShrinkage = log((xyGrowth-1)*growthRedistributionScale + 1)/(x+y)/growthRedistributionTime;
			shrunkRates[0] = scaleFactorShrinkage*x;
			shrunkRates[1] = scaleFactorShrinkage*y;
			shrunkRates[2] = z;
			x = shrunkRates[0];
			y = shrunkRates[1];
			z = shrunkRates[2];
		}
		else{
			double growthRedistributionScaleComplementary = 2-growthRedistributionScale;
			double scaleFactorExpansion = log((xyGrowth-1)*growthRedistributionScaleComplementary + 1)/(x+y)/growthRedistributionTime;
			increasedRates[0] = scaleFactorExpansion*x;
			increasedRates[1] = scaleFactorExpansion*y;
			increasedRates[2] = z;
			x = increasedRates[0];
			y = increasedRates[1];
			z = increasedRates[2];
		}
	}
}


double ShapeBase::getCurrentVolume(){
	double J = 0;
	for (int iter =0; iter<numberOfGaussPoints;++iter){
		J +=detFs[iter]*gaussWeights[iter];
	}
	double currentVolume =J * ReferenceShape->Volume;
	return currentVolume;
}

gsl_matrix* ShapeBase::getCurrentFe(){
	gsl_matrix* FeAvr = gsl_matrix_calloc(3,3);
	gsl_matrix* Fetmp = gsl_matrix_calloc(3,3);

	for (int iter =0; iter<numberOfGaussPoints;++iter){
		createMatrixCopy(Fetmp, FeMatrices[iter]);
		gsl_matrix_scale(Fetmp,gaussWeights[iter]);
		gsl_matrix_add(FeAvr, Fetmp);
	}
	gsl_matrix_free(Fetmp);
	return FeAvr;
}

void ShapeBase::relaxElasticForces(){
	gsl_matrix* tmpFgForInversion =gsl_matrix_calloc(nDim,nDim);
	createMatrixCopy(Fg,TriPointF);
	createMatrixCopy(tmpFgForInversion,Fg);
	bool inverted = InvertMatrix(tmpFgForInversion, InvFg);
	if (!inverted){
		cerr<<"Fg not inverted!!"<<endl;
	}
	double detFg = determinant3by3Matrix(Fg);
	GrownVolume = detFg*ReferenceShape->Volume;
	VolumePerNode = GrownVolume/nNodes;
	gsl_matrix_free(tmpFgForInversion);
}

void ShapeBase::calculateFgFromRates(double dt, double x, double y, double z, gsl_matrix* rotMat, gsl_matrix* increment, int sourceTissue, double zMin, double zMax){
	double tissueWeight;
	bool continueCalaculation = isGrowthRateApplicable(sourceTissue, tissueWeight, zMin, zMax);
	if (continueCalaculation){
		scaleGrowthForZRedistribution(x,y,z,sourceTissue);
		double gx = exp(x*tissueWeight*dt);
		double gy = exp(y*tissueWeight*dt);
		double gz = exp(z*tissueWeight*dt);
		//double growthMutationMultiplier = getGrowthMutationMultiplier();
		//gx = (gx - 1)*growthMutationMultiplier + 1;
		//gy = (gy - 1)*growthMutationMultiplier + 1;
		//gz = (gz - 1)*growthMutationMultiplier + 1;
		gsl_matrix_set(increment,0,0,gx);
		gsl_matrix_set(increment,1,1,gy);
		gsl_matrix_set(increment,2,2,gz);
		gsl_matrix* temp = gsl_matrix_calloc(3,3);
		//R * increment
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, rotMat, increment, 0.0, temp);
		//increment * R^T
		gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, temp, rotMat, 0.0, increment);
		gsl_matrix_free(temp);
	}
	else{
		gsl_matrix_set_identity(increment);
	}
}

void ShapeBase::calculateShapeChangeIncrementFromRates(double dt, double rx, double ry, double rz, gsl_matrix* increment){
    if (Id == 0){
    	cout<<" rate: "<<rx<<" "<<ry<<" "<<rz<<" increment: "<<exp(rx*dt)<<" "<<exp(ry*dt)<<" "<<exp(rz*dt)<<endl;
    }
    gsl_matrix_set(increment,0,0,exp(rx*dt));
	gsl_matrix_set(increment,1,1,exp(ry*dt));
	gsl_matrix_set(increment,2,2,exp(rz*dt));
}

void ShapeBase::calculateFgFromGridCorners(int gridGrowthsInterpolationType, double dt, GrowthFunctionBase* currGF, gsl_matrix* increment, int sourceTissue,  int IndexX, int IndexY, double FracX, double FracY){
	double tissueWeight;
	bool continueCalaculation = isGrowthRateApplicable(sourceTissue,tissueWeight,currGF->zMin,currGF->zMax);
	if (continueCalaculation){
		//taking growth data around 4 grid points
		//
		// Grid shape:
		//    [point 2] ------------- [point 3]
		//       |                        |
		//       |<--fracX----> (o)       |
		//       |               |        |
		//       |               |        |
		//       |             fracY      |
		//       |               |        |
		//    [point 0] ------------- [point 1]
		//
		double *growth0, *growth1, *growth2, *growth3;
		growth0 = new double[3];
		growth1 = new double[3];
		growth2 = new double[3];
		growth3 = new double[3];
		double *angles;
		angles = new double[4];
		bool *angleEliminated;
		angleEliminated = new bool[4];
		currGF->getGrowthProfileAt4Corners(IndexX, IndexY, growth0, growth1, growth2, growth3, angles, angleEliminated);
		double growth[3]= {0.0,0.0,0.0};
		double angle = 0;
		if (gridGrowthsInterpolationType == 0){
			//using the growth rate at the grid point:
			//if fraction is below 0,5, I will use the index availabe.
			//If it is above 0.5, it is within the range of next groid point, I will use index +1:
			if(FracX > 0.5){
				if (FracY > 0.5){
					for (int i=0; i<3; ++i){
						growth[i] = growth3[i];
					}
				}
				else{
					for (int i=0; i<3; ++i){
						growth[i] = growth1[i];
					}
				}
			}
			else{
				if (FracY > 0.5){
					for (int i=0; i<3; ++i){
						growth[i] = growth2[i];
					}
				}
				else{
					for (int i=0; i<3; ++i){
						growth[i] = growth0[i];
					}
				}
			}
		}
		else if (gridGrowthsInterpolationType == 1){
			//calculating the angle fraction eliminated, if any:
			double FracEliminated = 0.0;
			if (angleEliminated[0]){ FracEliminated += (1.0-FracX)*(1.0-FracY);	}
			if (angleEliminated[1]){ FracEliminated += FracX*(1.0-FracY);		}
			if (angleEliminated[2]){ FracEliminated += (1.0-FracX)*FracY;		}
			if (angleEliminated[3]){ FracEliminated += FracX*FracY;				}
			//taking the linear interpolation of 4 angles at 4 grid points:
			angle = angles[0]*(1.0-FracX)*(1.0-FracY)+angles[1]*FracX*(1.0-FracY)+angles[2]*(1.0-FracX)*FracY+angles[3]*FracX*FracY;
			if (FracEliminated>0){
				if (FracEliminated >= 0.9999999){
					angle = 0.0; //if all the angles should be eliminated because all corners have low aspect ratio, then angle is arbitrary, selected as zero
				}
				else{
					angle /= (1.0-FracEliminated); //normalising the sum to the eliminated averaging
				}
			}
			//taking the linear interpolation of 4 growth rates at 4 grid points
			for (int axis =0; axis<3; axis++){
				growth[axis]  = growth0[axis]*(1.0-FracX)*(1.0-FracY)+growth1[axis]*FracX*(1.0-FracY)+growth2[axis]*(1.0-FracX)*FracY+growth3[axis]*FracX*FracY;
				growth[axis] *= tissueWeight;
				//if (tissuePlacement == 0){
					//Prevented basal element z growth! Correct this later!!!
					//growth[2] = 0;
				//}
			}
		}
		//write the increment from obtained growth:
		if (isActinMimicing){
			growth[2] = 0; //no z-growth in actin mimicing apical surfaces.
		}
		scaleGrowthForZRedistribution(growth[0],growth[1],growth[2],sourceTissue);

		for (int axis =0; axis<3; axis++){
			double gAxis = exp(growth[axis]*dt);
			//double growthMutationMultiplier = getGrowthMutationMultiplier();
			//gAxis = (gAxis - 1) * growthMutationMultiplier + 1;
			gsl_matrix_set(increment,axis,axis,gAxis);
		}

		//Rotate the growth if the angel is not zero:
		if (angle != 0.0){
			gsl_matrix* rotMat  = gsl_matrix_calloc(3,3);
			double c = cos(angle);
			double s = sin(angle);
			gsl_matrix_set(rotMat,0,0,  c );
			gsl_matrix_set(rotMat,0,1, -1.0*s);
			gsl_matrix_set(rotMat,0,2,  0.0);
			gsl_matrix_set(rotMat,1,0,  s);
			gsl_matrix_set(rotMat,1,1,  c);
			gsl_matrix_set(rotMat,1,2,  0.0);
			gsl_matrix_set(rotMat,2,0,  0.0);
			gsl_matrix_set(rotMat,2,1,  0.0);
			gsl_matrix_set(rotMat,2,2,  1.0);
			gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
			//R * increment
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, rotMat, increment, 0.0, temp);
			//increment * R^T
			gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, temp, rotMat, 0.0, increment);
			gsl_matrix_free(temp);
			gsl_matrix_free(rotMat);
		}
		delete[] growth0;
		delete[] growth1;
		delete[] growth2;
		delete[] growth3;
		delete[] angles;
		delete[] angleEliminated;
	}
	else{
		gsl_matrix_set_identity(increment);
	}
}

gsl_matrix* ShapeBase::getGrowthIncrement(){
	return growthIncrement;
}

void ShapeBase::updateGrowthIncrement(gsl_matrix* columnar, gsl_matrix* peripodial ){
	gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
	if (tissueType == 0){//columnar layer element, no peripodial application necessary
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, columnar, growthIncrement, 0.0, temp);
		createMatrixCopy(growthIncrement, temp);
	}
	else if (tissueType == 1){//peripodial layer element, no columnar application necessary
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, peripodial, growthIncrement, 0.0, temp);
		createMatrixCopy(growthIncrement, temp);
	}
	else if (tissueType == 2){//linker between columnar and peripodial layer element, the growths are already weighted, need to apply both
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, columnar, growthIncrement, 0.0, temp);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, peripodial, temp, 0.0, growthIncrement);
	}
	for (int i=0; i<3; ++i){
		//this is used for display purposes of the simulation. As each new value is added to growthIncrement, I an update this directly, and the result is already cumulative of multiple growth functions
		GrowthRate[i] = gsl_matrix_get(growthIncrement,i,i);
	}
	gsl_matrix_free(temp);
}

void ShapeBase::updateShapeChangeIncrement(gsl_matrix* columnarShapeChangeIncrement){
	//cout<<" Element: "<<Id<<endl;
	//displayMatrix(columnarShapeChangeIncrement,"currIncrementalShapeChangeIncrement");
	gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, columnarShapeChangeIncrement, shapeChangeIncrement, 0.0, temp);
	gsl_matrix_memcpy(shapeChangeIncrement, temp);
	//displayMatrix(shapeChangeIncrement,"updatedShapeChangeIncrement");
}

void ShapeBase::setRelativePosInBoundingBox(double x, double y){
	relativePosInBoundingBox[0] = x;
	relativePosInBoundingBox[1] = y;
}

void ShapeBase::getRelativePosInBoundingBox(double* relativePos){
	relativePos[0] =  relativePosInBoundingBox[0];
	relativePos[1] =  relativePosInBoundingBox[1];
}

void ShapeBase::setInitialRelativePosInBoundingBox(){
	initialRelativePosInBoundingBox[0] = relativePosInBoundingBox[0];
	initialRelativePosInBoundingBox[1] = relativePosInBoundingBox[1];
}

void ShapeBase::setInitialZPosition(double zMin, double TissueHeight){
	//maximum z of the tissue and the tissue height given as input
	double* c = getCentre();
	initialRelativePositionInZ = 1.0 - ( (c[2] - zMin)/TissueHeight );
	//I cannot work with maximum tissue thickness, as there may be peripodial membrane on top
	//the tissue height is applicable to columnar layer only. I need to calculate from the bottom, then convert
	//such that apical surface will have a value of 0, and basal surface will have 1.0;
	//Apical surface is 0, basal surface is 1;
	if (initialRelativePositionInZ < 0){
		initialRelativePositionInZ = 0;
	}
	if (initialRelativePositionInZ > 1.0){
		initialRelativePositionInZ = 1.0;
	}
	delete[] c;
}
void ShapeBase::getInitialRelativePosInBoundingBox(double* relativePos){
	relativePos[0] =  initialRelativePosInBoundingBox[0];
	relativePos[1] =  initialRelativePosInBoundingBox[1];
}

void ShapeBase::convertRelativePosToGridIndex(double* relpos, int& indexX, int &indexY, double &fracX, double &fracY, int nGridX, int nGridY){
	//cout<<"relpos: "<<relpos[0]<<" "<<relpos[1]<<endl;
	relpos[0] *= (float) (nGridX-1);
	relpos[1] *= (float) (nGridY-1);
	indexX = floor(relpos[0]);
	fracX  = relpos[0] - indexX;
	indexY = floor(relpos[1]);
	fracY  = relpos[1] - indexY;
	//cout<<" indexX "<<indexX<<" fracX "<<fracX<<" indexY "<<indexY<<" fracY "<<fracY<<endl;
	if (indexX >= nGridX-1) { //this is for the point that is exactly the point determining the bounding box high end in X, or the side elements for columnar parameter generation (outside the bounding box by definition
		indexX = nGridX-2;
		fracX = 1.0;
		//cout<<" in if 1, indexX: "<<indexX<<" fracX: "<<fracX<<endl;
	}else if (indexX<0){
		indexX = 0;
		fracX = 0.0;
		//cout<<" in if 2, indexX: "<<indexX<<" fracX: "<<fracX<<endl;
	}
	if (indexY >= nGridY-1) {//this is for the point that is exactly the point determining the bounding box high end in X, or the side elements for columnar parameter generation (outside the bounding box by definition
		indexY = nGridY-2;
		fracY = 1.0;
		//cout<<" in if 3, indexY: "<<indexY<<" fracY: "<<fracY<<endl;
	}else if (indexY<0){
		indexY = 0;
		fracY = 0.0;
		//cout<<" in if 4, indexY: "<<indexY<<" fracY: "<<fracY<<endl;
	}
}

void 	ShapeBase::readNodeIds(int* inpNodeIds){
	/**
	 *  This function will take the input of an int pointer that contains the Ids of the nodes that
	 *  construct the element, and write them into the int array ShapeBase#NodeIds of size ShapeBase#nNodes.
	 *  There is no checkpoint to ensure the size of the array input pointer points to is equal to ShapeBase#nNodes.
	 *
	 */
	for (int i=0; i<nNodes; ++i){
		this->NodeIds[i] = inpNodeIds[i];
	}
}

void 	ShapeBase::updateNodeIdsForRefinement(int* inpNodeIds){
	readNodeIds(inpNodeIds);
}

void 	ShapeBase::displayName(){
	/**
	 *  This function will write the shape type and shape id to standard output.
	 */
	cout<<"Type: "<<this->ShapeType<<" Id: "<<this->Id<<endl;
}

void 	ShapeBase::setPositionMatrix(vector<Node*>& Nodes){
	/**
	 *  This function will take the address of a Nodes pointers vector (It should be the Simulation#Nodes vector)
	 *  It will go through the ShapeBase#NodeIds in a nested loop of ShapeBase#nNodes and ShapeBase#nDim,
	 *  read the positions of the corresponding node from the input vector storing the pointers to all the nodes,
	 *  and write the position information into the ShapeBase#Positions array. This is a double storing of the
	 *  same information, yet is practical for elasticity calculations.
	 *
	 */
	const int n = nNodes;
	const int dim = nDim;
	Positions = new double*[n];
	for (int i = 0; i<nNodes; ++i){
		Positions[i] = new double[dim];
		for (int j = 0; j<dim; ++j){
			Positions[i][j] = Nodes[NodeIds[i]]->Position[j];
		}
	}
}

void 	ShapeBase::setTissuePlacement(vector<Node*>& Nodes){
	/**
	 *  This function will take the address of a Nodes pointers vector (It should be the Simulation#Nodes vector)
	 *  Checking the nodes that construct the element, the function will decide the placement of the element
	 *  within the tissue. The placement of the element within the tissue is defined by the ShapeBase#tissuePlacement
	 *  integer such that:
	 *  	- Apical:  ShapeBase#tissuePlacement = 1
	 *  	- Basal:   ShapeBase#tissuePlacement = 0
	 *  	- Mid-line or spans the whole tissue: ShapeBase#tissuePlacement = 2
	 *  	- Lateral: ShapeBase#tissuePlacement = 3
	 *
	 *  The function will first decide which type of nodes the element is composed of.
	 *
	 */
	bool hasApicalNode = false;
	bool hasBasalNode = false;
	bool hasLateralNode = false;
	spansWholeTissue = false;
	for (int i = 0; i<nNodes; ++i){
		if (Nodes[NodeIds[i]]->tissuePlacement == 1){
			hasApicalNode = true;
		}
		else if (Nodes[NodeIds[i]]->tissuePlacement == 0){
			hasBasalNode = true;
		}
		else if (Nodes[NodeIds[i]]->tissuePlacement == 3){
			hasLateralNode = true;
		}
	}
	/**
	*  Lateral elements can be composed of any combination of lateral, mid-line, apical and basal nodes.
	*  Any element that contains at least one lateral node, must be a lateral element.
	*
	*/
	if (hasLateralNode){
		tissuePlacement = 3;
	}
	else{
		/**
		* Elements that does not contain any lateral nodes, but apical nodes must be apical.
		* The only exception to this when the whole tissue is spanned by a single layer of
		* elements. Then the element will have apical and basal nodes, and it will be defined as
		* a mid-line node, with ShapeBase#spansWholeTissue set to true.
		*
		*/
		if (hasApicalNode){
			if (hasBasalNode){
				//the element spans through the whole tissue, the mid-line value should be used
				tissuePlacement = 2;
				spansWholeTissue = true;
			}
			else{
				//the element has only apical and midline nodes, it is apical
				tissuePlacement = 1;
			}
		}
		else if (hasBasalNode){
			/**
			* Elements that does not contain any lateral or apical nodes, but
			* does contain basal nodes then must be basal.
			*
			*/
			//the element only has basal and mid-line nodes, it is basal
			tissuePlacement = 0;
		}
		else{
			/**
			* If the element does not contain any lateral, apical or basal nodes, it must be that the
			* element is composed wholly of mid-line nodes. Then, the element lies in the mid-layer of
			* tissue.
			*
			*/
			//the element has only mid-line nodes, it is mid-line
			tissuePlacement = 2;
		}
	}
}

void ShapeBase::setECMMimicing(bool IsECMMimicing){
	this->isECMMimicing = IsECMMimicing;
	//setting poisson ratio to be zero, so the ECM elements will not be thinning.
	this->v = 0;
	//calculateInitialThickness();
	//setECMMimicingElementThicknessGrowthAxis();
}

void ShapeBase::setActinMimicing(bool isActinMimicing){
	this->isActinMimicing = isActinMimicing;
}

void 	ShapeBase::setTissueType(vector<Node*>& Nodes){
	bool hasColumnarNode = false;
	bool hasPeripodialNode = false;
	bool hasLinkerNode = false;
	for (int i = 0; i<nNodes; ++i){
		if (Nodes[NodeIds[i]]->tissueType == 0){
			hasColumnarNode = true;
		}
		else if (Nodes[NodeIds[i]]->tissueType == 1){
			hasPeripodialNode = true;
		}
		else if (Nodes[NodeIds[i]]->tissueType == 2){
			hasLinkerNode = true;
		}
	}
	if (hasLinkerNode){
		tissueType = 2;
	}
	else if (hasPeripodialNode){
		//ASK LINKER ZONE BEFORE THIS, SOME LINKER ELEMENTS CAN HAVE LINKER NODES AND OTHER TISSUE NODES, NO COLUMNAR ELEMENT OR PERIPODIAL ELEMENT SHOULD HAVE A LINKER NODE
		tissueType = 1;
		setGrowthWeightsViaTissuePlacement( 1.0);//default is set to be columnar, I will not set this for linkers, as they are set in the initiation of peripodial membrane
	}
	else if (hasColumnarNode){
		//ASK PERIPODIAL MEMBRANE BEFORE THIS, SOME PERIPODIAL ELEMENTS CAN HAVE COLUMNAR NODES, AND SOME LINKER ELEMENTS CAN HAVE COLUMNAR NODES. NO COLUMNAR ELEMENT SHOULD HAVE A PERIPODIAL NODE
		tissueType = 0;
	}
	else {
		cerr<<"Element is not placed into tissue correctly, Id: "<<Id<<endl;
	}
	//cout<<"Element : "<<Id<<" hasColumnarNode: "<<hasColumnarNode<<" hasPeripodialmNode "<<hasPeripodialNode<<" tissueType: "<<tissueType<<endl;
}

void 	ShapeBase::setGrowthWeightsViaTissuePlacement (double periWeight){
	peripodialGrowthWeight = periWeight;
	columnarGrowthWeight = 1.0 - peripodialGrowthWeight;
	//cout<<" Element: "<<Id<<" peripodialness: "<<peripodialGrowthWeight<<" columnarness: "<<columnarGrowthWeight<<endl;
}

void 	ShapeBase::setReferencePositionMatrix(){
	/**
	*  This function will allocate the position array of the ReferenceShape (ReferenceShape#Positions)
	*  Then the reference will be equated to the current position of the element.
	*  It is essential this function is called at simulation initiation, and after any subsequent modifications
	*  made to the reference structure of the tissue. But it should not be called in cases where
	*  the current shape of the element has progressed to differ from its reference (such as after loading
	*  steps of a save file). In saves, the reference shapes are set at the initiation of the system,
	*  and the positions of saves steps are loaded afterwards.
	*
	*/
	const int n = nNodes;
	const int dim = nDim;
	ReferenceShape -> Positions = new double*[n];
	for (int i = 0; i<nNodes; ++i){
		ReferenceShape -> Positions[i] = new double[dim];
		for (int j = 0; j<dim; ++j){
			ReferenceShape -> Positions[i][j] = Positions[i][j];
		}
	}
}

void ShapeBase::setFg(gsl_matrix* currFg){
    gsl_matrix_memcpy (Fg, currFg);
    gsl_matrix* tmpFgForInversion =gsl_matrix_calloc(nDim,nDim);
    createMatrixCopy(tmpFgForInversion, Fg);
    bool inverted = InvertMatrix(tmpFgForInversion, InvFg);
    if (!inverted){
    	cerr<<"Fg cannot be inverted!"<<endl;
    }
    gsl_matrix_free(tmpFgForInversion);
}

void ShapeBase::setYoungsModulus(double E){
	this -> E = E;
}

void ShapeBase::setViscosity(double viscosity){
	this -> internalViscosity = viscosity;
}

void ShapeBase::setCellMigration(bool migratingBool){
	cellsMigrating = migratingBool;
}

bool ShapeBase::getCellMigration(){
	return cellsMigrating;
}

void ShapeBase::setViscosity(double viscosityApical,double viscosityBasal, double viscosityMid){
	this -> internalViscosity = viscosityMid;
	if (tissuePlacement == 0 ){
		this -> internalViscosity = viscosityBasal;
	}
	else if(tissuePlacement == 1 ){
		this -> internalViscosity = viscosityApical;
	}
	this -> originalInternalViscosity = internalViscosity;
}

void ShapeBase::setViscosity(double viscosityApical,double viscosityBasal){
	this -> internalViscosity = 0.5*(viscosityApical+viscosityBasal);
	if (tissuePlacement == 0 ){
		this -> internalViscosity = viscosityBasal;
	}
	else if(tissuePlacement == 1 ){
		this -> internalViscosity = viscosityApical;
	}
}

void 	ShapeBase::updateShapeFromSave(ifstream& file){
	file >> IsAblated;
	updateNodeIdsFromSave(file);
	updateReferencePositionMatrixFromSave(file);
	//displayName();
	//displayPositions();
	//displayReferencePositions();
	//displayNodeIds();
}

void 	ShapeBase::updateNodeIdsFromSave(ifstream& file){
	for (int i = 0; i<nNodes; ++i){
		int savedId;
		file >> savedId;
		NodeIds[i] = savedId;
	}
}

bool 	ShapeBase::readNodeIdData(ifstream& file){
	for (int i = 0; i<nNodes; ++i){
		int savedId;
		file >> savedId;

		if (NodeIds[i] != savedId){
			cout<<"NodeId "<<NodeIds[i]<<" savedId "<<savedId<<endl;
			return false;
		}
	}
	return true;
}

void 	ShapeBase::updateReferencePositionMatrixFromSave(ifstream& file){
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			double savedPos;
			file >> savedPos;
			ReferenceShape -> Positions[i][j] = savedPos;
			//cout<<"savedPos: "<<savedPos<<endl;
		}
	}
}

void 	ShapeBase::updateReferencePositionMatrixFromInput(double** inputPos){
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			ReferenceShape -> Positions[i][j] = inputPos[i][j];
		}
	}
}

bool	ShapeBase::readReferencePositionData(ifstream& file){
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			double savedPos;
			file >> savedPos;
			//cout<<Id<<" "<<i<<" "<<j<<" ReferenceShape->Positions: "<<ReferenceShape -> Positions[i][j]<<" savedPos "<<savedPos<<endl;
			if (ReferenceShape -> Positions[i][j] != savedPos){
				//the positions are not equal, it may be an issue of rounding, my satisfactory precision is 2%
				float percentError = (ReferenceShape -> Positions[i][j] - savedPos) / ReferenceShape -> Positions[i][j]*100.0;
				if (percentError>2.0 || percentError< -2.0){
					cout<<"ReferenceShape->Positions: "<<ReferenceShape -> Positions[i][j]<<" savedPos "<<savedPos<<" percent Error: "<<percentError<<endl;
					return false;
				}
			}
		}
	}
	return true;
}

void 	ShapeBase::updateReferencePositionMatrixFromMeshInput(ifstream& file){
	updateReferencePositionMatrixFromSave(file);
}

void ShapeBase::updateElementVolumesAndTissuePlacementsForSave(vector<Node*>& Nodes){
	calculateReferenceVolume();
	setTissuePlacement(Nodes);
	setTissueType(Nodes);
}

void 	ShapeBase::displayNodeIds(){
	for (int i=0; i<nNodes;++i){
			cout<<NodeIds[i]<<"  ";
		cout<<endl;
	}
}

void 	ShapeBase::displayPositions(){
	for (int i=0; i<nNodes;++i){
		for (int j =0; j<nDim; ++j){
			cout<<Positions[i][j]<<"  ";
		}
		cout<<endl;
	}
}

void 	ShapeBase::displayReferencePositions(){
	for (int i=0; i<nNodes;++i){
		for (int j =0; j<nDim; ++j){
			cout<<ReferenceShape ->Positions[i][j]<<"  ";
		}
		cout<<endl;
	}
}

int*	ShapeBase::getIdentifierColour(){
	return IdentifierColour;
}

void 	ShapeBase::getStrain(int type, float &StrainMag){
	StrainMag = 0.0;
	if (type == 0){
		//this is the average strain
        //for (int i=0; i<3; ++i){
        //   StrainMag += gsl_matrix_get(Strain,i,0) ;
        //}
		//StrainMag /= 3;
		//This is volumetric strain, the total volume change:
		gsl_matrix* Fe = getFe();
		StrainMag = determinant3by3Matrix(Fe)-1 ;
		gsl_matrix_free(Fe);
		//StrainMag = 1;
		//for (int i=0; i<3; ++i){
		//	StrainMag *= (1+gsl_matrix_get(Strain,i,0)) ;
		//}
		//StrainMag -= 1.0;
	}
	else if (type == 1){
		//DV
        StrainMag = gsl_matrix_get(Strain,0,0);
	}
	else if (type == 2){
		//AP
        StrainMag = gsl_matrix_get(Strain,1,0);
	}
	else if (type == 3){
		//AB
        StrainMag = gsl_matrix_get(Strain,2,0);
	}
	else if (type == 4){
		//xy
        StrainMag = gsl_matrix_get(Strain,3,0);
	}
	else if (type == 5){
		//yz
        StrainMag = gsl_matrix_get(Strain,4,0);
	}
	else if (type == 3){
		//xz
        StrainMag = gsl_matrix_get(Strain,5,0);
	}
	else{
		return;
	}
}

void 	ShapeBase::getNodeBasedPysProp(int type, int NodeNo, vector<Node*>& Nodes, float& PysPropMag){
	PysPropMag = 0.0;
	if (type == 0){
		PysPropMag = Nodes[NodeIds[NodeNo]] -> externalViscosity[2];
	}
}

void 	ShapeBase::getPysProp(int type, float &PysPropMag, double dt){
	if (type == 1){
		PysPropMag = getInternalViscosity();
	}
	else if (type == 2){
		PysPropMag = getYoungModulus();
	}
	else if (type == 3 ){
		PysPropMag = getPoissonRatio();
	}
	else if (type == 4){
		double* growth;
		growth = getGrowthRate();
        double timescale = 24*60.0*60.0; //reporting per 24 hours
        //for (int i =0 ; i< nDim-1 ; ++i){ //reporting only x & y
        //for (int i =2 ; i< nDim ; ++i){ //reporting only z
        for (int i =0; i< nDim ; ++i){//reporting x,y,z
			//growth is in form exp(r*dt), get r first, then adjust the time scale, and report the exponential form still:
			//And I want to rate of volume growth, that is x*y*z
			double value = exp(log(growth[i])/dt*timescale);
			PysPropMag *= value;
		}
        //converting to percentage increase per 24 hours:
        //PysPropMag -= 1.0;
        //PysPropMag *= 100.0;
	}
	else if (type == 5){
		PysPropMag = GrownVolume/ReferenceShape->Volume;
	}
	else if (type == 6){
		//calculate the emergent volume and shape
		PysPropMag = calculateEmergentShapeOrientation();
	}
	else if (type == 7){
		double* shapechange;
		shapechange = getShapeChangeRate();
		PysPropMag = shapechange[2];
	}
}

double	ShapeBase::calculateEmergentShapeOrientation(){
	//cout<<"calculating emergent shape orientation"<<endl;
	//I want to know in which direction the emergent shape is oriented.
	//I need to have the combination of growth gradient and deformation gradient.
	//It will reflect how the clones would "look":
	double currEmergentVolume = calculateCurrentGrownAndEmergentVolumes();
    gsl_matrix* E;
    gsl_matrix* C = calculateCauchyGreenDeformationTensor(TriPointF);
	E = calculateEForNodalForcesKirshoff(C);
	gsl_matrix* strainBackUp = gsl_matrix_calloc(6,1);
	for (int i=0;i<5;++i){
		gsl_matrix_set(strainBackUp,i,0,gsl_matrix_get(Strain,i,0));
	}
	gsl_matrix_set_zero(Strain);
	gsl_matrix_set(Strain,0,0, gsl_matrix_get(E,0,0));
	gsl_matrix_set(Strain,1,0, gsl_matrix_get(E,1,1));
	gsl_matrix_set(Strain,2,0, gsl_matrix_get(E,2,2));
	gsl_matrix_set(Strain,3,0, 2.0*gsl_matrix_get(E,0,1));
	gsl_matrix_set(Strain,4,0, 2.0*gsl_matrix_get(E,2,1));
	gsl_matrix_set(Strain,5,0, 2.0*gsl_matrix_get(E,0,2));
	double e1 = 0.0, e2 = 0.0, e3 = 0.0;
	gsl_matrix* eigenVec = gsl_matrix_calloc(3,3);
	calculatePrincipalStrains2D(e1,e2,e3,eigenVec);
	//The Eigen vector matrix is a 3 by 3 matrix, but only stores the 2D vectors in its upper corner 2x2 terms
	//The vectors are written in columns of the matrix.
	//the strain I have here is Green strain, I would like to convert it back to deformation
	//gradient terms. Since  E = 1/2 *(Fe^T*Fe-I):
	double F11 = pow(e1*2+1,0.5);
	double F22 = pow(e2*2+1,0.5);
	double AR = F11/F22;
	//cout<<"AR is : "<<AR<<endl;
	//cout<<" currEmergentVolume "<<currEmergentVolume<<endl;
	if (AR < 1.0){
		AR = 1.0/AR;
		emergentShapeLongAxis[0] = AR * gsl_matrix_get(eigenVec,0,1);
		emergentShapeLongAxis[1] = AR * gsl_matrix_get(eigenVec,1,1);
		emergentShapeShortAxis[0] = gsl_matrix_get(eigenVec,0,0);
		emergentShapeShortAxis[1] = gsl_matrix_get(eigenVec,1,0);
	}
	else{
		emergentShapeLongAxis[0] = AR * gsl_matrix_get(eigenVec,0,0);
		emergentShapeLongAxis[1] = AR * gsl_matrix_get(eigenVec,1,0);
		emergentShapeShortAxis[0] = gsl_matrix_get(eigenVec,0,1);
		emergentShapeShortAxis[1] = gsl_matrix_get(eigenVec,1,1);
	}
	//copy the real strains back on the strains matrix
	for (int i=0;i<5;++i){
		gsl_matrix_set(Strain,i,0,gsl_matrix_get(strainBackUp,i,0));
	}
	gsl_matrix_free(C);
	gsl_matrix_free(E);
	gsl_matrix_free(strainBackUp);
	gsl_matrix_free(eigenVec);
	return currEmergentVolume/ReferenceShape->Volume;


}

void 	ShapeBase::displayIdentifierColour(){
	cout <<" IdentifierColour:  "<<IdentifierColour[0]<<" "<<IdentifierColour[1]<<" "<<IdentifierColour[2]<<endl;
}
/*
void 	ShapeBase::resetCurrStepShapeChangeData(){
	for (int i=0;i<3;++i){
		CurrShapeChangeToAdd[i] = 0.0;
	}
	CurrShapeChangeStrainsUpToDate = false;
	IsChangingShape = false;
}
*/
void 	ShapeBase::changeShapeByFsc(double dt){
    gsl_matrix* FscIncrement = gsl_matrix_calloc(nDim,nDim); ///< The increment of shape change that will be induced this step
    if (rotatedGrowth){
    	double rTemp[3] = {0.0,0.0,0.0};
		for (int i = 0; i<3; ++i){
			rTemp[i] = gsl_matrix_get(GrowthStrainsRotMat,i,0)*ShapeChangeRate[0]+gsl_matrix_get(GrowthStrainsRotMat,i,1)*ShapeChangeRate[1]+gsl_matrix_get(GrowthStrainsRotMat,i,2)*ShapeChangeRate[2];
			if ( (ShapeChangeRate[i] <0 && rTemp[i] >0.0 ) || (ShapeChangeRate[i] >0 && rTemp[i] < 0.0) ){
				rTemp[i] *= -1.0;
			}
		}
		for (int i = 0; i<3; ++i){
			ShapeChangeRate[i] = rTemp[i];
		}
	}
    for (int i=0; i<3 ;++i){
    	gsl_matrix_set(FscIncrement,i,i, exp(ShapeChangeRate[i]*dt));
    }
    gsl_matrix* temp1 = gsl_matrix_calloc(nDim,nDim);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, FscIncrement, Fsc, 0.0, temp1);
	gsl_matrix_memcpy(Fsc, temp1);
	gsl_matrix* tmpFscForInversion = gsl_matrix_calloc(nDim,nDim);
	createMatrixCopy(tmpFscForInversion,Fsc);
	bool inverted = InvertMatrix(tmpFscForInversion, InvFsc);
	if (!inverted){
		cerr<<"Fsc not inverted!!"<<endl;
	}
	//double detFsc = determinant3by3Matrix(Fsc);
	//freeing matrices allocated in this function
	gsl_matrix_free(FscIncrement);
	gsl_matrix_free(temp1);
	gsl_matrix_free(tmpFscForInversion);
}

void ShapeBase::setPlasticDeformationIncrement(double xx, double yy, double zz){
	gsl_matrix_set(plasticDeformationIncrement,0,0,xx);
	gsl_matrix_set(plasticDeformationIncrement,1,1,yy);
	gsl_matrix_set(plasticDeformationIncrement,2,2,zz);
}

void 	ShapeBase::growShapeByFg(){
    if (rotatedGrowth){
        gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
        //R^T * growthIncrement
        gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0, GrowthStrainsRotMat, growthIncrement, 0.0, temp);
        //R^T * growthIncrement * R
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, GrowthStrainsRotMat, 0.0, growthIncrement);
    	//rotate shape change increment:
    	gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0, GrowthStrainsRotMat, shapeChangeIncrement, 0.0, temp);
    	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp, GrowthStrainsRotMat, 0.0, shapeChangeIncrement);
    	gsl_matrix_free(temp);
    }
    //incrementing Fg with current growth rate, plastic deformation rate, and shape changes:
    gsl_matrix* temp1 = gsl_matrix_calloc(nDim,nDim);
    gsl_matrix* temp2 = gsl_matrix_calloc(nDim,nDim);
    gsl_matrix* temp3 = gsl_matrix_calloc(nDim,nDim);
    //adding plastic deformation, this increment is in already in correct orientation:
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, plasticDeformationIncrement,growthIncrement, 0.0, temp1);
    //adding shape change, this increment is in already in correct orientation:
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, shapeChangeIncrement,temp1, 0.0, temp2);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temp2, Fg, 0.0, temp3);
    gsl_matrix_memcpy(Fg, temp3);
    gsl_matrix* tmpFgForInversion = gsl_matrix_calloc(nDim,nDim);
    createMatrixCopy(tmpFgForInversion,Fg);
    bool inverted = InvertMatrix(tmpFgForInversion, InvFg);
    if (!inverted){
        cerr<<"Fg not inverted!!"<<endl;
    }
    double detFg = determinant3by3Matrix(Fg);
    GrownVolume = detFg*ReferenceShape->Volume;
    VolumePerNode = GrownVolume/nNodes;
    //freeing matrices allocated in this function
    gsl_matrix_free(temp1);
    gsl_matrix_free(temp2);
    gsl_matrix_free(temp3);
    gsl_matrix_free(tmpFgForInversion);
}

void ShapeBase::displayDebuggingMatrices(){
	//double a = gsl_matrix_get(FeMatrices[0],0,0);
	//if (a>(1.0+10E-4) || a <(1-1E-5) ){
		cout<<" Prism "<<Id<<" rotatedGrowth: "<<rotatedGrowth<<endl;
		displayMatrix(Fg, " Fg");
		displayMatrix(FeMatrices[0], " FeMatrices[0]");
		displayMatrix(FeMatrices[1], " FeMatrices[1]");
		displayMatrix(FeMatrices[2], " FeMatrices[2]");
		displayMatrix(GrowthStrainsRotMat, " GrowthStrainsRotMat");
	//}
}

void ShapeBase::addMigrationIncrementToGrowthIncrement(gsl_matrix* migrationIncrement){
		gsl_matrix* temp1 = gsl_matrix_calloc(3,3);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, migrationIncrement,growthIncrement, 0.0, temp1);
		gsl_matrix_memcpy(growthIncrement, temp1);
		gsl_matrix_free( temp1 );
}


double 	ShapeBase::calculateCurrentGrownAndEmergentVolumes(){
	calculateReferenceVolume();
	double detFg = determinant3by3Matrix(Fg);
    GrownVolume = detFg*ReferenceShape->Volume;
    /*if (Id == 302){
    	cout<<"Element: "<<Id<<" For displaySave, detFg: "<<detFg<<endl;
    	displayMatrix(Fg,"FgForDebugging");
    }*/
	calculateTriPointFForRatation();
	double detF =determinant3by3Matrix(TriPointF);
	double emergentVolume = detF*ReferenceShape->Volume;
    return emergentVolume;

}

bool ShapeBase::isMyosinViaEllipsesAppliedToElement(bool isApical, bool isLateral, vector <int> & myosinEllipseBandIds, int numberOfMyosinAppliedEllipseBands){
	if (!isECMMimicing){
		if ( isLateral //all elements will feel the myosin regardless of placement
			 || (tissuePlacement == 1 && isApical ) //the myosin function is apical, and the element is apical
			 || (tissuePlacement == 0 && !isApical  ) //the myosin function is not apical, and the element is basal
			 || (atBasalBorderOfECM   && !isApical  ) //the myosin function is not apical (applied to basal) and element is at the layer above the "basal elements, but there is basal ECM, therefore the element is basal of the tissue (atBasalBorderOfECM)
			){
			if(insideEllipseBand){ //this covers the tissue type as well as the position, ellipses are applied to columnar and linker zones only.
				for (int myoRangeCounter =0; myoRangeCounter<numberOfMyosinAppliedEllipseBands; ++myoRangeCounter){
					if (coveringEllipseBandId == myosinEllipseBandIds[myoRangeCounter]){
						return true;
					}
				}
			}
		}
	}
	return false;
}

bool ShapeBase::isActinStiffnessChangeAppliedToElement(bool ThereIsWholeTissueStiffnessPerturbation, bool ThereIsApicalStiffnessPerturbation, bool ThereIsBasalStiffnessPerturbation, bool ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation, bool ThereIsBasolateralStiffnessPerturbation, vector <int> &stiffnessPerturbationEllipseBandIds, int numberOfStiffnessPerturbationAppliesEllipseBands ){
	if (!isECMMimicing){
		if( ThereIsWholeTissueStiffnessPerturbation  //the whole columnar tissue is perturbed
			|| ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation	// there is relaxation on the apical surface and stiffenning on the rest of the tissue, further checks needed while calculating the rate
			|| (ThereIsBasolateralStiffnessPerturbation && !tissuePlacement == 1) //there is only basolateral stiffness perturbation, without affecting apical sides
			|| (tissuePlacement == 0 && ThereIsBasalStiffnessPerturbation    ) //the basal surface is perturbed and element is basal
			|| (tissuePlacement == 1 && ThereIsApicalStiffnessPerturbation   ) // the apical surface is perturbed and element is apical
			|| (atBasalBorderOfECM   && ThereIsBasalStiffnessPerturbation    ) //the basal surface is perturbed and element is at the layer above the "basal elements, but there is basal ECM, therefore the element is basal of the tissue (atBasalBorderOfECM)
			){
			if(insideEllipseBand){ //this covers the tissue type as well as the position, ellipses are applied to columnar and linker zones only.
				for (int stiffnessPerturbationRangeCounter =0; stiffnessPerturbationRangeCounter<numberOfStiffnessPerturbationAppliesEllipseBands; ++stiffnessPerturbationRangeCounter){
					if (coveringEllipseBandId == stiffnessPerturbationEllipseBandIds[stiffnessPerturbationRangeCounter]){
						return true;
					}
				}
			}
		}
	}
	return false;
}

bool ShapeBase::isShapeChangeAppliedToElement(vector<int> &ellipseBandIds, bool applyBasalECM, bool applyToLateralECM, bool applyApically, bool applyBasally, bool applyMidLayer ){
	bool checkForEllipseId = false;
	if (tissueType == 1){ //element is on the peripodial membrane, I am not applyign this to peripodial
		return false;
	}
	if (isECMMimicing){
		if  (   (applyBasalECM  && tissuePlacement == 0 )
			 || (applyToLateralECM && isECMMimimcingAtCircumference && !tissuePlacement == 0) //do not grow the basal element twice
			){
			checkForEllipseId = true;
		}
	}
	else{
		if  (   (applyBasally  && (tissuePlacement == 0 || atBasalBorderOfECM) )
				 || (applyApically && tissuePlacement == 1)
				 || (applyMidLayer && tissuePlacement == 2)
			){
				checkForEllipseId = true;
		}
	}
	int nEllipseBands = ellipseBandIds.size();
	//The element qualifies for shape change in this function, is it inside the ellipse band?
	if(checkForEllipseId && insideEllipseBand){ //this covers the tissue type as well as the position, ellipses are applied to columnar and linker zones only.
		for (int shapeChangeEllipseBandCounter =0; shapeChangeEllipseBandCounter<nEllipseBands; ++shapeChangeEllipseBandCounter){
			if (coveringEllipseBandId == ellipseBandIds[shapeChangeEllipseBandCounter]){
				return true;
			}
		}
	}
	return false;
}

bool ShapeBase::isECMChangeAppliedToElement(bool changeApicalECM, bool changeBasalECM, vector<int> &ECMChangeEllipseBandIds, int numberOfECMChangeEllipseBands){
	if (isECMMimicing){
		if  (    (changeApicalECM && tissuePlacement == 1 )
			  || (changeBasalECM  && tissuePlacement == 0 )
			  || (isECMMimimcingAtCircumference)
			  //|| (changeStiffnessBasalECM  && tissuePlacement == 2 )
			){
			if(insideEllipseBand){
				for (int ECMReductionRangeCounter = 0; ECMReductionRangeCounter<numberOfECMChangeEllipseBands; ++ECMReductionRangeCounter){
					if (coveringEllipseBandId == ECMChangeEllipseBandIds[ECMReductionRangeCounter]){
						return true;
					}
				}
			}

		}
	}
	return false;
}


void ShapeBase::calculateStiffnessPerturbationRate(bool ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation, double stiffnessPerturbationBeginTimeInSec, double stiffnessPerturbationEndTimeInSec, double stiffnessChangedToFractionOfOriginal){
    double totalTimePerturbationWillBeAppliedInSec = stiffnessPerturbationEndTimeInSec-stiffnessPerturbationBeginTimeInSec;
    if (totalTimePerturbationWillBeAppliedInSec <0){
        stiffnessPerturbationRateInSec = 0;
        return;
    }
    if (ThereIsBasolateralWithApicalRelaxationStiffnessPerturbation){
    	//the used rate will be different for apical elements and all the remaining elements.
    	//I do not need to check for ECM, as this is called for only the elements that has
    	//applied stiffness perturbations, which already excluded ECM elements.
    	if(tissuePlacement == 1){ //element is apical.
    		//This will not be feasible for elements that span the whole disc. Then you cannot
    		//do a baso-lateral change for elements that cover the whole tissue!
    		//the element is apical, whatever I am applying to the basal side, I will apply the inverse to the apical side.
    		// If the baso-lateral side is doubling, apical surface will halve.
    		stiffnessChangedToFractionOfOriginal = 1.0/stiffnessChangedToFractionOfOriginal;
    	}
    }
    stiffnessPerturbationRateInSec =  (stiffnessChangedToFractionOfOriginal - 1.0)/totalTimePerturbationWillBeAppliedInSec;
    if (stiffnessPerturbationRateInSec<0){
    	minimumValueOfStiffnessMultiplier = stiffnessChangedToFractionOfOriginal;
    }
    else{
    	maximumValueOfStiffnessMultiplier = stiffnessChangedToFractionOfOriginal;
    }
	//cout<<"Element "<<Id<<" updated  rate: "<<stiffnessPerturbationRateInSec<<"compartment "<< compartmentType<<" identity frac "<<compartmentIdentityFraction<<endl;

}

void ShapeBase::updateStiffnessMultiplier(double dt){
    //stiffnessMultiplier *= ( 1.0 + stiffnessPerturbationRateInSec*dt);
	stiffnessMultiplier += stiffnessPerturbationRateInSec*dt;
	if (stiffnessMultiplier<minimumValueOfStiffnessMultiplier){
		stiffnessMultiplier = minimumValueOfStiffnessMultiplier;
	}
	if (stiffnessMultiplier>maximumValueOfStiffnessMultiplier){
		stiffnessMultiplier = maximumValueOfStiffnessMultiplier;
	}
	//cout<<" element:  "<<Id<<" stiffnessMultiplier: "<<stiffnessMultiplier<<endl;
}

void ShapeBase::calculateStiffnessFeedback(double dt){
	if (tissueType == 0 && tissuePlacement == 1){ //apical columnar layer element
		gsl_matrix* Fe = getCurrentFe();
		//taking into account only x&y deformation
		gsl_matrix_set(Fe,2,2,1.0);
		double detFe = determinant3by3Matrix(Fe);
		if (detFe > 0){
			double eqStiffnessMultiplier = 1 + (detFe -1)*10.0;
			cout<<"Id: "<<Id<<" determinant of Fe: "<<detFe<<" eq Actin: "<<eqStiffnessMultiplier<<" actin before update: "<<stiffnessMultiplier<<" ";
			double rate =1.0/3600 * dt;
			if (rate > 1.0) {rate =1.0;}
			stiffnessMultiplier += rate* (eqStiffnessMultiplier-stiffnessMultiplier);
			cout<<" actin after update: "<<stiffnessMultiplier<<endl;
			if (stiffnessMultiplier < 1.0){
				stiffnessMultiplier = 1;
			}
		}
		gsl_matrix_free(Fe);
	}
}

bool	ShapeBase::assignSoftHinge(double lowHingeLimit, double highHingeLimit,double softnessLevel){
	if (!isECMMimicing){
		if (relativePosInBoundingBox[0]>lowHingeLimit && relativePosInBoundingBox[0] < highHingeLimit){
			stiffnessMultiplier *= softnessLevel;
			updateElasticProperties();
		}
	}
}

bool	ShapeBase::checkZCappingInRemodelling(bool volumeConserved, double zRemodellingLowerThreshold, double zRemodellingUpperThreshold, gsl_matrix* increment, gsl_matrix* eigenVec){
	//I want to do some further checks on the z axis, I
	//need to find what the growth on z axis will be:
	bool zCapped = false;
	gsl_matrix* tempForZCap = gsl_matrix_calloc(nDim,nDim);
	gsl_matrix* incrementToTestZCapping = gsl_matrix_calloc(nDim,nDim);
	//eigenVec * increment
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, eigenVec, increment, 0.0, tempForZCap);
	//eigenVec * increment * eigenVec^T
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, tempForZCap, eigenVec, 0.0, incrementToTestZCapping);
	//The increment I have calculated is non-volume conserved, and non-scaled
	//for z. There are limitations to how much z axis can be remodelled.
	//I will check those now and cap the z deformation if necessary:
	double Fzz = gsl_matrix_get(incrementToTestZCapping,2,2);
	if (volumeConserved){
		double det = determinant3by3Matrix(incrementToTestZCapping);
		double scale = 1.0/pow (det,1.0/3.0);
		double zIncrement = Fzz*scale;
		if (zRemodellingSoFar*zIncrement >= zRemodellingUpperThreshold){
			zCapped = true;
		}
		if (zRemodellingSoFar*zIncrement <= zRemodellingLowerThreshold){
			zCapped = true;
		}
	}
	else{
		if (zRemodellingSoFar*Fzz >= zRemodellingUpperThreshold){
			zCapped = true;
		}
		if (zRemodellingSoFar*Fzz <= zRemodellingLowerThreshold){
			zCapped = true;
		}
	}
	gsl_matrix_free(tempForZCap);
	gsl_matrix_free(incrementToTestZCapping);
	return zCapped;
}

void	ShapeBase::checkIfInsideEllipseBands(int nMarkerEllipseRanges, vector<double> markerEllipseBandXCentres,vector<double> markerEllipseBandR1Ranges, vector<double> markerEllipseBandR2Ranges, vector<Node*>& Nodes){
	for (int i=0;i<nMarkerEllipseRanges; ++i){	
		double dx  = relativePosInBoundingBox[0] - markerEllipseBandXCentres [i];
		double dy = relativePosInBoundingBox[1];
		if ( (markerEllipseBandR1Ranges[2*i]> 0 && dx <0) || (markerEllipseBandR1Ranges[2*i]< 0 && dx >0)){
			double dxOverR1 = dx/markerEllipseBandR1Ranges[2*i];
			double dyOverR2 = dy/markerEllipseBandR2Ranges[2*i];
			double d_squareLower = dxOverR1*dxOverR1 + dyOverR2*dyOverR2;
			dxOverR1 = dx/markerEllipseBandR1Ranges[2*i+1];
			dyOverR2 = dy/markerEllipseBandR2Ranges[2*i+1];
			double d_squareUpper = dxOverR1*dxOverR1 + dyOverR2*dyOverR2;
			//For element to be inside the band, calculated distance square values
			//should be larger than 1 for the smaller ellipse, and smaller than 1 for the
			//larger ellipse, hence in between two ellipses.
			if (d_squareLower> 1 && d_squareUpper<1){
				insideEllipseBand = true;
				coveringEllipseBandId = i;
				for (int j =0 ; j<nNodes; ++j){
					int currNodeId = NodeIds[j];
					Nodes[currNodeId]->insideEllipseBand=true;
					Nodes[currNodeId]->coveringEllipseBandId  = i;
					//cout<<"Node "<<currNodeId<<" is inside ellipse"<<i<<endl;
				}
			}
		}		
	}
}

void	ShapeBase::calculatePlasticDeformation3D(bool volumeConserved, double dt, double plasticDeformationHalfLife, double zRemodellingLowerThreshold, double zRemodellingUpperThreshold){
	double e1 = 0.0, e2 = 0.0, e3 = 0.0;
	gsl_matrix* eigenVec = gsl_matrix_calloc(3,3);
	//by default, I will ignore z deformations for columnat tissue.
	//If the element is mimicing an explicit ECM, then it should be able to deform in z too,.
	//Linker zone elements are allowed to deform in z, as we have no information on what they are doing.
	//They do change their z height, and grow in weird patterns. So let them be...
	//bool ignoreZ = false;
	bool checkZCapping = true;
	if( isECMMimicing || tissueType == 2){
		//ignoreZ = false;
		checkZCapping = false;
		//if( Id == 42){
		//	cout<<" Id : "<<Id<<" isECMMimicing: "<<isECMMimicing<<" tissueType: "<<tissueType<<endl;
		//}
	}
	calculatePrincipalStrains3D(e1,e2,e3,eigenVec);
	if (Id ==0){cout<<"e1, e2, e3: "<<e1<<" "<<e2<<" "<<e3<<endl;}
	//displayMatrix(eigenVec,"eigenVec");
	//NowI have the green strain in principal direction in the orientation of the element internal coordinats.
	//I can simply grow the element in this axis, to obtain some form of plastic growth.
	//the strain I have here is Green strain, I would like to convert it back to deformation
	//gradient terms. Since  E = 1/2 *(Fe^T*Fe-I):
	double F11 = pow(e1*2+1,0.5);
	double F22 = pow(e2*2+1,0.5);
	double F33 = pow(e3*2+1,0.5);
	//half life of plastic deformation:
	//Maria's aspect ratio data shows, normalised to initial aspect ratio, if a tissue
	//is stretched to an aspect ratio of 2.0, and relaxed in 20 minutes, it relaxes back to original shape.
	//On the other hand, if it is stretched for 3 hr, it relaxed to an aspect ratio of 1.2.
	//Then the elastic deformation gradient, starting from 2.0, relaxes to a values such that
	//1.2 * Fe = 2.0 -> Fe = 1.6667. Then the deformation I am calculating decays from 1.0 to 0.66667
	//N(0) = 1.0, N(3hr) = 0.66667, then this gives me
	//a half life of 5.12 hr ( N(t) = N(0) * 2 ^ (-t/t_{1/2}) )
	//This is the value set into modelinput file
	double tau = plasticDeformationHalfLifeMultiplier * plasticDeformationHalfLife/(log(2)); // (mean lifetime tau is half life / ln(2))
	double F11t = (F11-1)*exp(-1.0*dt/tau) + 1;
	double F22t = (F22-1)*exp(-1.0*dt/tau) + 1;
	double F33t = (F33-1)*exp(-1.0*dt/tau) + 1;
	F11 = F11/F11t;
	F22 = F22/F22t;
	F33 = F33/F33t;
	//writing onto an incremental matrix:
	gsl_matrix*  increment = gsl_matrix_calloc(3,3);
	gsl_matrix_set(increment,0,0,F11);
	gsl_matrix_set(increment,1,1,F22);
	gsl_matrix_set(increment,2,2,F33);
	bool zCapped = false;
	if (checkZCapping){
		zCapped = checkZCappingInRemodelling(volumeConserved, zRemodellingLowerThreshold, zRemodellingUpperThreshold, increment, eigenVec);
	}
	//If the element is ECM mimicking but not lateral, then I do not want any z remodelling:
	//I need to write the specific conditions that cover the setup where there is explicit ECM, but no
	//peripodial, therefore, the circumference elements should not be zCapped, they should be treated as lateral:
	if ((isECMMimicing) ){
		//first condition will exempt the lateral elements.
		//second condition will exempt circumferential elements that are assigned to be ECM
		if (tissueType != 2) {
			if (!isECMMimimcingAtCircumference){
				zCapped = true;
			}
		}
	}
	if (zCapped){
		calculatePrincipalStrains2D(e1,e2,e3,eigenVec);
		F11 = pow(e1*2+1,0.5);
		F22 = pow(e2*2+1,0.5);
		F11t = (F11-1)*exp(-1.0*dt/tau) + 1;
		F22t = (F22-1)*exp(-1.0*dt/tau) + 1;
		F33t = 1.0;
		F11 = F11/F11t;
		F22 = F22/F22t;
		F33 = 1.0;
		gsl_matrix_set(increment,0,0,F11);
		gsl_matrix_set(increment,1,1,F22);
		gsl_matrix_set(increment,2,2,F33);
	}
	//if( Id == 42){
	//	cout<<" Id : "<<Id<<" volumeConserved: "<<volumeConserved<<" zCapped: "<<zCapped<<endl;
	//}
	//If I am conserving the volume, I need to scale:
	if (volumeConserved){
		double det = determinant3by3Matrix(increment);
		if (zCapped){
			double scale = 1.0/pow (det,1.0/2.0); //scale the size with the square root of the determinant to keep the volume conserved. Then set z to 1 again, we do not want to affect z growth, remodelling is in x & y
			gsl_matrix_scale(increment,scale);
			//I know the eigen vector of the 2nd (in 3 dim indexed from 0) column is z axis, if z is capped
			//Then the incremental value at 2,2 position (in 3 dim indexed from 0) is z
			gsl_matrix_set(increment,2,2,1);
		}
		else{
			double scale = 1.0/pow (det,1.0/3.0); //scaling in x,y, and z
			gsl_matrix_scale(increment,scale);
		}
	}
	//The element is ECM mimicing, and it is lateral. As a temporary solution, I simply will not
	//allow any lateral elements to shrink as part of their remodelling.
	if(isECMMimicing && tissueType == 2){
		for (int i= 0; i<3; ++i){
			double value = gsl_matrix_get(increment,i,i);
			if (value < 1.0){
				gsl_matrix_set(increment,i,i,1.0);
			}
		}
	}
	//The growth I would like to apply now is written on the increment.
	//I would like to rotate the calculated incremental "growth" to be aligned with the coordinate
	//system defined by the eigen vectors.
	//In a growth setup, I calculate the rotationa matrix to be rotation by a certain angle.
	//Here, the rotation matrix is the eigen vector matrix itself. The eigen vector matrix
	//will rotate the identity matrix upon itself.
	gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);
	//eigenVec*increment
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, eigenVec, increment, 0.0, temp);
	//eigenVec*increment *eigenVec^T
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, temp, eigenVec, 0.0, plasticDeformationIncrement);
	//update z remodelling so far to keep track and not update beyond limits
	zRemodellingSoFar *= gsl_matrix_get(plasticDeformationIncrement,2,2);
	//cout<<"plastic deformation of element "<<Id<<endl;
	//displayMatrix(plasticDeformationIncrement,"plasticDeformationIncrement-afternormalplasticdeformation");
	//if (Id == 104 || Id == 104 ){
		//cout<<" Prism "<<Id<<" e1, e2, e3: "<<e1<<" "<<e2<<" "<<e3<<endl;
		//displayMatrix(plasticDeformationIncrement, " plasticDeformationIncrement");
		//displayMatrix(eigenVec," eigenVec");
		//displayMatrix(temp," temp");
		//cout<<"    F11,  F22,  F33 : "<<F11<<" "<<F22<<" "<<F33<<endl;
		//cout<<"    F11t, F22t, F33t: "<<F11t<<" "<<F22t<<" "<<F33t<<endl;
		//cout<<"    plasticDeformationHalfLife: "<<plasticDeformationHalfLife<<" tau: "<<tau<<" exp(-1.0*dt/tau): "<<exp(-1.0*dt/tau)<<endl;
		//cout<<"    zRemodellingSoFar: "<<zRemodellingSoFar<<" zCapped: "<<zCapped<<" limits: "<<zRemodellingLowerThreshold<<" "<<zRemodellingUpperThreshold<<endl;
		//displayMatrix(Fg,"FgBeforeUpdateByPlasticity");
		//displayMatrix(Strain,"Strain");
	//}
	gsl_matrix_free(temp);
	gsl_matrix_free(eigenVec);
	gsl_matrix_free(increment);
}


void 	ShapeBase::CalculateGrowthRotationByF(){
    gsl_matrix* rotMat = gsl_matrix_alloc(3,3);
    gsl_matrix_set_identity(rotMat);
    //rotatedGrowth = false;
    //updating the F for the current shape positions
    //(not using leftovers from previous iteration)
    calculateTriPointFForRatation();
    rotatedGrowth = calculate3DRotMatFromF(rotMat);
    if (rotatedGrowth){
        rotatedGrowth = disassembleRotationMatrixForZ(rotMat);
        if (rotatedGrowth){
            gsl_matrix_transpose(rotMat);
            gsl_matrix_memcpy(GrowthStrainsRotMat,rotMat);
        }
    }
    //freeing matrices allocated in this function
    gsl_matrix_free(rotMat);
}

void 	ShapeBase::calculateTriPointFForRatation(){
	gsl_matrix_set_zero(TriPointF);
    gsl_matrix* currF = gsl_matrix_calloc(nDim,nDim);
    //The point order is established in shape function derivative calculation!
    //Make sure the weights fir in with the order - eta zeta nu:
    //double points[3][3]={{1.0/6.0,1.0/6.0,0.0},{2.0/3.0,1.0/6.0,0.0},{1.0/6.0,2.0/3.0,0.0}};
    for (int iter =0; iter<numberOfGaussPoints;++iter){
        //cout<<"Calculating gauss point: "<<eta<<" "<<nu<<" "<<zeta<<endl;
        calculateCurrTriPointFForRotation(currF,iter);
        gsl_matrix_scale(currF,gaussWeights[iter]);
        gsl_matrix_add(TriPointF, currF);
    }
    gsl_matrix_free(currF);
}

bool 	ShapeBase::disassembleRotationMatrixForZ(gsl_matrix* rotMat){
	//From Extracting euler angles from a rotation matrix by Mike Day of insomniac games:
	//To extract a rotation of R_x(tet_1) R_y(tet_2) R_z(tet_3) from a matrix M where:
	//  M = [ m00   m01   m02
	//        m10   m11   m12
	//        m20   m21   m22]
	// and c1 denote cos(tet_1) and s1 denote sin(tet_1):
	//
	//	tet_1 = atan2(m12,m22) = atan2(s1c2, c1c2)
	//     c2 = sqrt(m00*m00 + m01*m01)
	//  tet_2 = atan2(-m02,c2)
	//  tet_3 = atan2(m01,m00)
    double tethaZ = atan2(gsl_matrix_get(rotMat,0,1),gsl_matrix_get(rotMat,0,0));
    //rotatedGrowth_tethaZ = tethaZ;
    if (tethaZ > 0.017 || tethaZ < -0.017){ //0.017 rad is ~1 degrees
        //rotation is more than 1 degrees, element incremental growth should be rotated
        double c = cos(tethaZ);
        double s = sin(tethaZ);
        gsl_matrix_set(rotMat,0,0,  c );
        gsl_matrix_set(rotMat,0,1, -1.0*s);
        gsl_matrix_set(rotMat,0,2,  0.0);
        gsl_matrix_set(rotMat,1,0,  s);
        gsl_matrix_set(rotMat,1,1,  c);
        gsl_matrix_set(rotMat,1,2,  0.0);
        gsl_matrix_set(rotMat,2,0,  0.0);
        gsl_matrix_set(rotMat,2,1,  0.0);
        gsl_matrix_set(rotMat,2,2,  1.0);
        return true;
    }
    else{
        return false;   //rotation is less than 1 degrees;
    }
}

void 	ShapeBase::calculateEmergentRotationAngles(){
	calculateTriPointFForRatation();
	gsl_matrix* rotMat = gsl_matrix_calloc(3,3);
	calculate3DRotMatFromF(rotMat);
	gsl_matrix_free (rotMat);
}

bool 	ShapeBase::calculate3DRotMatFromF(gsl_matrix* rotMat){
    //if (Id == 0) {displayMatrix(TriPointF,"TriPointF");}
    gsl_matrix* Sgsl = gsl_matrix_alloc (3, 3);
    gsl_matrix* V = gsl_matrix_alloc (3, 3);
    gsl_matrix* R = gsl_matrix_alloc (3, 3);
    gsl_vector* Sig = gsl_vector_alloc (3);
    gsl_vector* workspace = gsl_vector_alloc (3);

    /*
    //Added to test polar decomposition:
    gsl_vector* eigenValues = gsl_vector_calloc(3);
    gsl_matrix* eigenVec = gsl_matrix_calloc(3,3);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(3);
	gsl_matrix* FTF = gsl_matrix_alloc (3, 3);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0, TriPointF, TriPointF,0.0, FTF);
	gsl_eigen_symmv(FTF, eigenValues, eigenVec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eigenValues, eigenVec, GSL_EIGEN_SORT_ABS_ASC);
	gsl_matrix* U = gsl_matrix_calloc(3,3);
	for (int i=0;i<3; ++i){
		double Sii = gsl_vector_get(eigenValues,i);
		Sii = pow(Sii,0.5);
		gsl_matrix_view U2 = gsl_matrix_submatrix (eigenVec, 0, i, 3, 1);
		gsl_matrix* Uincrement = gsl_matrix_calloc(3,3);
		//cout<<" i, "<<i<<" Sii "<<endl;
		gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, &U2.matrix, &U2.matrix,0.0, Uincrement);
		gsl_matrix_scale (Uincrement,Sii);
		gsl_matrix_add(U,Uincrement);
	}
	//inversing U
	gsl_matrix* invU = gsl_matrix_calloc(3,3);
	bool inverted = InvertMatrix(U, invU);
	if (!inverted){
		cerr<<"U not inverted!!"<<endl;
	}
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, TriPointF, invU,0.0, rotMat);
    double polarTethaZNew = atan2(gsl_matrix_get(rotMat,0,1),gsl_matrix_get(rotMat,0,0));

	gsl_vector_free (eigenValues);
	gsl_matrix_free (FTF);
	gsl_matrix_free (U);
	gsl_matrix_free (invU);
	gsl_matrix_free (eigenVec);
	//End OF Polar Decomposition
	*/
    //Singular Value Decomposition
    createMatrixCopy (Sgsl,TriPointF);
    (void) gsl_linalg_SV_decomp (Sgsl, V, Sig, workspace); //This function does have a return value, but I do not need to use it, hence cast it to void!

    //Sgsl ended up as U, I need the rotation matrix U V^T in the decomposition A = U S V^T (jose writes as F = V S U^T in emails, so be careful)
    gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, Sgsl, V,0.0, rotMat);

    double det = determinant3by3Matrix(rotMat);
    if (det<0){
        cout<<"Error! Flipped element, Id: "<<Id<<endl;
        isFlipped = true;
    }
    gsl_matrix_free (Sgsl);
    gsl_matrix_free (V);
    gsl_matrix_free (R);
    gsl_vector_free (Sig);
    gsl_vector_free (workspace);
    //added to compare polar vs single value decomposition
    //double SVDTethaZ = atan2(gsl_matrix_get(rotMat,1,0),gsl_matrix_get(rotMat,0,0));
    //double SVDTethaZNew = atan2(gsl_matrix_get(rotMat,0,1),gsl_matrix_get(rotMat,0,0));
    //cout<<" Element "<<Id<<" tetha from polar: "<<polarTethaZ<<" tetha from SVD: "<<SVDTethaZ<<" tetha from SVDCorrected: "<<SVDTethaZCorrected<<"new version each: "<<polarTethaZNew<<" "<<SVDTethaZNew<<" "<<SVDTethaZCorrectedNew<<endl;


    //Now I need to check if there is only numerical error accumulationg on rotMat, or there is an actual rotation (above 1 degrees):
    double threshold = 0.017; //this is sine 1 degrees
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            if(i != j){
                if (gsl_matrix_get(rotMat,i,j)>threshold || gsl_matrix_get(rotMat,i,j)< (-1.0*threshold)) {
                    return true;
                }
            }
        }
    }
    return false; //none of the off - diagonal terms of the matrix are above the threshold, the current rotation is only numerical error.
}

void ShapeBase::mutateElement(double growthFold, double growthRatePerHour){
	isMutated = true;
	//growth will be uniform in x and y
	mutationGrowthRatePerSec= growthRatePerHour/2.0/3600.0;
	mutationGrowthFold= growthFold;
}

void ShapeBase::updateGrowthByMutation(double dt){
	if (mutationGrowthFold>0){
		//the fold increase is not zero, which means I should be inducing fold increase in growth:
		//I will take the curretn increment, calculate the determinant (absolute growth).
		//Then I will redistribute this growth to x & y:
		double growthPerDt = determinant3by3Matrix(growthIncrement);
		double ratePerSec = log(growthPerDt)/dt;
		double growthPer24hrs = exp(ratePerSec*3600.0*24.0);
		double newGrowthPer24hrs = mutationGrowthFold*growthPer24hrs;
		double newGrowthRatePerSec = log(newGrowthPer24hrs)/24.0/3600.0/2.0;
		//cout<<"current growth: "<<growthPerDt<<" ratePerSec: "<<ratePerSec<<" GrowthPer24hrs:" <<growthPer24hrs<<endl;
		//cout<<" newGrowthPer24hrs: "<<newGrowthPer24hrs<<" newGrowthRatePerSec = "<<newGrowthRatePerSec<<endl;
		setGrowthRate(dt,newGrowthRatePerSec,newGrowthRatePerSec,0.0);
		updateGrowthIncrementFromRate();
	}
	else{
		//overwriting up any growth that might be there, with uniform growth in x & y:
		setGrowthRate(dt,mutationGrowthRatePerSec,mutationGrowthRatePerSec,0.0);
		updateGrowthIncrementFromRate();
	}
}

void 	ShapeBase::calculateRelativePosInBoundingBox(double boundingBoxXMin, double boundingBoxYMin, double boundingBoxLength, double boundingBoxWidth){
	relativePosInBoundingBox = getCentre();
	relativePosInBoundingBox[0] = (relativePosInBoundingBox[0] -boundingBoxXMin) / boundingBoxLength;
	relativePosInBoundingBox[1] = (relativePosInBoundingBox[1] - boundingBoxYMin) / boundingBoxWidth;
    if (Id == 9382){ cout<<" relativePosInBoundingBox: "<<relativePosInBoundingBox[0]<<" "<<relativePosInBoundingBox[1]<<endl;}

}
/*
void 	ShapeBase::calculateRelativePosInBoundingBox(double columnarBoundingBoxXMin, double columnarBoundingBoxYMin, double columnarBoundingBoxLength, double columnarBoundingBoxWidth, double peripodialBoundingBoxXMin, double peripodialBoundingBoxYMin, double peripodialBoundingBoxLength, double peripodialBoundingBoxWidth){
	columnarRelativePosInBoundingBox = getCentre();
	if (tissueType != 0){ //the tissue is not columnar, so there is peripodial membrane
		peripodialRelativePosInBoundingBox[0] = columnarRelativePosInBoundingBox[0];
		peripodialRelativePosInBoundingBox[1] = columnarRelativePosInBoundingBox[1];
		peripodialRelativePosInBoundingBox[0] = (peripodialRelativePosInBoundingBox[0] - peripodialBoundingBoxXMin) / peripodialBoundingBoxLength;
		peripodialRelativePosInBoundingBox[1] = (peripodialRelativePosInBoundingBox[1] - peripodialBoundingBoxYMin) / peripodialBoundingBoxWidth;
		columnarRelativePosInBoundingBox[0] = peripodialRelativePosInBoundingBox[0];
		columnarRelativePosInBoundingBox[1] = peripodialRelativePosInBoundingBox[1];
	}
	else{
		columnarRelativePosInBoundingBox[0] = (columnarRelativePosInBoundingBox[0] - columnarBoundingBoxXMin) / columnarBoundingBoxLength;
		columnarRelativePosInBoundingBox[1] = (columnarRelativePosInBoundingBox[1] - columnarBoundingBoxYMin) / columnarBoundingBoxWidth;
	}
	//cout<<"Element: "<<Id<<" RelPos: "<<columnarRelativePosInBoundingBox[0]<<" "<<columnarRelativePosInBoundingBox[1]<<" "<<peripodialRelativePosInBoundingBox[0]<<" "<<peripodialRelativePosInBoundingBox[1]<<endl;
	//double* a = new double[3];
	//a = getRelativePosInBoundingBox();
	//cout<<" a: "<< a[0]<<" "<<a[1]<<endl;
	//delete[] a;
}
*/
bool ShapeBase::isElementFlippedInPotentialNewShape(int nodeId, double newX, double newY, double newZ){
    const int n = nNodes;
    const int dim = nDim;
    gsl_matrix* currF = gsl_matrix_alloc(dim,dim);
    gsl_matrix* currFe = gsl_matrix_alloc(dim,dim);
    gsl_matrix* CurrShape = gsl_matrix_alloc(n,dim);
    getPos(CurrShape);
    //displayMatrix(CurrShape, "CurrShape_beforeUpdate");
    for (int i=0;i<n;++i){
    	if (NodeIds[i] == nodeId){
    		gsl_matrix_set(CurrShape,i,0,newX);
    		gsl_matrix_set(CurrShape,i,1,newY);
    		gsl_matrix_set(CurrShape,i,2,newZ);
    	}
    }
    //displayMatrix(CurrShape, "CurrShape_afterUpdate");
    //displayNodeIds();
    bool elementWillFlip = false;
    for (int pointNo =0; pointNo<numberOfGaussPoints;++pointNo){
		//calculating dx/de (Jacobian) and reading out dX/de, shape function derivaties:
		gsl_matrix* ShapeFuncDer = ShapeFuncDerivatives[pointNo];
		gsl_matrix* InvdXde = InvdXdes[pointNo];
		gsl_matrix* Jacobian = gsl_matrix_calloc(dim, dim);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShapeFuncDer, CurrShape, 0.0, Jacobian);
		gsl_matrix_transpose(Jacobian);
		//calculating F:
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Jacobian, InvdXde, 0.0, currF);
		//calculating Fe:
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, currF, InvFg, 0.0, currFe);	///< Removing growth
		double detFe = determinant3by3Matrix(currFe);
		if (detFe<0){
		  //element will be flipped if I collapse its node
			elementWillFlip = true;
		}
		gsl_matrix_free(Jacobian);
		if (elementWillFlip){
			break;
		}
    }
	gsl_matrix_free(CurrShape);
	gsl_matrix_free(currFe);
	gsl_matrix_free(currF);
	return elementWillFlip;
}


void ShapeBase::checkForCollapsedNodes(int TissueHeightDiscretisationLayers, vector<Node*>& Nodes, vector<ShapeBase*>& Elements){
	bool elementCollapsed = false;
		for (int j =0 ; j<nNodes; ++j){
		int nodeId = NodeIds[j];
		int collapsedNodeNumber = Nodes[nodeId]->collapsedWith.size();
		if (collapsedNodeNumber>0){elementCollapsed = true;break;}

		//is this node collapsed with another one of elements own nodes?
		/*for (int k=0;k<collapsedNodeNumber; ++k){
			for (int l = j+1 ; l<nNodes; ++l){
				if (Nodes[nodeId]->collapsedWith[k] == NodeIds[l]){
					elementCollapsed = true;
					break;
				}
			}
			if (elementCollapsed){
				break;
			}
		}
		if (elementCollapsed){
			break;
		}*/
	}
	if (elementCollapsed){
		//assignEllipseBandIdToNodes(Nodes);
		insideEllipseBand = true;
		int selectedEllipseBandId = 100;
		if (tissuePlacement == 1){ //apical collapse: ECM relaxation, cell shortening, volume redistribution to shrink top
			selectedEllipseBandId = 100;
		}
		else{ //basal collapse, volume redistribution to shrink bottom
			selectedEllipseBandId = 101;
		}
		coveringEllipseBandId = selectedEllipseBandId;
		cout<<"Id : "<<Id<<" assigned "<<selectedEllipseBandId<<" in function checkForCollapsedNodes "<<endl;
		assignEllipseBandIdToWholeTissueColumn(TissueHeightDiscretisationLayers,Nodes,Elements);
	}
}

bool ShapeBase::hasEnoughNodesOnCurve(vector<Node*>&Nodes){
	int threshold = 3;
	int nodesOnCurveCounter = 0;
	for (int i=0;i<nNodes; i++){
		if (Nodes[NodeIds[i]]->onFoldInitiation){
			nodesOnCurveCounter++;
		}
		if (nodesOnCurveCounter>=threshold){
			return true;
		}
	}
	return false;
}

void ShapeBase::assignEllipseBandIdToWholeTissueColumn(int TissueHeightDiscretisationLayers, vector<Node*>& Nodes, vector<ShapeBase*>& Elements){
	for (int i=0; i < TissueHeightDiscretisationLayers;++i){
		int idOfElementOnSameColumn = elementsIdsOnSameColumn[i];
		Elements[idOfElementOnSameColumn]->assignEllipseBandId(Nodes,coveringEllipseBandId);
	}
}

void ShapeBase::assignEllipseBandId(vector<Node*>& Nodes, int selectedEllipseBandId){
	insideEllipseBand = true;
	coveringEllipseBandId = selectedEllipseBandId;
	assignEllipseBandIdToNodes(Nodes);
}

void ShapeBase::assignEllipseBandIdToNodes(vector<Node*>& Nodes){
	for (int j =0 ; j<nNodes; ++j){
		int nodeId = NodeIds[j];
		Nodes[nodeId]->insideEllipseBand = true;
		Nodes[nodeId]->coveringEllipseBandId = coveringEllipseBandId;
		Nodes[nodeId]->checkOwnersforEllipseAsignment = true;
	}
}

void 	ShapeBase::displayRelativePosInBoundingBox(){
		cout<<"Element: "<<Id<<"  relative position in the tissue bounding box: "<<relativePosInBoundingBox[0]<<" "<<relativePosInBoundingBox[1]<<endl;
}

bool 	ShapeBase::checkPackingToThisNodeViaState(int ColumnarLayerDiscretisationLayers, Node* NodePointer){
	if(IsAblated){
		//if the element is ablated, do not pack against it
		return false;
	}
	if (tissueType == 2){
		//the element is on the lateral section of the tissue, linking peripodial to columnar, no packing on this element
		return false;
	}
	if(ColumnarLayerDiscretisationLayers>1){
		//If the columnar layer is discretised into multiple layers, the apical elements should be checked against apical nodes,
		// and basal nodes should be checked against basal elements. The midline elements should not have packing, BUT on  a single layer tissue, all is midline, therefore
		// this check would not be valid.
		if ( tissuePlacement == 2 ){	//tissue placement of the element is midline in a multi-layered columnar layer, it should not pack to anything
			return false;
		}
		if (NodePointer->tissuePlacement != tissuePlacement){
			//apical nodes pack to apical elements and basal nodes pack to basal elements only
			return false;
		}
	}
	//The node and element are positioned correctly to be able to pack, then does the element belong to the node?
	bool pointBelongsToElement = DoesPointBelogToMe(NodePointer->Id);
	if (pointBelongsToElement){
		return false;
	}
	return true;
}

bool 	ShapeBase::DoesPointBelogToMe(int IdNode){
	for (int i = 0; i<nNodes; ++i){
		if (NodeIds[i] == IdNode){
			return true;
		}
	}
	return false;
}

double 	ShapeBase::determinant3by3Matrix(double* rotMat){
	double det =0.0;
	det  =  rotMat[0]*(rotMat[4]*rotMat[8]-rotMat[5]*rotMat[7]);
	det -= rotMat[1]*(rotMat[3]*rotMat[8]-rotMat[5]*rotMat[6]);
	det += rotMat[2]*(rotMat[3]*rotMat[7]-rotMat[4]*rotMat[6]);
	return det;
}

double 	ShapeBase::determinant3by3Matrix(boost::numeric::ublas::matrix<double>& Mat){
	double det =0.0;
	det  =  Mat(0,0)*(Mat(1,1)*Mat(2,2)-Mat(1,2)*Mat(2,1));
	det -= Mat(0,1)*(Mat(1,0)*Mat(2,2)-Mat(1,2)*Mat(2,0));
	det += Mat(0,2)*(Mat(1,0)*Mat(2,1)-Mat(1,1)*Mat(2,0));
	return det;
}

double 	ShapeBase::determinant3by3Matrix(gsl_matrix* Mat){
    double det =0.0;

    det  = gsl_matrix_get(Mat,0,0)* ( gsl_matrix_get(Mat,1,1)*gsl_matrix_get(Mat,2,2)-gsl_matrix_get(Mat,1,2)*gsl_matrix_get(Mat,2,1) );
    det -= gsl_matrix_get(Mat,0,1)* ( gsl_matrix_get(Mat,1,0)*gsl_matrix_get(Mat,2,2)-gsl_matrix_get(Mat,1,2)*gsl_matrix_get(Mat,2,0) );
    det += gsl_matrix_get(Mat,0,2)* ( gsl_matrix_get(Mat,1,0)*gsl_matrix_get(Mat,2,1)-gsl_matrix_get(Mat,1,1)*gsl_matrix_get(Mat,2,0) );
    return det;
}
double 	ShapeBase::determinant2by2Matrix(boost::numeric::ublas::matrix<double>& Mat){
	double det = Mat(0,0) * Mat(1,1) - Mat(0,1) * Mat(1,0);
	return det;
}

void	ShapeBase::calculateRotationAngleSinCos(double* u, double* v, double& c, double& s){
	//aligning u onto v:
	c = dotProduct3D(u,v);
	if (c > 1.0){
		c = 1.0;
		s = 0.0;

	}
	else if( c<-1.0){
		c = -1.0;
		s = 0.0;
	}
	else{
		double tet = acos(c);
		s = sin(tet);
	}
}

void	ShapeBase::calculateRotationAxis(double* u, double* v,double* rotAx, double c){
	//aligning u onto v:
	if (c>-0.99998){
		crossProduct3D(u,v,rotAx);
		double dummy = normaliseVector3D(rotAx);
	}
	else{
		//the angle is 180 degree, the standard rotation axis calculation will be wrong, I am rotating over x axis at all times;
		rotAx[0]= 1;rotAx[1]= 0;rotAx[2]= 0;
	}
}

void	ShapeBase::constructRotationMatrix(double c, double s, double* rotAx, double* rotMat){
	rotMat[0] = c + rotAx[0]*rotAx[0]*(1 - c);
	rotMat[1] = rotAx[0]*rotAx[1]*(1 - c) - rotAx[2]*s;
	rotMat[2] = rotAx[0]*rotAx[2]*(1 - c) + rotAx[1]*s;

	rotMat[3] = rotAx[1]*rotAx[0]*(1 - c) + rotAx[2]*s;
	rotMat[4] = c + rotAx[1]*rotAx[1]*(1 - c);
	rotMat[5] = rotAx[1]*rotAx[2]*(1 - c) - rotAx[0]*s;

	rotMat[6] = rotAx[2]*rotAx[0]*(1 - c) - rotAx[1]*s;
	rotMat[7] = rotAx[2]*rotAx[1]*(1 - c) + rotAx[0]*s;
	rotMat[8] = c + rotAx[2]*rotAx[2]*(1 - c);
}

void	ShapeBase::rotateVectorByRotationMatrix(double* u,double* rotMat){
	double x = rotMat[0]*u[0]+rotMat[1]*u[1]+rotMat[2]*u[2];
	double y = rotMat[3]*u[0]+rotMat[4]*u[1]+rotMat[5]*u[2];
	double z = rotMat[6]*u[0]+rotMat[7]*u[1]+rotMat[8]*u[2];
	u[0] = x;
	u[1] = y;
	u[2] = z;
}

void	ShapeBase::rotateVectorByRotationMatrix(double* u,gsl_matrix* rotMat){
	double x = gsl_matrix_get(rotMat,0,0)*u[0]+gsl_matrix_get(rotMat,0,1)*u[1]+gsl_matrix_get(rotMat,0,2)*u[2];
	double y = gsl_matrix_get(rotMat,1,0)*u[0]+gsl_matrix_get(rotMat,1,1)*u[1]+gsl_matrix_get(rotMat,1,2)*u[2];
	double z = gsl_matrix_get(rotMat,2,0)*u[0]+gsl_matrix_get(rotMat,2,1)*u[1]+gsl_matrix_get(rotMat,2,2)*u[2];
	u[0] = x;
	u[1] = y;
	u[2] = z;
}


void  ShapeBase::rotateReferenceElementByRotationMatrix(double* rotMat){
	//cout<<"rotating the reference matrix of element: "<<Id<<endl;
	for (int i=0; i<nNodes; ++i){
		double * u;
		u = new double[3];
		for (int j=0; j<nDim; ++j){
			u[j] = ReferenceShape->Positions[i][j];
		}
		rotateVectorByRotationMatrix(u,rotMat);
		for (int j=0; j<nDim; ++j){
			ReferenceShape->Positions[i][j] = u[j];
		}
		delete[] u;
	}
}

void	ShapeBase::calculateForces(vector <Node*>& Nodes, gsl_matrix* displacementPerDt, bool recordForcesOnFixedNodes, double **FixedNodeForces){
	if (ShapeDim == 3){		//3D element
        calculateForces3D(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces);
    }
}

void ShapeBase::writeInternalForcesTogeAndgv(gsl_matrix* ge, gsl_matrix* gvInternal, double** SystemForces, vector <Node*>& Nodes){
	//now all the forces are written on SysyemForces
    //Now I will add the forces into ge, this step can be made faster by separating calculate forces function into two,
    //and filling up either ge or System forces depending on the solution method:
	for (int i = 0; i< nNodes; ++i){
        for ( int j=0; j<nDim; j++){
        	int indexI = nDim*NodeIds[i]+j;
        	double elementalvalue = gsl_matrix_get(ElementalElasticSystemForces,i,j);
        	double matrixValue = gsl_matrix_get(ge,indexI,0);
        	//if (isnan(elementalvalue)){
        	//	cout<<" element: "<<Id<<" ge dimention: "<<j<<" is NaN: "<<elementalvalue<<endl;
        	//}
            gsl_matrix_set(ge, indexI,0,matrixValue + elementalvalue);
        	elementalvalue = gsl_matrix_get(ElementalInternalViscousSystemForces,i,j);
        	matrixValue = gsl_matrix_get(gvInternal,indexI,0);
            gsl_matrix_set(gvInternal, indexI,0,matrixValue + elementalvalue);
        }
    }
    int counter = 0;
    for (int i = 0; i<nNodes; ++i){
        for (int j = 0; j<nDim; ++j){
            if (!Nodes[NodeIds[i]]->FixedPos[j]){
               // cout<<"element: "<<Id<<" writing to system forces for NodeID: "<<NodeIds[i]<<" dim: "<<j<<" initial value "<<SystemForces[NodeIds[i]][j]<<" Felastic: "<<gsl_matrix_get(ElementalElasticSystemForces,i,j)<<" Fvisc: "<<gsl_matrix_get(ElementalInternalViscousSystemForces,i,j)<<endl;
                SystemForces[NodeIds[i]][j] = SystemForces[NodeIds[i]][j] + gsl_matrix_get(ElementalElasticSystemForces,i,j) + gsl_matrix_get(ElementalInternalViscousSystemForces,i,j);

            }
            /*else if(recordForcesOnFixedNodes){
                FixedNodeForces[NodeIds[i]][j] = FixedNodeForces[NodeIds[i]][j] - gsl_matrix_get(TriPointg,counter,0);
            }*/
            counter++;
        }
    }
}



void	ShapeBase::calculateForces3D(vector <Node*>& Nodes,  gsl_matrix* displacementPerDt, bool recordForcesOnFixedNodes, double **FixedNodeForces){
    int dim = nDim;
    int n = nNodes;
    //calculating F and B in a 3 point gaussian:
    gsl_matrix* TriPointge  = gsl_matrix_calloc(dim*n,1);
    gsl_matrix* TriPointgv  = gsl_matrix_calloc(dim*n,1);
    gsl_matrix_set_zero(TriPointF);
    gsl_matrix_set_zero(ElementalElasticSystemForces);
    gsl_matrix_set_zero(ElementalInternalViscousSystemForces);
    gsl_matrix* currge = gsl_matrix_calloc(dim*n,1);
    gsl_matrix* currgv = gsl_matrix_calloc(dim*n,1);
    gsl_matrix* currF = gsl_matrix_calloc(dim,dim);
    //The point order is established in shape function derivative calculation!
    //Make sure the weights fir in with the order - eta zeta nu:
    for (int iter =0; iter<numberOfGaussPoints;++iter){
    	calculateCurrNodalForces(currge, currgv, currF, displacementPerDt, iter);
        gsl_matrix_scale(currge,gaussWeights[iter]);
        gsl_matrix_add(TriPointge, currge);
        gsl_matrix_scale(currgv,gaussWeights[iter]);
        gsl_matrix_add(TriPointgv, currgv);
        gsl_matrix_scale(currF,gaussWeights[iter]);
        gsl_matrix_add(TriPointF, currF);
    }
    int counter = 0;
    for (int i = 0; i<nNodes; ++i){
            for (int j = 0; j<nDim; ++j){
            	if (!Nodes[NodeIds[i]]->FixedPos[j]){
            		double value = gsl_matrix_get(ElementalElasticSystemForces,i,j);
            		value -= gsl_matrix_get(TriPointge,counter,0);
            		gsl_matrix_set(ElementalElasticSystemForces,i,j,value);
            		if(isnan(value)){
						cout<<"force from ElementalElasticSystemForces for element "<<Id<<" dim: "<<j<<" is nan"<<endl;
					}
            		value = gsl_matrix_get(ElementalInternalViscousSystemForces,i,j);
            		value -= gsl_matrix_get(TriPointgv,counter,0);
            		gsl_matrix_set(ElementalInternalViscousSystemForces,i,j,value);
            		if(isnan(value)){
            			cout<<"force from ElementalInternalViscousSystemForces for element "<<Id<<" dim: "<<j<<" is nan"<<endl;
            		}
				}
				counter++;
            }
    }
    //freeing matrices allocated in this function
    gsl_matrix_free(TriPointge);
    gsl_matrix_free(TriPointgv);
    gsl_matrix_free(currge);
    gsl_matrix_free(currgv);
    gsl_matrix_free(currF);
    //cout<<"Element: "<<Id<<endl;
    //displayMatrix(Fg, "Fg");
}


gsl_matrix* ShapeBase::calculateCauchyGreenDeformationTensor(gsl_matrix* Fe){
	//calculating C (C = (Fe^T*Fe):
	gsl_matrix* C =  gsl_matrix_alloc(nDim, nDim);
	gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0, Fe, Fe,0.0, C);
	return C;
}

gsl_matrix* ShapeBase::calculateEForNodalForcesKirshoff(gsl_matrix* C){
    //calculating E ( E = 1/2 *(Fe^T*Fe-I) ; E = 1/2 *(C-I):):
    gsl_matrix* E =  gsl_matrix_alloc(nDim, nDim);
	createMatrixCopy(E,C);
    gsl_matrix* I = gsl_matrix_alloc(nDim, nDim);
    gsl_matrix_set_identity(I);
    gsl_matrix_sub(E,I);
    gsl_matrix_scale(E, 0.5);
    gsl_matrix_free(I);
    return E;
}

gsl_matrix* ShapeBase::calculateSForNodalForcesKirshoff(gsl_matrix* E){
    //calculating S: (S = D:E)
    gsl_matrix_set_zero(Strain);
    gsl_matrix* compactS = gsl_matrix_calloc(6,1);
    gsl_matrix_set(Strain,0,0, gsl_matrix_get(E,0,0));
    gsl_matrix_set(Strain,1,0, gsl_matrix_get(E,1,1));
    gsl_matrix_set(Strain,2,0, gsl_matrix_get(E,2,2));
    gsl_matrix_set(Strain,3,0, 2.0*gsl_matrix_get(E,0,1));
    gsl_matrix_set(Strain,4,0, 2.0*gsl_matrix_get(E,2,1));
    gsl_matrix_set(Strain,5,0, 2.0*gsl_matrix_get(E,0,2));
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, D, Strain,0.0, compactS);

    gsl_matrix* S =  gsl_matrix_alloc(nDim, nDim);
    gsl_matrix_set(S,0,0,gsl_matrix_get(compactS,0,0));
    gsl_matrix_set(S,1,1,gsl_matrix_get(compactS,1,0));
    gsl_matrix_set(S,2,2,gsl_matrix_get(compactS,2,0));
    gsl_matrix_set(S,1,0,gsl_matrix_get(compactS,3,0));
    gsl_matrix_set(S,1,2,gsl_matrix_get(compactS,4,0));
    gsl_matrix_set(S,2,0,gsl_matrix_get(compactS,5,0));
    gsl_matrix_set(S,0,1,gsl_matrix_get(compactS,3,0));
    gsl_matrix_set(S,2,1,gsl_matrix_get(compactS,4,0));
    gsl_matrix_set(S,0,2,gsl_matrix_get(compactS,5,0));

    gsl_matrix_free(compactS);
    return S;
}

gsl_matrix* ShapeBase::calculateSForNodalForcesNeoHookean(gsl_matrix* invC, double lnJ){
	//S = mu (I - C^-1) + lambda (lnJ) C^-1
	gsl_matrix* S =  gsl_matrix_alloc(nDim, nDim);
	createMatrixCopy(S,invC);
	//displayMatrix(S,"S1");
	gsl_matrix* I = gsl_matrix_alloc(nDim, nDim);
	gsl_matrix_set_identity(I);
    gsl_matrix_sub(I,invC);  //(I - C^-1)
    gsl_matrix_scale(I, mu); // mu (I - C^-1)
    gsl_matrix_scale(S, lambda*lnJ); //lambda (lnJ) C^-1
    gsl_matrix_add(S,I); // mu (I - C^-1) + lambda (lnJ) C^-1
    gsl_matrix_free(I);
	return S;
}

void ShapeBase::updateLagrangianElasticityTensorNeoHookean(gsl_matrix* invC, double lnJ, int pointNo){
    //calculating 4th order tensor C, for convenience the matrix D81 is used for both Kirshoff materials and neo-Hookean materials in the code.
    //The documentation lists  Lagrangian Elasticity Tensor with C for neo-Hookean, and with D for Kirshoff materials.
    //lambda is Lame s first parameter and mu is the shear modulus .
	double multiplier = 2*(mu - lambda*lnJ);
	for (int I = 0; I<nDim; ++I){
        for (int J = 0; J<nDim; ++J){
            for (int K = 0; K<nDim; ++K){
                for (int L = 0; L<nDim; ++L){
                    double Iijkl = 0.5* (gsl_matrix_get(invC,I,K)*gsl_matrix_get(invC,J,L) + gsl_matrix_get(invC,I,L)*gsl_matrix_get(invC,J,K));
                    D81[pointNo][I][J][K][L] = lambda*gsl_matrix_get(invC,I,J)*gsl_matrix_get(invC,K,L) + multiplier * Iijkl;
                    //cout<<"D81["<<pointNo<<"]["<<I<<"]["<<J<<"]["<<K<<"]["<<L<<"]"<<D81[pointNo][I][J][K][L]<<endl;
                    /*if (isnan(D81[pointNo][I][J][K][L])){
                    	double invCIJ = gsl_matrix_get(invC,I,J);
                    	double invCKL = gsl_matrix_get(invC,K,L);
                    	cout<<" element: "<<Id<<" D81[pointNo][I][J][K][L] is nan in initial calculation. inverseC[I,J]: "<<invCIJ<<" invCKL: "<<invCKL<<" Iijkl "<<Iijkl<<" multiplier: "<<multiplier<<" lambda: "<<lambda<<" mu: "<<mu<<" lnJ: "<<lnJ<<endl;
                    }*/
                    //D81[I][J][K][L] = lambda*gsl_matrix_get(invC,I,J)*gsl_matrix_get(invC,K,L) + multiplier * gsl_matrix_get(invC,I,K)*gsl_matrix_get(invC,J,L);
                }
            }
        }
    }
}

gsl_matrix* ShapeBase::calculateCompactStressForNodalForces(double detFe, gsl_matrix* Fe, gsl_matrix* S, gsl_matrix* Stress){
    //calculating stress (stress = (detFe)^-1 Fe S Fe^T):
    //double detFe = determinant3by3Matrix(Fe);
    gsl_matrix* tmpMat1 =  gsl_matrix_calloc(nDim, nDim);
    //cout<<"detFe: "<<detFe<<endl;
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Fe, S,0.0, tmpMat1);
    gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, tmpMat1, Fe,0.0, Stress);
    gsl_matrix_scale(Stress, 1.0/detFe);
    gsl_matrix* compactStress =  gsl_matrix_calloc(6,1);
    gsl_matrix_set(compactStress,0,0,gsl_matrix_get(Stress,0,0));
    gsl_matrix_set(compactStress,1,0,gsl_matrix_get(Stress,1,1));
    gsl_matrix_set(compactStress,2,0,gsl_matrix_get(Stress,2,2));
    gsl_matrix_set(compactStress,3,0,gsl_matrix_get(Stress,0,1));
    gsl_matrix_set(compactStress,4,0,gsl_matrix_get(Stress,2,1));
    gsl_matrix_set(compactStress,5,0,gsl_matrix_get(Stress,0,2));

    gsl_matrix_free(tmpMat1);
    return compactStress;
}

void ShapeBase::calculateInvJShFuncDerSWithFe(gsl_matrix* currFe, gsl_matrix* InvDXde, gsl_matrix* ShapeFuncDerStack, gsl_matrix *invJShFuncDerSWithFe){
	//I want InvJe, normally J InvDXde = F, I can get Je from
	// Je InvDXde = Fe
	// but I can also get InvJe directly from:
	// InvJe Je InvdXde = InvJe Fe => I InvdXde = InvJe Fe => InvdXde InvFe = InvJe I => InvJe = InvdXde InvFe
	gsl_matrix* tmpFeforInversion =  gsl_matrix_calloc(nDim,nDim);
	gsl_matrix* InvFe = gsl_matrix_calloc(nDim,nDim);
	gsl_matrix* InvJe = gsl_matrix_calloc(nDim,nDim);
	createMatrixCopy(tmpFeforInversion,currFe);
	bool inverted = InvertMatrix(tmpFeforInversion, InvFe);
	if (!inverted){
		cerr<<"Fe not inverted!!"<<endl;
	}
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, InvDXde, InvFe,0.0, InvJe);

	int dim2 = nDim*nDim;
	//Generating the inverse Jacobian(elastic) stack:
	gsl_matrix* InvJacobianElasticStack =  gsl_matrix_calloc(dim2,dim2);
	for (int i =0; i<nDim; i++){
		for (int m=0; m<nDim; ++m){
			for (int n=0; n<3; ++n){
				gsl_matrix_set(InvJacobianElasticStack,i*nDim+m,i*nDim+n,gsl_matrix_get(InvJe,n,m));
			}
		}
	}

	//I am calculating this for k calculation, in case there is growth. Under conditions that there is no growth, this function is not necessary,
	//the values of invJShFuncDerSWithF and  invJShFuncDerS will be equal
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, InvJacobianElasticStack, ShapeFuncDerStack,0.0, invJShFuncDerSWithFe);
	gsl_matrix_free(tmpFeforInversion);
	gsl_matrix_free(InvFe);
	gsl_matrix_free(InvJe);
	gsl_matrix_free(InvJacobianElasticStack);
}

gsl_matrix* ShapeBase::calculateBTforNodalForces(gsl_matrix* InvJacobianStack, gsl_matrix* ShapeFuncDerStack, gsl_matrix* B, gsl_matrix *invJShFuncDerS){
     //calculating B:
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, InvJacobianStack, ShapeFuncDerStack,0.0, invJShFuncDerS);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, CoeffMat, invJShFuncDerS,0.0, B);
    //generating B^T:
    gsl_matrix* BT = gsl_matrix_alloc(nNodes*nDim,6);
    gsl_matrix_transpose_memcpy(BT,B);
    return BT;
}

gsl_matrix* ShapeBase::calculateInverseJacobianStackForNodalForces(gsl_matrix* Jacobian){
    int dim2 = nDim*nDim;
    //invrting the Jacobian:
    gsl_matrix* tmpJacobianForInversion =  gsl_matrix_calloc(nDim,nDim);
    gsl_matrix* InvJacobian = gsl_matrix_calloc(nDim,nDim);
    createMatrixCopy(tmpJacobianForInversion,Jacobian);
    bool inverted = InvertMatrix(tmpJacobianForInversion, InvJacobian);
    if (!inverted){
        cerr<<"Jacobian not inverted!!"<<endl;
    }
    //displayMatrix(Jacobian,"Jacobian");
    //Generating the inverse Jacobian stack:
    gsl_matrix* InvJacobianStack =  gsl_matrix_calloc(dim2,dim2);
    for (int i =0; i<nDim; i++){
        for (int m=0; m<nDim; ++m){
            for (int n=0; n<3; ++n){
                gsl_matrix_set(InvJacobianStack,i*nDim+m,i*nDim+n,gsl_matrix_get(InvJacobian,n,m));
            }
        }
    }
    gsl_matrix_free(tmpJacobianForInversion);
    gsl_matrix_free(InvJacobian);
    return InvJacobianStack;
}

gsl_matrix* ShapeBase::calculateVelocityGradientTensor(gsl_matrix* B, gsl_matrix* displacementPerDt){
	/**
	 * Inputs:
	 * -# The elemental B matrix (6 , ShapeBase#nDim x ShapeBase#nNodes).
	 * -# The displacement of all nodes of the system, divided by the time
	 * step (ShapeBase#nDim x Simulation#nNodes, 1).
	 *
	 * Output:
	 * -# Velocity gradient tensor in Voigt notation (6, 1).
	 *
	 * This function calculates the velocity gradient tensor from elemental B matrix and elemental displacement.
	 * The elemental B matrix is composed of a stack of B matrices for each node of the element:
	 * 	\f{eqnarray*}{
        	\textbf{B}  &=& \left[ \left[ \textbf{B}_{0} \right] \left[ \textbf{B}_{1} \right] ... \left[ \textbf{B}_{nNode}\right] \right]\\
						&=& \left[
		\begin{bmatrix}
			\partial_x N_0 	& 0 				& 0	\\
			0 				& \partial_y N_0  	& 0	\\
			0 				& 0 				& \partial_z N_0 \\
			\partial_y N_0 	& \partial_x N_0 	& 0\\
			\partial_z N_0 	& 0 				& \partial_x N_0 \\
			0 				& \partial_z N_0 	& \partial_y N_0
		\end{bmatrix}

		\begin{bmatrix}
			\partial_x N_1 	& 0 				& 0	\\
			0 				& \partial_y N_1  	& 0	\\
			0 				& 0 				& \partial_z N_1 \\
			\partial_y N_1 	& \partial_x N_1 	& 0\\
			\partial_z N_1 	& 0 				& \partial_x N_1 \\
			0 				& \partial_z N_1 	& \partial_y N_1
		\end{bmatrix}

		...

		\begin{bmatrix}
			\partial_x N_{nNode} 	& 0 					& 0	\\
			0 						& \partial_y N_{nNode}  & 0	\\
			0 						& 0 					& \partial_z N_{nNode} \\
			\partial_y N_{nNode} 	& \partial_x N_{nNode} 	& 0\\
			\partial_z N_{nNode} 	& 0 					& \partial_x N_{nNode} \\
			0 						& \partial_z N_{nNode} 	& \partial_y N_{nNode}
		\end{bmatrix}
		\right]
		\f}
	 * The elemental displacement matrix is extracted from the system displacement matrix via
	 * the function ShapeBase#constructElementalDisplacementMatrix. The displacement is calculated
	 * as the displacement of a node from its position at the end of last time step, \f$ u_{n}\f$ to the position
	 * at the current Newton-Raphson iteration \f$ u_{k}\f$. With the velocities (displacement per time step),
	 * and the \f$\textbf{B}\f$ matrix, velocity gradient tensor can be calculated through:
	 * \f{eqnarray*}{
	 	 	 	 \boldsymbol{l} & = \boldsymbol{B}  \boldsymbol{v_{n+1}}\nonumber \\
								& = \boldsymbol{B} \frac{{u_{n+1}^{k} - u_{n}}} {\delta t}.
		\f}
	 *
	 * Procedure:
	 * - construct the ElementalDisplacementMatrix.
	 */
	gsl_matrix* elementalDisplacementMatrix = constructElementalDisplacementMatrix(displacementPerDt);
	/**
	 * - Allocate the velocity gradient tensor in Voigt notation.
	 */
	gsl_matrix* l =  gsl_matrix_calloc(6,1);
	/**
	 * - calculate velocity gradient tensor.
	 */
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, B, elementalDisplacementMatrix,0.0, l);
    /**
	 * - free allocated memory.
	 */
	gsl_matrix_free(elementalDisplacementMatrix);
	//displayMatrix(l,"l_forForceCalc");
    /**
	 * - return velocity gradient tensor.
	 */
	return l;
}

gsl_matrix* ShapeBase::calculateRateOfDeformationTensor(gsl_matrix* l){
	/**
	 * Inputs:
	 * -# The velocity gradient tensor given in Voigt notation (6 x 1).
	 *
	 * Output:
	 * -# Rate of deformation tensor (3 x 3).
	 *
	 * This function primarily rearranges the elements of the
	 * velocity gradient tensor given in Voigt notation, to from the rate of deformation tensor
	 *
	 * Procedure:
	 * - Allocate the memory for rate of deformation tensor
	 */
	gsl_matrix* d =  gsl_matrix_calloc(3,3);
	/**
	 * - Write the terms of velocity gradient tensor into rate of deformation tensor
	 */
	gsl_matrix_set(d,0,0,gsl_matrix_get(l,0,0) );
	gsl_matrix_set(d,1,1,gsl_matrix_get(l,1,0) );
	gsl_matrix_set(d,2,2,gsl_matrix_get(l,2,0) );
	gsl_matrix_set(d,0,1,0.5*gsl_matrix_get(l,3,0));
	gsl_matrix_set(d,2,1,0.5*gsl_matrix_get(l,4,0));
	gsl_matrix_set(d,0,2,0.5*gsl_matrix_get(l,5,0));
	gsl_matrix_set(d,1,0,0.5*gsl_matrix_get(l,3,0));
	gsl_matrix_set(d,1,2,0.5*gsl_matrix_get(l,4,0));
	gsl_matrix_set(d,2,0,0.5*gsl_matrix_get(l,5,0));
	/**
	 * - Return rate of deformation tensor
	 */
	return d;
}

void ShapeBase::calculateViscousStress(gsl_matrix* d, gsl_matrix* viscousStress){
	/**
	 * Inputs:
	 * -# The rate of deformation matrix (ShapeBase#nDim x ShapeBase#nDim).
	 * -# the viscous stress of current Gauss point, the result will be written on this matrix
	 *
	 * This function will calculate the internal viscous stress of the element using rate of deformation matrix
	 * and ShapeBase#internalViscosity, \f$\eta\f$ , via:
	 * \f[\sigma^{v} = \eta  \textbf{d} \f]
	 *
	 * Procedure:
	 * - Copy rate of deformation tensor over to viscous stress tensor.
	 *
	 */
	createMatrixCopy(viscousStress, d);
	/**
	 *  - Scale with the internal viscosity to obtain viscous stress tensor
	 *
	 */
	gsl_matrix_scale(viscousStress,internalViscosity);
}

gsl_matrix* ShapeBase::constructElementalDisplacementMatrix(gsl_matrix* displacement){
	/**
	 * Inputs:
	 * -# The displacement matrix (ShapeBase#nDim x Simulation#nNodes, 1).
	 *
	 * This function calculates the elemental displacement matrix from the
	 * displacement matrix of the whole system, given as input. In
	 * current usage, under normal circumstances, the input matrix is displacement
	 * divided by time step. The displacement is calculated by Simulation#calculateDisplacementMatrix
	 * Both matrices, the displacement matrix of the whole system and the
	 * elemental displacement matrix are in vector form:
	 *
   	    \f$ displacement =
   	    \begin{bmatrix}
			\Delta x_{0}\\
			\Delta y_{0}\\
			\Delta z_{0}\\
			... ,\\
			\Delta x_{N}\\
			\Delta y_{N}\\
			\Delta z_{N}
			\end{bmatrix}
   	    \f$

	 */
	gsl_matrix* elementalDisplacementMatrix = gsl_matrix_calloc(nDim*nNodes,1);
	for (int i=0; i<nNodes; ++i){
		int index = NodeIds[i];
		for (int j=0; j<nDim; ++j){
			double value = gsl_matrix_get(displacement,index*nDim+j,0);
			gsl_matrix_set(elementalDisplacementMatrix,i*nDim+j,0,value);
		}
	}
	return elementalDisplacementMatrix;
}

void ShapeBase::calculateViscousForces(gsl_matrix*  gv, gsl_matrix*  BTdetFdetdXde, gsl_matrix* viscousStress){
	/**
	 * Inputs:
	 * -# Elemental matrix for internal viscous forces (ShapeBase#nDim x ShapeBase#nNodes, 1)
	 * the resulting forces will be written on this matrix
	 * -# Transpose of elemental B matrix, multiplied by the determinant
	 * of the deformation gradient, \f$\textbf{F}\f$, and the determinant of \f$ \delta \textbf{X}/\delta \boldsymbol{\xi}\f$
	 * -#  The viscous stresses calculated in ShapeBase#calculateViscousStress
	 *
	 * This function will calculate the elemental viscous forces from viscous stress, via:
	 * 	\f{eqnarray*}{
        	\textbf{g}^v &=& \int_{V} \textbf{B}^{T} \sigma^{v} dV \\
          	  	  	  	 &=& det(\textbf{F})det\left( \frac{\delta \textbf{X} }{\delta \boldsymbol{\xi}} \right) \textbf{B}^{T} \sigma^{v}
		\f}
	 * Procedure:
	 * - Allocate the memory for stress in Voigt notation
	 * */
	gsl_matrix* compactStress =  gsl_matrix_calloc(6,1);
	/**
	 * - Write stress in Voigt notation
	 */
    gsl_matrix_set(compactStress,0,0,gsl_matrix_get(viscousStress,0,0));
    gsl_matrix_set(compactStress,1,0,gsl_matrix_get(viscousStress,1,1));
    gsl_matrix_set(compactStress,2,0,gsl_matrix_get(viscousStress,2,2));
    gsl_matrix_set(compactStress,3,0,gsl_matrix_get(viscousStress,0,1));
    gsl_matrix_set(compactStress,4,0,gsl_matrix_get(viscousStress,2,1));
    gsl_matrix_set(compactStress,5,0,gsl_matrix_get(viscousStress,0,2));
	/**
	 * - Calculate nodal forces
	 */
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BTdetFdetdXde, compactStress,0.0, gv);
	/**
	 * - Free memory
	 */
    gsl_matrix_free(compactStress);
}

void ShapeBase::calculateImplicitKElastic(){
    //cout<<"calculating implicit K elastic for element: "<<Id<<endl;
    int dim = nDim;
    int n = nNodes;
    if (IsAblated){
    	gsl_matrix_set_zero(TriPointKe);
    }
    else{
    	//calculating Kelastic in a 6 point gaussian:
		gsl_matrix* currK = gsl_matrix_calloc(dim*n,dim*n);
		gsl_matrix_set_zero(TriPointKe);
		for (int iter =0; iter<numberOfGaussPoints;++iter){
			gsl_matrix_set_zero(currK);
			calculateElasticKIntegral1(currK,iter);
			calculateElasticKIntegral2(currK,iter);
			gsl_matrix_scale(currK,gaussWeights[iter]);
			gsl_matrix_add(TriPointKe, currK);
		}
	    gsl_matrix_free(currK);
    }
}

void ShapeBase::calculateImplicitKViscous(gsl_matrix* displacementPerDt, double dt){
	//This is a function called over each element
	//First function is the function to calculate the sum of all integrals,
	//each integral calculation will be listed below, individually.
	//The inputs are:
	//1) the displacement per dt for all nodes (uk-un)/dt,
	//2) the time step
	//the object has access to all the necessary matrices to carry out the calculation.
	int dim = nDim;
	int n = nNodes;
	if (IsAblated || internalViscosity == 0){
	    //This is for efficiency, I do not calculate if the
	    //current element is laser ablated
		gsl_matrix_set_zero(TriPointKv);
	}
	else{
		//calculating Kviscous in a 3 point gaussian:
		//assign a temporary matix, 18 x 18 for a prism (6 nodes, 3D).
		gsl_matrix* currK = gsl_matrix_calloc(dim*n,dim*n);
		//set the temporary matrix to zero.
		gsl_matrix_set_zero(TriPointKv);
		//the weights of the gauss points

		//define the matrices for velocity gradient, and the
		//term in parentheses in calculation of the first integral.
		//These will be calculated once per Gauss point for the element.
		gsl_matrix* velocityGradient = gsl_matrix_calloc(dim,dim);
		gsl_matrix* paranthesisTermForKv1 =  gsl_matrix_calloc(dim,dim);
	    //loop over Gauss points
		for (int iter =0; iter<numberOfGaussPoints;++iter){
			//set the temporary matrix to zero.
			gsl_matrix_set_zero(currK);
			//set the velocity gradient to zero.
			gsl_matrix_set_zero(velocityGradient);
			//set the parentheses term for first integral to zero.
			gsl_matrix_set_zero(paranthesisTermForKv1);
			//calculate the velocity gradient:
			calculateVelocityGradient(velocityGradient, displacementPerDt, iter);
			//calculate ( I / dt - velocityGradient) for first term:
			//set the term to identity:
			gsl_matrix_set_identity(paranthesisTermForKv1);
			//divide the term by dt to obtain I/dt:
			gsl_matrix_scale(paranthesisTermForKv1,1.0/dt);
			//substract velocity gradient to obtain ( I / dt - velocityGradient):
			gsl_matrix_sub(paranthesisTermForKv1,velocityGradient);
			/*if (Id == 0){
				displayMatrix(velocityGradient,"velocityGradient");
				displayMatrix(paranthesisTermForKv1,"paranthesisTermForKv1");
				displayMatrix(viscousStress[iter],"viscousStress");
			}*/
			//calculate the first integral:
			calculateViscousKIntegral1(currK, paranthesisTermForKv1, iter);
			//calculate the second integral:
			calculateViscousKIntegral2(currK, iter);
		    //scaling the resulting temporary matrix with Gauss point weight.
			gsl_matrix_scale(currK,gaussWeights[iter]);
		    //Adding the temporary matrix to the elemental Kviscous.
			gsl_matrix_add(TriPointKv, currK);
		}
		//free the memory allocated in this function
		gsl_matrix_free(currK);
		gsl_matrix_free(velocityGradient);
		gsl_matrix_free(paranthesisTermForKv1);
	}
};

void ShapeBase::writeKelasticToMainKatrix(gsl_matrix* K){
	//if (Id == 0) {
	//	displayMatrix(TriPointKe,"TriPointKeWhileWriting");
	//}

    for (int a=0; a<nNodes; ++a){
        for (int b=0; b<nNodes; ++b){
            int NodeId1 = NodeIds[a];
            int NodeId2 = NodeIds[b];
            NodeId1 *= nDim;
            NodeId2 *= nDim;
            for (int i=0; i<nDim; ++i){
                for (int j=0; j<nDim; ++j){
                    double valueij = gsl_matrix_get(K,NodeId1+i,NodeId2+j);
					valueij	+= gsl_matrix_get(TriPointKe,a*nDim+i,b*nDim+j);
                    gsl_matrix_set(K,NodeId1+i,NodeId2+j,valueij);
                	if (isnan(valueij)){
                		double tripointvalue = gsl_matrix_get(TriPointKe,a*nDim+i,b*nDim+j);
                		cout<<" element: "<<Id<<" K elastic dimention: "<<i<<" "<<j<<" is NaN after addition: "<<valueij<<" tri point value: "<<tripointvalue<<endl;
                	}
                }
            }
        }
    }
}

void ShapeBase::writeKviscousToMainKatrix(gsl_matrix* K){
	//if (Id == 0) {
	//	displayMatrix(TriPointKv,"TriPointKvWhileWriting");
	//}
	if (internalViscosity != 0){
		for (int a=0; a<nNodes; ++a){
			for (int b=0; b<nNodes; ++b){
				int NodeId1 = NodeIds[a];
				int NodeId2 = NodeIds[b];
				NodeId1 *= nDim;
				NodeId2 *= nDim;
				for (int i=0; i<nDim; ++i){
					for (int j=0; j<nDim; ++j){
						double valueij = gsl_matrix_get(K,NodeId1+i,NodeId2+j);
						valueij	+= gsl_matrix_get(TriPointKv,a*nDim+i,b*nDim+j);
						gsl_matrix_set(K,NodeId1+i,NodeId2+j,valueij);
	                	if (isnan(valueij)){
	                		double tripointvalue = gsl_matrix_get(TriPointKv,a*nDim+i,b*nDim+j);
	                		cout<<" element: "<<Id<<" K viscous dimention: "<<i<<" "<<j<<" is NaN after addition: "<<valueij<<" tri point value: "<<tripointvalue<<endl;
	                	}
					}
				}
			}
		}
    }
}

void	ShapeBase::calculateForceFromStress(int nodeId, gsl_matrix* Externalstress, gsl_matrix *ExternalNodalForces){
    gsl_matrix_set_zero(ExternalNodalForces);
    int nodeNo = 0;
    for (int i=0; i<nNodes; i++){
        if (NodeIds[i] == nodeId){
            nodeNo = i;
            break;
        }
    }
    for (int pointNo = 0; pointNo<3; pointNo++){
        gsl_matrix* BaT = gsl_matrix_calloc(nDim,6);
        gsl_matrix* Bb = gsl_matrix_calloc(6,nDim);
        gsl_matrix* B = Bmatrices[pointNo];
        consturctBaTBb(B, BaT,Bb,nodeNo,0);
        gsl_matrix* NodeForces = gsl_matrix_calloc(3,1);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BaT, Externalstress,0.0, NodeForces);
        gsl_matrix_scale(NodeForces,1.0/3.0);
        gsl_matrix_scale(NodeForces,detFs[pointNo]);
        gsl_matrix_scale(NodeForces,detdXdes[pointNo]);
        gsl_matrix_add(ExternalNodalForces,NodeForces);
        gsl_matrix_free(BaT);
        gsl_matrix_free(Bb);
        gsl_matrix_free(NodeForces);
    }
    //displayMatrix(ExternalNodalForces,"ExternalNodalForces");
}
void	ShapeBase::calculateElasticKIntegral1(gsl_matrix* currElementalK,int pointNo){
    gsl_matrix* invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    double detF = detFs[pointNo];
    double detdXde = detdXdes[pointNo];
    gsl_matrix* Fe = FeMatrices[pointNo];
	double detFe = determinant3by3Matrix(Fe);
    //finished calculating 4th order tensor D
	/*if(isnan(detF) || detF == 0){
		cout<<"element "<<Id<<" detF at 1st integration is nan/zero "<<detF<<endl;
	}
	if(isnan(detFe) || detFe == 0){
		cout<<"element "<<Id<<" detFe at 1st integration is nan/zero "<<detFe<<endl;
	}
	for (int i=0;i<3;++i){
		for (int j=0;j<3;++j){
			double value = gsl_matrix_get(Fe,i,j);
			if(isnan(value)){
				cout<<"element "<<Id<<" Fe at 1st integration is nan at dim: "<<i<<" "<<j<<endl;
			}
		}
	}
	 */
    for (int a =0; a<nNodes; ++a){
        for (int b=0; b<nNodes; ++b){
            gsl_matrix* Keab = gsl_matrix_calloc(3,3);
            double DNa[3] = {0.0,0.0,0.0};
            double DNb[3] = {0.0,0.0,0.0};

            for (int i=0;i<nDim;++i){
                // original version: DNa[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*a);
                // original version: DNb[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*b);
            	DNa[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*a);
                DNb[i] = gsl_matrix_get(invJShFuncDerS,i,nDim*b);
            }
            //cout<<" DNb from Fe: "<<DNb[0]<<" "<<DNb[1]<<" "<<DNb[2]<<" DNb from F: "<<DNbold[0]<<" "<<DNbold[1]<<" "<<DNbold[2]<<endl;
            //writing Kab:
            for (int i = 0 ; i<nDim; ++i){
                for (int k=0; k<nDim; ++k){
                    double value = 0;
                    //the sum over j,l,I,J,K,L, to get Kab(i,k):
                    for (int j = 0; j<nDim; ++j){
                        for (int l=0; l<nDim; ++l){
                            for (int I=0; I<nDim; ++I){
                                for (int J=0; J<nDim; ++J){
                                    for (int K=0; K<nDim; ++K){
                                        for (int L=0; L<nDim; ++L){
                                            value += (gsl_matrix_get(Fe,i,I)*gsl_matrix_get(Fe,j,J)*gsl_matrix_get(Fe,k,K)*gsl_matrix_get(Fe,l,L)*D81[pointNo][I][J][K][L]*DNb[l]*DNa[j]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    value *= detF*detdXde;
                    value /= detFe;
                    value += gsl_matrix_get(Keab,i,k);
                    gsl_matrix_set(Keab,i,k,value);
                }
            }
            //now I have Kab for current gauss point, I need to write in into currK:
            for (int i=0; i<nDim; ++i){
                for (int j=0; j<nDim; ++j){
                    double value = gsl_matrix_get(currElementalK,a*nDim+i, b*nDim+j);
                    value += gsl_matrix_get(Keab,i, j);
                    gsl_matrix_set(currElementalK,a*nDim+i, b*nDim+j,value);
                }
            }
            gsl_matrix_free(Keab);
        }
    }
}


void	ShapeBase::calculateElasticKIntegral2(gsl_matrix* currElementalK,int pointNo){
    gsl_matrix* invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    gsl_matrix* Stress = elasticStress[pointNo];
    double detF = detFs[pointNo];
    double detdXde = detdXdes[pointNo];

    gsl_matrix* DNaT = gsl_matrix_calloc(1,nDim);
    gsl_matrix* DNb = gsl_matrix_calloc(nDim,1);
    gsl_matrix* Keab2 = gsl_matrix_calloc(1,1);
    for (int a =0; a<nNodes; ++a){
        for (int b=0; b<nNodes; ++b){
            for (int i=0;i<nDim;++i){
                gsl_matrix_set(DNaT,0,i,gsl_matrix_get(invJShFuncDerS,i,nDim*a));
                gsl_matrix_set(DNb,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*b));
            }
            gsl_matrix* tmp1 = gsl_matrix_calloc(1,nDim);
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, DNaT, Stress,0.0, tmp1);
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, tmp1, DNb,0.0, Keab2);
            double value = gsl_matrix_get(Keab2,0,0)*detF*detdXde;
            for (int i=0; i<nDim; ++i){
                int index1 = a*nDim+i;
                int index2 = b*nDim+i;
                double addedValue = gsl_matrix_get(currElementalK,index1,index2) + value;//adding the calculated value to current K matirx
                gsl_matrix_set(currElementalK,index1,index2,addedValue);
            }
            gsl_matrix_free(tmp1);
        }
    }
    gsl_matrix_free(DNaT);
    gsl_matrix_free(DNb);
    gsl_matrix_free(Keab2);
}

void ShapeBase::calculateVelocityGradient( gsl_matrix* velocityGradient, gsl_matrix* displacementPerDt, int pointNo){
	gsl_matrix* invJShFuncDerS = invJShapeFuncDerStack[pointNo];
	for (int c=0; c<nNodes; ++c){
		gsl_matrix* delVc = gsl_matrix_calloc(nDim,nDim);
		//get \DelNc^T
		gsl_matrix* DNc = gsl_matrix_calloc(nDim,1);
		for (int i=0;i<nDim;++i){
			gsl_matrix_set(DNc,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*c));
		}
		//calculate velocity of node c:
		gsl_matrix* vc = gsl_matrix_calloc(nDim,1);
		int id = NodeIds[c];
		gsl_matrix_set(vc,0,0,gsl_matrix_get(displacementPerDt,id*nDim,0));
		gsl_matrix_set(vc,1,0,gsl_matrix_get(displacementPerDt,id*nDim+1,0));
		gsl_matrix_set(vc,2,0,gsl_matrix_get(displacementPerDt,id*nDim+2,0));
		//calculate vc *DNc^T:
		calculateOuterProduct(vc,DNc,delVc);
		//add the nodal calculation to velocity gradient:
		gsl_matrix_add(velocityGradient, delVc);
		//free memory allocated in this loop:
		gsl_matrix_free(delVc);
		gsl_matrix_free(DNc);
		gsl_matrix_free(vc);
	}
	/*if (Id == 0 ){
		displayMatrix(velocityGradient,"DelV");
	}*/
}

void ShapeBase::calculateViscousKIntegral1(gsl_matrix* currElementalK, gsl_matrix* paranthesisTermForKv1, int pointNo){
	gsl_matrix* BaT = gsl_matrix_calloc(nDim,6);
    gsl_matrix* Bb = gsl_matrix_calloc(6,nDim);
    gsl_matrix* BaTBb = gsl_matrix_calloc(nDim,nDim);
    gsl_matrix* KvabInt1 = gsl_matrix_calloc(nDim,nDim); //First integral in calculation of Kv for nodes a and b
    gsl_matrix* B = Bmatrices[pointNo];
    for (int a=0;a<nNodes;++a){
    	for (int b=0; b<nNodes; ++b){
    		consturctBaTBb(B, BaT,Bb,a,b);
    		//Bb matrix should have the last three rows as 0.5, this matrix stems from the
    		//"d" matrix, as opposed to viscous stresses, therefore the definition
    		//should have 0.5 on off diagonal terms, therefore the last three rows of Bb.
    		for (int i=3; i<6; ++i){
    			for (int j=0; j<nDim; ++j){
    				double value = gsl_matrix_get(Bb, i,j);
    				gsl_matrix_set(Bb,i,j,0.5*value);
    			}
    		}
    		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BaT, Bb,0.0, BaTBb);
    		//if (a == 3 && b ==3){
    		//	displayMatrix(BaTBb,"BaTBb_33");
    		//}
    		//the paranthesis term is: ( I / dt - velocityGradient)
    		//calculate all multiplication: BaT*Bb* (I/dt - \Del v):
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, BaTBb, paranthesisTermForKv1,0.0, KvabInt1);
			//if (Id == 0){// && a == 3 && b ==3){
			//	//displayMatrix(KvabInt1,"KvabInt1_beforeVolumeIntegraiton_33");
			//	cout<<"a,b: "<<a<<", "<<b<<" detF: "<<detFs[pointNo]<<" detdXde: "<<detdXdes[pointNo]<<" dV: "<<detFs[pointNo]*detdXdes[pointNo]<<" element volume "<<ReferenceShape->Volume<<endl;
			//	displayMatrix(BaTBb,"BaTBb");
			//}
			//volume integration:
    	    gsl_matrix_scale(KvabInt1,detFs[pointNo]);
    	    gsl_matrix_scale(KvabInt1,detdXdes[pointNo]);
    	    //scaling by viscosity:
    	    gsl_matrix_scale(KvabInt1,internalViscosity);
    	    //if (Id == 0 && a == 0 && b == 0){
    	    //	displayMatrix(KvabInt1,"KvabInt1_afterVolumeIntegraiton_00Element0");
    	    //	cout<<"dV: "<<detFs[pointNo]*detdXdes[pointNo];
    	    //	displayMatrix(BaTBb,"BaTBb");
    	    //}
    	    //Now KabvInt1 is a 3x3 (dim x dim) matrix, while currK is 18 x 18 (dim*n_node x dim*n_nodes)
    	    //currK is composed of 3x3 Kab matrices placed into the blocks (a,b), with indexing from zero,
    	    //for a = 1 and b=2, Kab will be placed in the 3x3 block covering indices (3,4,5) x (6,7,8) on K matrix.
			/*if (a == 3 && b == 3){
				cout<<"a: "<<a<<" b: "<<b<<endl;
				displayMatrix(KvabInt1,"KvabInt1");
			}*/
    	    for (int i=0; i<nDim; ++i){
        	    for (int j=0; j<nDim; ++j){
        	    	int index2 = a*nDim+i;
        	    	int index1 = b*nDim+j;
        	    	double value = gsl_matrix_get(KvabInt1,i,j);
        	    	double addedValue = gsl_matrix_get(currElementalK,index1,index2) + value; //adding the calculated value to current K matrix
        	    	gsl_matrix_set(currElementalK,index1,index2,addedValue);
        	    }
			}
    	}
    }
    //if (Id ==0 ){displayMatrix(KvabInt1,"KvabInt1");}
    gsl_matrix_free(BaT);
    gsl_matrix_free(Bb);
    gsl_matrix_free(BaTBb);
    gsl_matrix_free(KvabInt1);
}

void ShapeBase::calculateViscousKIntegral2(gsl_matrix* currElementalK,int pointNo){
    gsl_matrix* invJShFuncDerS = invJShapeFuncDerStack[pointNo];
    gsl_matrix* Stress = viscousStress[pointNo];
    gsl_matrix* DNa = gsl_matrix_calloc(nDim,1);
    gsl_matrix* DNb = gsl_matrix_calloc(nDim,1);
    gsl_matrix* KvabInt2 = gsl_matrix_calloc(nDim,nDim);
    for (int a =0; a<nNodes; ++a){
        for (int b=0; b<nNodes; ++b){
            for (int i=0;i<nDim;++i){
                gsl_matrix_set(DNa,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*a));
                gsl_matrix_set(DNb,i,0,gsl_matrix_get(invJShFuncDerS,i,nDim*b));
            }
            gsl_matrix* DNaDNbOuterProduct = gsl_matrix_calloc(nDim, nDim);
            gsl_matrix* DNbDNaOuterProduct = gsl_matrix_calloc(nDim, nDim);
            calculateOuterProduct(DNa, DNb, DNaDNbOuterProduct);
            calculateOuterProduct(DNb, DNa, DNbDNaOuterProduct);
            gsl_matrix* paranthesisTerm = gsl_matrix_calloc(nDim, nDim);
            gsl_matrix_memcpy(paranthesisTerm,DNaDNbOuterProduct);
            gsl_matrix_sub(paranthesisTerm,DNbDNaOuterProduct);
            //gsl_matrix_add(paranthesisTerm,DNbDNaOuterProduct);
            //gsl_matrix_scale(paranthesisTerm, 0.5);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Stress, paranthesisTerm,0.0, KvabInt2);
			/*if (Id == 0){
				//displayMatrix(paranthesisTerm,"paranthesisTerm");
				displayMatrix(Stress,"Stress");
				//displayMatrix(KvabInt2,"KvabInt2-beforeVolumeIntegration");
			}*/
			gsl_matrix_scale(KvabInt2, detFs[pointNo]);
			gsl_matrix_scale(KvabInt2, detdXdes[pointNo]);
			/*if (a == 3 && b == 3){
				cout<<"a: "<<a<<" b: "<<b<<endl;
				displayMatrix(DNaDNbOuterProduct,"DNaDNbOuterProduct");
				displayMatrix(DNbDNaOuterProduct,"DNbDNaOuterProduct");
				displayMatrix(paranthesisTerm,"paranthesisTerm");
				displayMatrix(Stress,"Stress");
				cout<<"detFs[pointNo]: "<<detFs[pointNo]<<endl;
				cout<<"detdXdes[pointNo]: "<<detdXdes[pointNo]<<endl;
				displayMatrix(KvabInt2,"KvabInt2");
			}*/
			for (int i=0; i<nDim; ++i){
				for (int j=0; j<nDim; ++j){
					int index2 = a*nDim+i;
					int index1 = b*nDim+j;
					double value = gsl_matrix_get(KvabInt2,i,j);
					double addedValue = gsl_matrix_get(currElementalK,index1,index2) + value; //adding the calculated value to current K matrix
					gsl_matrix_set(currElementalK,index1,index2,addedValue);
				}
			}
            gsl_matrix_free(DNaDNbOuterProduct);
            gsl_matrix_free(DNbDNaOuterProduct);
            gsl_matrix_free(paranthesisTerm);
        }
    }
    gsl_matrix_free(DNa);
    gsl_matrix_free(DNb);
    gsl_matrix_free(KvabInt2);
	/*if (Id == 0){
		displayMatrix(DNa,"DNa");
		displayMatrix(DNaT,"DNaT");
		displayMatrix(DNb,"DNb");
		displayMatrix(DNbT,"DNbT");
		displayMatrix(DNaDNbT,"DNaDNbT");
		displayMatrix(DNbDNaT,"DNbDNaT");
		//displayMatrix(paranthesisTerm,"paranthesisTerm");
		displayMatrix(Stress,"Stress");
		displayMatrix(KvabInt2,"KvabInt2-beforeVolumeIntegration");
	}*/
}

void	ShapeBase::calculateOuterProduct(gsl_matrix* a, gsl_matrix* b, gsl_matrix* outerProduct){
	//cout<<"inside outer product"<<endl;
	int size1 = a->size1;
	int size2 = a->size2;
	if ((int) b->size2 != size2){
		cerr<<"matrix dimension mismatch in outer product calculation"<<endl;
	}
	//a outer b =  a b^T
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, a, b,0.0, outerProduct);
	//cout<<"finalised outer product"<<endl;
}

gsl_matrix*	ShapeBase::calculateSymmetricisedTensorProduct(gsl_matrix* a, gsl_matrix* b){
	int size1 = a->size1;
	int size2 = a->size2;
	if ((int) b->size1 != size1){
		cerr<<"matrix dimension mismatch in symmetricised outer product calculation"<<endl;
	}
	if ((int) b->size2 != size2){
		cerr<<"matrix dimension mismatch in symmetricised outer product calculation"<<endl;
	}
	//calculating individual outer products a x b = a bT
	gsl_matrix* abOuterProduct = gsl_matrix_calloc(size1,size1);
	calculateOuterProduct(a,b,abOuterProduct);
	//gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, a, bT,0.0, abOuterProduct);
	gsl_matrix* baOuterProduct = gsl_matrix_calloc(size1,size1);
	calculateOuterProduct(b,a,abOuterProduct);
	//gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, b, aT,0.0, baOuterProduct);
	//calculating the averaged contraction term:
	gsl_matrix* averagedContraction = gsl_matrix_calloc(size1,size1);
	gsl_matrix_add(averagedContraction,abOuterProduct);
	gsl_matrix_add(averagedContraction,baOuterProduct);
	gsl_matrix_scale(averagedContraction,0.5);
	gsl_matrix_free(abOuterProduct);
	gsl_matrix_free(baOuterProduct);
	return averagedContraction;
}

void ShapeBase::consturctBaTBb(gsl_matrix* B, gsl_matrix* BaT, gsl_matrix* Bb, int a, int b){
    for (int i=0; i<6; ++i){
        for (int j=0; j<nDim; ++j){
            //double value = gsl_matrix_get(B,i,a*dim+j);
            gsl_matrix_set(BaT,j,i,gsl_matrix_get(B,i,a*nDim+j)); //transpose of Ba
            //value  = gsl_matrix_get(B,i,b*dim+j);
            gsl_matrix_set(Bb,i,j,gsl_matrix_get(B,i,b*nDim+j)); //Bb
            //displayMatrix(Bb,"Bb");
        }
    }
}

void	ShapeBase::fillNodeNeighbourhood(vector<Node*>& Nodes){
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nNodes; ++j){
			if ( i !=j ){
				int n = Nodes[NodeIds[i]]->immediateNeigs.size();
				bool alreadyOnList = false;
				for (int k=0; k<n; ++k){
					if (NodeIds[j] == Nodes[NodeIds[i]]->immediateNeigs[k]){
						alreadyOnList = true;
						break;
					}
				}
				if (!alreadyOnList){
					Nodes[NodeIds[i]]->immediateNeigs.push_back(NodeIds[j]);
				}
			}
		}
	}
}

void	ShapeBase::updatePositions(vector<Node*>& Nodes){
	for (int i = 0; i<nNodes; ++i){
		for (int j = 0; j<nDim; ++j){
			Positions[i][j] = Nodes[NodeIds[i]]->Position[j];
		}
	}
}

void	ShapeBase::updateReferencePositionsToCurentShape(){
	for (int i=0; i<nNodes; ++i){
		for (int j=0; j<nDim; ++j){
			ReferenceShape ->Positions[i][j] = Positions[i][j];
		}
	}
}

void 	ShapeBase::setGrowthRate(double dt, double rx, double ry, double rz){
	GrowthRate[0] = exp(rx*dt);
	GrowthRate[1] = exp(ry*dt);
	GrowthRate[2] = exp(rz*dt);
	//if (Id ==0) {displayMatrix(growthIncrement, "Element0growthIncrement_initialSetting");}
}

void 	ShapeBase::setGrowthRateExpFromInput(double x, double y, double z){
	GrowthRate[0] = x;
	GrowthRate[1] = y;
	GrowthRate[2] = z;
}

void 	ShapeBase::updateGrowthIncrementFromRate(){
	gsl_matrix_set_identity(growthIncrement);
	gsl_matrix_set(growthIncrement,0,0,GrowthRate[0]);
	gsl_matrix_set(growthIncrement,1,1,GrowthRate[1]);
	gsl_matrix_set(growthIncrement,2,2,GrowthRate[2]);
}

void 	ShapeBase::cleanMyosinForce(){
	for (int i=0; i<nNodes; ++i){
		MyoForce[i][0] = 0;
		MyoForce[i][1] = 0;
		MyoForce[i][2] = 0;
	}
}

bool ShapeBase::calculateIfInsideActiveStripe(double initialPoint,double endPoint, double stripeSize1, double stripeSize2){
	//All nodes and the centre of the element should be inside the active zone:
	//starting from the centre:
	double* c;
	c = new double[3];
	c = getCentre();
	for (int i=0; i<nNodes+1; ++i){
		//getting the node position
		double x;
		double y;
		double bufferFrac;
		if (i == 0 ){
			x = c[0];
			y = c[1];
			bufferFrac = 0.0;
		}
		else{
			x = Positions[i-1][0];
			y = Positions[i-1][1];
			bufferFrac = 0.1;
		}
		//is the node inside the active region:
		bool xInActiveZone = false;
		bool yInActiveZone = false;
		bool loopComplete = false;
		double lowEnd = initialPoint - bufferFrac*stripeSize1;
		double highEnd = lowEnd + stripeSize1 + 2.0*bufferFrac*stripeSize1;
		while (!xInActiveZone && !loopComplete){
			if(stripeSize1 == 0 || highEnd>endPoint){
				highEnd = endPoint + bufferFrac*stripeSize1;
				loopComplete = true;
			}
			//cout<<" x: "<<x<<" low: "<<lowEnd<<" high: "<<highEnd<<endl;
			if (x>= lowEnd && x<=highEnd){
				xInActiveZone = true;
			}
			lowEnd += stripeSize1*2.0;
			highEnd = lowEnd + stripeSize1 + 2.0*bufferFrac*stripeSize1;
		}
		if (xInActiveZone){
			//if x is in active zone, I will move on to check y:
			loopComplete = false;
			lowEnd = initialPoint;
			highEnd = lowEnd + stripeSize2 + 2.0*bufferFrac*stripeSize2;
			while (!yInActiveZone && !loopComplete){
				if(stripeSize2 == 0 || highEnd>endPoint){
					highEnd = endPoint + bufferFrac*stripeSize2;
					loopComplete = true;
				}
				if (y>= lowEnd && y<=highEnd){
					yInActiveZone = true;
				}
				lowEnd += stripeSize2*2.0;
				highEnd = lowEnd + stripeSize2 + 2.0*bufferFrac*stripeSize2;
			}
		}
		//cout<<"Element: "<<Id<<" point ( 0 for centre, i-1 for node): "<<i<<" pos: "<<x<<" "<<y<<" is inside: "<<xInActiveZone<<" "<<yInActiveZone<<endl;
		if (!xInActiveZone || !yInActiveZone ){
			//if this node is not in the active zone, then the element is not in the active zone
			delete[] c;
			return false;
		}
	}
	//I did not return the function in any of the nodes, then all nodes must be inside the active zone:
	delete[] c;
	return true;
};

double ShapeBase::getCmyosinUniformForNode (int TissuePlacement){
	if(TissuePlacement == 1) {	//apical node
		return cMyoUniform[0];
	}
	if(TissuePlacement == 0) {	//basal node
		return cMyoUniform[1];
	}
	return 0.0;
}

double ShapeBase::getCmyosinUnipolarForNode (int TissuePlacement){
	if(TissuePlacement == 1) {
		return cMyoUnipolar[0];
	}
	if(TissuePlacement == 0) {
		return cMyoUnipolar[1];
	}
	return 0.0;
}

gsl_matrix* ShapeBase::getMyosinDirection(){
    gsl_matrix* tmpMyoDir =gsl_matrix_calloc(2, nDim);
    createMatrixCopy(tmpMyoDir,myoPolarityDir);
    return tmpMyoDir;
}

void ShapeBase::getMyosinLevels (double *cMyo){
	cMyo[0] = cMyoUniform[0];
	cMyo[1] = cMyoUniform[1];
	cMyo[2] = cMyoUnipolar[0];
	cMyo[3] = cMyoUnipolar[1];
}

void ShapeBase::getEquilibriumMyosinLevels (double* cMyoEq){
	cMyoEq[0] = cMyoUniformEq[0];
	cMyoEq[1] = cMyoUniformEq[1];
	cMyoEq[2] = cMyoUnipolarEq[0];
	cMyoEq[3] = cMyoUnipolarEq[1];
}

void ShapeBase::setMyosinLevels (double cUni0, double cUni1, double cPol0, double cPol1){
	cMyoUniform[0] = cUni0;
	cMyoUniform[1] = cUni1;
	cMyoUnipolar[0] = cPol0;
	cMyoUnipolar[1] = cPol1;
}

void ShapeBase::setEquilibriumMyosinLevels (double cUni0, double cUni1, double cPol0, double cPol1){
	cMyoUniformEq[0] = cUni0;
	cMyoUniformEq[1] = cUni1;
	cMyoUnipolarEq[0] = cPol0;
	cMyoUnipolarEq[1] = cPol1;
}

void	ShapeBase::updateUniformEquilibriumMyosinConcentration(bool isApical, double cEqUniform){
	if (isApical){
		cMyoUniformEq[0] = cEqUniform;
	}
	else{
		cMyoUniformEq[1] = cEqUniform;
	}
	//cout<<"Element: "<<Id<<" eq upon signal: "<<cMyoUniformEq[0]<<" "<<cMyoUniformEq[1]<<endl;
}

void	ShapeBase::updateUnipolarEquilibriumMyosinConcentration(bool isApical, double cEqUnipolar, double orientationX, double orientationY){
	int indice = 1;
	if (isApical){
		indice = 0;
	}
	cMyoUnipolarEq[indice] = cEqUnipolar;
	gsl_matrix_set(myoPolarityDir,indice,0,orientationX);
	gsl_matrix_set(myoPolarityDir,indice,1,orientationY);
	gsl_matrix_set(myoPolarityDir,indice,2,0.0);
}


void 	ShapeBase::calculatePrincipalStrains2D(double& e1, double &e2,  double &e3, gsl_matrix* eigenVec){
	gsl_matrix* Strain2D = gsl_matrix_calloc(2,2);
	gsl_matrix_set(Strain2D,0,0, gsl_matrix_get(Strain,0,0));
	gsl_matrix_set(Strain2D,1,1, gsl_matrix_get(Strain,1,0));
	gsl_matrix_set(Strain2D,0,1, 0.5 * gsl_matrix_get(Strain,3,0));
	gsl_matrix_set(Strain2D,1,0, 0.5 * gsl_matrix_get(Strain,3,0));
	gsl_vector* eigenValues = gsl_vector_calloc(2);
	gsl_matrix* eigenVec2D = gsl_matrix_calloc(2,2);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(2);
	gsl_eigen_symmv(Strain2D, eigenValues, eigenVec2D, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eigenValues, eigenVec2D, GSL_EIGEN_SORT_ABS_ASC);
	e1 = gsl_vector_get(eigenValues,0);
	e2 = gsl_vector_get(eigenValues,1);
	e3 = 0;
	gsl_matrix_set_identity(eigenVec);
	for (int i=0; i<2; ++i){
		for (int j=0; j<2; ++j){
			gsl_matrix_set(eigenVec,i,j,gsl_matrix_get(eigenVec2D,i,j));
		}
	}
	gsl_vector_free(eigenValues);
	gsl_matrix_free(eigenVec2D);
	gsl_matrix_free(Strain2D);
}

void 	ShapeBase::calculatePrincipalStrains3D(double& e1, double &e2,  double &e3, gsl_matrix* eigenVec){
	gsl_matrix* Strain3D = gsl_matrix_calloc(3,3);
	gsl_matrix_set(Strain3D,0,0, gsl_matrix_get(Strain,0,0));
	gsl_matrix_set(Strain3D,1,1, gsl_matrix_get(Strain,1,0));
	gsl_matrix_set(Strain3D,0,1, 0.5 * gsl_matrix_get(Strain,3,0));
	gsl_matrix_set(Strain3D,1,0, 0.5 * gsl_matrix_get(Strain,3,0));
	gsl_matrix_set(Strain3D,2,2, gsl_matrix_get(Strain,2,0));
	gsl_matrix_set(Strain3D,2,1, 0.5 * gsl_matrix_get(Strain,4,0));
	gsl_matrix_set(Strain3D,1,2, 0.5 * gsl_matrix_get(Strain,4,0));
	gsl_matrix_set(Strain3D,0,2, 0.5 * gsl_matrix_get(Strain,5,0));
	gsl_matrix_set(Strain3D,2,0, 0.5 * gsl_matrix_get(Strain,5,0));
	gsl_vector* eigenValues = gsl_vector_calloc(3);
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(3);
	gsl_eigen_symmv(Strain3D, eigenValues, eigenVec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eigenValues, eigenVec, GSL_EIGEN_SORT_ABS_ASC);
	e1 = gsl_vector_get(eigenValues,0);
	e2 = gsl_vector_get(eigenValues,1);
	e3 = gsl_vector_get(eigenValues,2);
	gsl_vector_free(eigenValues);
	gsl_matrix_free(Strain3D);
}

double  ShapeBase::calculateVolumeForInputShapeStructure(double** shapePositions, int nTriangularFaces, int** triangularFaces, double* midPoint ){
	double totalVolume = 0;
	for (int i = 0; i<nTriangularFaces; ++i){
		//calculateTetrahedraVolume:
		int node0 = triangularFaces[i][0];
		int node1 = triangularFaces[i][1];
		int node2 = triangularFaces[i][2];
		double* vec1 = new double [(const int) nDim];
		double* vec2 = new double [(const int) nDim];
		double* vecMid = new double [(const int) nDim];
		for (int j=0 ;j < nDim; ++j){
				vec1   [j] = shapePositions[node1][j] - shapePositions[node0][j];
				vec2   [j] = shapePositions[node2][j] - shapePositions[node0][j];
				vecMid [j] = midPoint[j] - shapePositions[node0][j];
		}
		double* baseVec = new double [(const int) nDim];
		crossProduct3D(vec1,vec2,baseVec);
		double normBaseVec= calculateMagnitudeVector3D (baseVec);
		double baseArea= normBaseVec/2;
		if (i == 0){
			ReferenceShape->BasalArea = baseArea;
		}
		double height = dotProduct3D(vecMid,baseVec) / normBaseVec;
		if (height <0){
			height *= (-1.0);
		}
		ReferenceShape->height = height;
		double currVolume = 1.0/3.0 *(height *baseArea);
		totalVolume += currVolume;
		delete[] vec1;
		delete[] vec2;
		delete[] vecMid;
		delete[] baseVec;
	}
	return totalVolume;
}

void 	ShapeBase::calculatePrincipalStrainAxesOnXYPlane(double& e1, double &e2, double& tet){
	//principal strains:
	//e1,e2 = (exx + eyy) /2  +- sqrt( ( (exx - eyy)/2 ) ^2 + exy ^2)
	//extension is taken to be positive, therefore the most extended axis will be e1.
	//angle of the strains (direction of e1):
	// tan (2*tetha) = (2 exy ) / ( exx - eyy )


	//remodellingPlaneRotationMatrix
	double exx, eyy, exy;
	if (tissueType == 2 && ShapeType == 1){ //the matrix is only calculated for prisms of lateral tissue type
		gsl_matrix* rotatedStrain = gsl_matrix_calloc(3,3);
		gsl_matrix_set(rotatedStrain,0,0, gsl_matrix_get(Strain,0,0));
		gsl_matrix_set(rotatedStrain,1,1, gsl_matrix_get(Strain,1,0));
		gsl_matrix_set(rotatedStrain,2,2, gsl_matrix_get(Strain,2,0));
		gsl_matrix_set(rotatedStrain,0,1, 0.5 * gsl_matrix_get(Strain,3,0));
		gsl_matrix_set(rotatedStrain,1,0, 0.5 * gsl_matrix_get(Strain,3,0));
		gsl_matrix_set(rotatedStrain,2,1, 0.5 * gsl_matrix_get(Strain,4,0));
		gsl_matrix_set(rotatedStrain,1,2, 0.5 * gsl_matrix_get(Strain,4,0));
		gsl_matrix_set(rotatedStrain,0,2, 0.5 * gsl_matrix_get(Strain,5,0));
		gsl_matrix_set(rotatedStrain,2,0, 0.5 * gsl_matrix_get(Strain,5,0));

		gsl_matrix* tmp = gsl_matrix_calloc(3,3);
		//R^T * rotatedStrain
		gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0, remodellingPlaneRotationMatrix, rotatedStrain, 0.0, tmp);
		//R^T * rotatedStrain * R
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, tmp, remodellingPlaneRotationMatrix, 0.0, rotatedStrain);
		exx = gsl_matrix_get(rotatedStrain,0,0);
		eyy = gsl_matrix_get(rotatedStrain,1,1);
		exy = gsl_matrix_get(rotatedStrain,	0,1);
		gsl_matrix_free(rotatedStrain);
		gsl_matrix_free(tmp);

	}
	else{
		exx = gsl_matrix_get(Strain,0,0);
		eyy = gsl_matrix_get(Strain,1,0);
		exy = gsl_matrix_get(Strain,3,0)/2.0;
	}

	double difference = (exx - eyy)/2.0;
	//double sumTerm = (exx + eyy) /2.0 ;
	//double sqrootTerm = pow ( difference*difference +  exy*exy, 0.5);
	//e1 = sumTerm + sqrootTerm;
	//e2 = sumTerm - sqrootTerm;
	//double tan2Tet = exy/difference;
	tet = atan2(exy,difference)/2.0;
	//tet = atan(tan2Tet)/2;
	//I can calculate e1 and e2, but I will not know which direction they
	//correspond to, I should instead, calculate the rotation from a quaternian,
	//and therefore obtain the main strain direction:
	double c = cos(tet);
	double s = sin(tet);
	gsl_matrix* Rot = gsl_matrix_calloc(2,2);
	gsl_matrix* RotT = gsl_matrix_calloc(2,2);
	gsl_matrix* tmp = gsl_matrix_calloc(2,2);
	gsl_matrix* newStrain = gsl_matrix_calloc(2,2);
	gsl_matrix_set(Rot,0,0,c);
	gsl_matrix_set(Rot,0,1,s);
	gsl_matrix_set(Rot,1,0,(-1.0)*s);
	gsl_matrix_set(Rot,1,1,c);
	gsl_matrix_set(RotT,0,0,c);
	gsl_matrix_set(RotT,1,0,s);
	gsl_matrix_set(RotT,0,1,(-1.0)*s);
	gsl_matrix_set(RotT,1,1,c);
	gsl_matrix_set(newStrain,0,0,exx);
	gsl_matrix_set(newStrain,0,1,exy);
	gsl_matrix_set(newStrain,1,0,exy);
	gsl_matrix_set(newStrain,1,1,eyy);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, RotT, newStrain, 0.0, tmp);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, tmp, Rot, 0.0, newStrain);
	e1 = gsl_matrix_get(newStrain,0,0);
	e2 = gsl_matrix_get(newStrain,1,1);
	//cout<<"Id: "<<Id<<" tan2Tet "<<tan2Tet<<" 2Tet "<<atan2(exy,difference)<<" tet: "<<tet;
	if (e2>e1){
		//the main strain is in the direction perpendicular to the calculated angle.
		double tmp = e1;
		e1 = e2;
		e2 = tmp;
		tet = tet + M_PI/2.0;
	}
	//cout<<"after correction tet: "<<tet<<endl;
	//displayMatrix(Strain, "strain");
	//displayMatrix(newStrain, "newStrain");
	gsl_matrix_free(Rot);
	gsl_matrix_free(RotT);
	gsl_matrix_free(newStrain);
	gsl_matrix_free(tmp);

}

bool	ShapeBase::checkIfXYPlaneStrainAboveThreshold(double thres){
	double exx = gsl_matrix_get(Strain,0,0);
	double eyy = gsl_matrix_get(Strain,1,0);
	double exy = gsl_matrix_get(Strain,3,0)/2.0;
	double thresSquared = thres*thres;
	double diff= exx - eyy;
	if (diff*diff > thresSquared){
		//cout<<"Element: "<<Id<<" exx, eyy, exy: "<<exx<<" "<<eyy<<" "<<exy<<endl;
		return true;
	}
	if (exy*exy > thresSquared){
		//cout<<"Element: "<<Id<<" exx, eyy, exy: "<<exx<<" "<<eyy<<" "<<exy<<endl;
		return true;
	}
	return false;
}

void	ShapeBase::updateEquilibriumMyoWithFeedbackFromZero(double MyosinFeedbackCap){
	if (tissueType == 0 ){//|| tissueType == 1){
		//the feedback is applied only in peripodial membrane or the disc proper.
		//Linker zones are not affected.
		if (spansWholeTissue || tissuePlacement == 1){
			//The feedback is only on the apical surface of the tissue:
			//This feedback does not focus on the aspect ratio of the tissue.
			//Only stretch in one direction can result in myosin feedback.
			//If the tissue is under compression, with a high aspect ratio, there will be no feedback.
			double e1 = 0.0, e2 = 0.0, tet = 0.0;
			bool calculatePrincipalStrainDirection = checkIfXYPlaneStrainAboveThreshold(1E-5);
			if (calculatePrincipalStrainDirection){
				calculatePrincipalStrainAxesOnXYPlane(e1, e2, tet);
				double eEffectiveGreen = e1 - e2; //eEffective (eEffective = e1 - e2), is the difference between strains of major and minor axes.
				double eEffectiveEngineering = pow(2*eEffectiveGreen +1 , 0.5) - 1;

				double lowThres = 0.0;
				double upThres = 0.2;
				double diffThres = upThres-lowThres;
				//if(e1<lowThres){
				if(eEffectiveEngineering<lowThres){
					cMyoUnipolarEq[0] = 0.0;
				}else {
					//calculate concentration:
					//if (e1>upThres){
					if (eEffectiveEngineering>upThres){
						cMyoUnipolarEq[0] = MyosinFeedbackCap;
					}
					else{
						//cMyoUnipolarEq[0] = MyosinFeedbackCap*(e1 - lowThres)/diffThres;
						cMyoUnipolarEq[0] = MyosinFeedbackCap*(eEffectiveEngineering - lowThres)/diffThres;
					}
					//give the current direction:
					gsl_matrix_set(myoPolarityDir,0,0,cos(tet));
					gsl_matrix_set(myoPolarityDir,0,1,sin(tet));
					gsl_matrix_set(myoPolarityDir,0,2,0.0);
					//cout<<"Element: "<<Id<<" e1: "<<e1<<" e2 "<<e2 <<" tet: "<<tet<<" cmyo: "<<cMyoUnipolarEq[0]<<endl;
					//displayMatrix(Strain, "strain");
					//displayMatrix(myoPolarityDir, "myoPolarityDir");
				}
			}
		}
		/*double exx = gsl_matrix_get(Strain,0,0);
		double eyy = gsl_matrix_get(Strain,1,0);
		double feedbackStrain = 0.0;
		int directionIndex = 0; //by default the myosin response is in x
		if (eyy > exx){
			//if the stretch is higher in yy, then change the direction
			//and assign the strain to be used in feedback calculation.
			directionIndex = 1;
			feedbackStrain = eyy;
		}
		else{
			feedbackStrain = exx;
		}

		if(feedbackStrain<lowThres){
			cMyoUnipolarEq[0] = 0.0;
		}else if (feedbackStrain>upThres){
			cMyoUnipolarEq[0] = MyosinFeedbackCap;
			//reset direction:
			gsl_matrix_set(myoPolarityDir,0,0,0);
			gsl_matrix_set(myoPolarityDir,0,1,0.0);
			gsl_matrix_set(myoPolarityDir,0,2,0.0);
			//give the current direction:
			gsl_matrix_set(myoPolarityDir,0,directionIndex,1.0);
		}
		else{
			cMyoUnipolarEq[0] = feedbackStrain*1000.0;
			//reset direction:
			gsl_matrix_set(myoPolarityDir,0,0,0);
			gsl_matrix_set(myoPolarityDir,0,1,0.0);
			gsl_matrix_set(myoPolarityDir,0,2,0.0);
			//give the current direction:
			gsl_matrix_set(myoPolarityDir,0,directionIndex,1.0);
		}*/

	}
}

void	ShapeBase::updateEquilibriumMyoWithFeedbackFromFixedTotal(double totalMyosinLevel){
	if (tissueType == 0 ){//|| tissueType == 1){
		//the feedback is applied only in the disc proper.
		//Linker zones are not affected.
		if (spansWholeTissue || tissuePlacement == 1){
			//The feedback is only on the apical surface of the tissue:
			//This feedback does not focus on the aspect ratio of the tissue.
			//Only stretch in one direction can result in myosin feedback.
			//If the tissue is under compression, with a high aspect ratio, there will be no feedback.
			double e1 = 0.0, e2 = 0.0, tet = 0.0;
			bool calculatePrincipalStrainDirection = checkIfXYPlaneStrainAboveThreshold(1E-5);
			if (calculatePrincipalStrainDirection){
				calculatePrincipalStrainAxesOnXYPlane(e1, e2, tet);
				//The experimental curve between strain and the ratio of polarised to non-polar myosin:
				//Y = 0.6164* X + 1.038
				//X = strain ( 50% strain is 0.5, defined in engineering strain terms),
				//Y = ratio of polar/uniform myo
				//The strain in the model is Green strain, I need to do the conversion:
				//Exx = ((1+exx)^2 - 1);
				//Then the formulation is:
				// Y = 0.6164 * ( (2*Exx +1)^0.5 - 1 ) + 1.038
				//In the experiments, the strain is applied by the stretcher, in one direction only.
				//In the model, this should be the residual strain, as is the difference of the strains
				//of the major and minor axis, as long as, e1 is positive (there is stretch)
				if (e1>0){
					//cout<<"Element : "<<Id<<" has positive strain: [e1, e2, tet]: "<<e1<<" "<<e2<<" "<<tet<<endl;
					//there is stretch
					double eEffectiveGreen = e1 - e2; //eEffective (eEffective = e1 - e2), is the difference between strains of major and minor axes.
					double eEffectiveEngineering = pow(2*eEffectiveGreen +1 , 0.5) - 1;
					double ratioOfPolarToUniformMyo = 0.6164*eEffectiveEngineering + 1.038;
					cMyoUniformEq[0] = totalMyosinLevel / (1+ratioOfPolarToUniformMyo);
					cMyoUnipolarEq[0] = totalMyosinLevel - cMyoUniformEq[0];
					//give the current direction:
					gsl_matrix_set(myoPolarityDir,0,0,cos(tet));
					gsl_matrix_set(myoPolarityDir,0,1,sin(tet));
					gsl_matrix_set(myoPolarityDir,0,2,0.0);
					//cout<<"Element: "<<Id<<" e1: "<<e1<<" e2 "<<e2 <<" tet: "<<tet<<" cmyo: "<<cMyoUnipolarEq[0]<<endl;
					//displayMatrix(Strain, "strain");
					//displayMatrix(myoPolarityDir, "myoPolarityDir");
				}
			}
			else{
				cMyoUnipolarEq[0] = 0.0;
				cMyoUniformEq[0] = totalMyosinLevel;
			}
			//cout<<" element: "<<Id<< " e1, e2: "<<e1<<" "<<e2<<" uniform: "<<cMyoUniformEq[0]<<" polar: "<<cMyoUnipolarEq[0]<<endl;
		}

	}
}

void	ShapeBase::updateMyosinConcentration(double dt, double kMyo, bool thereIsMyosinFeedback, double MyosinFeedbackCap){
	double thresholdValue = 1E-8, thresholdFraction= 0.01;
	//the value of kMyo is taken form my thesis
	double currMyoDt[3] = {dt,dt*2.0,dt/2.0};
	double cFinal[3];
	//First value is the final value with the current time step,
	//second is with currTimeStep*2 and
	//third is with 0.5 currTimeStep;
	double cInitial, cEq;
	//0 for any polarity below
	if (thereIsMyosinFeedback){
		updateEquilibriumMyoWithFeedbackFromZero(MyosinFeedbackCap); // this function adds in polarised myosin to the tissue, independent of the uniform myosin levels
		//updateEquilibriumMyoWithFeedbackFromFixedTotal(MyosinFeedbackCap); //this function redistributes a total pool of myosin in between uniform and polarised myosin levels
	}
	bool useDiffusion = false;
	if (useDiffusion){
		for (int myoIter =0; myoIter<4; myoIter++){
			if (myoIter == 0){ //apical uniform
				cInitial = cMyoUniform[0];
				cEq = cMyoUniformEq[0];
			}
			else if (myoIter == 1){//basal uniform
				cInitial = cMyoUniform[1];
				cEq = cMyoUniformEq[1];
			}
			else if (myoIter == 2){//apical polar
				cInitial = cMyoUnipolar[0];
				cEq = cMyoUnipolarEq[0];
			}
			else if (myoIter == 3){//basal polar
				cInitial = cMyoUnipolar[1];
				cEq = cMyoUnipolarEq[1];
			}
			bool converged = false;
			while (!converged){
				int steps[3] = {dt/currMyoDt[0],dt/currMyoDt[1],dt/currMyoDt[2]};
				for (int j=0; j<3; ++j){
					cFinal[j] = cInitial;
					for (int i =0 ;i<steps[j]; ++i){
						cFinal[j] += (cEq - cFinal[j])*kMyo*currMyoDt[j];
					}
				}
				//check if the value of the current dt and half current dt are below the threshold:
				double diff = fabs((cFinal[0] - cFinal[2]));
				if ( diff < thresholdValue ){
					converged = true;
				}
				else if( fabs(diff / cFinal[2]) < thresholdFraction){
					converged = true;
				}
				else{
					currMyoDt[1] = currMyoDt[0];
					currMyoDt[0] = currMyoDt[2];
					currMyoDt[2] *= 0.5;
				}
				//need to implement increasing this
			}
			if (myoIter == 0){
				cMyoUniform[0] = cFinal[2];
			}
			else if (myoIter == 1){
				cMyoUniform[1] = cFinal[2];
			}
			else if (myoIter == 2){
				cMyoUnipolar[0] = cFinal[2];
			}
			else if (myoIter == 3){
				cMyoUnipolar[1] = cFinal[2];
			}
		}
	}
	else{
		cMyoUniform[0] = cMyoUniformEq[0];
		cMyoUniform[1] = cMyoUniformEq[1];
		cMyoUnipolar[0] = cMyoUnipolarEq[0];
		cMyoUnipolar[1] = cMyoUnipolarEq[1];
	}
	//if (cMyoUnipolar[0] > 0 || cMyoUnipolar[1]>0){
		//cout<<" Element: "<<Id<<" unipolar myo:  "<<cMyoUnipolar[0]<<" "<<cMyoUnipolar[1]<<" eq: "<<cMyoUnipolarEq[0]<<" "<<cMyoUnipolarEq[1]<<" polarity dir: "<<endl;
		//displayMatrix(myoPolarityDir,"myoPolarityDir");
		//cout<<" Element: "<<Id<<" uniform myo:  "<<cMyoUniform[0]<<" "<<cMyoUniform[1]<<" eq: "<<cMyoUniformEq[0]<<" "<<cMyoUniformEq[1]<<endl;
	//}
}

void ShapeBase::adjustCMyosinFromSave(){
	//if (cMyoUnipolarEq[0] < 1E-5){cMyoUnipolar[0] = 0.0;}
	//if (cMyoUnipolarEq[1] < 1E-5){cMyoUnipolar[1] = 0.0;}
	cMyoUnipolar[0] = 0.0;
	cMyoUnipolar[1] = 0.0;
}


void 	ShapeBase::setShapeChangeInrementToIdentity(){
	gsl_matrix_set_identity(shapeChangeIncrement);
}

void 	ShapeBase::setShapeChangeRate(double x, double y, double z, double xy, double yz, double xz){
	ShapeChangeRate[0] = x;
	ShapeChangeRate[1] = y;
	ShapeChangeRate[2] = z;
	ShapeChangeRate[3] = xy;
	ShapeChangeRate[4] = yz;
	ShapeChangeRate[5] = xz;
}

//R Fginc R^T version:
void ShapeBase::calculateCurrentGrowthIncrement(gsl_matrix* resultingGrowthIncrement, double dt, double growthx, double growthy, double growthz, gsl_matrix* ShearAngleRotationMatrix){
	gsl_matrix_set_identity(resultingGrowthIncrement);
	gsl_matrix_set(resultingGrowthIncrement,0,0,exp(growthx*dt));
	gsl_matrix_set(resultingGrowthIncrement,1,1,exp(growthy*dt));
	gsl_matrix_set(resultingGrowthIncrement,2,2,exp(growthz*dt));
	gsl_matrix* temp = gsl_matrix_calloc(nDim,nDim);;
	//R * resultingGrowthIncrement
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, ShearAngleRotationMatrix, resultingGrowthIncrement, 0.0, temp);
	//R * resultingGrowthIncrement * R^T
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,1.0, temp, ShearAngleRotationMatrix, 0.0, resultingGrowthIncrement);
	gsl_matrix_free(temp);

}

void 	ShapeBase::updateShapeChangeRate(double x, double y, double z, double xy, double yz, double xz){
	ShapeChangeRate[0] += x;
	ShapeChangeRate[1] += y;
	ShapeChangeRate[2] += z;
	ShapeChangeRate[3] += xy;
	ShapeChangeRate[4] += yz;
	ShapeChangeRate[5] += xz;
}

bool 	ShapeBase::InvertMatrix(gsl_matrix* input, gsl_matrix* inverse){
    // Define the dimension n of the matrix
    // and the signum s (for LU decomposition)
    int s;

    // Define all the used matrices
    gsl_permutation * perm = gsl_permutation_alloc (input->size1);

    // Make LU decomposition of matrix m
    gsl_linalg_LU_decomp (input, perm, &s);

    // Invert the matrix m
    gsl_linalg_LU_invert (input, perm, inverse);

    return true;
}

bool 	ShapeBase::InvertMatrix(boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse){
	//Matrix inversion routine.
	//Uses lu_factorize and lu_substitute in uBLAS to invert a matrix
	using namespace boost::numeric::ublas;
	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix<double> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<double> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

void	ShapeBase::crossProduct3D(double* u, double* v, double* cross){
	cross[0] = u[1]*v[2] - u[2]*v[1];
	cross[1] = u[2]*v[0] - u[0]*v[2];
	cross[2] = u[0]*v[1] - u[1]*v[0];
}

void	ShapeBase::crossProduct3D(gsl_vector* u, gsl_vector* v, gsl_vector* cross){
	gsl_vector_set(cross,0, ( gsl_vector_get(u,1)*gsl_vector_get(v,2) - gsl_vector_get(u,2)*gsl_vector_get(v,1) ) );
	gsl_vector_set(cross,1, ( gsl_vector_get(u,2)*gsl_vector_get(v,0) - gsl_vector_get(u,0)*gsl_vector_get(v,2) ) );
	gsl_vector_set(cross,2, ( gsl_vector_get(u,0)*gsl_vector_get(v,1) - gsl_vector_get(u,1)*gsl_vector_get(v,0) ) );
}

double	ShapeBase::calculateMagnitudeVector3D(double* v){
	double mag = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
	mag = pow(mag,0.5);
	return mag;
}

double	ShapeBase::normaliseVector3D(double* v){
	double mag2 = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
	if (fabs(mag2) > 1E-14 && fabs(mag2 - 1.0f) > 1E-14) {
		double mag = pow(mag2,0.5);
		v[0] /= mag;
		v[1] /= mag;
		v[2] /= mag;
		return mag;
	}
	else{
		return 0;
	}
}
void	ShapeBase::normaliseVector3D(gsl_vector* v){
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	double z = gsl_vector_get(v,2);
	double mag2 = x*x + y*y + z*z;
	if (fabs(mag2) > 1E-14 && fabs(mag2 - 1.0f) > 1E-14) {
		double mag = pow(mag2,0.5);
		gsl_vector_scale(v,1.0/mag);
	}
}

double	ShapeBase::getNormVector3D(gsl_vector* v){
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	double z = gsl_vector_get(v,2);
	double mag2 = x*x + y*y + z*z;
	return pow(mag2,0.5);
}

double 	ShapeBase::dotProduct3D(double* u, double* v){
	double dot = 0;
	dot = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
	return dot;
}


void 	ShapeBase::displayMatrix(boost::numeric::ublas::matrix<double>& mat, string matname){
	int m = mat.size1();
	int n = mat.size2();
	cout<<matname<<": "<<endl;

	for (int i =0; i<m; i++){
		for (int j =0; j<n; j++){
			cout.precision(4);
			cout.width(6);
			cout<<mat(i,j)<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
}

void 	ShapeBase::displayMatrix(gsl_matrix* mat, string matname){
    int m = mat->size1;
    int n = mat->size2;
    cout<<matname<<": "<<endl;

    for (int i =0; i<m; i++){
        for (int j =0; j<n; j++){
            cout.precision(4);
            cout.width(6);
            cout<<gsl_matrix_get(mat,i,j)<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}

void 	ShapeBase::displayMatrix(gsl_vector* mat, string matname){
    int m = mat->size;
    cout<<matname<<": "<<endl;

    for (int i =0; i<m; i++){
        cout.precision(4);
        cout.width(6);
        cout<<gsl_vector_get(mat,i)<<endl;
    }
    cout<<endl;
}

void 	ShapeBase::displayMatrix(boost::numeric::ublas::matrix<int>& mat, string matname){
	int m = mat.size1();
	int n = mat.size2();
	cout<<matname<<": "<<endl;

	for (int i =0; i<m; i++){
		for (int j =0; j<n; j++){
			cout.precision(4);
			cout.width(6);
			cout<<mat(i,j)<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
}

void	ShapeBase::displayMatrix(boost::numeric::ublas::vector<double>& vec, string matname){
	int m = vec.size();
	cout<<matname<<": "<<endl;
	for (int i =0; i<m; i++){
		cout.precision(4);
		cout.width(6);
		cout<<vec(i)<<" ";
	}
	cout<<endl;
}

void 	ShapeBase::assignVolumesToNodes(vector <Node*>& Nodes){
	for (int i=0; i<nNodes; i++){
        Nodes[NodeIds[i]]->mass += VolumePerNode;
	}
}

void 	ShapeBase::calculateExposedLateralAreaBasalSide(){
	double Threshold = 1E-5;
	exposedLateralAreaBasalSide = 0;
	int id0 = exposedLateralAreaBasalSideNodeIds[0];
	int id1 = exposedLateralAreaBasalSideNodeIds[1];
	int id2 = exposedLateralAreaBasalSideNodeIds[2];
	int id3 = exposedLateralAreaBasalSideNodeIds[3];

	double sideVec1[3];
	double sideVec2[3];
	double Side1 = 0.0;
	double Side2 = 0.0;
	double costet = 0.0;
	double Area = 0.0;
	for (int i = 0; i<3; ++i){
		sideVec1[i]= Positions[id1][i] - Positions[id0][i];
		sideVec2[i]= Positions[id2][i] - Positions[id0][i];
		costet += sideVec1[i] * sideVec2[i];
		Side1  += sideVec1[i] * sideVec1[i];
		Side2  += sideVec2[i] * sideVec2[i];
	}
	if (Side1 > Threshold && Side2 > Threshold){
		Side1 = pow(Side1,0.5);
		Side2 = pow(Side2,0.5);
		costet /= (Side1*Side2);
		double sintet = pow((1-costet*costet),0.5);
		Area = Side1* Side2 * sintet / 2.0;
	}
	exposedLateralAreaBasalSide += Area;
	for (int i = 0; i<3; ++i){
		sideVec1[i]= Positions[id3][i] - Positions[id2][i];
		sideVec2[i]= Positions[id0][i] - Positions[id2][i];
		costet += sideVec1[i] * sideVec2[i];
		Side1  += sideVec1[i] * sideVec1[i];
		Side2  += sideVec2[i] * sideVec2[i];
	}
	if (Side1 > Threshold && Side2 > Threshold){
		Side1 = pow(Side1,0.5);
		Side2 = pow(Side2,0.5);
		costet /= (Side1*Side2);
		double sintet = pow((1-costet*costet),0.5);
		Area = Side1* Side2 * sintet / 2.0;
	}
	exposedLateralAreaBasalSide  += Area;
	//cout<<" Element "<<Id<<" exposedLateralAreaBasalSide: "<<exposedLateralAreaBasalSide<<endl;

}

void 	ShapeBase::calculateExposedLateralAreaApicalSide(){
	double Threshold = 1E-5;
	exposedLateralAreaApicalSide = 0;
	int id0 = exposedLateralAreaApicalSideNodeIds[0];
	int id1 = exposedLateralAreaApicalSideNodeIds[1];
	int id2 = exposedLateralAreaApicalSideNodeIds[2];
	int id3 = exposedLateralAreaApicalSideNodeIds[3];

	double sideVec1[3];
	double sideVec2[3];
	double Side1 = 0.0;
	double Side2 = 0.0;
	double costet = 0.0;
	double Area = 0.0;
	for (int i = 0; i<3; ++i){
		sideVec1[i]= Positions[id1][i] - Positions[id0][i];
		sideVec2[i]= Positions[id2][i] - Positions[id0][i];
		costet += sideVec1[i] * sideVec2[i];
		Side1  += sideVec1[i] * sideVec1[i];
		Side2  += sideVec2[i] * sideVec2[i];
	}
	if (Side1 > Threshold && Side2 > Threshold){
		Side1 = pow(Side1,0.5);
		Side2 = pow(Side2,0.5);
		costet /= (Side1*Side2);
		double sintet = pow((1-costet*costet),0.5);
		Area = Side1* Side2 * sintet / 2.0;
	}
	exposedLateralAreaApicalSide += Area;
	for (int i = 0; i<3; ++i){
		sideVec1[i]= Positions[id3][i] - Positions[id2][i];
		sideVec2[i]= Positions[id0][i] - Positions[id2][i];
		costet += sideVec1[i] * sideVec2[i];
		Side1  += sideVec1[i] * sideVec1[i];
		Side2  += sideVec2[i] * sideVec2[i];
	}
	if (Side1 > Threshold && Side2 > Threshold){
		Side1 = pow(Side1,0.5);
		Side2 = pow(Side2,0.5);
		costet /= (Side1*Side2);
		double sintet = pow((1-costet*costet),0.5);
		Area = Side1* Side2 * sintet / 2.0;
	}
	exposedLateralAreaApicalSide  += Area;
	//cout<<" Element "<<Id<<" exposedLateralAreaApicalSide: "<<exposedLateralAreaApicalSide<<endl;

}


void 	ShapeBase::calculateViscositySurfaces(){
	if (elementHasExposedLateralApicalSurface){
		calculateExposedLateralAreaApicalSide();
	}
	if (elementHasExposedLateralBasalSurface){
		calculateExposedLateralAreaBasalSide();
	}
	if (elementHasExposedApicalSurface){
		calculateApicalArea();
	}
	if (elementHasExposedBasalSurface){
		calculateBasalArea();
	}
}

void ShapeBase::assignViscositySurfaceAreaToNodes(vector <Node*>& Nodes){
	if (elementHasExposedLateralApicalSurface){
		for(int i=0;i<nLateralSurfaceAreaNodeNumber; ++i){
			Nodes[NodeIds[exposedLateralAreaApicalSideNodeIds[i]]]->viscositySurface+=exposedLateralAreaApicalSide/nLateralSurfaceAreaNodeNumber;
		}
	}
	if (elementHasExposedLateralBasalSurface){
		for(int i=0;i<nLateralSurfaceAreaNodeNumber; ++i){
			Nodes[NodeIds[exposedLateralAreaBasalSideNodeIds[i]]]->viscositySurface+=exposedLateralAreaBasalSide/nLateralSurfaceAreaNodeNumber;
		}
	}
	if (elementHasExposedApicalSurface){
		for(int i=0;i<nSurfaceAreaNodeNumber; ++i){
			Nodes[NodeIds[exposedApicalSurfaceNodeIds[i]]]->viscositySurface+=ApicalArea/nSurfaceAreaNodeNumber;
		}
	}
	if (elementHasExposedBasalSurface){
		for(int i=0;i<nSurfaceAreaNodeNumber; ++i){
			Nodes[NodeIds[exposedBasalSurfaceNodeIds[i]]]->viscositySurface+=BasalArea/nSurfaceAreaNodeNumber;
		}
	}
}

/*
void 	ShapeBase:: assignSurfaceAreaToNodes(vector <Node*>& Nodes){
    double multiplier = 1.0;
    if (ShapeType ==1 ){ multiplier = 0.5;}
    for (int i=0; i<nNodes; i++){
        Nodes[NodeIds[i]]->surface +=ReferenceShape->BasalArea/(multiplier*nNodes);
	}
}*/

void 	ShapeBase::calculateZProjectedAreas(){
    double Threshold = 1E-5;
    int id0 = 0, id1 = 1, id2 = 2; // this is correct for basal side, I will change it for apical calculation
    for (int tissueSide = 0; tissueSide<2; tissueSide++){
        if ( tissueSide == 1){
            //I am calculating basal area,
            id0 = 3;
            id1 = 4;
            id2 = 5;
        }
        double sideVec1[2];
        double sideVec2[2];
        double Side1 = 0.0;
        double Side2 = 0.0;
        double costet = 0.0;
        double Area = 0.0;
        for (int i = 0; i<2; ++i){
            sideVec1[i]= Positions[id1][i] - Positions[id0][i];
            sideVec2[i]= Positions[id2][i] - Positions[id0][i];
            costet += sideVec1[i] * sideVec2[i];
            Side1  += sideVec1[i] * sideVec1[i];
            Side2  += sideVec2[i] * sideVec2[i];
        }
        if (Side1 > Threshold && Side2 > Threshold){
            Side1 = pow(Side1,0.5);
            Side2 = pow(Side2,0.5);
            costet /= (Side1*Side2);
            double sintet = pow((1-costet*costet),0.5);
            Area = Side1* Side2 * sintet / 2.0;
        }
        if(tissueSide == 0){
            ZProjectedBasalArea = Area;
        }
        else{
            ZProjectedApicalArea = Area;
        }
    }
}




void 	ShapeBase::assignZProjectedAreas(vector <Node*> Nodes){
    if (ShapeType == 1 ){ //only written for prisms
        for (int i=0; i<3; i++){
            Nodes[NodeIds[i]]->zProjectedArea +=ZProjectedBasalArea/3.0;
        }
        for (int i=3; i<6; i++){
            Nodes[NodeIds[i]]->zProjectedArea +=ZProjectedApicalArea/3.0;
        }
    }
}

void 	ShapeBase:: assignElementToConnectedNodes(vector <Node*>& Nodes){
	for (int i=0; i<nNodes; i++){
		Nodes[NodeIds[i]]->connectedElementIds.push_back(Id);
		double weightfFraction = (ReferenceShape->Volume/nNodes)/Nodes[NodeIds[i]]->mass;
		Nodes[NodeIds[i]]->connectedElementWeights.push_back(weightfFraction);
		if(tissueType == 2){
			Nodes[NodeIds[i]]->hasLateralElementOwner = true;
		}
	}
}

void 	ShapeBase::removeMassFromNodes(vector <Node*>& Nodes){
	for (int i=0; i<nNodes; i++){
			Nodes[NodeIds[i]]->mass -= VolumePerNode;
			//updating the weight fractions of the elements on the node due to elimination of the ablated element:
			int n = Nodes[NodeIds[i]]->connectedElementIds.size();
			double scaler = 1.0;
			for (int j=0;j<n;++j){
				if (Nodes[NodeIds[i]]->connectedElementIds[j]==Id){
					scaler = 1.0 - Nodes[NodeIds[i]]->connectedElementWeights[j];
					Nodes[NodeIds[i]]->connectedElementWeights[j]  = 0.0;
					break;
				}
			}
			for (int j=0;j<n;++j){
				Nodes[NodeIds[i]]->connectedElementWeights[j] /= scaler;
			}
			//All weights are normlised as the sum will make 1.0. Now I do not want this element in the weighing,
			//it does not have a mass anymore, therefore I will multiply all the remaining weights with (1-w_ablated);
		}
}

void 	ShapeBase::checkDisplayClipping(double xClip, double yClip, double zClip){
	IsClippedInDisplay=false;
	IsXSymmetricClippedInDisplay=false;
	IsYSymmetricClippedInDisplay=false;
	for (int j=0; j<nNodes; ++j){
		 if((-1.0)*Positions[j][0]>xClip){
			 IsXSymmetricClippedInDisplay = true;
		 }
		 if((-1.0)*Positions[j][1]<yClip){
			 IsYSymmetricClippedInDisplay = true;
		 }
		 if(Positions[j][0]>xClip){
			 IsClippedInDisplay = true;
			 return;
		 }
		 if(Positions[j][1]<yClip){
			 IsClippedInDisplay = true;
			 return;
		 }
		 if(Positions[j][2]>zClip){
			 IsClippedInDisplay = true;
			 return;
		 }
	 }
}

void 	ShapeBase::doesElementNeedRefinement(double areaThreshold, int surfacedentifier){
	if (surfacedentifier == 0){
		//checking for basal surface
		if ( BasalArea > areaThreshold){
			willBeRefined = true;
		}
	}
	else if (surfacedentifier == 1){
		//checking for apical surface
		if ( ApicalArea > areaThreshold){
			willBeRefined = true;
		}
	}
	if (willBeRefined){
		cout<<" Element "<<Id<<" will be  refined: "<<willBeRefined<<" apical area: "<<ApicalArea<<" basal area: "<<BasalArea<<endl;
	}
}

/*
double ShapeBase::calculateECMThickness(vector <Node*>& Nodes){
	double d;
	if (isECMMimicing){
		if(tissueType == 2){
			int basalIndex = -1;
			int apicalIndex = -2;
			int thirdIndex = -3;
			for (int i =0; i<3; ++i){
				if (Nodes[NodeIds[i]]->tissuePlacement == 0 ){
					basalIndex = i;
					break;
				}
			}
			for (int i =0; i<3; ++i){
				if (Nodes[NodeIds[i]]->tissuePlacement != 0 ){
					apicalIndex = i;
					break;
				}
			}
			for (int i =0; i<3; ++i){
				if (i != apicalIndex && i!=basalIndex  ){
					thirdIndex = i;
					break;
				}
			}
			bibap
			double * u = new double [3];
			double * v = new double [3];
			u[0] = Positions[apicalIndex][0] - Positions[basalIndex][0];
			u[1] = Positions[apicalIndex][1] - Positions[basalIndex][1];
			u[2] = Positions[apicalIndex][2] - Positions[basalIndex][2];
			v[0] = Positions[thirdIndex][0] - Positions[basalIndex][0];
			v[1] = Positions[thirdIndex][1] - Positions[basalIndex][1];
			v[2] = Positions[thirdIndex][2] - Positions[basalIndex][2];
		}
		else{
			gsl_matrix_set_identity(ECMThicknessPlaneRotationalMatrix);
		}
	}
	return d;
}

void ShapeBase::calculateInitialECMThickness(){
	if (isECMMimicing){
		initialECMThickness = calculateECMThickness();

	}
}*/

void ShapeBase::scaleGrowthIncrement(double multiplier){
	gsl_matrix_scale(growthIncrement,multiplier);
}

double ShapeBase::getGrowthMutationMultiplier(){
	return 1.0;
	/*float pouchScale= 1.0;
	float hingeScale = 0.5;
	float notumScale = 1.0;
	double scale = 1.0;
	if (compartmentType == 0){
		scale = pouchScale;
	}
	else if (compartmentType == 1){
		scale = hingeScale;
	}
	else if (compartmentType == 2){
		scale = notumScale;
	}
	else{
		return 1.0;
	}
	scale = (scale - 1)*compartmentIdentityFraction +1;
	return scale;*/
}

void ShapeBase::setECMMimicingElementThicknessGrowthAxis(){
	if (isECMMimicing){
		if(tissueType == 2){
			double* rotAx = new double[3];
			rotAx[0] = Positions[3][0] - Positions[0][0];
			rotAx[1] = Positions[3][1] - Positions[0][1];
			rotAx[2] = Positions[3][2] - Positions[0][2];
			double dummy = normaliseVector3D(rotAx);
			//I rotation angle  will come from peripodialness. If peripodialness is 1, the angle is pi, if it is 0, the angle is 0.
			//I will rotate the coordinates via y axis in psi degrees:
			double psi = (-1.0) * peripodialGrowthWeight*M_PI;
			// I want to give growth to ECM elements such that each growth will happen in this axis. I will give the growth in z axis.
			// Therefore the rotation should be at -phi, so that the vector will look at (-z)
			psi *= (-1.0);
			double cosPsi = cos(psi);
			double sinPsi = sin(psi);
			double * rotMat = new double [9];
			constructRotationMatrix(cosPsi, sinPsi, rotAx, rotMat);
			for (int i=0; i<3; ++i){
				for (int j=0; j<3;++j){
					gsl_matrix_set(ECMThicknessPlaneRotationalMatrix,i,j,rotMat[i*3+j]);
				}
			}
			delete[] rotAx;
			delete[] rotMat;
		}
		else{
			gsl_matrix_set_identity(ECMThicknessPlaneRotationalMatrix);
		}
	}
	else{
		gsl_matrix_set_identity(ECMThicknessPlaneRotationalMatrix);
	}

}

bool ShapeBase::areanyOfMyNodesAtCircumference(vector<Node*>& Nodes){
	bool thereIsNodeAtCircumference = false;
	for (int i=0; i< nNodes; ++i){
		if (Nodes[NodeIds[i]]->atCircumference){
			thereIsNodeAtCircumference = true;
			break;
		}
	}
	return thereIsNodeAtCircumference;
}

void ShapeBase::setLateralElementsRemodellingPlaneRotationMatrix(double systemCentreX, double systemCentreY){
	if(tissueType == 2){
		if (ShapeType == 1){
			//calculating the z vector I am interested in:
			//rotation axis will be the axis in the xy plane, for a prism, this is the
			//vector from node 0 to 3
			double* rotAx = new double[3];
			rotAx[0] = Positions[3][0] - Positions[0][0];
			rotAx[1] = Positions[3][1] - Positions[0][1];
			rotAx[2] = Positions[3][2] - Positions[0][2];
			double dummy = normaliseVector3D(rotAx);
			//I rotation angle  will come from peripodialness. If peripodialness is 1, the angle is pi, if it is 0, the angle is 0.
			//I will rotate the coordinates via y axis in psi degrees:
			double psi = (-1.0) * peripodialGrowthWeight*M_PI;
			double cosPsi = cos(psi);
			double sinPsi = sin(psi);
			double * rotMat = new double [9];
			constructRotationMatrix(cosPsi, sinPsi, rotAx, rotMat);
			for (int i=0; i<3; ++i){
				for (int j=0; j<3;++j){
					gsl_matrix_set(remodellingPlaneRotationMatrix,i,j,rotMat[i*3+j]);
				}
			}
			//Now I will align the y vector with the rotation axis, so that all the elemetns of the system will be consistent:
			//rotate the y unit vector with the matrix to see where it will be:
			gsl_matrix* y = gsl_matrix_calloc(3,1);
			gsl_matrix* tmp = gsl_matrix_calloc(3,1);
			gsl_matrix_set(tmp,1,0,1);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, remodellingPlaneRotationMatrix, tmp, 0.0, y);
			//now find the rotation you need to bring this y vector on top of rotation axis:
			double c,s;
			double* u = new double [3];
			double* rotAx2 = new double [3];
			u[0] = gsl_matrix_get(y,0,0);
			u[1] = gsl_matrix_get(y,1,0);
			u[2] = gsl_matrix_get(y,2,0);
			gsl_matrix_free(y);
			//calculating the rotation angle between my current y vector and the rotation axis:
			calculateRotationAngleSinCos(u, rotAx, c, s);
			//calculating the corresponding rotation axis:
			calculateRotationAxis(u, rotAx, rotAx2, c);
			//The rotation axis will be equal to my new z axis, or 180 degrees rotated version
			gsl_matrix* z = gsl_matrix_calloc(3,1);
			gsl_matrix_set(tmp,0,0,0);
			gsl_matrix_set(tmp,1,0,0);
			gsl_matrix_set(tmp,2,0,1);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, remodellingPlaneRotationMatrix, tmp, 0.0, z);
			double* v = new double[3];
			v[0] = gsl_matrix_get(z,0,0);
			v[1] = gsl_matrix_get(z,1,0);
			v[2] = gsl_matrix_get(z,2,0);
			gsl_matrix_free(z);
			double cRotAxes = dotProduct3D(rotAx2,v);
			//If cosince is negative, then the rotation angle I am using should be altered as well tet -> (-tet)
			//then cos tet = costet, sintet = -sintet
			if (cRotAxes <0){
				s *= (-1.0);
			}
			//now I can construct the rotation matrix:
			constructRotationMatrix(c, s, v, rotMat);
			//then write it on a gsl matrix and ad to current rotation matrix:
			gsl_matrix* temRotMat = gsl_matrix_calloc(3,3);
			for (int i=0; i<3; ++i){
				for (int j=0; j<3;++j){
					gsl_matrix_set(temRotMat,i,j,rotMat[i*3+j]);
				}
			}
			gsl_matrix* tmp2 = gsl_matrix_calloc(3,3);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, temRotMat, remodellingPlaneRotationMatrix, 0.0, tmp2);
		    createMatrixCopy(remodellingPlaneRotationMatrix,tmp2);
			delete[] u;
			delete[] v;
			delete[] rotAx;
			delete[] rotMat;
			delete[] rotAx2;
			gsl_matrix_free(tmp);
			gsl_matrix_free(temRotMat);
			gsl_matrix_free(tmp2);
		}
		else{
			//this axis calculation will only work for prisms!!!
			cout<<"Error! Trying to calculate remodellingPlaneRotationMatrix for non-prism element, update code, I will delte the matrix to make code break here!"<<endl;
			cerr<<"Error! Trying to calculate remodellingPlaneRotationMatrix for non-prism element, update code, I will delte the matrix to make code break here!"<<endl;
			gsl_matrix_free(remodellingPlaneRotationMatrix);
		}
	}
	//This function should not have been called if the element is not a linker!
}
