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
    	double PeriR1[4];
    	double PeriR2[2];
	double condensation;
	double sideLen;
	double sideLenInit; 
	vector <double> posx,posy;//,posz;
    	vector <int> tissueType; //0: columnar layer, 1: peripodium, 2 for linker
    	vector <double> Peripodialposx, Peripodialposy;
    	vector <int> PeripodialtissueType; //0: columnar layer, 1: peripodium, 2 for linker
	vector <int> atBorder;	//1 id the node is at the circumference of the setup, 0 otherwise
	double tetha;
	double dtet[4];
    	double Peripodialdtet[4];
	double Circumference[4];
    	double PeripodialCircumference[4];
	int Vborder[4];
    	int PeripodialVborder[4];
	bool tethaAtZero;
	bool r1failed;	
	vector <int> Links0,Links1;	
	vector <int*>triangles;
	double dForNodeGeneration;
	bool symmetricY, symmetricX;
    	int peripodialLayers;

	void calculateCircumference();
    	void calculatePeripodialCircumference();
	void calculateCurrentBorderNumber();
    	void calculatePeripodialBorderNumber();
	void calculatedtet();
    	void calculatePeripodialdtet();
	bool updateRadia();
	void addMidLine();
    	void addPeripodialMembrane(double ABHeight, double ABLayers, double lumenHeightFrac, double peripodialHeightFrac, double peripodialSideCurveFrac);
	void Tesselate2D();
	void readInTesselation2D();
	void writeVectors2D(ofstream &vectorsForGnuplot);
	void writeNodes2D(ofstream &nodesForGnuplot);
    void writePeripodialNodes2D(ofstream &nodesForGnuplot);
	void writeMeshFileForSimulation(double zHeight, int zLayers);
	void writeTriangleMeshFileForSimulation(double zHeight);
	void addEquidistantRingMin();
    void addPeripodialRings();
    void addPeripodialMidlines(double ABHeight, double ABLayers, double lumenHeightFrac, double peripodialHeightFrac);
    void minimiseCurrentPositions(bool addingPeripodial, int n , double d, int iterLen, int iterWid, vector<double> &CurrPosX, vector <double> &CurrPosY);
	void cleanUpCurrentPositions(double d, vector<double> &CurrPosX, vector <double> &CurrPosY);
	void cleanUpAgainsAlreadyExisting(double d, vector<double> &CurrPosX, vector <double> &CurrPosY);
	void addRectangle();
	void calculateAverageSideLength();
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
	this->symmetricX = false;

    PeriR1[0] = 0.0;
    PeriR1[1] = 0.0;
    PeriR2[0] = 0.0;
    PeriR2[1] = 0.0;
	tethaAtZero=true;
	r1failed = false;
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

void EllipseLayoutGenerator::calculatePeripodialCircumference(){
    for (int iterSide = 0 ;iterSide <4; iterSide++){
        int iterLen,iterWid;
        if (iterSide == 0){iterLen = 0; iterWid = 0;}		//iterSide == 0 is putting in (+)ve x - (+)ve y outer layer (peripodial side)
        else if (iterSide == 1){iterLen = 1; iterWid = 0;}	//iterSide == 1 is putting in (-)ve x - (+)ve y outer layer (columnar side)
        else if (iterSide == 2){iterLen = 2; iterWid = 1;}	//iterSide == 2 is putting in (+)ve x - (+)ve y inner layer (peripodial side)
        else if (iterSide == 3){iterLen = 3; iterWid = 1;}	//iterSide == 3 is putting in (-)ve x - (+)ve y inner layer (columnar side)
        double h = ((PeriR1[iterLen]-PeriR2[iterWid])*(PeriR1[iterLen]-PeriR2[iterWid]))/((PeriR1[iterLen]+PeriR2[iterWid])*(PeriR1[iterLen]+PeriR2[iterWid]));
        PeripodialCircumference[iterSide] = pi * (PeriR1[iterLen] + PeriR2[iterWid]) * (1.0 + (3*h / (10.0 + pow(4.0 - 3*h,0.5) )));
    }

}
void EllipseLayoutGenerator::calculateCurrentBorderNumber(){
	for (int i =0; i<4; ++i){
		Vborder[i] = Circumference[i] / (sideLen*1.8);
	}
}

void EllipseLayoutGenerator::calculatePeripodialBorderNumber(){
    for (int i =0; i<4; ++i){
        PeripodialVborder[i] = PeripodialCircumference[i] / (sideLen*1.8);
    }
}

void EllipseLayoutGenerator::calculatedtet(){
	for (int i =0; i<4; ++i){
		int n  = Vborder[i]/2;	
		dtet[i] = pi/n;
	}
}

void EllipseLayoutGenerator::calculatePeripodialdtet(){
	for (int i =0; i<4; ++i){
		int n  = PeripodialVborder[i]/2;
		Peripodialdtet[i] = pi/n;
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
                ellipseSum = (CurrPosX[i]/PeriR1[iterLen]) * (CurrPosX[i]/PeriR1[iterLen]) +  (CurrPosY[i]/PeriR2[iterWid]) * (CurrPosY[i]/PeriR2[iterWid]);
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

void EllipseLayoutGenerator::cleanUpAgainsAlreadyExisting(double d, vector<double> &CurrPosX, vector <double> &CurrPosY){
	//now doing the same check with already added points:
	double threshold = d*0.7;
	threshold *= threshold;		
	int n = CurrPosX.size();	
	int i=1;	
	while (i < n-1){
		for (int j=0; j< posx.size(); ++j){
			double dx = CurrPosX[i] - posx[j];
			double dy = CurrPosY[i] - posy[j];
			double d2 = dx*dx + dy*dy;
			//cerr<<"i: "<<i<<" j "<<j<<" d: "<<d2 <<endl;
			if (d2 < threshold){
				vector<double>::iterator it;
				it = CurrPosX.begin();
				it +=i;
				CurrPosX.erase(it);
				it = CurrPosY.begin();	
				it +=i;
				CurrPosY.erase(it);
				n = CurrPosX.size();
				i =0;			
				break;				
			}		
		}
		i++;	
	}
}


void EllipseLayoutGenerator::addPeripodialRings(){
    int iterLen, iterWid;
    for (int i = 0; i<4; i++){ //this is for left and right sides of outer(i=0,1) or inner (i = 2,3) loop
        if (i == 0){iterLen = 0; iterWid = 0;}		//iterSide == 0 is putting in (+)ve x - (+)ve y (outer loop)
        else if (i == 1){iterLen = 1; iterWid = 0;}	//iterSide == 1 is putting in (-)ve x - (+)ve y (outer loop)
        else if (i == 2){iterLen = 2; iterWid = 1;} //iterSide == 0 is putting in (+)ve x - (+)ve y (inner loop)
        else if (i == 3){iterLen = 3; iterWid = 1;}	//iterSide == 1 is putting in (-)ve x - (+)ve y (inner loop)
        double CQ = PeripodialCircumference[i]/4.0;
        int n= CQ / sideLen;
        Peripodialdtet[i]= pi/2.0/n;
        double d = CQ / n;
        cout<<" circumference of ring: "<<CQ<<" number of nodes: "<<n<<" dtet: "<<Peripodialdtet[i]<<" d: "<<d<<endl;

        dForNodeGeneration  = d;
        double x0 = PeriR1[i];
        double y0 = 0;
        vector <double> CurrPosX, CurrPosY;
        CurrPosX.push_back(x0);
        CurrPosY.push_back(y0);
        tetha = Peripodialdtet[i];
        for (int j= 0; j<n-1; ++j){
            double x = PeriR1[iterLen]*cos(tetha);
            double y = PeriR1[iterWid]*sin(tetha);
            CurrPosX.push_back(x);
            CurrPosY.push_back(y);
            tetha += Peripodialdtet[i];
        }
        CurrPosX.push_back(0);
        CurrPosY.push_back(PeriR2[iterWid]);
        cout<<" number of nodes after initial addition: "<<CurrPosX.size()<<endl;
        minimiseCurrentPositions(true , n , d, iterLen, iterWid, CurrPosX,CurrPosY);
        cout<<" number of nodes after minimisation : "<<CurrPosX.size()<<endl;
        cleanUpCurrentPositions(d, CurrPosX,CurrPosY);
        cout<<" number of nodes after clean up : "<<CurrPosX.size()<<endl;
        for (int j = 0; j<CurrPosX.size(); ++j){
            cout<<"node "<<j<<" pos: "<<CurrPosX[j]<<" "<<CurrPosY[j]<<endl;
        }
        n = CurrPosX.size();
        int startindex = 0, endIndex = n;
        //rotate the values for columnar nodes
        if (i == 1 || i == 3 ){endIndex = n-1;for (int j = 0; j<n; ++j){CurrPosX[j] *= (-1.0);}}
        for (int j = startindex; j<endIndex; ++j){
            Peripodialposx.push_back(CurrPosX[j]);
            Peripodialposy.push_back(CurrPosY[j]);
            PeripodialtissueType.push_back(2); //all nodes are linker at this stage
        }
    }
}


void EllipseLayoutGenerator::addPeripodialMidlines(double ABHeight, double ABLayers, double lumenHeightFrac, double peripodialHeightFrac){
    //Add the base on the edge of columnar layer:
    double d = ABHeight / ABLayers;
    double midPoint = ABHeight* (1+ 0.5*lumenHeightFrac);
    for (int j=1; j< ABLayers; j++){
        Peripodialposx.push_back(j*d-midPoint);
        Peripodialposy.push_back(0);
        PeripodialtissueType.push_back(2); //all nodes are linker at this stage

    }
    double hPeripodial = peripodialHeightFrac*ABHeight;
    peripodialLayers = ceil(hPeripodial/d );
    d = hPeripodial/peripodialLayers;
    if (peripodialLayers>1){
        //peripodial membrane is also divided into layers
        for (int j=1; j< peripodialLayers; j++){
            cout<<" distance for peripodial: "<<d<<" layer no: "<<j<<" midpoint: "<<midPoint<<endl;
            Peripodialposx.push_back(j*d+PeriR1[2]);
            Peripodialposy.push_back(0);
            PeripodialtissueType.push_back(2); //all nodes are linker at this stage
        }
    }

}

void EllipseLayoutGenerator::addEquidistantRingMin(){	
	int iterSide = 0; 		
	int iterLen, iterWid;
	int loopend = 4;
	if (symmetricX && symmetricY) {loopend = 1;}
	else if (symmetricY) {loopend = 2;}
	for (iterSide = 0; iterSide<loopend; iterSide++){
		if (iterSide == 0){iterLen = 0; iterWid = 0;}		//iterSide == 0 is putting in (+)ve x - (+)ve y 
		else if (iterSide == 1){iterLen = 1; iterWid = 0;}	//iterSide == 1 is putting in (-)ve x - (+)ve y 
		else if (iterSide == 2){iterLen = 1; iterWid = 1;}	//iterSide == 2 is putting in (-)ve x - (-)ve y 
		else if (iterSide == 3){iterLen = 0; iterWid = 1;}	//iterSide == 3 is putting in (+)ve x - (-)ve y 
		double CQ = Circumference[iterSide]/4.0;
		int n= CQ / sideLen;
		dtet[iterSide]= pi/2.0/n;
		double d = CQ / n;
		dForNodeGeneration  = d;
		double x0 = r1[iterLen];
		double y0 = 0;
		vector <double> CurrPosX, CurrPosY;
		CurrPosX.push_back(x0);
		CurrPosY.push_back(y0);
		tetha = dtet[iterSide];
		for (int i= 0; i<n-1; ++i){
			double x = r1[iterLen]*cos(tetha);
			double y = r2[iterWid]*sin(tetha);
			CurrPosX.push_back(x);
			CurrPosY.push_back(y);
			tetha += dtet[iterSide];
		}
		CurrPosX.push_back(0);
		CurrPosY.push_back(r2[iterWid]);
		//cout<<" iterSide: "<<iterSide<<" before minimiser - current size: "<<CurrPosX.size()<<endl;
        minimiseCurrentPositions(false, n , d, iterLen, iterWid, CurrPosX,CurrPosY);
		//cout<<" iterSide: "<<iterSide<<" before clean-up - current size: "<<CurrPosX.size()<<endl;
		cleanUpCurrentPositions(d, CurrPosX,CurrPosY);
		//cout<<" iterSide: "<<iterSide<<" after clean-up against current - current size: "<<CurrPosX.size()<<endl;
		n = CurrPosX.size();
		int startindex = 0, endIndex = n;
		if (iterSide == 1){endIndex = n-1;for (int i = 0; i<n; ++i){CurrPosX[i] *= (-1.0);}}
		else if (iterSide == 2){startindex = 1;for (int i = 0; i<n; ++i){CurrPosX[i] *= (-1.0);CurrPosY[i] *= (-1.0);}}
		else if (iterSide == 3){startindex = 1; endIndex = n-1;for (int i = 0; i<n; ++i){CurrPosY[i] *= (-1.0);}}
		cleanUpAgainsAlreadyExisting(d, CurrPosX,CurrPosY);
		//cout<<" iterSide: "<<iterSide<<" after clean-up against all - current size: "<<CurrPosX.size()<<endl;
		n = CurrPosX.size();	
		for (int i = startindex; i<endIndex; ++i){
			posx.push_back(CurrPosX[i]);
			posy.push_back(CurrPosY[i]);		
			tissueType.push_back(0); //all nodes are columnar at this stage
		}
	}	
}

void EllipseLayoutGenerator::addRectangle(){	
	int n = r2[0] / sideLen;
	double d = r2[0] / n;
	dForNodeGeneration  = d;
	double x0 = r1[0];
	double y0 = 0;
	vector <double> CurrPosX, CurrPosY;
	CurrPosX.push_back(x0);
	CurrPosY.push_back(y0);
	tetha = dtet[0];
	for (int i= 0; i<n-1; ++i){
		double x = r1[0];
		double y = (i+1)*d;
		CurrPosX.push_back(x);
		CurrPosY.push_back(y);
	}
	CurrPosX.push_back(r1[0]);
	CurrPosY.push_back(r2[0]);
	n = r1[0] / sideLen;
	d = r1[0] / n;
	for (int i= 0; i<n-1; ++i){
		double x = r1[0] - (i+1)*d;
		double y = r2[0];
		CurrPosX.push_back(x);
		CurrPosY.push_back(y);
	}
	CurrPosX.push_back(0);
	CurrPosY.push_back(r2[0]);
	
	//now rotating to add the remaining points:
	n = CurrPosX.size();
	for (int i = n-2; i>-1; --i){
		CurrPosX.push_back((-1.0)*CurrPosX[i]);
		CurrPosY.push_back(CurrPosY[i]);	
	}
	n = CurrPosX.size();
	for (int i = n-2; i>0; --i){
		CurrPosX.push_back(CurrPosX[i]);
		CurrPosY.push_back((-1.0)*CurrPosY[i]);		
	}
	n = CurrPosX.size();
	for (int i = 0; i<n; ++i){
		posx.push_back(CurrPosX[i]);
		posy.push_back(CurrPosY[i]);
		tissueType.push_back(0); //all nodes are columnar at this stage
	}
}

bool EllipseLayoutGenerator::updateRadia(){	
	if (r1[0]>sideLen*2.0 && r1[1]>sideLen*2.0 && r2[0] > sideLen*2.0 && r2[1] > sideLen*2.0){
		r1[0] -=sideLen;
		r1[1] -=sideLen;			
		r2[0] -=sideLen;
		r2[1] -=sideLen;
		return true; 	
	}
	else{
		if (r1[0]<=sideLen*1.5 || r1[1]<=sideLen*1.5){
			r1failed = true;
		}
		return false;	
	}
}
void EllipseLayoutGenerator::addMidLine(){
	vector <double> CurrPosX, CurrPosY;
	CurrPosX.push_back(0.0);	
	CurrPosY.push_back(0.0);
	if (r1failed){
		//r1 was the shorter axe, I will add the midline along r2 (vertical)	:
		double y = sideLen*0.8; //I want this line dense, so that the elements will be symmetric on this axis
		while (y<(r2[0]-0.5*sideLen)){
			CurrPosX.push_back(0.0);
			CurrPosY.push_back(y);
			y +=sideLen;	
		}
		y = sideLen;
		while (y<(r2[1]-0.5*sideLen)){
			CurrPosX.push_back(0.0);
			CurrPosY.push_back((-1.0)*y);
			y +=sideLen;
		}	
	}
	else{
		double x = sideLen*0.8;
		while (x<(r1[0]-0.5*sideLen)){
			CurrPosX.push_back(x);
			CurrPosY.push_back(0.0);
			x +=sideLen;
		}
		x = sideLen;
		while (x<(r1[1]-0.5*sideLen)){
			CurrPosX.push_back((-1.0)*x);
			CurrPosY.push_back(0.0);
			x +=sideLen;	
		}
	}
	//now eliminating any node that may be too close to already existing nodes:
	double threshold = dForNodeGeneration*0.7;
	threshold *= threshold;	
	int n = CurrPosX.size();	
	int i=1;	
	while (i < n-1){
		for (int j=0; j< posx.size(); ++j){
			double dx = CurrPosX[i] - posx[j];
			double dy = CurrPosY[i] - posy[j];
			double d2 = dx*dx + dy*dy;
			if (d2 < threshold){
				vector<double>::iterator it;
				it = CurrPosX.begin();
				it +=i;
				CurrPosX.erase(it);
				it = CurrPosY.begin();	
				it +=i;
				CurrPosY.erase(it);
				n = CurrPosX.size();
				i--;			
				break;				
			}		
		}
		i++;	
	}
	n = CurrPosX.size();
	for (int i=0;i<n; ++i){
		posx.push_back(CurrPosX[i]);
		posy.push_back(CurrPosY[i]);
		tissueType.push_back(0);
	}
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

    //starting the peripodial membrane
    ofstream peripodialNodesForGnuplot;
    peripodialNodesForGnuplot.open("./PeripodialNodesPreTesselation.out",ofstream::trunc);
    writePeripodialNodes2D(peripodialNodesForGnuplot);

    ofstream peripodialPointsForTesselation;
    peripodialPointsForTesselation.open("./PeripodialPoints.node",ofstream::trunc);
    peripodialPointsForTesselation<<Peripodialposx.size()<<" 2 0 1"<<endl; //dim, attribute number , border markers on or of
    for (int i =0; i< Peripodialposx.size(); ++i){
        peripodialPointsForTesselation<<i<<" "<<Peripodialposx[i]<<" "<<Peripodialposy[i]<<endl;
    }
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


void EllipseLayoutGenerator::writePeripodialNodes2D(ofstream &nodesForGnuplot){
    cout<<" in writing, number of nodes: "<<Peripodialposx.size()<<endl;
    for (int i =0; i< Peripodialposx.size(); ++i){

        double x = Peripodialposx[i];
        double y = Peripodialposy[i];
        cout<<" i: "<<i<<" x "<<x<<" y "<<y<<endl;
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

void EllipseLayoutGenerator::writeMeshFileForSimulation(double zHeight, int zLayers){
	cout<<" Writing triangle file, size of points: "<<posx.size()<<" size of tissue shape: "<<tissueType.size()<<endl;
	ofstream MeshFile;
	MeshFile.open("./MeshFile.out",ofstream::trunc);	
	MeshFile<<posx.size()*(zLayers+1)<<endl;
	cerr<<"posx.size(): "<<posx.size()<<" "<<posx.size()*(zLayers+1)<<endl;
	int n = posx.size();
	for (int i=0;i<n;++i){
		MeshFile<<posx[i]<<"	"<<posy[i]<<"	"<<0<<"    0	"<<tissueType[i]<<"	"<<atBorder[i]<<endl;	//x-coord   y-coord   z-coord   basal-node identifier(0) tissueType(columnar-0, peripodium-2) flag for border nodes
		
	}
	double dzHeight = zHeight/zLayers;
	double currzHeight = dzHeight;
	for (int layers=0; layers<zLayers; ++layers){
		int nodePositionIdentifier = 2; //0=basal. 1=apical, 2=midline
		if (layers == zLayers-1){nodePositionIdentifier=1;}
		for (int i=0;i<n;++i){
			MeshFile<<posx[i]<<"	"<<posy[i]<<"	"<<currzHeight<<"    "<<nodePositionIdentifier<<"	"<<tissueType[i]<<"	"<<atBorder[i]<<endl;	
		}
		currzHeight += dzHeight;
	}
	int nTri = triangles.size();
	MeshFile<<nTri*zLayers<<endl;
	int currOffset =0;
	for (int layers=0; layers<zLayers; ++layers){
		for (int i =0; i<nTri; ++i){
			double refpos[6][3];
			int nodes[6];
			for (int j=0;j<3;++j){			
				nodes[j]=triangles[i][j]+currOffset;
			}
			for (int j=3;j<6;++j){			
				nodes[j]=triangles[i][j-3]+currOffset+n;
			}	
			//writing Shapetype
			MeshFile<<"1	";
			//writing nodes of prism
			for (int j=0;j<6;++j){
				MeshFile<<nodes[j]<<"     ";
			}
			//writing positions for reference prism:
			//basal nodes:
			double currzHeight = dzHeight*layers;
			for (int j=0;j<3;++j){
				MeshFile<<posx[triangles[i][j]]<<"   "<<posy[triangles[i][j]]<<"   "<<currzHeight<<"    ";		
			}
			currzHeight = dzHeight*(layers+1);
			//apical nodes:
			for (int j=0;j<3;++j){
				MeshFile<<posx[triangles[i][j]]<<"   "<<posy[triangles[i][j]]<<"   "<<currzHeight<<"   ";		
			}
			MeshFile<<endl;
		}
		currOffset += n;
	}
}

void EllipseLayoutGenerator::writeTriangleMeshFileForSimulation(double zHeight){
	ofstream MeshFile;
	MeshFile.open("./TriangleFile.out",ofstream::trunc);	
	MeshFile<<posx.size()<<endl;
	int n = posx.size();
	for (int i=0;i<n;++i){
		MeshFile<<posx[i]<<"	"<<posy[i]<<"	0	    2	"<<tissueType[i]<<"	"<<atBorder[i]<<endl;	//x-coord   y-coord   z-coord   mid-node identifier(0) tissueType(columnar-0,peripodium-2) flag for border nodes
		
	}
	int nTri = triangles.size();
	MeshFile<<nTri<<endl;
	for (int i =0; i<nTri; ++i){
		double refpos[3][3];
		int nodes[3];
		for (int j=0;j<3;++j){			
			nodes[j]=triangles[i][j];
		}
			
		//writing Shapetype as truiangle (4)
		MeshFile<<"4	";
		//writing nodes of triangle
		for (int j=0;j<3;++j){
			MeshFile<<nodes[j]<<"     ";
		}
		//writing the slab height for the triangle
		MeshFile<<zHeight<<"     ";
		//writing positions for reference triangle:
		for (int j=0;j<3;++j){
			MeshFile<<posx[triangles[i][j]]<<"   "<<posy[triangles[i][j]]<<"   "<<0<<"    ";		
		}
		MeshFile<<endl;		
	}
}

void EllipseLayoutGenerator::addPeripodialMembrane(double ABHeight, double ABLayers, double lumenHeightFrac, double peripodialHeightFrac, double peripodialSideCurveFrac){
    double tissueTop = ABHeight*(1+ lumenHeightFrac + peripodialHeightFrac);
    double midPoint = ABHeight* (1+ 0.5*lumenHeightFrac);
    PeriR1[0] = tissueTop-midPoint;
    PeriR1[1] = midPoint;
    PeriR1[2] = ABHeight * 0.5 * lumenHeightFrac;
    PeriR1[3] = ABHeight * 0.5 * lumenHeightFrac;

    PeriR2[0] = ABHeight * peripodialSideCurveFrac;
    PeriR2[1] = ABHeight * 0.5 * lumenHeightFrac;

    calculatePeripodialCircumference();
    cout<<"PeriR1: "<<PeriR1[0]<<" "<<PeriR1[1]<<" "<<PeriR1[2]<<" "<<PeriR1[3]<<" PeriR2 "<<PeriR2[0]<<" "<<PeriR2[1]<<endl;
    cout<<"circumferences: "<<PeripodialCircumference[0]<<" "<<PeripodialCircumference[1]<<" "<<PeripodialCircumference[2]<<" "<<PeripodialCircumference[3]<<endl;

    calculatePeripodialdtet();
    calculatePeripodialBorderNumber();
    addPeripodialRings();
    addPeripodialMidlines(ABHeight, ABLayers, lumenHeightFrac, peripodialHeightFrac);
}



bool readinputs (int argc, char **argv, double* parameters){
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
	else if (parameters[0] == 2){
		if(argc<4){
			cerr<<"Please give lenght1, length2 of rectangle (in this order)"<<endl;
			return 0;
		}
		else{
			const char* inpstring1 = argv[2];
			parameters[1] = atof(inpstring1);
			parameters[2]  = parameters[1] ;
			const char* inpstring2 = argv[3];
			parameters[3] = atof(inpstring2);
			parameters[4]  = parameters[3];			
		}
		offset = 4;
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

int main(int argc, char **argv)
{	
	double * parameters;
	parameters = new double[9];
	bool success = readinputs (argc, argv, parameters);
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

	bool calculateNextRing = true;
	int counter =0;
	while (calculateNextRing) {
		Lay01.calculateCircumference();
		Lay01.calculateCurrentBorderNumber();
		Lay01.calculatedtet();
		if (GlobalShape == 1) {	//Generating ellipse tissue	
			Lay01.addEquidistantRingMin();	
		}
		else if (GlobalShape == 2) {	//Generating rectangular tissue
			Lay01.addRectangle();
		}
		calculateNextRing = Lay01.updateRadia();
	}    
	//outside the ring calculation, now I will add points in the middle, in horizontal or vertical,
	//depending on which axis failed first, r1 or r2:
	Lay01.addMidLine();
    double lumenHeightFrac = 0.25;
    double peripodialHeightFrac = 0.45;
    double peripodialSideCurveFrac = 0.5;
    Lay01.addPeripodialMembrane(ABHeight, ABLayers, lumenHeightFrac, peripodialHeightFrac, peripodialSideCurveFrac);
	//cout<<"before Tesselate2D, n =  "<<	Lay01.posx.size()<<endl;
	Lay01.Tesselate2D();
	Lay01.readInTesselation2D();
	Lay01.calculateAverageSideLength();

	Lay01.writeMeshFileForSimulation(ABHeight,ABLayers);	
	Lay01.writeTriangleMeshFileForSimulation(ABHeight);
	//output the points for plotting:
	int n=Lay01.posx.size();	
	//cerr<<"r1: "<<Lay01.r1[0]<<" "<<Lay01.r1[1]<<" r2: "<<Lay01.r2[0]<<" "<<Lay01.r2[1]<<endl;
	Lay01.calculateAverageSideLength();
}


