/*
 * MuscleFibre.h
 *
 *  Created on: 9 May 2017
 *      Author: melda
 */

#ifndef MUSCLEFIBRE_H_
#define MUSCLEFIBRE_H_

//#include "../TissueFolding/SourceCode/Simulation.h"

#include <vector>
using namespace std;

class MuscleFibres{
private:
	int nFibres;
	vector <double*> fibreForces;
	vector <double> initiationTimes;
	vector <double> maturationTimes;
	vector <double*> applicationAngles;
	vector <double> forceFractions;


	vector <int*> nodeIdsFirstQuadrant;
	vector <int*> nodeIdsSecondQuadrant;
	vector <int*> nodeIdsThirdQuadrant;
	vector <int*> nodeIdsFourthQuadrant;
	double totalMassFirstQuadrant;
	double totalMassSecondQuadrant;
	double totalMassThirdQuadrant;
	double totalMassFourthQuadrant;


	void calculateCurrentForceFractions(int simTimeInHr){
		//cout<<" inside calculateCurrentForceFractions: "<<endl;
		for (int fibreIndex =0; fibreIndex< nFibres; ++fibreIndex){
			//cout<<" calculating fibre "<<fibreIndex<<endl;
			//cout<<" time range: "<<initiationTimes[fibreIndex]<<" " << maturationTimes[fibreIndex]<<endl;
			if (simTimeInHr < initiationTimes[fibreIndex]){
				forceFractions[fibreIndex] = 0.0;
			}
			else if (simTimeInHr > maturationTimes[fibreIndex]){
				forceFractions[fibreIndex] = 1.0;
			}
			else{
				forceFractions[fibreIndex] = (simTimeInHr-initiationTimes[fibreIndex])/(maturationTimes[fibreIndex]-initiationTimes[fibreIndex]);
			}
		}
		//cout<<" finished  calculateCurrentForceFractions"<<endl;
	}

	void calculateTotalWeightOfQuadrants(vector<Node*> &Nodes){
		cout<<" inside calculateTotalWeightOfQuadrants: "<<endl;
		totalMassFirstQuadrant  = 0;
		totalMassSecondQuadrant = 0;
		totalMassThirdQuadrant  = 0;
		totalMassFourthQuadrant = 0;
		int n = nodeIdsFirstQuadrant.size();
		for (int i=0; i<n; ++i){
			totalMassFirstQuadrant += Nodes[i]->mass;
		}
		n = nodeIdsSecondQuadrant.size();
		for (int i=0; i<n; ++i){
			totalMassSecondQuadrant += Nodes[i]->mass;
		}
		n = nodeIdsThirdQuadrant.size();
		for (int i=0; i<n; ++i){
			totalMassThirdQuadrant += Nodes[i]->mass;
		}
		n = nodeIdsFourthQuadrant.size();
		for (int i=0; i<n; ++i){
			totalMassFourthQuadrant += Nodes[i]->mass;
		}
	}

	void freeVector(vector <int*> &vectorToDelete){
		while(!vectorToDelete.empty()){
			int* tmp_IntPtr;
			tmp_IntPtr = vectorToDelete.back();
			vectorToDelete.pop_back();
			delete tmp_IntPtr;
		}
	}

	void freeVector(vector <double*> &vectorToDelete){
		while(!vectorToDelete.empty()){
			double* tmp_DoublePtr;
			tmp_DoublePtr = vectorToDelete.back();
			vectorToDelete.pop_back();
			delete tmp_DoublePtr;
		}
	}
public:
	MuscleFibres(int nFibres, vector<double> orientationAngles, vector <double> fibreForceMagnitudes, vector <double> initTimesInHr, vector <double > matureTimesinHr, vector<double> initAngle, vector<double> endAngle ){
		this->nFibres = nFibres;
		for (int i=0; i<nFibres; ++i){
			//calculate the force for this fibre:
			//convert angle to radians:
			double angle = orientationAngles[i]/180.0*M_PI;
			double ForceMag = fibreForceMagnitudes[i];
			double * F = new double[2];
			F[0] = ForceMag*cos(angle);
			F[1] = ForceMag*sin(angle);
			this->fibreForces.push_back(F);
			this->initiationTimes.push_back(initTimesInHr[i]);
			this->maturationTimes.push_back(matureTimesinHr[i]);
			double * angleRange = new double[2];
			//convert to radian:
			angleRange[0] = initAngle[i]/180.0*M_PI;
			angleRange[1] = endAngle[i]/180.0*M_PI;
			applicationAngles.push_back( angleRange);
			forceFractions.push_back(0.0);
		}
		totalMassFirstQuadrant  = 0;
		totalMassSecondQuadrant = 0;
		totalMassThirdQuadrant  = 0;
		totalMassFourthQuadrant = 0;
	}

	~MuscleFibres(){
		freeVector(fibreForces);
		freeVector(applicationAngles);
		freeVector(nodeIdsFirstQuadrant);
		freeVector(nodeIdsSecondQuadrant);
		freeVector(nodeIdsThirdQuadrant);
		freeVector(nodeIdsFourthQuadrant);
	}

	void assignNodesForFibreAttachement(vector<Node*> &Nodes, bool symmetricX, bool symmetricY, double systemX, double systemY){
		for (vector<Node*>::iterator itNode=Nodes.begin(); itNode<Nodes.end(); ++itNode){
			if ((*itNode)->atCircumference){
				//the node is at circumference, what is its position at the circumference:
				double dx = (*itNode)->Position[0] - systemX;
				double dy = (*itNode)->Position[1] - systemY;
				double x1 = 1.0, y1 = 0; //get the angle the vector is making with x axis
				double dot = x1*dx + y1*dy;     // dot product between [x1, y1] and [x2, y2]
				double det = x1*dy - y1*dx;      // determinant
				double angle = atan2(det, dot);  // atan2(y, x) or atan2(sin, cos)
				for (int i =0; i< nFibres; ++i){
					if (angle >= applicationAngles[i][0] && angle <= applicationAngles[i][1]){
						int * fibreAndNodeId = new int[2];
						fibreAndNodeId[0] = i;
						fibreAndNodeId[1] = (*itNode)->Id;
						nodeIdsFirstQuadrant.push_back(fibreAndNodeId);
					}
					if (!symmetricX){
						//there is no x symmetricity, I need to check the angles for this range too
						double newBoundary[2] = {M_PI - applicationAngles[i][1], M_PI - applicationAngles[i][0]};
						if (angle >= newBoundary[0] && angle <= newBoundary[1]){
							int * fibreAndNodeId = new int[2];
							fibreAndNodeId[0] = i;
							fibreAndNodeId[1] = (*itNode)->Id;
							nodeIdsSecondQuadrant.push_back(fibreAndNodeId);
						}
					}
					if (!symmetricY){
						double newBoundary[2] = {2.0*M_PI - applicationAngles[i][1], 2.0*M_PI - applicationAngles[i][0]};
						if (angle >= newBoundary[0] && angle <= newBoundary[1]){
							int * fibreAndNodeId = new int[2];
							fibreAndNodeId[0] = i;
							fibreAndNodeId[1] = (*itNode)->Id;
							nodeIdsFourthQuadrant.push_back(fibreAndNodeId);
						}
					}

					if (!symmetricX && !symmetricY){
						double newBoundary[2] = {M_PI + applicationAngles[i][0], M_PI + applicationAngles[i][1]};
						if (angle >= newBoundary[0] && angle <= newBoundary[1]){
							int * fibreAndNodeId = new int[2];
							fibreAndNodeId[0] = i;
							fibreAndNodeId[1] = (*itNode)->Id;
							nodeIdsThirdQuadrant.push_back(fibreAndNodeId);
						}
					}
				}
			}
		}
	}

	void addFibreForces(vector<Node*> &Nodes, gsl_matrix* gExt, double simTimeInHr){
		calculateTotalWeightOfQuadrants(Nodes);
		calculateCurrentForceFractions(simTimeInHr);
		int n = nodeIdsFirstQuadrant.size();
		//cout<<" n 1st quadrant: "<<n<<endl;
		for (int i=0; i<n; ++i){
			int fibreId = nodeIdsFirstQuadrant[i][0];
			if (forceFractions[fibreId] == 0.0 ){
				continue;
			}
			int nodeId =  nodeIdsFirstQuadrant[i][1];
			double nodeMass = Nodes[nodeId]->mass;
			double Fx = forceFractions[fibreId] * nodeMass * fibreForces[fibreId][0]/totalMassFirstQuadrant;
			double Fy = forceFractions[fibreId] * nodeMass * fibreForces[fibreId][1]/totalMassFirstQuadrant;
			cout<<" First quadrant, node "<<nodeId<<" force : "<<Fx <<" "<<Fy<<endl;
			int indice =nodeId*3;
			Fx += gsl_matrix_get(gExt,indice,0);
			gsl_matrix_set(gExt,indice,0,Fx);
			Fy += gsl_matrix_get(gExt,indice+1,0);
			gsl_matrix_set(gExt,indice+1,0,Fy);
		}
		n = nodeIdsSecondQuadrant.size();
		//cout<<" n 2nd quadrant: "<<n<<endl;
		for (int i=0; i<n; ++i){
			int fibreId = nodeIdsSecondQuadrant[i][0];
			if (forceFractions[fibreId] == 0.0 ){
				continue;
			}
			int nodeId =  nodeIdsSecondQuadrant[i][1];
			double nodeMass = Nodes[nodeId]->mass;
			//in the second quadrant, the force component in x is flipped:
			double Fx = forceFractions[fibreId] * nodeMass * (-1.0) * fibreForces[fibreId][0]/totalMassSecondQuadrant;
			double Fy = forceFractions[fibreId] * nodeMass * fibreForces[fibreId][1]/totalMassSecondQuadrant;
			int indice =nodeId*3;
			Fx += gsl_matrix_get(gExt,indice,0);
			gsl_matrix_set(gExt,indice,0,Fx);
			Fy += gsl_matrix_get(gExt,indice+1,0);
			gsl_matrix_set(gExt,indice+1,0,Fy);
		}
		n = nodeIdsThirdQuadrant.size();
		//cout<<" n 3rd quadrant: "<<n<<endl;
		for (int i=0; i<n; ++i){
			int fibreId = nodeIdsThirdQuadrant[i][0];
			if (forceFractions[fibreId] == 0.0 ){
				continue;
			}
			int nodeId =  nodeIdsThirdQuadrant[i][1];
			double nodeMass = Nodes[nodeId]->mass;
			//in the second quadrant, the force component in y is flipped:
			double Fx = forceFractions[fibreId] * nodeMass * fibreForces[fibreId][0]/totalMassThirdQuadrant;
			double Fy = forceFractions[fibreId] * nodeMass * (-1.0) * fibreForces[fibreId][1]/totalMassThirdQuadrant;
			int indice =nodeId*3;
			Fx += gsl_matrix_get(gExt,indice,0);
			gsl_matrix_set(gExt,indice,0,Fx);
			Fy += gsl_matrix_get(gExt,indice+1,0);
			gsl_matrix_set(gExt,indice+1,0,Fy);
		}
		n = nodeIdsFourthQuadrant.size();
		//cout<<" n 4th quadrant: "<<n<<endl;
		for (int i=0; i<n; ++i){
			int fibreId = nodeIdsFourthQuadrant[i][0];
			if (forceFractions[fibreId] == 0.0 ){
				continue;
			}
			int nodeId =  nodeIdsFourthQuadrant[i][1];
			double nodeMass = Nodes[nodeId]->mass;
			//in the second quadrant, the force components in x&y are flipped:
			double Fx = forceFractions[fibreId] * nodeMass * (-1.0) * fibreForces[fibreId][0]/totalMassFourthQuadrant;
			double Fy = forceFractions[fibreId] * nodeMass * (-1.0) * fibreForces[fibreId][1]/totalMassFourthQuadrant;
			int indice =nodeId*3;
			Fx += gsl_matrix_get(gExt,indice,0);
			gsl_matrix_set(gExt,indice,0,Fx);
			Fy += gsl_matrix_get(gExt,indice+1,0);
			gsl_matrix_set(gExt,indice+1,0,Fy);
		}
	}
};



#endif /* MUSCLEFIBRE_H_ */
