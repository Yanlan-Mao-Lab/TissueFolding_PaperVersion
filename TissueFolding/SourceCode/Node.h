#ifndef Node_H
#define Node_H

#include <stdio.h>
#include <iostream>
#include <vector>

using namespace std;
class ShapeBase;
/**
 *  The node class
 *  */
class Node{
private:

public:
	Node(int id, int dim, double* pos, int tissuePos, int tissueType);	//<Constructer of the node
	~Node();
	bool 			*FixedPos;							///< The boolean array stating if the node's position is fixed in any direction, format: [x y z], true = fixed
	int 			Id;									///< The unique identification number of the node
	int 			nDim;								///< The number of dimensions of the node, (2 or 3)
	double 			*Position;							///< The pointer to the position array of the node. The array itself is declared within the constructor, depending on nDim
	double 			*initialPosition;					///< The pointer to the initial position array of the node. The array itself is declared within the constructor, depending on nDim.
	double 			*RKPosition;						///< The pointer to the position array for position during a Runge-Kutta step array of the node. The array itself is declared within the constructor, depending on nDim
	bool			externalViscositySetInFixing[3];	///< The boolean array stating if the external viscosity of any axis has been set in node fixing options. The node fixing is carried out before the physical parameter settings in most cases. The boolean check is carried out not to overwrite the existing set viscosity in normal viscosity assignment.
	double 			externalViscosity[3];				///< External viscosity of the node, defined by its placement within the tissue. This can be defined as an external adhesion, ECM remodelling, or any other form of viscosity.
	double 			initialExternalViscosity[3];		///< External viscosity of the node before any modifications by perturbation (for example ECM stiffness/ viscosity change)
	double 			maximumExternalViscosity[3];
	double 			minimumExternalViscosity[3];
	double 			ECMViscosityChangePerHour[3];		///< The change in ECM viscosity per one hour. The double array of size (1,3) stores the viscosity change in x, y, and z directions, respectively.
	double			displacement;						///< the displacement of the node from previous time step;
	int 			tissuePlacement;					///< The tissue placement is 0 for basal nodes, 1 for apical nodes, and 2 for middle range
	int 			tissueType;		 					///< The tissue type is 0 for columnar layer, 1 for peripodial membrane, and 2 for linker zone
	bool 			atCircumference;					///< Boolean defining if the node is at the circumference of the columnar layer of the tissue.
	double 			mass;								///< The mass of the node, calculated via the elements that use the node as a vertex
	//double 			surface;						///< The surface of the node, calculated via apical or basal elements. Lateral surfaces are not included
	double			viscositySurface;					///< The surface of the node, calculated for application of external viscosity via surface. It is positive for a surface that is to feel external viscosity, zero otherwise.
	double  		zProjectedArea;         			///< The surface of the node, as projected in Z, calculated from apical or pasal surfces of elements, lateral surfaces are not included.
   	vector <int> 	immediateNeigs;						///< The list of Id's for immediate neighbours of the node, i.e. the nodes that are shared by the elements that utilise the owner of this node.
	vector <int> 	connectedElementIds;				///< The list of Id's for elements that are utilising this node.
	vector <double>	connectedElementWeights;			///< The list of weights (normalised mass) for elements that are utilising this node, order is linked to Node#connectedElementIds.
	int 			symmetricEquivalentId;				///< The id of the node that this node will move symmetrically in y, if the tissue is symmetric and only half is simulated.
	bool			hasLateralElementOwner;				///< The boolean stating if any lateral element uses this node
	bool			atSymmetricityBorder;				///< The boolean stating if the node is at the border of symmetricity
	bool			insideEllipseBand;					///< The boolean stating if the node is inside a marker ellipse
	int				coveringEllipseBandId;				///< The id of the marker ellipse that the node is covered by. If hte node is not covered by any marker ellipse, the value is -1.
	bool			allOwnersECMMimicing;				///< The boolean stating all the elements making use of the node are ECM mimicking elements.
	bool			isMaster[3];
	int				slaveTo[3];
	vector<int>		collapsedWith;
	int 			adheredTo;
	bool			attachedToPeripodial;
	bool			onFoldInitiation;
	bool			checkOwnersforEllipseAsignment;    ///< When the ellipse ids are assigned to nodes nby owner elemetns, you can end up with elements that have all their nodes engulfed in an ellipse, but the element itself is not assigned into an ellipse. This flag will check for that.
	bool			positionUpdateOngoing;
	bool 			positionUpdateCounter;
	bool haveCommonOwner(Node* nodeSlave);
	int  getCommonOwnerId(Node* nodeSlave);
	void setExternalViscosity(double ApicalVisc,double BasalVisc, bool extendExternalViscosityToInnerTissue);///< The function to set the viscosity of the node.
	bool checkIfNeighbour(int IdToCheck); 				///< The function to check if the node with input Id (IdToCheck) is an immediate neighbour of the node
	bool checkIfNodeHasPacking();						///< The function to check if the node is eligible for packing.
	void getCurrentPosition(double* pos);				///< return the current position of the node
	void displayConnectedElementIds();					///< This function will print out a list of connected element Id's
	void displayConnectedElementWeights();				///< This function will print out the weights of the connected elements, in the order of  Id s given in connectedElementIds
	void displayPosition();
	void addToImmediateNeigs(int newNodeId);			///< This function adds the input node Id (int) to the list of neighbours of the node (Node#immediateNeigs)
	void addToConnectedElements(int newElementId, double volumePerNode);	///< This function adds the input newElementId (int) to the list of elements connected by this node, updating the mass, and weights of mass per connected element in the process.
	void removeFromConnectedElements(int newElementId, double volumePerNode);///< This function removes the input newElementId (int) from the list of elements connected by this node, updating the mass, and weights of mass per connected element in the process.
	bool isMyNeig(int nodeId);
	bool isNeigWithMyCollapsedNodes(int NodeId, vector<Node*>& Nodes);
	void getNewCollapseListAndAveragePos(vector<int> &newCollapseList, double* avrPos, bool* fix, vector<Node*>& Nodes, int masterNodeId);
	void clearDuplicatesFromCollapseList();
	void collapseOnNode(vector<int> &newCollapseList, double* avrPos, bool* fix, vector<Node*>& Nodes, int masterNodeId);
	void collapseOnNode(vector<Node*>& Nodes, int masterNodeId);
	void collapseOnNodeInStages( vector<int> &newCollapseList, double* avrPos, bool* fix, vector<Node*>& Nodes, int masterNodeId);
	void updatePositionTowardsPoint(double* avrPos,bool* fix);
	bool isECMChangeAppliedToNode(bool changeApicalECM, bool changeBasalECM, vector<int> &ECMChangeEllipseBandIds, int numberOfECMChangeEllipseBands);
	int getId();										///< The function returns the Id (Node#Id) of the node
};
#endif
