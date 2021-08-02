#ifndef ReferenceShapeBase_H
#define ReferenceShapeBase_H

#include <iostream>
#include <stdio.h>
#include <string>


using namespace std;

class ReferenceShapeBase{
protected:
	int ShapeType;		///<The integer defining the type of the shape, it is a negative value at the same basis of ShapeBase#ShapeType: prisms ShapeType = -1
	int Id;				///<The unique identifier of the reference shape, it is identical to the owner shape, it is set inside the constructor.
	int nNodes;			///<The number of nodes of the reference element, it is based on ReferenceShapeBase#ShapeType, through function ReferenceShapeBase#setNodeNumber.


	void setShapeType(string TypeName); //<The function setting the shape type of the reference element
	void setNodeNumber();				//<The function setting the number of nodes of the reference element, selection based on ReferenceShapeBase#ShapeType
public:
	double** 	Positions; 	///< The pointer to the position matrix of the reference element. The array itself is declared within the constructor, depending on ShapeBase#nDim and ShapeBase#nNodes, its size being [nNodes][nDim].
	double		Volume;		///< The volume of reference element, calculated assuming regular shapes (e.g. Perpendicular prism of equal top and bottom surfaces).
	double 		BasalArea;  ///< The basal area of the reference shape
	double 		height; 	///< Now being used as Reference element height in reference element area calculation for myosin forces. Slab height for 2D elements, value is -100 for 3D elements. Is obsolete now, must be deleted together with any 2D element option. Mesh file reading should be cleared in parallel.
	ReferenceShapeBase(string SyapeType, int Id); ///<Constructer of the ReferenceShapeBase class
	~ReferenceShapeBase();	 ///<Destructor of the ReferenceShapeBase class
	//int getShapeType();
};

#endif
