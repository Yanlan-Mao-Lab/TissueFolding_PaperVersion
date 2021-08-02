/*
 * NewtonRaphsonSolver.cpp
 *
 *  Created on: 26 Apr 2016
 *      Author: melda
 */

#include "NewtonRaphsonSolver.h"
#include <math.h>

//#include "Node.h"
//#include <gsl/gsl_linalg.h>

using namespace std;

NewtonRaphsonSolver::NewtonRaphsonSolver(int dim, int n){
    threshold = 1E-8;

	nDim = dim;
	nNodes = n;
	externalViscosityVolumeBased = false; //using external surface
	numericalParametersSet = false;

	un = gsl_matrix_calloc(nDim*nNodes,1);
	ge = gsl_matrix_calloc(nDim*nNodes,1);
	gvInternal = gsl_matrix_calloc(nDim*nNodes,1);
	gvExternal = gsl_matrix_calloc(nDim*nNodes,1);
	gExt = gsl_matrix_calloc(nDim*nNodes,1);
	gSum = gsl_vector_calloc(nDim*nNodes);
	uk = gsl_matrix_calloc(nDim*nNodes,1);
	displacementPerDt = gsl_matrix_calloc(nDim*nNodes,1);
	deltaU = gsl_vector_calloc(nDim*nNodes);
	K = gsl_matrix_calloc(nDim*nNodes,nDim*nNodes);
	Knumerical = 0;
	boundNodesWithSlaveMasterDefinition = false;
}

NewtonRaphsonSolver::~NewtonRaphsonSolver(){
	cout<<"inside deleting NR solver"<<endl;

	gsl_matrix_free(uk);
	cout<<" deleted uk"<<endl;
    gsl_matrix_free(un);
	cout<<" deleted un"<<endl;
	gsl_matrix_free(ge);
	cout<<" deleted ge"<<endl;
	gsl_matrix_free(gvExternal);
	cout<<" deleted gvExternal"<<endl;
	gsl_matrix_free(gvInternal);
	cout<<" deleted gvInternal"<<endl;
	gsl_matrix_free(gExt);
	cout<<" deleted gExt"<<endl;
	gsl_vector_free(gSum);
	cout<<" deleted gSum"<<endl;
	gsl_matrix_free(displacementPerDt);
	cout<<" deleted displacementPerDt"<<endl;
	gsl_vector_free(deltaU);
	cout<<" deleted deltaU"<<endl;
	//cout<<"K(0,0): "<<gsl_matrix_get(K,0,0)<<endl;
	//cout<<"(K->size1,K->size2): "<<K->size1<<" "<<K->size2<<endl;
	gsl_matrix_free(K);
	cout<<" deleted K"<<endl;
	if (Knumerical != 0){
		gsl_matrix_free(Knumerical);
	}
	else{
		free(Knumerical);
	}
	cout<<" deleted Knumerical"<<endl;
}

void NewtonRaphsonSolver::reInitiateMatricesAfterRefinement(int n){
	//I have added nodes to the system after refinement.
	//I need to change the number of nodes recorded in NR solver, and
	//re-allocate the matrices with new node size.
	nNodes = n;

	gsl_matrix_free(un);
	gsl_matrix_free(ge);
	gsl_matrix_free(gvInternal);
	gsl_matrix_free(gvExternal);
	gsl_matrix_free(gExt);
	gsl_vector_free(gSum);
	gsl_matrix_free(uk);
	gsl_matrix_free(displacementPerDt);
	gsl_vector_free(deltaU);
	gsl_matrix_free(K);

	un = gsl_matrix_calloc(nDim*nNodes,1);
	ge = gsl_matrix_calloc(nDim*nNodes,1);
	gvInternal = gsl_matrix_calloc(nDim*nNodes,1);
	gvExternal = gsl_matrix_calloc(nDim*nNodes,1);
	gExt = gsl_matrix_calloc(nDim*nNodes,1);
	gSum = gsl_vector_calloc(nDim*nNodes);
	uk = gsl_matrix_calloc(nDim*nNodes,1);
	displacementPerDt = gsl_matrix_calloc(nDim*nNodes,1);
	deltaU = gsl_vector_calloc(nDim*nNodes);
	K = gsl_matrix_calloc(nDim*nNodes,nDim*nNodes);
}


void NewtonRaphsonSolver::setMatricesToZeroAtTheBeginningOfIteration(bool thereIsNumericalCalculation){
	gsl_matrix_set_zero(un);
	gsl_matrix_set_zero(ge);
	gsl_matrix_set_zero(gvInternal);
	gsl_matrix_set_zero(gvExternal);
	gsl_matrix_set_zero(gExt);
	gsl_vector_set_zero(gSum);
	gsl_matrix_set_zero(uk);
	gsl_matrix_set_zero(displacementPerDt);
	gsl_vector_set_zero(deltaU);
	gsl_matrix_set_zero(K);
	if (thereIsNumericalCalculation){
		if (!numericalParametersSet){
			Knumerical = gsl_matrix_calloc(K->size1,K->size2);
			numericalParametersSet = true;
		}
		else{
			gsl_matrix_set_zero(Knumerical);
		}
	}
}

void NewtonRaphsonSolver::setMatricesToZeroInsideIteration(){
	gsl_matrix_set_zero(ge);
	gsl_matrix_set_zero(gvInternal);
	gsl_matrix_set_zero(gvExternal);
	gsl_vector_set_zero(gSum);
	gsl_matrix_set_zero(gExt);
	gsl_matrix_set_zero(K);


}

void NewtonRaphsonSolver::constructUnMatrix(vector <Node*>& Nodes){
    for (int i = 0; i<nNodes; ++i ){
        for (int j=0; j<nDim; ++j){
            gsl_matrix_set(un,3*i+j,0,Nodes[i]->Position[j]);
        }
    }
}

void NewtonRaphsonSolver::initialteUkMatrix(){
    gsl_matrix_memcpy(uk,un);
}


void NewtonRaphsonSolver::calculateDisplacementMatrix(double dt){
	gsl_matrix_memcpy(displacementPerDt,uk);
	gsl_matrix_sub(displacementPerDt,un);
	gsl_matrix_scale(displacementPerDt,1.0/dt);
}

void NewtonRaphsonSolver::calculateBoundKWithSlavesMasterDoF(){
	if (boundNodesWithSlaveMasterDefinition){

		/*int n =slaveMasterList.size();
		cout<<" number of master/slave couples: "<<n<<endl;
		for (vector< vector<int> >::iterator itSlaveMasterCouple=slaveMasterList.begin(); itSlaveMasterCouple<slaveMasterList.end(); ++itSlaveMasterCouple){
			//add all forces on slave to master, set slave forces to zero (g_bound = N^T g)
			int slaveIndex = (*itSlaveMasterCouple)[0];
			int masterIndex = (*itSlaveMasterCouple)[1];
			cout<<" slave/master: "<<slaveIndex<<" 	"<<masterIndex<<endl;
		}*/

		//cout<<" binding slaves to masters"<<endl;
		//Original calculation necessary is carried out by matrices N and I, where:
		//N = matrix of size (nDim*nNodes, nDim*nNodes)
		//I = matrix of size (nDim*nNodes, nDim*nNodes)
		//N is identity except:
		//	For degrees of freedom  i that are slaves, the diagonal of N should be zero;
		//	For degrees of freedom  i that are slaves to degrees of freedom  j, N(i,j) should be 1.
		//I is zero matrix except:
		//	For degrees of freedom  i that are slave, I(i,i) should be 1
		//Then,
		//	K_bound = N^T K N + I
		//	g_bound = N^T g
		//Due to efficient memory usage, these operations will be carried out on a row/column basis, rather than
		//using full matrices
		//	In other words, I do not wish to allocate two more (nDim*nNodes, nDim*nNodes) matrices in definition, and
		//	another (nDim*nNodes, nDim*nNodes) matrix as temporary during  N^T K N calculation

		//int n= slaveMasterList.size();
		int totalNumberOfDoF = K->size1; //nDim * nNodes
		//for( int slaveIterator =0; slaveIterator <n; ++slaveIterator){
		for (vector< vector<int> >::iterator itSlaveMasterCouple=slaveMasterList.begin(); itSlaveMasterCouple<slaveMasterList.end(); ++itSlaveMasterCouple){
			//add all forces on slave to master, set slave forces to zero (g_bound = N^T g)
			int slaveIndex = (*itSlaveMasterCouple)[0];
			int masterIndex = (*itSlaveMasterCouple)[1];
			//if (slaveIndex == 17631 || slaveIndex == 17633){
			//	continue;
			//}
			double slaveForce = gsl_vector_get(gSum,slaveIndex);
			double masterForce = gsl_vector_get(gSum,masterIndex);
			masterForce +=  slaveForce;
			gsl_vector_set(gSum,slaveIndex,0);
			gsl_vector_set(gSum,masterIndex,masterForce);
			//Start K manipulation - K_bound = N^T K N + I
			//add all elements of the slave row to master row - N^T K
			// format of view: [ matrix*, origin i, origin j, rows, columns ]
			gsl_matrix_view KSlaveRow = gsl_matrix_submatrix (K, slaveIndex, 0, 1, totalNumberOfDoF);
			gsl_matrix_view KMasterRow = gsl_matrix_submatrix (K, masterIndex, 0, 1, totalNumberOfDoF);
			gsl_matrix_add(&(KMasterRow.matrix),&(KSlaveRow.matrix));
			gsl_matrix_set_zero(&(KSlaveRow.matrix));
			//add all elements of the slave column to master column - N^T K N
			gsl_matrix_view KSlaveColumn = gsl_matrix_submatrix (K, 0, slaveIndex, totalNumberOfDoF,1);
			gsl_matrix_view KMasterColumn = gsl_matrix_submatrix (K,0, masterIndex, totalNumberOfDoF,1);
			gsl_matrix_add(&(KMasterColumn.matrix),&(KSlaveColumn.matrix));
			gsl_matrix_set_zero(&(KSlaveColumn.matrix));
			//make K(slave,slave) equal to 1 - N^T K N + I
			gsl_matrix_set(K,slaveIndex,slaveIndex,1);
			//double forceOnSlave = gsl_vector_get(gSum,slaveIndex);
			//cout<<"slave index: "<<slaveIndex<<" force: "<<forceOnSlave<<endl;
			//views do not need to be freed, they are addresses to original matrices
		}
	}
}

void NewtonRaphsonSolver::equateSlaveDisplacementsToMasters(){
	int n = slaveMasterList.size();
	for( int slaveIterator = 0; slaveIterator <n; ++slaveIterator){
		int slaveIndex = slaveMasterList[slaveIterator][0];
		int masterIndex = slaveMasterList[slaveIterator][1];
		//if (slaveIndex == 17631 || slaveIndex == 17633){
		//	continue;
		//}
		double masterDisplacement = gsl_vector_get(deltaU,masterIndex);
		//double slaveDisplacement =  gsl_vector_get(deltaU,slaveIndex);
		//cout<<" slave-master ("<<slaveIndex<<" "<<masterIndex<<") displacement before correction: "<< slaveDisplacement<<" "<<masterDisplacement<<endl;
		gsl_vector_set(deltaU,slaveIndex,masterDisplacement);
		//masterDisplacement = gsl_vector_get(deltaU,masterIndex);
		//slaveDisplacement =  gsl_vector_get(deltaU,slaveIndex);
		//cout<<" slave-master ("<<slaveIndex<<" "<<masterIndex<<")displacement after correction: "<< slaveDisplacement<<" "<<masterDisplacement<<endl;
	}
}

void NewtonRaphsonSolver::calcutateFixedK(vector <Node*>& Nodes){
    int dim = 3;
    int Ksize = K->size1;
    for(int i=0; i<nNodes; i++){
        for (int j=0; j<dim; ++j){
            //if ( (Xfixed && j == 0) || (Yfixed && j == 1) || (Zfixed && j == 2) ){
            if (Nodes[i]->FixedPos[j]){
                int index1 = i*dim+j;
                gsl_vector_set(gSum,index1,0.0); // making the forces zero
                for (int k =0; k<Ksize; ++k){
                    double value =0.0;
                    if (index1 == k ){value =1.0;}
                    gsl_matrix_set(K, index1, k, value);
                    gsl_matrix_set(K, k, index1, value); //K is symmetric;
                }
            }
        }
    }
}

void NewtonRaphsonSolver::calculateForcesAndJacobianMatrixNR(vector <Node*>& Nodes, vector <ShapeBase*>& Elements, double dt, bool recordForcesOnFixedNodes, double **FixedNodeForces){
	const int maxThreads = omp_get_max_threads();
	omp_set_num_threads(maxThreads);
	#pragma omp parallel for //private(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces, outputFile, dt)
		for( vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
			if (!(*itElement)->IsAblated){
				(*itElement)->calculateForces(Nodes, displacementPerDt, recordForcesOnFixedNodes, FixedNodeForces);
			}
			//cout<<"finished calculating forces in NR"<<endl;
			(*itElement)->calculateImplicitKElastic(); //This is the stiffness matrix, elastic part of the jacobian matrix
			//cout<<"finished calculating ImplicitK elastic in NR"<<endl;
			(*itElement)->calculateImplicitKViscous(displacementPerDt, dt); //This is the viscous part of jacobian matrix
			//cout<<"finished calculating ImplicitK viscous in NR"<<endl;
		}
}

void NewtonRaphsonSolver::writeForcesTogeAndgvInternal(vector <Node*>& Nodes, vector <ShapeBase*>& Elements, double** SystemForces){
    for(vector<ShapeBase*>::iterator  itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
        if (!(*itElement)->IsAblated){
        	(*itElement)->writeInternalForcesTogeAndgv(ge,gvInternal,SystemForces,Nodes);
        }
    }
}

void NewtonRaphsonSolver::writeImplicitElementalKToJacobian(vector <ShapeBase*>& Elements){
    //writing all elements K values into big K matrix:
    for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	//if element is ablated, current elemental K matrix will be identity
    	//(*itElement)->writeKelasticToMainKatrix(K);
    	//activate this for internal viscosity
    	(*itElement)->writeKviscousToMainKatrix(K);
    }
    //Elements[0]->displayMatrix(K,"normalKonlyViscous");
    for(vector<ShapeBase*>::iterator itElement=Elements.begin(); itElement<Elements.end(); ++itElement){
    	//if element is ablated, current elemental K matrix will be identity
    	(*itElement)->writeKelasticToMainKatrix(K);
    	//activate this for internal viscosity
    }
}

void NewtonRaphsonSolver::calculateExternalViscousForcesForNR(vector <Node*>& Nodes){
    //the mass is already updated for symmetricity boundary nodes, and the viscous forces will be calculated correctly,
	//as mvisdt is already accounting for the doubling of mass
    //gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, mvisc, displacementPerDt,0.0, gvExternal);
	for (int i = 0; i<nNodes; ++i ){
		for (int j=0; j<nDim; ++j){
			//double matrixValue = gsl_matrix_get(mvisc,0,3*i+j);
			if (externalViscosityVolumeBased){
				double massTimesViscosity = Nodes[i]->mass*Nodes[i]->externalViscosity[j];
				double displacementValue = gsl_matrix_get(displacementPerDt,3*i+j,0);
				gsl_matrix_set(gvExternal,3*i+j,0,massTimesViscosity*displacementValue);
			}
			else{
				double surfaceAreaTimesViscosity = Nodes[i]->viscositySurface*Nodes[i]->externalViscosity[j];
				double displacementValue = gsl_matrix_get(displacementPerDt,3*i+j,0);
				gsl_matrix_set(gvExternal,3*i+j,0,surfaceAreaTimesViscosity*displacementValue);
				if (isnan(surfaceAreaTimesViscosity)){
					cout<<" node: "<<i<<" dimention: "<<j<<" surfaceAreaTimesViscosity is nan: "<<surfaceAreaTimesViscosity<<endl;
				}
				if (isnan(displacementValue)){
					cout<<" node: "<<i<<" dimention: "<<j<<" displacementValue is nan: "<<displacementValue<<endl;
				}
			}
		}
	}
	//added this line with sign change
    gsl_matrix_scale(gvExternal,-1.0);
}

void NewtonRaphsonSolver::addImplicitKViscousExternalToJacobian(vector <Node*>& Nodes, double dt){
    //gsl_matrix_add(K,mviscPerDt);
	//cout<<"Jacobian due to explicit viscosity" <<endl;
	for (int i = 0; i<nNodes; ++i ){
		if (externalViscosityVolumeBased){
			double massPerDt = Nodes[i]->mass/dt;
			for (int j=0; j<nDim; ++j){
				double curKValue = gsl_matrix_get(K,3*i+j,3*i+j);
				double massTimesViscosityPerDt = massPerDt*Nodes[i]->externalViscosity[j];
				//double matrixValuePerDt = gsl_matrix_get(mvisc,0,3*i+j)/dt;
				gsl_matrix_set(K,3*i+j,3*i+j,massTimesViscosityPerDt+curKValue);
			}
		}
		else{
			double surfaceAreaPerDt = Nodes[i]->viscositySurface/dt;
			for (int j=0; j<nDim; ++j){
				double curKValue = gsl_matrix_get(K,3*i+j,3*i+j);
				double surfaceAreaTimesViscosityPerDt = surfaceAreaPerDt*Nodes[i]->externalViscosity[j];
				gsl_matrix_set(K,3*i+j,3*i+j,surfaceAreaTimesViscosityPerDt+curKValue);
			}
		}
	}
}

void NewtonRaphsonSolver::checkJacobianForAblatedNodes(vector <int> & AblatedNodes){
	int nAblatedNode = AblatedNodes.size();
	for (int a = 0; a<nAblatedNode; ++a){
		int NodeId = AblatedNodes[a]*3;
		for (int aa= 0; aa<3; ++aa){
			double Kdiagonal = gsl_matrix_get(K,NodeId+aa,NodeId+aa);
			if (Kdiagonal == 0){
				gsl_matrix_set(K,NodeId+aa,NodeId+aa,1);
			}
		}
	}
}

void NewtonRaphsonSolver::calculateSumOfInternalForces(){
	for (int i=0; i<nDim*nNodes; ++i){
		gsl_vector_set(gSum,i,gsl_matrix_get(ge,i,0)+gsl_matrix_get(gvInternal,i,0)+gsl_matrix_get(gvExternal,i,0));
	}
}

void NewtonRaphsonSolver::addExernalForces(){
	//cout<<"after random forces addition"<<endl;
	for (int i=0; i<nDim*nNodes; ++i){
		gsl_vector_set(gSum,i,gsl_vector_get(gSum,i)+gsl_matrix_get(gExt,i,0));
		double value = gsl_vector_get(gSum,i);
		if (isnan(value)){
		      cout<<" gSUM is nan at matrix point: "<<i<<endl;
		}
	}
}



void NewtonRaphsonSolver::solveForDeltaU(){
    const int nmult  = nDim*nNodes;
    int *ia = new int[nmult+1];
    double *b = new double[nmult];
    vector <int> ja_vec;
    vector <double> a_vec;
    constructiaForPardiso(ia, nmult, ja_vec, a_vec);
    const int nNonzero = ja_vec.size();
    int* ja = new int[nNonzero];
    double* a = new double [nNonzero];
    writeKinPardisoFormat(nNonzero, ja_vec, a_vec, ja, a);
    writeginPardisoFormat(b,nmult);
    int error = solveWithPardiso(a, b, ia, ja, nmult);
    if (error != 0){cerr<<"Pardiso solver did not return success!!"<<endl;}

    if (boundNodesWithSlaveMasterDefinition){
    	equateSlaveDisplacementsToMasters();
    }
    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;
}





// PARDISO prototype. //
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, double *, int    *,    int *, int *,   int *, int *,   int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);


int NewtonRaphsonSolver::solveWithPardiso(double* a, double*b, int* ia, int* ja, const int n_variables){

    // I am copying my libraries to a different location for this to work:
    // On MAC:
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libgfortran.3.dylib /usr/local/lib/
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libgomp.1.dylib /usr/local/lib/
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libquadmath.0.dylib /usr/local/lib/
    // cp libpardiso500-MACOS-X86-64.dylib usr/local/lib
    //
    // compilation:
    // g++ pardiso_sym.cpp -o pardiso_sym  -L./ -L/usr/local/lib -L/usr/lib/  -lpardiso500-MACOS-X86-64 -llapack


    // On ubuntu,
    // cp libpardiso500-GNU461-X86-64.so /usr/lib/
    //
    // sometimes linux cannot recognise liblapack.so.3gf or liblapack.so.3.0.1 or others like this, are essentially liblapack.so
    // on ubuntu you can get this solved by installing liblapack-dev:
    // sudo apt-get install liblapack-dev
    //
    // compilation:
    // gcc test.cpp -o testexe  -L/usr/lib/  -lpardiso500-GNU461-X86-64  -fopenmp  -llapack

    //
    // also for each terminal run:
    // export OMP_NUM_THREADS=1
    // For mkl this is :
    // export MKL_PARDISO_OOC_MAX_CORE_SIZE=10000
    // export MKL_PARDISO_OOC_MAX_SWAP_SIZE=2000
    // fo
    // MSGLVL: the level of verbal output, 0 is no output.

    int    n = n_variables;
    int    nnz = ia[n];
    int    mtype = 11;        /* Real unsymmetric matrix */

    /* RHS and solution vectors. */
    int      nrhs = 1;          /* Number of right hand sides. */
    double   x[n_variables];//, diag[n_variables];
    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */
    void    *pt[64];

    /* Pardiso control parameters. */
    int      iparm[64];
    double   dparm[64];
    int      maxfct, mnum, phase, error, msglvl, solver;

    iparm[60] = 1; //use in-core version when there is enough memory, use out of core version when not.

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    char    *var;
    int      i;// k;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */


/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters.                                */
/* -------------------------------------------------------------------- */

    error = 0;
    solver = 0; /* use sparse direct solver */
    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         return 1;
    }
    else
        //printf("[PARDISO]: License check was successful ... \n");

    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    iparm[2]  = num_procs;

    maxfct = 1;		    /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */

    iparm[10] = 0; /* no scaling  */
    iparm[12] = 0; /* no matching */

    msglvl = 0;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */

/* -------------------------------------------------------------------- */
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
    for (i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] += 1;
    }

/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */
    bool carryOutDebuggingChecks = false;
    if (carryOutDebuggingChecks){
        pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
        if (error != 0) {
            printf("\nERROR in consistency of matrix: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* ..  pardiso_chkvec(...)                                              */
/*     Checks the given vectors for infinite and NaN values             */
/*     Input parameters (see PARDISO user manual for a description):    */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

    if (carryOutDebuggingChecks){
        pardiso_chkvec (&n, &nrhs, b, &error);
        if (error != 0) {
            printf("\nERROR  in right hand side: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* .. pardiso_printstats(...)                                           */
/*    prints information on the matrix to STDOUT.                       */
/*    Use this functionality only for debugging purposes                */
/* -------------------------------------------------------------------- */
    if (carryOutDebuggingChecks){
        pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
        if (error != 0) {
            printf("\nERROR right hand side: %d", error);
            exit(1);
        }
    }
/* -------------------------------------------------------------------- */
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */
    phase = 11;
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
//cout<<"symbolic factorisation"<<endl;
    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    //printf("\nReordering completed ... ");
    //printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
    //printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

/* -------------------------------------------------------------------- */
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */
    phase = 22;
//    iparm[32] = 1; /* compute determinant */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
//cout<<"numerical factorisation"<<endl;
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    //printf("\nFactorization completed ...\n ");

/* -------------------------------------------------------------------- */
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */
   /* phase = 33;

    iparm[7] = 1;       // Max numbers of iterative refinement steps.

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);

    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }
    bool displayResult = false;
    if (displayResult){
        printf("\nSolve completed ... ");
        printf("\nThe solution of the system is: ");
        for (i = 0; i < n; i++) {
            printf("\n x [%d] = % f", i, x[i] );
        }
        printf ("\n");
    }
    //Write x into deltaU:
    for (int i=0; i<n_variables; ++i){
        gsl_vector_set(deltaU,i,x[i]);
    }
    */
/* -------------------------------------------------------------------- */
/* ..  Back substitution with tranposed matrix A^t x=b                  */
/* -------------------------------------------------------------------- */

	phase = 33;
	//iparm[4]  = 61;	 /*changing the precision of convergence with pre-conditioning, not sure what it does, I added as trial, but did not change anything */
	iparm[7]  = 1;       /* Max numbers of iterative refinement steps. */
	iparm[11] = 1;       /* Solving with transpose matrix. */

	pardiso (pt, &maxfct, &mnum, &mtype, &phase,
			 &n, a, ia, ja, &idum, &nrhs,
			 iparm, &msglvl, b, x, &error,  dparm);

	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	bool displayResult = false;
	if (displayResult){
		printf("\nSolve completed ... ");
		printf("\nThe solution of the system is: ");
		for (i = 0; i < n; i++) {
			printf("\n x [%d] = % f", i, x[i] );
		}
		printf ("\n");
	}
    //Write x into deltaU:
    for (int i=0; i<n_variables; ++i){
        gsl_vector_set(deltaU,i,x[i]);
    }

/* -------------------------------------------------------------------- */
/* ..  Convert matrix back to 0-based C-notation.                       */
/* -------------------------------------------------------------------- */
    for (i = 0; i < n+1; i++) {
        ia[i] -= 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] -= 1;
    }

/* -------------------------------------------------------------------- */
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */
    phase = -1;                 /* Release internal memory. */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
    return 0;
}


void NewtonRaphsonSolver::constructiaForPardiso(int* ia, const int nmult, vector<int> &ja_vec, vector<double> &a_vec){
    double negThreshold = -1E-13, posThreshold = 1E-13;
    //count how many elements there are on K matrix and fill up ia:
    int counter = 0;
    for (int i =0; i<nmult; ++i){
        bool wroteiaForThisRow = false;
        for (int j=0; j<nmult; ++j){
            double Kvalue = gsl_matrix_get(K,i,j);
            if (Kvalue>posThreshold || Kvalue<negThreshold){
                ja_vec.push_back(j);
                a_vec.push_back(Kvalue);
                if (!wroteiaForThisRow){
                    //cout<<"writing is for row "<<i<<" column is: "<<j<<endl;
                    ia[i] = counter;
                    wroteiaForThisRow = true;
                }
                counter++;
            }
        }
    }
    ia[nmult] = counter;
}




void NewtonRaphsonSolver::writeKinPardisoFormat(const int nNonzero, vector<int> &ja_vec, vector<double> &a_vec, int* ja, double* a){
    //now filling up the int & double arrays for ja, a
    for (int i=0 ; i<nNonzero; ++i){
        ja[i] = ja_vec[i];
        a[i]  = a_vec [i];
    }
}

void NewtonRaphsonSolver::writeginPardisoFormat(double* b, const int n){
    for (int i=0; i<n; ++i){
        b[i] = gsl_vector_get(gSum,i);
    }
}


bool NewtonRaphsonSolver::checkConvergenceViaDeltaU(){
    bool converged = true;

    double d = gsl_blas_dnrm2 (deltaU);
    //displayMatrix(deltaU,"deltaUInConvergence");
    if (d>threshold){
        converged = false;
        cout<<" not  yet converged via du: norm "<<d<<endl;
    }
    else{
        cout<<"converged with displacement: norm"<<d<<endl;
    }
    return converged;
}

bool NewtonRaphsonSolver::checkConvergenceViaForce(){
    bool converged = true;
    double d = gsl_blas_dnrm2 (gSum);
    if (d>threshold){
        converged = false;
        cout<<" not  yet converged via forces: norm "<<d<<endl;
    }
    else{
        cout<<"converged with forces: norm"<<d<<endl;
    }
    return converged;
}

void NewtonRaphsonSolver::updateUkInIteration(){
    int n = uk->size1;
    for (int i=0; i<n;++i){
    	double newValue = gsl_matrix_get(uk,i,0)+gsl_vector_get(deltaU,i);
        gsl_matrix_set(uk,i,0,newValue);
    }
}

void NewtonRaphsonSolver::calculateDifferenceBetweenNumericalAndAnalyticalJacobian(vector <Node*>& Nodes, bool displayMatricesDuringNumericalCalculation){
	//The normal K still includes the values for fixed nodes. Should correct the fixed nodes,
	//then calculate the difference
	calcutateFixedK(Nodes);
	//calculate the difference:
	gsl_matrix* Kdiff = gsl_matrix_calloc(nDim*nNodes,nDim*nNodes);
	double d = 0;
	for (int i=0; i<nDim*nNodes; ++i){
		for (int j=0; j<nDim*nNodes; ++j){
			double value = gsl_matrix_get(Knumerical,i,j) - gsl_matrix_get(K,i,j);
			gsl_matrix_set(Kdiff,i,j,value);
			d += value*value;
		}
	}
	d = pow(d,0.5);
	if(displayMatricesDuringNumericalCalculation){
		displayMatrix(K,"normalK");
		displayMatrix(Kdiff,"differenceKMatrix");
	}
	cout<<"norm of difference between numerical and analytical K: "<<d<<endl;
	gsl_matrix_free(Kdiff);
}

void NewtonRaphsonSolver::useNumericalJacobianInIteration(){
	gsl_matrix_memcpy(K,Knumerical);
}

void NewtonRaphsonSolver::displayMatrix(gsl_matrix* mat, string matname){
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

void NewtonRaphsonSolver::displayMatrix(gsl_vector* mat, string matname){
    int m = mat->size;
    //int n = mat->size2;
    cout<<matname<<": "<<endl;
    for (int i =0; i<m; i++){
		cout.precision(4);
		cout.width(6);
		cout<<gsl_vector_get(mat,i)<<endl;
    }
}

bool NewtonRaphsonSolver::checkIfCombinationExists(int dofSlave, int dofMaster){
	int n= slaveMasterList.size();
	for(int i=0; i<n;++i){
		if(slaveMasterList[i][0] == dofSlave && slaveMasterList[i][1] == dofMaster){
			return false; //continue addition? false
		}
		if(slaveMasterList[i][1] == dofSlave && slaveMasterList[i][0] == dofMaster){
			return false; //continue addition? false
		}
	}
	return true;
}

void NewtonRaphsonSolver::checkMasterUpdate(int& dofMaster, int& masterId){
	//order is [slave][master]
	int n= slaveMasterList.size();
	for(int i=0; i<n;++i){
		if(slaveMasterList[i][0] == dofMaster){
			dofMaster = slaveMasterList[i][1];
			int dim = dofMaster % 3;
			masterId = (dofMaster-dim)/3;
			break;
		}
	}
}

void NewtonRaphsonSolver::cleanPeripodialBindingFromMaster(int masterDoF, vector<Node*>& Nodes){
	//order is [slave][master]
	//cout<<" searching for master: "<<masterDoF<<endl;
	for (vector< vector<int> >::iterator iter = slaveMasterList.begin(); iter != slaveMasterList.end(); ) {
	    //cout<<"checking slaveMasterList of size "<<slaveMasterList.size()<<" slave/master "<< (*iter)[0]<<" / "<<(*iter)[1]<<endl;
		bool deletedItem = false;
		if ((*iter)[1] == masterDoF){
	    	//I have found a couple where this is a master, is the slave a peripodiaal node?
	    	int dofSlave = (*iter)[0];
	    	int dim = dofSlave % 3;
	    	int slaveId = (dofSlave-dim)/3;
	    	if (Nodes[slaveId]->tissueType == 1){
	    		//clearing the peripodial slave
	    		//cout<<"deleting binding from peripodial node "<<Nodes[slaveId]<<" dim: "<<dim<<endl;
	    		Nodes[slaveId]->slaveTo[dim] = -1;
	    		iter = slaveMasterList.erase(iter);
	    		deletedItem = true;
	    	}
	    }
	    if(!deletedItem){
	        ++iter;
	    }
	}

}

bool NewtonRaphsonSolver::checkIfSlaveIsAlreadyMasterOfOthers(int dofSlave, int dofMaster){
	int n= slaveMasterList.size();
	bool madeChange = false;
	for(int i=0; i<n;++i){
		if(slaveMasterList[i][1] == dofSlave){
			cout<<"making change, slave was master of "<<slaveMasterList[i][0]<<endl;
			//proposed slave is already a master, update the master to the proposed master
			slaveMasterList[i][1] = dofMaster;
			madeChange = true;
		}
	}
	return madeChange;
}
/*
extern int mkl_get_max_threads();

void NewtonRaphsonSolver::solveForDeltaUMKL(){

    const int nmult  = nDim*nNodes;
    MKL_INT *ia = new MKL_INT[nmult+1];
    double *b = new double[nmult];
    vector <int> ja_vec;
    vector <double> a_vec;
    constructiaForPardisoMKL(ia, nmult, ja_vec, a_vec);
    const int nNonzero = ja_vec.size();
    MKL_INT* ja = new MKL_INT[nNonzero];
    double* a = new double [nNonzero];
    writeKinPardisoFormatMKL(nNonzero, ja_vec, a_vec, ja, a);
    writeginPardisoFormat(b,nmult);
    int error = solveWithPardisoMKL(a, b, ia, ja, nmult);
    if (error != 0){cerr<<"Pardiso solver did not return success!!"<<endl;}

    if (boundNodesWithSlaveMasterDefinition){
    	equateSlaveDisplacementsToMasters();
    }
    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;
}

void NewtonRaphsonSolver::setupControlParametersMKL(MKL_INT* iparm, MKL_INT& maxfct, MKL_INT& mnum, MKL_INT& msglvl,  MKL_INT& error){
	for (int i = 0; i < 64; i++) {
	    iparm[i] = 0;
	}
	iparm[0] = 1; // No solver default
	iparm[1] = 2; // Fill-in reordering from METIS
	iparm[2]  = mkl_get_max_threads(); // Numbers of processors, value of MKL_NUM_THREADS
	iparm[3]  =  0;    	// No iterative-direct algorithm
	iparm[4]  =  0;    	// No user fill-in reducing permutation
	iparm[5]  =  0;    	// Write solution into x
	iparm[7]  =  2;    	// Max numbers of iterative refinement steps
 	iparm[9]  =  13;   	// Perturb the pivot elements with 1E-13
    iparm[10] =  0;   	//  no scaling  -  //THis was 1 in the original independent license version!!
    iparm[11] =  0;   	// Use nonsymmetric permutation and scaling MPS
    iparm[12] =  0;   	//  no matching
	iparm[17] = -1;		// Output: Number of nonzeros in the factor LU
	iparm[18] = -1;  	// Output: Mflops for LU factorization
	iparm[19] =  0;   	// Output: Numbers of CG Iterations

	maxfct =1 ; // Maximum number of numerical factorizations.
	mnum= 1; // Which factorization to use.
	msglvl = 0; // Print statistical information in file
   	error = 0;  // Initialize error flag
}

int NewtonRaphsonSolver::solveWithPardisoMKL(double* a, double*b, int* ia, int* ja, const int n_variables){
	//This is the mkl version

    // I am copying my libraries to a different location for this to work:
    // On MAC:
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libgfortran.3.dylib /usr/local/lib/
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libgomp.1.dylib /usr/local/lib/
    // cp /usr/local/lib/gcc/x86_64-apple-darwin14.4.0/4.7.4/libquadmath.0.dylib /usr/local/lib/
    // cp libpardiso500-MACOS-X86-64.dylib usr/local/lib
    //
    // compilation:
    // g++ pardiso_sym.cpp -o pardiso_sym  -L./ -L/usr/local/lib -L/usr/lib/  -lpardiso500-MACOS-X86-64 -llapack


    // On ubuntu,
    // cp libpardiso500-GNU461-X86-64.so /usr/lib/
    //
    // sometimes linux cannot recognise liblapack.so.3gf or liblapack.so.3.0.1 or others like this, are essentially liblapack.so
    // on ubuntu you can get this solved by installing liblapack-dev:
    // sudo apt-get install liblapack-dev
    //
	//MKL path:
	//	export MKLROOT=/opt/intel/compilers_and_libraries/linux/mkl/
	//compilation - using g++ and openMP:
	//	g++  -DMKL_ILP64 -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed  -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl  matrixsolver.cpp -o matrixsolver


    //
    // also for each terminal run:
    // export OMP_NUM_THREADS=1
    // For mkl this is :
    // export MKL_PARDISO_OOC_MAX_CORE_SIZE=10000
    // export MKL_PARDISO_OOC_MAX_SWAP_SIZE=2000
    // MSGLVL: the level of verbal output, 0 is no output.

    MKL_INT    n = n_variables;
    int    nnz = ia[n];
    MKL_INT    mtype = 11;        // Real unsymmetric matrix //

    // RHS and solution vectors.
    MKL_INT      nrhs = 1;          // Number of right hand sides.
    double   x[n_variables];//, diag[n_variables];
    // Internal solver memory pointer pt,                  //
    // 32-bit: int pt[64]; 64-bit: long int pt[64]         //
    // or void *pt[64] should be OK on both architectures  //
    void    *pt[64];

    // Pardiso control parameters. //
    MKL_INT*  iparm = new MKL_INT[64];
    double   dparm[64];
    MKL_INT  maxfct, mnum, phase, error, msglvl, solver;

    iparm[60] = 1; //use in-core version when there is enough memory, use out of core version when not.

    // Number of processors. //
    int      num_procs;

    // Auxiliary variables. //
    char    *var;
    int      i;
    double   ddum;              // Double dummy
    MKL_INT  idum;              // Integer dummy.


// --------------------------------------------------------------------
// ..  Setup Pardiso control parameters.
// --------------------------------------------------------------------
    setupControlParametersMKL(iparm, maxfct, mnum, msglvl, error);d

// --------------------------------------------------------------------
// .. Initialize the internal solver memory pointer. This is only
// necessary for the FIRST call of the PARDISO solver.
// --------------------------------------------------------------------
	for (i = 0; i < 64; i++) {
	    pt[i] = 0;
	}

// -------------------------------------------------------------------- //
// ..  Convert matrix from 0-based C-notation to Fortran 1-based        //
//     notation.                                                        //
// -------------------------------------------------------------------- //
    for (i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] += 1;
    }



// --------------------------------------------------------------------
// .. Reordering and Symbolic Factorization. This step also allocates
// all memory that is necessary for the factorization.
// --------------------------------------------------------------------
	//cout<<"Reordering and symbolic factorisation"<<endl;
	phase = 11;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
    //printf("\nReordering completed ... ");
    //printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
    //printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

// --------------------------------------------------------------------
// .. Numerical factorization
// --------------------------------------------------------------------
//cout<<"Numerical factorisation"<<endl;
	phase = 22;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	//printf("\nFactorization completed ... ");


// --------------------------------------------------------------------
// .. Back substitution and iterative refinement.
// --------------------------------------------------------------------
	phase = 33;
	iparm[7] = 2;
	// Max numbers of iterative refinement steps.
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
 		exit(3);
	}

	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	bool displayResult = false;
	if (displayResult){
		printf("\nSolve completed ... ");
		printf("\nThe solution of the system is: ");
		for (i = 0; i < n; i++) {
			printf("\n x [%d] = % f", i, x[i] );
		}
		printf ("\n");
	}
    //Write x into deltaU:
    for (int i=0; i<n_variables; ++i){
        gsl_vector_set(deltaU,i,x[i]);
    }

// -------------------------------------------------------------------- //
// ..  Convert matrix back to 0-based C-notation.                       //
// -------------------------------------------------------------------- //
    for (i = 0; i < n+1; i++) {
        ia[i] -= 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] -= 1;
    }

// -------------------------------------------------------------------- //
// ..  Termination and release of memory.                               //
// -------------------------------------------------------------------- //
    phase = -1;                 // Release internal memory.

// Release internal memory.
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    return 0;
}

void NewtonRaphsonSolver::constructiaForPardisoMKL(MKL_INT* ia, const int nmult, vector<int> &ja_vec, vector<double> &a_vec){
    double negThreshold = -1E-13, posThreshold = 1E-13;
    //count how many elements there are on K matrix and fill up ia:
    int counter = 0;
    for (int i =0; i<nmult; ++i){
        bool wroteiaForThisRow = false;
        for (int j=0; j<nmult; ++j){
            double Kvalue = gsl_matrix_get(K,i,j);
            if (Kvalue>posThreshold || Kvalue<negThreshold){
                ja_vec.push_back(j);
                a_vec.push_back(Kvalue);
                if (!wroteiaForThisRow){
                    //cout<<"writing is for row "<<i<<" column is: "<<j<<endl;
                    ia[i] = counter;
                    wroteiaForThisRow = true;
                }
                counter++;
            }
        }
    }
    ia[nmult] = counter;
}

void NewtonRaphsonSolver::writeKinPardisoFormatMKL(const int nNonzero, vector<int> &ja_vec, vector<double> &a_vec, MKL_INT* ja, double* a){
    //now filling up the int & double arrays for ja, a
    for (int i=0 ; i<nNonzero; ++i){
        ja[i] = ja_vec[i];
        a[i]  = a_vec [i];
    }
}
*/
