/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2015   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "../sm/EngineeringModels/pod2.h"
#include "../sm/EngineeringModels/nlinearstatic.h"
#include "../sm/Elements/structuralelement.h"
#include "nummet.h"
#include "timestep.h"
#include "metastep.h"
#include "error.h"
#include "verbose.h"
#include "sparsenonlinsystemnm.h"
#include "nrsolver.h"
#include "calmls.h"
#include "outputmanager.h"
#include "datastream.h"
#include "classfactory.h"
#include "timer.h"
#include "contextioerr.h"
#include "sparsemtrx.h"
#include "errorestimator.h"
#include "mathfem.h"
#include "dofmanager.h"
#include "dof.h"
#include "unknownnumberingscheme.h"
#include "util.h"


#include "engngm.h"
#include "oofemtxtdatareader.h"
#include "skylineu.h"
#include "node.h"
#include "sparsemtrx.h"
#include "exportmodulemanager.h"

namespace oofem {
REGISTER_EngngModel(POD2);

POD2 :: POD2(int i, EngngModel *_master) : NonLinearStatic(i, _master)
{

  hyperReductionFlag = true;
  epsilonPOD = 1.e-8;
}


POD2 :: ~POD2()
{
}


int
POD2 :: instanciateYourself(DataReader *dr, InputRecord *ir, const char *dataOutputFileName, const char *desc)
{
    int result;
    result = NonLinearStatic :: instanciateYourself(dr, ir, dataOutputFileName, desc);
    ir->finish();
    // instanciate slave problems
    result &= this->instanciateSlaveProblems();
    return result;
}



IRResultType
POD2 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    NonLinearStatic :: initializeFrom(ir);
    IR_GIVE_OPTIONAL_FIELD(ir, this->hyperReductionFlag, _IFT_POD2_hr);
    IR_GIVE_OPTIONAL_FIELD(ir, this->epsilonPOD, _IFT_POD2_rtolortho);

    IR_GIVE_FIELD(ir, inputStreamNames, _IFT_POD2_slaveprob);

    /*    inputStreamNames.resize(2);
    IR_GIVE_FIELD(ir, inputStreamNames [ 0 ], _IFT_POD_prob1);
    IR_GIVE_FIELD(ir, inputStreamNames [ 1 ], _IFT_POD_prob2);*/
    
    return IRRT_OK;
}
  



  
void 
POD2 :: solveYourself()
{
  // compute POD reduced basis
  this->computeReducedBasis();
  
  StructuralEngngModel :: solveYourself();
}
  
  
void
POD2 :: solveYourselfAt(TimeStep *tStep)
{
    proceedStep(1, tStep);
}
  
  
void
POD2 :: terminate(TimeStep *tStep)
{

    this->doStepOutput(tStep);
    this->printReactionForces(tStep, 1);
    // update load vectors before storing context
    fflush( this->giveOutputStream() );
    this->updateLoadVectors(tStep);
    this->saveStepContext(tStep);

}





void
POD2 :: proceedStep(int di, TimeStep *tStep)
{



 if ( initFlag ) {
        //
        // first step  create space for stiffness Matrix
        //
        if ( !stiffnessMatrix ) {
            stiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        }

        if ( !stiffnessMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        if ( nonlocalStiffnessFlag ) {
            if ( !stiffnessMatrix->isAsymmetric() ) {
                OOFEM_ERROR("stiffnessMatrix does not support asymmetric storage");
            }
        }

        stiffnessMatrix->buildInternalStructure( this, di, EModelDefaultEquationNumbering() );
    }

#if 0
    if ( ( mstep->giveFirstStepNumber() == tStep->giveNumber() ) ) {
 #ifdef VERBOSE
        OOFEM_LOG_INFO("Resetting load level\n");
 #endif
        if ( mstepCumulateLoadLevelFlag ) {
            cumulatedLoadLevel += loadLevel;
        } else {
            cumulatedLoadLevel = 0.0;
        }
        this->loadLevel = 0.0;
    }
#endif

    if ( loadInitFlag || controlMode == nls_directControl ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling reference load\n");
#endif
        //
        // assemble the incremental reference load vector
        //

	FloatArray iLV, iLVP;

        this->assembleIncrementalReferenceLoadVectors(iLV,incrementalLoadVectorOfPrescribed,
                                                      refLoadInputMode, this->giveDomain(di), tStep);

	incrementalLoadVector.beTProductOf(hyperReducedBasisMatrix, iLV);

        loadInitFlag = 0;
    }

    if ( tStep->giveNumber() == 1 ) {
        int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
        totalDisplacement.resize(neq);
        totalDisplacement.zero();
        incrementOfDisplacement.resize(neq);
        incrementOfDisplacement.zero();

	totalReducedCoordinate.resize(hyperReducedBasisMatrix.giveNumberOfColumns());
	totalReducedCoordinate.zero();
	incrementOfReducedCoordinate.resize(hyperReducedBasisMatrix.giveNumberOfColumns());
	incrementOfReducedCoordinate.zero();
    }

    //
    //    ->   BEGINNING OF LOAD (OR DISPLACEMENT) STEP  <-
    //

    //
    // set-up numerical model
    //
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
    //
    // call numerical model to solve arise problem
    //
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "\n\nSolving       [step number %5d.%d, time = %e]\n\n", tStep->giveNumber(), tStep->giveVersion(), tStep->giveIntrinsicTime() );
#endif

    FloatArray extrapolatedForces;
    FloatArray *extrapolatedForcesPtr = &extrapolatedForces;
    if ( this->initialGuessType == IG_Original ) {
      incrementOfDisplacement.zero();
      this->assembleExtrapolatedForces( extrapolatedForces, tStep, ElasticStiffnessMatrix, this->giveDomain(di) );
      extrapolatedForces.negated();
    } else {
      incrementOfDisplacement.zero();
      extrapolatedForcesPtr = NULL;
    }


    
    



    //totalDisplacement.printYourself();
    if ( initialLoadVector.isNotEmpty() ) {
      numMetStatus = nMethod->solveHR(*stiffnessMatrix, incrementalLoadVector, & initialLoadVector, extrapolatedForcesPtr,
                                      totalReducedCoordinate, incrementOfReducedCoordinate, internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep, hyperReducedBasisMatrix);
    } else {
      numMetStatus = nMethod->solveHR(*stiffnessMatrix, incrementalLoadVector, NULL, extrapolatedForcesPtr,
                                      totalReducedCoordinate, incrementOfReducedCoordinate, internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep, hyperReducedBasisMatrix);
    }
    ///@todo Martin: ta bort!!!
    //this->updateComponent(tStep, NonLinearLhs, this->giveDomain(di));
    ///@todo Use temporary variables. updateYourself() should set the final values, while proceedStep should be callable multiple times for each step (if necessary). / Mikael
    OOFEM_LOG_RELEVANT("Equilibrium reached at load level = %f in %d iterations\n", cumulatedLoadLevel + loadLevel, currentIterations);
    prevStepLength =  currentStepLength;
    

}


void
POD2 :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
//
// updates some component, which is used by numerical method
// to newly reached state. used mainly by numerical method
// when new tangent stiffness is needed during finding
// of new equilibrium stage.
//
{
    switch ( cmpn ) {
    case NonLinearLhs:
      {
	if ( stiffMode == nls_tangentStiffness ) {
	  stiffnessMatrix->zero(); // zero stiffness matrix
#ifdef VERBOSE
	  OOFEM_LOG_DEBUG("Assembling tangent stiffness matrix\n");
#endif
	  this->assemble(*stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
			 EModelDefaultEquationNumbering(), d);
	  
	} else if ( ( stiffMode == nls_secantStiffness ) || ( stiffMode == nls_secantInitialStiffness && initFlag ) ) {
#ifdef VERBOSE
	  OOFEM_LOG_DEBUG("Assembling secant stiffness matrix\n");
#endif
	  stiffnessMatrix->zero(); // zero stiffness matrix
	  this->assemble(*stiffnessMatrix, tStep, TangentAssembler(SecantStiffness),
			 EModelDefaultEquationNumbering(), d);
	  initFlag = 0;
	} else if ( ( stiffMode == nls_elasticStiffness ) && ( initFlag ||
							       ( this->giveMetaStep( tStep->giveMetaStepNumber() )->giveFirstStepNumber() == tStep->giveNumber() ) || (updateElasticStiffnessFlag) ) ) {
#ifdef VERBOSE
	  OOFEM_LOG_DEBUG("Assembling elastic stiffness matrix\n");
#endif
	  stiffnessMatrix->zero(); // zero stiffness matrix
	  this->assemble(*stiffnessMatrix, tStep, TangentAssembler(ElasticStiffness),
			 EModelDefaultEquationNumbering(), d);
	  initFlag = 0;
	} else {
	// currently no action , this method is mainly intended to
	// assemble new tangent stiffness after each iteration
	// when secantStiffMode is on, we use the same stiffness
	// during iteration process
	}
	
	FloatMatrix KA, ATKA;
	SkylineUnsym *skyReducedStiffnessMatrix;
	std :: unique_ptr<SkylineUnsym> srStiffnessMatrix;
	std :: unique_ptr<SkylineUnsym> testMtrx;
	std :: unique_ptr< SparseMtrx > reducedStiffnessMatrix;

	srStiffnessMatrix =  std::unique_ptr<SkylineUnsym>(static_cast<SkylineUnsym*> (stiffnessMatrix.release()));

	
	srStiffnessMatrix->times(hyperReducedBasisMatrix, KA);
	for(int i = 1; i <= interfaceDofs.giveSize(); i++) {
	  for(int j = 1; j < hyperReducedBasisMatrix.giveNumberOfColumns(); j++) {
	    KA.at(interfaceDofs.at(i),j) = 0;
	  }
	}
	ATKA.beTProductOf(hyperReducedBasisMatrix,KA);
	//transform to skyline ant to unique_pntr
	srStiffnessMatrix->initializeFromFloatMatrix(ATKA);
	stiffnessMatrix =  std::unique_ptr<SparseMtrx>(static_cast<SparseMtrx*> (srStiffnessMatrix.release()));


	/*
	  stiffnessMatrix->times(hyperReducedBasisMatrix, KA);
	for(int i = 1; i <= interfaceDofs.giveSize(); i++) {
	  for(int j = 1; j < hyperReducedBasisMatrix.giveNumberOfColumns(); j++) {
	    KA.at(interfaceDofs.at(i),j) = 0;
	  }
	}
	ATKA.beTProductOf(hyperReducedBasisMatrix,KA);
      //transform to skyline ant to unique_pntr
	skyReducedStiffnessMatrix->initializeFromFloatMatrix(ATKA);
	skyReducedStiffnessMatrix = new SkylineUnsym(ATKA);
	reducedStiffnessMatrix.reset(skyReducedStiffnessMatrix);
      
      

      *stiffnessMatrix = *reducedStiffnessMatrix;
      //  testMtrx = static_cast<SkylineUnsym*> (stiffnessMatrix);
      testMtrx = std::unique_ptr<SkylineUnsym>(static_cast<SkylineUnsym*>(stiffnessMatrix.release()));
      FloatArray x;
      (*testMtrx).factorized()->backSubstitutionWith(x);
	*/

      //stiffnessMatrix = std::move(reducedStiffnessMatrix);
      // stiffnessMatrix.release();
      //      stiffnessMatrix.reset(skyReducedStiffnessMatrix);
      }
      break;
    case InternalRhs:
      {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Updating internal forces\n");
#endif
	FloatArray iF;
	totalDisplacement.beProductOf(hyperReducedBasisMatrix,totalReducedCoordinate );
	incrementOfDisplacement.beProductOf(hyperReducedBasisMatrix,incrementOfReducedCoordinate );
        // update internalForces and internalForcesEBENorm concurrently
        this->giveInternalForces(iF, true, d->giveNumber(), tStep);
	for(int i = 1; i <= interfaceDofs.giveSize(); i++) {
	  iF.at(interfaceDofs.at(i)) = 0;
	}
	internalForces.beTProductOf(hyperReducedBasisMatrix,iF);
        break;
      }
    default:
        OOFEM_ERROR("Unknown Type of component.");
    }

}




void
POD2 :: computeReducedBasis()
{

  // Solve slave problems and save their solution as snapshots 
  int nSnapshots = 0;
  int smstep = 1, sjstep = 1;
  EngngModel *sp;
  for ( int i = 1; i <= this->giveNumberOfSlaveProblems(); i++ ) {
    sp = this -> giveSlaveProblem(i);
    FILE *out = sp->giveOutputStream();

    //	sp -> solveYourself();
    sp->giveTimer()->startTimer(EngngModelTimer :: EMTT_AnalysisTimer);
    if ( sp->giveCurrentStep() ) {
      smstep = sp->giveCurrentStep()->giveMetaStepNumber();
      sjstep = sp->giveMetaStep(smstep)->giveStepRelativeNumber( sp->giveCurrentStep()->giveNumber() ) + 1;
    } 
    // solve slave problems to get snapshots 
    for ( int imstep = smstep; imstep <= nMetaSteps; imstep++, sjstep = 1 ) { //loop over meta steps
      MetaStep *activeMStep = sp->giveMetaStep(imstep);
      // update state according to new meta step
      sp->initMetaStepAttributes(activeMStep);
      int nTimeSteps = activeMStep->giveNumberOfSteps();
      for ( int jstep = sjstep; jstep <= nTimeSteps; jstep++ ) { //loop over time steps
	nSnapshots++;
	sp->giveTimer()->startTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
	sp->giveTimer()->initTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
	sp->giveNextStep();
	// renumber equations if necessary. Ensure to call forceEquationNumbering() for staggered problems
	if ( sp->requiresEquationRenumbering( sp->giveCurrentStep() ) ) {
	  sp->forceEquationNumbering();
	}
	OOFEM_LOG_DEBUG("Number of equations %d\n", sp->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering()) );

	sp->solveYourselfAt( sp->giveCurrentStep() );
	sp->updateYourself( sp->giveCurrentStep() );

	sp->giveTimer()->stopTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
	sp->terminate( sp->giveCurrentStep() );
	
	double _steptime = sp->giveSolutionStepTime();
	OOFEM_LOG_INFO("EngngModel info: user time consumed by solution step %d: %.2fs\n", sp->giveCurrentStep()->giveNumber(), _steptime);
	fprintf(out, "\nUser time consumed by solution step %d: %.3f [s]\n\n",sp->giveCurrentStep()->giveNumber(), _steptime);
	
	// take the snapshot
	Domain *d = sp->giveDomain(1);
	IntArray dofIDArray;
	FloatArray vec;
	FloatArray FEM_vars(0);

	for (int iNode = 1; iNode <= d->giveNumberOfDofManagers(); iNode ++ ) {	
	  DofManager* dofManager = d->giveDofManager(iNode);
	  dofManager->giveCompleteMasterDofIDArray(dofIDArray);
	  dofManager->giveUnknownVector(vec, dofIDArray, VM_Total, sp->giveCurrentStep(), true);
	  // @todo delate following row
	  //	  vec.at(3) = vec.at(3)*4;
	  FEM_vars.append(vec);
	
	    //	  }		  
	}

	// add snapshot to dof basis
	nDofs = dofIDArray.giveSize();
	this->subspaceExpansion(FEM_vars, 1);

      }
    }
  }
  
  bool ret = this->subspaceSelection();
  //@todo uncomment these two lines
  if(hyperReductionFlag) {
    this->buildReducedDomain(giveDomain(1), sp->giveCurrentStep());
  } else {
    hyperReducedBasisMatrix = reducedBasisMatrix;

  }

}



EngngModel *
POD2 :: giveSlaveProblem(int i)
{
    if ( ( i > 0 ) && ( i <= this->giveNumberOfSlaveProblems() ) ) {
        return this->emodelList[i-1].get();
    } else {
        OOFEM_ERROR("Undefined problem");
    }

    return NULL;
}


int
POD2 :: instanciateSlaveProblems()
{
    //first instantiate master problem if defined
    EngngModel *timeDefProb = NULL;
    emodelList.resize(inputStreamNames.size());
    if ( timeDefinedByProb ) {
        OOFEMTXTDataReader dr( inputStreamNames [ timeDefinedByProb - 1 ] );
        std :: unique_ptr< EngngModel >prob( InstanciateProblem(& dr, this->pMode, this->contextOutputMode, NULL) );
        timeDefProb = prob.get();
        emodelList[timeDefinedByProb-1] = std::move(prob);
    }

    for ( int i = 1; i <= (int)inputStreamNames.size(); i++ ) {
        if ( i == timeDefinedByProb ) {
            continue;
        }

        OOFEMTXTDataReader dr( inputStreamNames [ i - 1 ] );
        //the slave problem dictating time needs to have attribute master=NULL, other problems point to the dictating slave
        std :: unique_ptr< EngngModel >prob( InstanciateProblem(& dr, this->pMode, this->contextOutputMode, NULL) );
        emodelList[i-1] = std::move(prob);
    }

    return 1;
}

 


void 
POD2 :: subspaceExpansion(FloatArray& FEM_vars, double weight) 
{

   FloatArray residuum;
   residuum = FEM_vars;
   int nVectors = reducedBasisMatrix.giveNumberOfColumns();
   //orthogonalize FEM_var against columns of reducedBasisMatrix   
   this->orthogonalize(residuum, reducedBasisMatrix);   
   // if projection error, normalize projection residual and add it to base
   //@todo check the norms and define eps
   if ( norm(residuum) > 1.e-6*norm(FEM_vars) ) {
     residuum.times(1./norm(residuum));
     reducedBasisMatrix.resizeWithData(FEM_vars.giveSize(),nVectors+1);
     reducedBasisMatrix.setColumn(residuum,nVectors+1);
   }

   // add related reduced variables to rbCoords, a = A'*FEM_var
   FloatArray a;
   a.beTProductOf(reducedBasisMatrix, FEM_vars);

   rbCoords.resizeWithData(a.giveSize(),rbCoords.giveNumberOfColumns()+1);
   rbCoords.setColumn(a,rbCoords.giveNumberOfColumns());
     
    //   rbCoords.add(a);
   // update covariance:  Covariance += weight*(a*a')
   if(a.giveSize()>covarianceMatrix.giveNumberOfRows()) {
     covarianceMatrix.resizeWithData(a.giveSize(),a.giveSize());
   }
   FloatMatrix aa;
   aa.beDyadicProductOf(a,a);
   covarianceMatrix.add(aa);
   
}


void 
POD2 :: orthogonalize(FloatArray &answer, const FloatMatrix& A)
{
  //int one = 1;

  for (int k=1;k<=A.giveNumberOfColumns();k++){
    // compute d{k} = -1.*A(:,k)'*R{k}
    double d;
    FloatArray v;
    v.beColumnOf(A,k);
    d = - dot(v, answer);
    // compute R{k+1} = R{k}+d{k}*A(:,k)
    answer.add(d,v);
  }
}


bool 
POD2 ::subspaceSelection()
{
   FloatArray W;
   FloatMatrix U, Ucut;

   // compute [U,W] = eigs(U); 
   // on exit U hold eigenvectors, W eigenvalues sorted in descending order
   int info;
   info = covarianceMatrix.jaco_(W, U, 1.e-12);
   if (info !=0) {
     OOFEM_ERROR("Error in computation of eigenvaluse");
   }
   // find the most relevant eigenvalues
   double treshold = epsilonPOD * W.at(W.giveIndexMaxElem());
   int cutoff;
   for (cutoff = 0; cutoff < W.giveSize(); cutoff++) {
     if (W(cutoff) < treshold) { 
	break;
     }
   }
   // restriction of eigenvectors to selected ones, U = U(:,1:cutoff);
   //@todo check whether beSubMatrixOf is numbered form 1 to size or form 0 to sie-1
   Ucut.beSubMatrixOf(U,1,U.giveNumberOfRows(), 1, cutoff);
    FloatMatrix tmpMatrix;
   // reducedBasisMatrix update : reducedBasisMatrix = reducedBasisMatrix*U;
   tmpMatrix.beProductOf(reducedBasisMatrix,Ucut);
   reducedBasisMatrix = tmpMatrix;
   // Covariance update : Covariance =  U'*Covariance*U;
   tmpMatrix.beProductOf(covarianceMatrix, Ucut);

   //@todo check whether this product should be symmetric
   /*   covarianceMatrix.plusProductSymmUpper(Ucut,tmpMatrix,1);
	covarianceMatrix.symmetrized();*/
   covarianceMatrix.beTProductOf(Ucut, tmpMatrix);
   // rbCoords update : rbCoords = U'*rbCoords;
    tmpMatrix.beTProductOf(U, rbCoords);
   //@todo add or set to tmpVector???
   //rbCoords.add(tmpMatrix);
    rbCoords = tmpMatrix;
  
   return true;
}



void 
POD2 :: buildReducedDomain(Domain *domain, TimeStep* tStep) 
{


  // initialize selected nodes and elements 
  selectedNodes.resize(domain->giveNumberOfDofManagers());
  selectedNodes.zero();
  IntArray nodesMultiplicity(domain->giveNumberOfDofManagers());
  nodesMultiplicity.zero();
  selectedElements.resize(domain->giveNumberOfElements());
  selectedElements.zero();


  IntArray tmp;
  // perform EIM for dof   
  bool ret = this->EIM(tmp);




  if (ret == false) {
    OOFEM_ERROR("solver failed in dof EIM");
  }
  // add Dof_EIM selected nodes to RID
  for (int i = 1; i <= tmp.giveSize(); i++) {
    //@todo this is working if the u is stored first for all nodes then v, and so on
    //    int index = tmp.at(i)%domain->giveNumberOfDofManagers();
    // @todo this is working if all nodes have the same number of dofs
    // @todo should there be +1 or not???
    int index = tmp.at(i)/nDofs;//+1;
    selectedNodes.at(index) = 1;
  }


  // add nodes with nonzero dirichlet boundary condition
  for(int iNode = 1; iNode <= domain->giveNumberOfDofManagers();iNode++) {
    DofManager *dofManager = domain->giveDofManager(iNode);
    FloatArray vec;
    IntArray dofIDArray;
    dofManager->giveCompleteMasterDofIDArray(dofIDArray);
    dofManager->givePrescribedUnknownVector(vec, dofIDArray,VM_Total, tStep);
    if(vec.giveSize()) {
      for(int i = 1; i<= vec.giveSize(); i++) {
	if(vec.at(i) != 0) {
	  selectedNodes.at(iNode) = 1;
	  break;
	}
      }
    }

  }




  //@todo introduce paramter something like
  //    extendReducedDomain(domain, extensionParam);
  extendReducedDomain(domain, 1);

  for(int i = 1; i <= domain->giveNumberOfElements(); i++) {
    if(selectedElements.at(i) == 1)
      domain->giveElement(i)->deActivateYourself();

  }




  IntArray selectedDofs(reducedBasisMatrix.giveNumberOfRows());
  int numberOfReducedDomainDofs = 0;
  int numberOfReducedDomainNodes = 0;
  int numberOfReducedDomainElements = 0;

  int limit = max(selectedNodes.giveSize(), selectedElements.giveSize());



  reducedDomain = domain->Clone();


  

  for( int i = 1; i <= limit; i++) {
    if(i <= selectedNodes.giveSize()) {
      if(selectedNodes.at(i) == 1) {
	numberOfReducedDomainNodes++;
	reducedDomain->setDofManager(numberOfReducedDomainNodes, domain->giveDofManager(i));
	//@todo odkomentit tohle	
	//selectedNodes.at(numberOfReducedDomainNodes) = i;

	//@todo predelat s pouzitim max poctu dofu v domene?
	IntArray dofIDArray;
	domain->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
	for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
	  if(!reducedDomain->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)) {
	    numberOfReducedDomainDofs++;
	    selectedDofs.at(numberOfReducedDomainDofs) = (i-1)*dofIDArray.giveSize() + k;
	  }
	}
      }
    }
    
    if(i <= selectedElements.giveSize()) {
      if(selectedElements.at(i) == 1) {
	numberOfReducedDomainElements++;
	reducedDomain->setElement(numberOfReducedDomainElements, domain->giveElement(i));
      }      
    }
  }
  




  selectedDofs.resizeWithValues(numberOfReducedDomainDofs);
  //@todo odkomentit tohle  
//selectedNodes.resizeWithValues(numberOfReducedDomainNodes);

  // @todo reduced domain???
  reducedDomain->resizeDofManagers(numberOfReducedDomainNodes);
  reducedDomain->resizeElements(numberOfReducedDomainElements);
  
  //  reducedDomain->giveElement(1)->giveLocationArray(testingLocationArray3, EModelDefaultEquationNumbering());




  


  IntArray allColumns(reducedBasisMatrix.giveNumberOfColumns());
  for(int i = 1; i <= reducedBasisMatrix.giveNumberOfColumns(); i++) {
    allColumns.at(i) = i;
  }
  hyperReducedBasisMatrix.resize(selectedDofs.giveSize(),allColumns.giveSize());
  hyperReducedBasisMatrix.beSubMatrixOf(reducedBasisMatrix,selectedDofs,allColumns);
  

  //domainList.resize(2);
  //domainNeqs.resize(2);
  //domainNeqs.at(2) = hyperReducedBasisMatrix.giveNumberOfRows();
  
  //domainPrescribedNeqs.at(2) = 0;

  // this is needed just for post-processing in paraview!!!!!!!!!!!!!!!!!!!!!!
  // change internal component references from labels to assigned local numbers
  /*int num;
  std :: map< int, int >dofManLabelMap, elemLabelMap;
  for ( int i = 1; i <= reducedDomain->giveNumberOfDofManagers(); i++ ) {
    num = reducedDomain->giveDofManager(i)->giveNumber();
    if ( dofManLabelMap.find(num) == dofManLabelMap.end() ) {
      // label does not exist yet
      dofManLabelMap [ num ] = i;
    } else {
      OOFEM_ERROR("iDofmanager entry already exist (label=%d)", num);
    }
  }
  for ( int i = 1; i <= reducedDomain->giveNumberOfElements(); i++ ) {
    num = reducedDomain->giveElement(i)->giveNumber();
    if ( elemLabelMap.find(num) == elemLabelMap.end() ) {
      // label does not exist yet
      elemLabelMap [ num ] = i;
    } else {
      OOFEM_ERROR("Element entry already exist (label=%d)", num);
    }
  }
  reducedDomain->BuildElementPlaceInArrayMap();
  reducedDomain->BuildDofManPlaceInArrayMap();
  MapBasedEntityRenumberingFunctor labelToLocNumFunctor(dofManLabelMap, elemLabelMap);
  for ( int i = 1; i < reducedDomain->giveNumberOfDofManagers(); i++ ) {
    reducedDomain->giveDofManager(i)-> updateLocalNumbering(labelToLocNumFunctor);
  }
  
  for ( int i = 1; i < reducedDomain->giveNumberOfElements(); i++ ) {
    reducedDomain->giveElement(i)->updateLocalNumbering(labelToLocNumFunctor);
  }

  */
  Domain *dm = domain->createReducedDomain(numberOfReducedDomainNodes,selectedNodes,numberOfReducedDomainElements,selectedElements);

  if(hyperReductionFlag) {
    reducedDomain = dm;    
    //reducedDomain = domain;
  }

  this-> setDomain(1, reducedDomain, 0);
  if(hyperReductionFlag) {
    exportModuleManager->reInitialize();
  }
  //this-> setDomain(1, dm, 0);
  

  //?????
  reducedDomain = this->giveDomain(1);
  //this->domainNeqs.resizeWithValues(2);
  //domainPrescribedNeqs.resizeWithValues(2);
  //domainNeqs.at(2) = hyperReducedBasisMatrix.giveNumberOfRows();
  if(!hyperReductionFlag)
    this->forceEquationNumbering(1);

  // create reducedDomain




  

  // testing nodes are : 
  // - selected nodes connected only to selected elements
  // - not Dirichlet nodes


  IntArray nonDomainNodes(domain->giveNumberOfDofManagers());
  // this should be moved to extendRID perharps?
  for (int i = 1; i <= domain->giveNumberOfElements(); i++) {
    if(selectedElements.at(i) == 0) {
      for (int j = 1; j <= domain->giveElement(i)->giveNumberOfDofManagers(); j++) {
	nonDomainNodes.at(domain->giveElement(i)->giveNode(j)->giveNumber()) = 1;
      }
    }    
  }

  interfaceDofs.resize(numberOfReducedDomainDofs);
  interfaceDofs.zero();



  int nIDof = 0;
  int nBC = 0;
  int index = 0;
if(hyperReductionFlag) {
    for(int i = 1; i <= reducedDomain->giveNumberOfDofManagers(); i++) {
      IntArray dofIDArray;
      if(nonDomainNodes.at(reducedDomain->giveDofManager(i)->giveGlobalNumber()) == 1) {
	reducedDomain->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
	for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
	  if(!reducedDomain->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)) {
	    nIDof++;
	    interfaceDofs.at(nIDof) = (i-1)*dofIDArray.giveSize() + k - nBC;	
	  } else {
	    nBC++;
	  }
	}
      } else {
	reducedDomain->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
	for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
	  if(reducedDomain->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)){
	    nBC++;	  
	  }	
	}
      }
    }
 } else {
  for(int i = 1; i <= numberOfReducedDomainNodes; i++) {
    IntArray dofIDArray;
    if(nonDomainNodes.at(reducedDomain->giveDofManager(i)->giveNumber()) == 1) {
      this->giveDomain(1)->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
      for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
	if(!this->giveDomain(1)->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)) {
	  nIDof++;
	  interfaceDofs.at(nIDof) = (i-1)*dofIDArray.giveSize() + k - nBC;	
	} else {
	  nBC++;
	}
      }
    } else {
      this->giveDomain(1)->giveDofManager(i)->giveCompleteMasterDofIDArray(dofIDArray);
      for (int k = 1; k <= dofIDArray.giveSize(); k++ ) {
	if(this->giveDomain(1)->giveDofManager(i)->giveDofWithID(dofIDArray.at(k))->hasBc(tStep)){
	  nBC++;	  
	}	
      }
    }
  }
 }
  interfaceDofs.resizeWithValues(nIDof);

    // create set of interface nodes
  // int nSets = this->giveDomain(1)->giveNumberOfSets();
  //this->giveDomain(1)->resizeSets(nSets+1);
  
  Domain *reducedDomainTest = reducedDomain->Clone();

    
   

}


void 
POD2 :: extendReducedDomain(Domain *domain, int times)
{


  // add connected elements to the selected ones
  for (int count=0; count<times; count++) {
    // select the elements attached to selected nodes
    for (int i = 1; i <= domain->giveNumberOfElements(); i++) {
      for (int j = 1; j<= domain->giveElement(i)->giveNumberOfDofManagers(); j++) {
	//	int rk = mesh.elements[i]->nodes[j]->give_rank();
	//	if (selectedNodes.contains(domain->giveElement(i)->giveNode(j)->giveNumber())) {
	if (selectedNodes.at(domain->giveElement(i)->giveNode(j)->giveNumber())) {
	  selectedElements.at(i) = 1; 
	  break;
	}
      }
    }
    // select the nodes attached to selected elements
    for (int i = 1; i <= domain->giveNumberOfElements(); i++) {
      if(! selectedElements.at(i))
	continue;
      for (int j = 1; j <= domain->giveElement(i)->giveNumberOfDofManagers(); j++) {
	//@todo rozmyslet tohle
	selectedNodes.at(domain->giveElement(i)->giveNode(j)->giveNumber()) = 1;
      }
    }
  }
}


bool 
POD2 :: EIM(IntArray &RID_ddl) 
{
  int rk_max; 
  FloatArray Resid, tmp;
  FloatMatrix V;
  // V = A(:,1)/A(imax,1); RID_ddl = [RID_ddl,imax];  
  Resid.beColumnOf(reducedBasisMatrix,1);
  V.resize(Resid.giveSize(), 1);
  V.setColumn(Resid, 1);
  //  Resid.beColumnOf(V,1);
  //@todo max or max(abs)?
  rk_max = Resid.giveIndexMaxAbsElem(); 
  Resid.times(1/fabs(Resid.at(rk_max)));
  // put index of maximal value to the RID_ddl
  RID_ddl.resizeWithValues(RID_ddl.giveSize()+1);
  RID_ddl.at(RID_ddl.giveSize()) = rk_max;
  
  FloatArray gamma, b;
  FloatMatrix VDDL, tmpMtrx;
  IntArray kArray;
  for (int k=2; k<= reducedBasisMatrix.giveNumberOfColumns(); k++) {
    // gamma = V(RID_ddl,:)\A(RID_ddl,k);
    kArray.resize(1);
    kArray.at(1) = k;
    tmpMtrx.beSubMatrixOf(reducedBasisMatrix, RID_ddl, kArray);
    b.beColumnOf(tmpMtrx, 1);

    IntArray lArray, rArray;
    lArray.resize(V.giveNumberOfColumns());
    rArray.resize(V.giveNumberOfColumns());
    for(int i = 1; i <= V.giveNumberOfColumns(); i ++) {
      lArray.at(i) = RID_ddl.at(i);
      rArray.at(i) = i;
    }
    VDDL.resize(V.giveNumberOfColumns(),V.giveNumberOfColumns());
    VDDL.beSubMatrixOf(V,lArray,rArray);
    // VDDL.beSubMatrixOf(V, RID_ddl.at(1),RID_ddl.at(RID_ddl.giveSize()) , 1, V.giveNumberOfColumns());
    // solve gamma = b/VDDL;
    VDDL.solveForRhs(b, gamma);    
    // R = A(:,k) - V*gamma;
    Resid.beColumnOf(reducedBasisMatrix,k);
    tmp.beProductOf(V,gamma);
    tmp.times(-1.);
    Resid.add(tmp);

    // [val,imax] = max(abs(R)); RID_ddl = [RID_ddl,imax]; V = [V, R/R(imax)];
    rk_max = Resid.giveIndexMaxAbsElem(); 
    Resid.times(1./fabs(Resid.at(rk_max)));    
    // put index of maximal value to the RID_ddl
    RID_ddl.resizeWithValues(RID_ddl.giveSize()+1);
    RID_ddl.at(RID_ddl.giveSize()) = rk_max;
    // @todo maybe doesnt change size, addSubVectorCol may be better???
    V.resizeWithData(V.giveNumberOfRows(), V.giveNumberOfColumns()+1);
    V.setColumn(Resid, V.giveNumberOfColumns());
    
  }
  
  return true;
}



} // end namespace oofem
