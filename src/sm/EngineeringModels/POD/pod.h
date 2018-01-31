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

#ifndef pod_h
#define pod_h

#include "../sm/EngineeringModels/nlinearstatic.h"
#include "../sm/EngineeringModels/POD/reducedstate.h"
#include "../sm/EngineeringModels/POD/hyperreduction.h"
#include "../sm/EngineeringModels/POD/podvtkxmlexportmodule.h"
#include "../sm/EngineeringModels/POD/reduceddomainnumberingscheme.h"

#include "inputrecord.h"
#include "floatmatrix.h"


#include <string>

#define _IFT_POD_Name "pod"
#define _IFT_POD_hr "hr"
#define _IFT_POD_epsortho "epsortho"
#define _IFT_POD_plotPodBasis "plotpodbasis"
#define _IFT_POD_evaluateError "evalerr"
#define _IFT_POD_numberOfReducedModes "nrmodes"

#define _IFT_POD_performAnalysis "performanalysis"

#define _IFT_POD_performSnapshots "performsnapshots"
#define _IFT_POD_slaveprob "slaveprob"


#define _IFT_POD_saveReducedStateToFile "saverstofile"
#define _IFT_POD_initReducedStateFromFile "initrsfromfile"


#define _IFT_POD_dofsReducedStateInputFileName "dofsrsinfile"
//#define _IFT_POD_dofsInRcFname "dofsinrcfname"
//#define _IFT_POD_dofsInCmFname "dofsincmfname"

#define _IFT_POD_stressReducedStateInputFileName "stressrsinfile"
//#define _IFT_POD_stressInRcFname "stressinrcfname"
//#define _IFT_POD_stressInCmFname "stressincmfname"

#define _IFT_POD_dofsReducedStateOutputFileName "dofsrsoutfile"
//#define _IFT_POD_dofsOutRcFname "dofsoutrcFname"
//#define  _IFT_POD_dofsOutCmFname "dofsoutcmfname"

#define _IFT_POD_stressReducedStateOutputFileName "stressrsoutfile"
//#define _IFT_POD_stressOutRcFname "stressoutrcfname"
//#define _IFT_POD_stressOutCmFname "stressoutcmfname"

#define _IFT_POD_separateBasis "separatebasis"
#define _IFT_POD_nDofGroups "numberdofgroups"
#define _IFT_POD_dofGroups "dofgroups"

#define _IFT_POD_dofWeightsFlag "dofweightsflag"
#define _IFT_POD_dofWeights "dofweights"



#define _IFT_POD_testFlag "testflag"

#define _IFT_POD_computeResultsOutsideRIDFlag "computeresultsoutsiderid"
#define _IFT_POD_outOfRIDElementSet "outofridelementset"

///@name Input fields for NonLinearStatic
//@{
//@}

namespace oofem {

/**
 * huhu popis POD
**/
class POD : public NonLinearStatic
{
protected:

  bool testFlag;


  ///  hyper reduction, 0 - no hyper reduction other hyper reduction is applied
  bool hyperReductionFlag;
  /// print basis flag 
  bool plotPodBasisFlag;
  /// evaluate ||u_fem - u_rom|| error flag
  bool evaluateErrorFlag;
  /// initialization of reduced state form file
  bool initReducedStateFromFileFlag;
  /// save reduced state to file
  bool saveReducedStateToFileFlag;
  /// perform snaphots
  bool performSnapshotsFlag;
  /// perform analysis
  bool performAnalysisFlag;
  /// compute results outside the reduced integration domain
  bool computeResultsOutsideRIDFlag;
  /// compute basis separately for different dofs
  bool separateBasisFlag;
  /// use different weights for different dofs
  bool dofWeightsFlag;

  /// number of reduced modes taken into account
  int nReducedModes;

  /// tolerance in orthogonalization of snapshots
  double epsOrtho;
  /// L^inf norm of residuum
  double rInfNorm;
  /// testing - to be removed
  double rInfNormXi, rInfNormU;
  /// L^2 norm of residuum
  double r2Norm;
  /// constant in error estimator
  double beta;
  /// reduced domainequation numbering
  ReducedDomainNumberingScheme *equationNumbering;
  /// number of stress components
  int nStressComponents;

  /// number of dof groups
  int numberOfDofGroups;
  /// number of dofs in dof groups
  IntArray nDofGroups;
  /// dof groups
  IntArray dofGroups;
  /// array of dof groups
  IntArray* dofGroupArray;
  /// array location arrays for different dof groups
  IntArray* dofGroupIDs;
  /// array containing weight for each dof type
  FloatArray dofWeights;
  /// long array containing weights for all dofs in the mesh
  FloatArray dofWeightsArray;
  /// number of dofs
  int numberOfDomainDofs;
  

  /// vector of slave problems 
  std :: vector< std :: string >inputStreamNames;
  /// List of engineering models to solve sequentially to obtain snapshots
  std :: vector< std :: unique_ptr< EngngModel > >emodelList;
 
  /// input file of reduced state for dofs 
  std :: string dofsReducedStateInputFileName;
  /// input file of reduced state for stress
  std :: string stressReducedStateInputFileName;
  /// output file of reduced state for dofs
  std :: string dofsReducedStateOutputFileName;
  /// output file of reduced state for stress
  std :: string stressReducedStateOutputFileName;

  /// Hyper-reduced reduced basis
  FloatMatrix hrReducedBasisMatrix;
  /// reduced basis for postprocessing
  FloatMatrix rbPostprocess;

  /// reduce state for dofs and stress, can be possibly extended to arbitrary internal variable, i.e. damage ...
  ReducedState *reducedState_dofs;
  ReducedState *reducedState_stress;
  /// hyper reduction
  HyperReduction *hyperReduction;
  /// Export module for printing basis
  PODVTKXMLExportModule *podVTKXMLExportModule;

  /// total reduced coordinates and its increment
  FloatArray totalReducedCoordinate,  incrementOfReducedCoordinate;
  /// dofIDMatrix
  std :: vector<IntArray> dofIDMatrix;
  /// set of elements outside RID  where the results will be computed in postprocessing
  IntArray outOfRIDElementSet;
public:
    POD(int i, EngngModel * _master = NULL);
    virtual ~POD();

    virtual const char *giveClassName() const { return "POD"; }
    virtual const char *giveInputRecordName() const { return _IFT_POD_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);  
    virtual int instanciateYourself(DataReader *dr, InputRecord *ir, const char *outFileName, const char *desc);


    virtual EngngModel *giveSlaveProblem(int i);
    virtual int giveNumberOfSlaveProblems() { return (int)inputStreamNames.size(); }

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void assembleIncrementalReferenceLoadVectors2(FloatArray &_incrementalLoadVector, FloatArray &_incrementalLoadVectorOfPrescribed, SparseNonLinearSystemNM :: referenceLoadInputModeType _refMode, Domain *sourceDomain, TimeStep *tStep);

    virtual void updateYourself(TimeStep *tStep);
    virtual void terminate(TimeStep *tStep);
    double giveReducedBasisValue(TimeStep *tStep, DofManager *dm, int dofID, int iSnapshot);

    virtual void printOutputAt(FILE *file, TimeStep *tStep);



protected:

    void proceedStep(int di, TimeStep *tStep);
    
    void proceedStep2(int di, TimeStep *tStep);
    int forceEquationNumbering2(int id);


    int instanciateSlaveProblems();  
    virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);
    void computeReducedBasis();
    void takeSnapshot_dofs(FloatArray &answer, TimeStep *tStep, Domain *d, std:: vector<IntArray> &dofIDMatrix);
    void takeSnapshot_stress(FloatArray &answer, TimeStep *tStep, Domain *d, int &stressSize);
    void buildReducedDomain();
    void computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di);

    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof);

    void postProcessResults(TimeStep *tStep);
    int forceEquationNumberingPostProcessing(int id);

    FloatMatrix giveReducedBasisMatrix_dofs();
    void assembleDofGroupIDs();
    void assembleDofWeights();

    double computeBeta(TimeStep *tStep);


    void printError(TimeStep *tStep);
    void printPodBasisMatrix();
    void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep);

};
} // end namespace oofem
#endif // pod_h
