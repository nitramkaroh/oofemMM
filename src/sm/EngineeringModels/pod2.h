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

#ifndef pod2_h
#define pod2_h

#include "../sm/EngineeringModels/nlinearstatic.h"
#include "inputrecord.h"
#include "floatmatrix.h"

#define _IFT_POD2_Name "pod2"
#define _IFT_POD2_hr "hr"
#define _IFT_POD2_rtolortho "rtolortho"

#define _IFT_POD2_slaveprob "slaveprob"

/*#define _IFT_POD_prob1 "prob1"
#define _IFT_POD_prob2 "prob2"
*/



///@name Input fields for NonLinearStatic
//@{
//@}

namespace oofem {
  class FloatMatrix;
/**
 * huhu popis POD
**/
class POD2 : public NonLinearStatic
{
protected:
  ///  hyper reduction, 0 - no hyper reduction other hyper reduction is applied
  bool hyperReductionFlag;
  /// Optional parameter which specify problems to define load time functions
  int timeDefinedByProb;
  /// relative error for Gramm-Schmidt orthonormalization
  double epsilonPOD;
  /// reduced basis coordinates
  FloatMatrix rbCoords;
  /// matrix containing reducedBasis
  FloatMatrix reducedBasisMatrix;
  /// covariance matrix
  FloatMatrix covarianceMatrix;
  /// vector of slave problems 
  std :: vector< std :: string >inputStreamNames;
  /// List of engineering models to solve sequentially to obtain snapshots
  std :: vector< std :: unique_ptr< EngngModel > >emodelList;

  /// Selected nodes for hyper-reduction
  IntArray selectedNodes;
  /// Selected elements for hyper-reduction
  IntArray selectedElements;
  /// hyper-reduced basis matrix
  FloatMatrix hyperReducedBasisMatrix;
  /// hyper-reduced domain
  Domain* reducedDomain;
  /// interface dofs
  IntArray interfaceDofs;
  /// total reduced coordinates and its increment
  FloatArray totalReducedCoordinate,  incrementOfReducedCoordinate;



  /// @todo remove this
  int nDofs;
public:
    POD2(int i, EngngModel * _master = NULL);
    virtual ~POD2();

    virtual IRResultType initializeFrom(InputRecord *ir);  

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);

    virtual void terminate(TimeStep *tStep);




    // identification
  virtual const char *giveClassName() const { return "POD2"; }
  virtual const char *giveInputRecordName() const { return _IFT_POD2_Name; }
  
  virtual int instanciateYourself(DataReader *dr, InputRecord *ir, const char *outFileName, const char *desc);

   virtual EngngModel *giveSlaveProblem(int i);
   virtual int giveNumberOfSlaveProblems() { return (int)inputStreamNames.size(); }

protected:
    void proceedStep(int di, TimeStep *tStep);
    int instanciateSlaveProblems();

    ///compute reduced basis functions
    void computeReducedBasis();  
    void subspaceExpansion(FloatArray& FEM_vars, double weight) ;
    void orthogonalize(FloatArray &answer, const FloatMatrix& A);
    bool subspaceSelection();

    //hyper-reduction methods
    void buildReducedDomain(Domain *domain, TimeStep *tStep);
    void extendReducedDomain(Domain *domain, int times);

    // empirical interpolation method
    bool EIM(IntArray &);

    virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);


};
} // end namespace oofem
#endif // pod2_h
