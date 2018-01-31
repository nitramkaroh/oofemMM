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

#ifndef hyperreduction_h
#define hyperreduction_h

#include "inputrecord.h"

#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "domain.h"
#include "timestep.h"


#define _IFT_HyperReduction_Name "hyperreduction"
#define _IFT_HyperReduction_ridExtension "ridextension"
#define _IFT_HyperReduction_ridNodeSet "ridnodeset"
#define _IFT_HyperReduction_igaMode "igamode"
#define _IFT_HyperReduction_EIMDofsSelection "eimdofs"
#define _IFT_HyperReduction_EIMStressSelection "eimstress"
#define _IFT_HyperReduction_ridElement2ReduceSet "ridelem2reduceset"




namespace oofem {
/**
 * huhu popis HyperReduction class
**/
class HyperReduction 
{
protected:
  /// isogeometric mode flag
  bool igaModeFlag;
  /// EIM dof selection flag
  bool eimDofsSelectionFlag;
  /// EIM stress selection flag
  bool eimStressSelectionFlag;
 /// Selected nodes for hyper-reduction
  IntArray selectedNodes;
  /// Selected elements for hyper-reduction
  IntArray selectedElements;
  /// Selected iga(integration rule) elements
  std :: map <int, IntArray> *selectedIGAElements;
  /// Selected dofs
  IntArray selectedDofs;
  /// hyper-reduced basis matrix
  FloatMatrix hyperReducedBasisMatrix;
  /// Domain to reduce
  Domain *domain;
  /// hyper-reduced domain
  Domain* reducedDomain;
  /// interface dofs
  IntArray interfaceDofs;
  /// extension of reduced domain
  int ridExtension;
  /// set of nodes which will be included in reduced integration domain
  IntArray ridNodeSetNumbers;
  /// RID set of elements to be reduced
  IntArray ridElement2ReduceSetNumbers;
  /// this is used to evaluate error otherwise unuseful?
  IntArray selectedDofsIDMask;
  /// map from dof number to node number
  std :: map< int, int> *dof2NodeMap;

public:
  HyperReduction(){;}
  HyperReduction(Domain *domain);

  virtual ~HyperReduction();

  virtual IRResultType initializeFrom(InputRecord *ir);  
// identification
  virtual const char *giveClassName() const { return "HyperReduction"; }
  virtual const char *giveInputRecordName() const { return _IFT_HyperReduction_Name; }


   void computeReducedDomainNodes_dofs(TimeStep *tStep, const FloatMatrix &reducedBasisMatrix);
   void computeReducedDomainNodes_stress(TimeStep *tStep, const FloatMatrix &reducedBasisMatrix, int nStressComponents);

   void buildReducedDomain(TimeStep *tStep, const FloatMatrix &reducedBasisMatrix);
   void buildReducedDomain2(TimeStep *tStep, const FloatMatrix &reducedBasisMatrix);


  
   const FloatMatrix &giveHyperReducedBasisMatrix() const {return hyperReducedBasisMatrix;}
   const IntArray &giveInterfaceDofs() const{return interfaceDofs;}
   const IntArray &giveSelectedNodes() const{return selectedNodes;}
   const IntArray &giveSelectedDofsIDMask() const{return selectedDofsIDMask;}
   void initializeYourself(Domain *d);
   Domain *giveReducedDomain(){return reducedDomain;}
protected:

   // empirical interpolation method
   bool EIM(IntArray &answer, const FloatMatrix &rbM);
   void extendReducedDomain(int nTimes, IntArray &nonDomainNodes);
   void extendIGAReducedDomain(int nTimes, IntArray &nonDomainNodes);

};
} // end namespace oofem
#endif // hyperreduction_h
