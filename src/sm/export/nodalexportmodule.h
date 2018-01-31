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

#ifndef nodalexportmodule_h_
#define nodalexportmodule_h_

#include <vector>

#include "exportmodule.h"
#include "internalstatevaluetype.h"
#include "nodalrecoverymodel.h"

///@name Input fields for NodalExportModule
//@{
#define _IFT_NodalExportModule_Name "nodalexportmodule"
/// list of internal vars to export
#define _IFT_NodalExportModule_internalVarsToExport "vars"
/// nodal recovery model (smoother) type
#define _IFT_NodalExportModule_stype "stype"
/// nodes to export
#define _IFT_NodalExportModule_nodes "nodes"
#define _IFT_NodalExportModule_elements "elements"
#define _IFT_NodalExportModule_localNumbers "localnumbers"

//@}

namespace oofem {
/**
 * (Under development) The Nodal export module enables oofem to export the results to a textfile containing the nodal values of variables extrapolated from 
 * gauss points
 * 
 * @author Martin Horak
 */
class OOFEM_EXPORT NodalExportModule : public ExportModule
{
protected:
    /// list of InternalStateType values, identifying the selected vars for export
    IntArray internalVarsToExport;
    IntArray nodeList;
    IntArray elementList;
    IntArray localNumbers;
    /// Smoother type.
    NodalRecoveryModel :: NodalRecoveryModelType stype;    
    NodalRecoveryModel *smoother;
    /// Default region set
    Set defaultElementSet;

public:
     NodalExportModule(int n, EngngModel * e);
     virtual ~NodalExportModule();
     virtual IRResultType initializeFrom(InputRecord *ir); 
     virtual void initialize();

     // virtual void terminate();
     virtual void doOutput(TimeStep *tStep, bool forcedOutput = false);
     virtual void doOutputNode(TimeStep *tStep, FILE *FID, int iNode, int iElement, int localNumber);
     virtual void getNodalVariableFromIS(FloatArray &answer, int nodeNumber, int elementNumber,int localNumber,TimeStep *tStep, InternalStateType type);

     virtual Set *giveRegionSet(int i);

    
     virtual NodalRecoveryModel *giveSmoother();

     FILE *giveOutputStream(TimeStep *tStep);
     std::string giveOutputFileName(TimeStep *tStep);


     virtual const char *giveClassName() const { return "NodalExportModule"; }
     virtual const char *giveInputRecordName() const { return _IFT_NodalExportModule_Name; }

};
} // end namespace oofem
#endif // nodalexportmodule_h_
