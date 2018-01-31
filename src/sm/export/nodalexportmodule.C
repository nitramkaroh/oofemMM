
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

#include "nodalexportmodule.h"
#include "gausspoint.h"
#include "timestep.h"
#include "engngm.h"
#include "node.h"
#include "mathfem.h"
#include "cltypes.h"
#include "material.h"
#include "classfactory.h"
#include "crosssection.h"
#include "set.h"


#include "nodalaveragingrecoverymodel.h"
#include "interfacetype.h"



#include <string>
#include <sstream>
#include <fstream>
#include <ctime>


namespace oofem {
REGISTER_ExportModule(NodalExportModule)

NodalExportModule :: NodalExportModule(int n, EngngModel *e) : ExportModule(n, e), internalVarsToExport(), nodeList(), defaultElementSet( 0, e->giveDomain(1))
{
    smoother = NULL;
    //redToFull = {1, 5, 9, 6, 3, 2}; //position of xx, yy, zz, yz, xz, xy in tensor
}

NodalExportModule :: ~NodalExportModule()
{
    if ( this->smoother ) {
        delete this->smoother;
    }
}


IRResultType
NodalExportModule :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    int val;

    ExportModule :: initializeFrom(ir);

    nodeList.resize(0);
    IR_GIVE_FIELD(ir, nodeList, _IFT_NodalExportModule_nodes); // Macro
    IR_GIVE_FIELD(ir, elementList, _IFT_NodalExportModule_elements);
    IR_GIVE_FIELD(ir, localNumbers, _IFT_NodalExportModule_localNumbers);
    IR_GIVE_FIELD(ir, internalVarsToExport, _IFT_NodalExportModule_internalVarsToExport); // Macro - see internalstatetype.h
 
    val = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_NodalExportModule_stype); // Macro
    stype = ( NodalRecoveryModel :: NodalRecoveryModelType ) val;
       
    return IRRT_OK;
}


void
NodalExportModule :: initialize()
{
    if ( this->smoother ) {
        delete this->smoother;
        this->smoother = NULL;
    }

    defaultElementSet.addAllElements();


}



NodalRecoveryModel *
NodalExportModule :: giveSmoother()
{
    Domain *d = emodel->giveDomain(1);

    if ( this->smoother == NULL ) {
        this->smoother = classFactory.createNodalRecoveryModel(this->stype, d);
    }

    return this->smoother;
}






void
NodalExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    FILE *FID;
    FID = giveOutputStream(tStep);
    //    Domain *domain  = emodel->giveDomain(1);
    //ndim=domain->giveNumberOfSpatialDimensions();
    // Output header
    // fprintf( FID, "%%%% OOFEM generated export file \n");
    //fprintf( FID, "%% Output for time %f\n", tStep->giveTargetTime() );
    

    for ( int iNode = 1; iNode <= nodeList.giveSize(); iNode++ ) {
      doOutputNode(tStep, FID, nodeList.at(iNode), elementList.at(iNode), localNumbers.at(iNode));
    }
    fclose(FID);

}
  

void
NodalExportModule :: doOutputNode(TimeStep *tStep, FILE *FID, int nodeNumber, int elementNumber, int localNumber)
{

    InternalStateType type;
    FloatArray val;
    //   fprintf(FID, "Node %d\n", nodeNumber);
    for ( int field = 1; field <= internalVarsToExport.giveSize(); field++ ) {
        type = ( InternalStateType ) internalVarsToExport.at(field);
	this->getNodalVariableFromIS(val, nodeNumber, elementNumber, localNumber,tStep, type);
	const char *name = __InternalStateTypeToString(type);
	//	fprintf(FID, "%s [", name);
	  for ( int component = 1 ; component <= val.giveSize(); component++ ) {
            fprintf( FID, "%e ", val.at(component) );
	  }

	  //fprintf(FID, "];\n");
	  fprintf(FID,"\n");
    }


    
    


}


FILE *
NodalExportModule :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;
    std :: string fileName = giveOutputFileName(tStep);
    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str());
    }

    return answer;
}

std :: string
NodalExportModule :: giveOutputFileName(TimeStep *tStep)
{
    return this->giveOutputBaseFileName(tStep) + ".node";
}




void
NodalExportModule :: getNodalVariableFromIS(FloatArray &answer, int nodeNumber, int elementNumber, int localNumber,TimeStep *tStep, InternalStateType type)
{
    // Recovers nodal values from Internal States defined in the integration points.
    // Should return an array with proper size supported by VTK (1, 3 or 9)
    // Domain *d = emodel->giveDomain(1
    this->giveSmoother();
    IntArray redIndx;

    Domain *d = emodel->giveDomain(1);
    Element *e = d->giveElement(elementNumber);
    Interface *interface;
    interface  = e->giveInterface(NodalAveragingRecoveryModelInterfaceType);
    if(interface) {      
      NodalAveragingRecoveryModelInterface *i = static_cast<NodalAveragingRecoveryModelInterface*>(interface);
      i->NodalAveragingRecoveryMI_computeNodalValue(answer, localNumber,  type, tStep);
    } else {
      this->smoother->recoverValues(* this->giveRegionSet(1), type, tStep);
      const FloatArray *val = NULL;
      FloatArray valueArray;
      //    InternalStateValueType valType = giveInternalStateValueType(type);
      
      int found = this->smoother->giveNodalVector( val, nodeNumber );
      if ( !found ) {
        valueArray.resize( redIndx.giveSize() );
        val = & valueArray;
      }    
      
      int isize = val->giveSize();
      answer.resize(isize);
      for ( int i = 1; i <= isize; i++ ) {
	answer.at(i) = val->at(i);
      }
    }
    
}



Set *NodalExportModule :: giveRegionSet(int i)
{
  /*    int setid = regionSets.at(i);
    if ( setid > 0 ) {
        return emodel->giveDomain(1)->giveSet(setid);
	} else {*/
        return & this->defaultElementSet;
	//}
}




} // end namespace oofem










