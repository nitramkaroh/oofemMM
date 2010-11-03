/* New class created by Milan Jirasek on 1 Feb 2010 */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2010   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

//
// class DofManExportModule
//

#ifndef dmexportmodule_h
#define dmexportmodule_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "exportmodule.h"
#include "domain.h"
#include "engngm.h"

namespace oofem{

/**
 Represents DofManager export module. 
 This module writes the coordinates of all dof managers
 along with the values of displacements
 for further processing.
*/
class DofManExportModule : public ExportModule
{
protected:

public:
 /// Constructor
 DofManExportModule (EngngModel* e);

 /// Destructor
 ~DofManExportModule();

 /// Initializes receiver acording to object description stored in input record
 virtual IRResultType initializeFrom (InputRecord* ir);
 /**
  Writes the output.
  @param tStep time step.
  */
 void doOutput (TimeStep* tStep);
 /// Returns class name of the receiver.
 virtual const char* giveClassName () const { return "DofManExportModuleClass" ;}

protected:
 /// returns the output stream for given solution step
 FILE* giveOutputStream (TimeStep*) ;
 
};
} // end namespace
#endif


