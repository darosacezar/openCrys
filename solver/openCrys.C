/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    OpenCrys

Description
    Solver for combined antisolvent + cooling crystallization 
    including micromixing and full crystal size distribution 
    Authors:
            Cezar A. da Rosa & Richard D. Braatz
    
    Please refer to the following articles:
    da Rosa, C. A.; Braatz, R. D. Multiscale Modeling and Simulation of Macromixing, 
    Micromixing, and Crystal Size Distribution in Radial Mixers/Crystallizers.
    DOI: 10.1021/acs.iecr.8b00359
    

    Turbulence modelling is based on k-epsilon models. Any k-epsilon model can be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "subCycle.H"
#include "twoPhaseMixtureThermoTransportModel/twoPhaseMixtureThermoTransportModel.H"
#include "PDFModel/PDFModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "itoa.H"
#include "PBESource.H"
#include <stdio.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    

    pimpleControl pimple(mesh);
    
    #include "readGravitationalAcceleration.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "CourantNo.H"
    #include "createPDF.H"
    #include "createPBE.H"
    #include "createTimeControls.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;



        // --- Pressure-velocity PIMPLE corrector loop
        
        while (pimple.loop())
        {            
	    thermoMixture.correct();
            #include "UEqn.H"
	    #include "EEqn.H"
            #include "CEqn.H"  
           

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
            if (pimple.turbCorr())
            {
                turbulence->correct();   
                pdf.correct();                          
                #include "PBEEqn.H"
            }
        }
        runTime.write();
	if(runTime.outputTime())
	{
            
	}

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
