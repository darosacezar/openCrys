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

Class
    Foam::PDFModel

Description
    A Finite mode 3 envinroment PDF model (micromixing)

SourceFiles
    PDFModel.C

\*---------------------------------------------------------------------------*/

#ifndef PDFModel_H
#define PDFModel_H

#include "IOdictionary.H"
#include "primitiveFieldsFwd.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "geometricOneField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class PDFModel Declaration
\*---------------------------------------------------------------------------*/

class PDFModel
:
    public IOdictionary
{
protected:

    // Protected data
	
	
	const volVectorField& U_;
	const volScalarField& rho_;
        const surfaceScalarField& rhoPhi_;
	const volScalarField& k_;
	const volScalarField& epsilon_;
	
	const dictionary& PDFdict_;
	
        dimensionedScalar C_phi_;
	dimensionedScalar Sc_t_;
	dimensionedScalar P_tol_;
	dimensionedScalar Cmu_;
	dimensionedScalar lamDiff_;
	dimensionedScalar S1_;
	dimensionedScalar S2_;
	
	
	volScalarField P1_;
        volScalarField P2_;
	volScalarField P3_;
	volScalarField S3_;
	volScalarField gamma_;
	volScalarField gammaS_;
	volScalarField mixFrac_;
	volScalarField mixVar_;
	
	volScalarField magSqrMixFrac_;
	
	volScalarField effDiff_;
	
	
	
	
	


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcParameters();
	


public:

    TypeName("PDFModel");


    // Constructors

        //- Construct from components
        PDFModel
        (
            const volVectorField& U,
	    const volScalarField& rho,
            const surfaceScalarField& rhoPhi,
	    const volScalarField& k,
	    const volScalarField& epsilon
        );


    //- Destructor
    virtual ~PDFModel()
    {}


    // Member Functions

        
        const volScalarField& P1() const
        {
            return P1_;
        }
        volScalarField& P1()
        {
            return P1_;
        }

        
	const volScalarField& P2() const
        {
            return P2_;
        }
        volScalarField& P2()
        {
            return P2_;
        }
        
        const volScalarField& P3() const
        {
            return P3_;
        }
        volScalarField& P3()
        {
            return P3_;
        }
        const volScalarField& effDiff() const
        {
            return effDiff_;
        }
        volScalarField& effDiff()
        {
            return effDiff_;
        }
        const dimensionedScalar& P_tol() const
        {
            return P_tol_;
        }
        dimensionedScalar& P_tol()
        {
            return P_tol_;
        }
        
         
        virtual tmp<volScalarField> specieSourceRx(dimensionedScalar& SE1, dimensionedScalar& SE2, volScalarField& sSource) const;

        
        virtual tmp<volScalarField> meanSpecie(dimensionedScalar& SE1, dimensionedScalar& SE2, volScalarField& SE3) const;
	
	virtual tmp<volScalarField> specieSourceEnv1(volScalarField& SE1) const;
	
	virtual tmp<volScalarField> specieSourceEnv2(volScalarField& SE2) const;
	
	virtual tmp<volScalarField> specieVolSourceEnv12Rx(volScalarField& SE1, volScalarField& SE2, volScalarField& dHrx) const;
	
        
        virtual void correct();
	virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

//#ifdef NoRepository
//#   include "PDFModel.C"
//#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
