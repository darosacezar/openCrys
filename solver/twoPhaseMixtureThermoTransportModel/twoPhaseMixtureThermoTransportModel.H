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
    Foam::twoPhaseMixtureThermoTransportModel

Description
    A two-phase incompressible transportModel

SourceFiles
    twoPhaseMixtureThermoTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseMixtureThermoTransportModel_H
#define twoPhaseMixtureThermoTransportModel_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixtureThermoTransportModel Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseMixtureThermoTransportModel
:
    public IOdictionary,
    public transportModel,
    public twoPhaseMixture
{
protected:

    // Protected data

        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;
	
	const volVectorField& U_;
        const surfaceScalarField& phi_;
	
	dimensionedScalar Tref_;
	dimensionedScalar Cp1_;
	dimensionedScalar Cp2_;
	dimensionedScalar k1_;
	dimensionedScalar k2_;
	dimensionedScalar rho1_0_;
	dimensionedScalar rho2_0_;

        volScalarField rho1_;
        volScalarField rho2_;
	volScalarField rho_;
	volScalarField Cp_;
	volScalarField he_;
	volScalarField T_;
	volScalarField alphaEff_;
	volScalarField K_;
	volScalarField nu_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcT();  
	void calcRho1();
	void calcRho2();
	void calcRho();
	void calcNu();
	void calcCp();
	void calcAlphaEff();
	void calcK();
	


public:

    TypeName("twoPhaseMixtureThermoTransportModel");


    // Constructors

        //- Construct from components
        twoPhaseMixtureThermoTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~twoPhaseMixtureThermoTransportModel()
    {}


    // Member Functions

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase1 density
        const volScalarField& rho1() const
        {
            return rho1_;
        }
        volScalarField& rho1()
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
	const volScalarField& rho2() const
        {
            return rho2_;
        }
        volScalarField& rho2()
        {
            return rho2_;
        }
        
        const volScalarField& rho() const
        {
            return rho_;
        }
        volScalarField& rho()
        {
            return rho_;
        }
        const volScalarField& Cp() const
        {
            return Cp_;
        }
        volScalarField& Cp()
        {
            return Cp_;
        }
        const volScalarField& he() const
        {
            return he_;
        }
        volScalarField& he()
        {
            return he_;
        }
        const volScalarField& T() const
        {
            return T_;
        }
        volScalarField& T()
        {
            return T_;
        }
        const volScalarField& alphaEff() const
        {
            return alphaEff_;
        }
        volScalarField& alphaEff()
        {
            return alphaEff_;
        }
        const volScalarField& K() const
        {
            return K_;
        }
        volScalarField& K()
        {
            return K_;
        }

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcT();  
	    calcRho1();
	    calcRho2();
	    calcRho();
	    calcNu();
	    calcCp();
	    calcAlphaEff();
	    calcK();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
