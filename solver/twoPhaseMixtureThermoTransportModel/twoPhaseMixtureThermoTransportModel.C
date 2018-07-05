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
    
Description
    This class implements all the necessary functions to deal with the macromixing between the solution and antisolvent

\*---------------------------------------------------------------------------*/

#include "twoPhaseMixtureThermoTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseMixtureThermoTransportModel, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void Foam::twoPhaseMixtureThermoTransportModel::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}
void Foam::twoPhaseMixtureThermoTransportModel::calcRho1()
{
    // Average kinematic viscosity calculated from dynamic viscosity
    const dimensionedScalar c1("c1", dimensionSet (1,-3,0,0,0,0,0), 10);
    rho1_ = rho1_0_;//this function can be used to input different functions for the density of the solution
}
void Foam::twoPhaseMixtureThermoTransportModel::calcRho2()
{
    // Average kinematic viscosity calculated from dynamic viscosity
    const dimensionedScalar c2("c2", dimensionSet (1,-3,0,0,0,0,0), 10);
    rho2_ = rho2_0_;//this function can be used to input different functions for the density of the antisolvent
}
void Foam::twoPhaseMixtureThermoTransportModel::calcRho()
{
    // Calculates mixture density
    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );
    rho_ = limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_;
    //Info << rho1_ << "\t" << rho2_ << "\n" << endl;
}
void Foam::twoPhaseMixtureThermoTransportModel::calcCp()
{
    // Calculates mixture specific heat
    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );
    Cp_ = (rho1_*limitedAlpha1*Cp1_+rho2_*(scalar(1)-limitedAlpha1)*Cp2_)/rho_;
}
void Foam::twoPhaseMixtureThermoTransportModel::calcAlphaEff()
{
    // Calculates mixture alphaEff
    const volScalarField& nut_ = U_.db().objectRegistry::lookupObject<volScalarField>("nut");
    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );
    alphaEff_ = (rho1_*limitedAlpha1*k1_+rho2_*(scalar(1)-limitedAlpha1)*k2_)/rho_/Cp_+rho_*nut_;
    //Info << rho1_ << "\t" << rho2_ << "\n" << endl;
}
void Foam::twoPhaseMixtureThermoTransportModel::calcT()
{
    // Calculates the Temperature based on the enthalpy value considering Cp=cte
    T_ = he_/Cp_+Tref_;
}
void Foam::twoPhaseMixtureThermoTransportModel::calcK()
{
    // Calculates the knitic energy
    K_ = scalar(0.5)*magSqr(U_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseMixtureThermoTransportModel::twoPhaseMixtureThermoTransportModel
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    
    U_(U),
    phi_(phi),
    
    Tref_("Tref", dimTemperature, nuModel1_->viscosityProperties().lookup("Tref")),
    Cp1_("Cp1", dimensionSet (0, 2, -2, -1, 0, 0, 0), nuModel1_->viscosityProperties().lookup("Cp")),
    Cp2_("Cp2", dimensionSet (0, 2, -2, -1, 0, 0, 0), nuModel2_->viscosityProperties().lookup("Cp")),	
    k1_("k1", dimensionSet (1, 1, -3, -1, 0, 0, 0), nuModel1_->viscosityProperties().lookup("k")),
    k2_("k2", dimensionSet (1, 1, -3, -1, 0, 0, 0), nuModel2_->viscosityProperties().lookup("k")),
    rho1_0_("rho", dimDensity, nuModel1_->viscosityProperties().lookup("rho")),
    rho2_0_("rho", dimDensity, nuModel2_->viscosityProperties().lookup("rho")),

    
    
    rho1_
    (
        IOobject
        (
            "rho1",
            U.time().timeName(),
            U.db()
        ),
        U_.mesh(),
        dimensionedScalar("rho1",rho1_0_),
        calculatedFvPatchScalarField::typeName
    ),
    rho2_
    (
        IOobject
        (
            "rho2",
            U.time().timeName(),
            U.db()
        ),
        U_.mesh(),
        dimensionedScalar("rho2",rho2_0_),
        calculatedFvPatchScalarField::typeName
    ),
    rho_
    (
        IOobject
        (
	    "rho",
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1_*rho1_+alpha2_*rho2_
    ),
    Cp_
    (
        IOobject
        (
	    "Cp",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_*Cp1_+alpha2_*Cp2_
    ),
    he_
    (
        IOobject
        (
            "he",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    
    T_
    (
        IOobject
        (
	    "T",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        he_/Cp_+Tref_
    ),
    alphaEff_
    (
        IOobject
        (
	    "alphaEff",
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        (alpha1_*rho1_*k1_ + alpha2_*rho2_*k2_)/rho_/Cp_
    ),
    K_
    (
        IOobject
        (
	    "K",
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        0.5*magSqr(U_)
    ) ,

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    calcT();  
    calcRho1();
    calcRho2();
    calcRho();
    calcNu();
    calcCp();
    //calcAlphaEff();
    calcK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::twoPhaseMixtureThermoTransportModel::mu() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            limitedAlpha1*rho1_*nuModel1_->nu()
          + (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseMixtureThermoTransportModel::muf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "muf",
            alpha1f*fvc::interpolate(nuModel1_->nu()*rho1_)
          + (scalar(1) - alpha1f)*fvc::interpolate(nuModel2_->nu()*rho2_)
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseMixtureThermoTransportModel::nuf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nuf",
            (
                alpha1f*fvc::interpolate(nuModel1_->nu()*rho1_)
              + (scalar(1) - alpha1f)*fvc::interpolate(nuModel2_->nu()*rho2_)
            )/(alpha1f*fvc::interpolate(rho1_) + (scalar(1) - alpha1f)*fvc::interpolate(rho2_))
        )
    );
}


bool Foam::twoPhaseMixtureThermoTransportModel::read()
{
    if (regIOobject::read())
    {
        if
        (
            nuModel1_().read
            (
                subDict(phase1Name_ == "1" ? "phase1": phase1Name_)
            )
         && nuModel2_().read
            (
                subDict(phase2Name_ == "2" ? "phase2": phase2Name_)
            )
        )
        {
           // nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
           // nuModel2_->viscosityProperties().lookup("rho") >> rho2_;
	  nuModel1_->viscosityProperties().lookup("Tref") >> Tref_;
	  nuModel1_->viscosityProperties().lookup("Cp") >> Cp1_;
	  nuModel2_->viscosityProperties().lookup("Cp") >> Cp2_;	
	  nuModel1_->viscosityProperties().lookup("k") >> k1_;
	  nuModel2_->viscosityProperties().lookup("k") >> k2_;
	  nuModel1_->viscosityProperties().lookup("rho") >> rho1_0_;
	  nuModel2_->viscosityProperties().lookup("rho") >> rho2_0_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
