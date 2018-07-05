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
    Creates the PDFModel class and all the necessary member functions
    Solve the PDFModel equations

\*---------------------------------------------------------------------------*/

#include "PDFModel.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "fvc.H"
#include "fvm.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PDFModel, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


void Foam::PDFModel::calcParameters()
{
  effDiff_ = rho_*(Cmu_.value()*k_*k_/Sc_t_/epsilon_+lamDiff_);
  scalar term1, term2;
  
  //const fvMesh& mesh = P1_.mesh();
  forAll(P1_.primitiveFieldRef(),celli)
  {
    //test P3 value to avoid undterminate values
    if(P3_[celli] >= P_tol_.value())
    {
      mixFrac_[celli] = std::min(std::max(S3_[celli]/P3_[celli],0.0),1.0);
    }
    else
    {
      mixFrac_[celli] = P1_[celli]*S1_.value()+P2_[celli]*S2_.value();
    }
  }
  forAll(P1_.boundaryFieldRef(), patchi)
  {
    if(P3_[patchi] >= P_tol_.value())
    {
      mixFrac_[patchi] = std::min(std::max(S3_[patchi]/P3_[patchi],0.0),1.0);
    }
    else
    {
      mixFrac_[patchi] = P1_[patchi]*S1_.value()+P2_[patchi]*S2_.value();
    }    
  }
  magSqrMixFrac_ = magSqr(fvc::grad(mixFrac_));
  forAll(P1_.primitiveFieldRef(),celli)
  {
  //calculate mixture variance in the 3rd environment
    mixVar_[celli] = P1_[celli]*(1.0-P1_[celli])-2.0*P1_[celli]*P3_[celli]*mixFrac_[celli]+P3_[celli]*(1.0-P3_[celli])*mixFrac_[celli]*mixFrac_[celli];
  
  //Conditionals to avoid undeterminate values (calculated using L'Hopital's rule)
    if((mixFrac_[celli] <= P_tol_.value()) || (mixFrac_[celli] >= (1.0-P_tol_.value())))
    {
      gamma_[celli] = C_phi_.value()*epsilon_[celli]/k_[celli];
    }
    else
    {
      if(P1_[celli] >= (1.0-P_tol_.value()))
      {
	gamma_[celli] = (C_phi_.value()*epsilon_[celli]/k_[celli])/((1.0-mixFrac_[celli])*(1.0-mixFrac_[celli]));
      }
      else if(P2_[celli] >= (1.0-P_tol_.value()))
      {
	gamma_[celli] = 0.0;
      }
      else if(P3_[celli] >= (1.0-P_tol_.value()))
      {
	gamma_[celli] = C_phi_.value()*epsilon_[celli]/k_[celli];
      }
      else
      {
	gamma_[celli] = (C_phi_.value()*epsilon_[celli]/k_[celli])*mixVar_[celli]/(P1_[celli]*(1.0-P1_[celli])*(1.0-mixFrac_[celli])*(1.0-mixFrac_[celli])+P2_[celli]*(1.0-P2_[celli])*mixFrac_[celli]*mixFrac_[celli]);
      }
    }   
    if(P3_[celli] > P_tol_.value())
    {
      gammaS_[celli] = 2.0*effDiff_[celli]*magSqrMixFrac_[celli]/((1.0-mixFrac_[celli])*(1.0-mixFrac_[celli])+mixFrac_[celli]*mixFrac_[celli]);
      //conditionals to improve stability of the numerical solution
      term1 = gamma_[celli]*P1_[celli]*(1.0-P1_[celli])-gammaS_[celli]*P3_[celli];
      term2 = gamma_[celli]*P2_[celli]*(1.0-P2_[celli])-gammaS_[celli]*P3_[celli];
      if((term1 < 0.0) || (term2 < 0.0))
      {
	gammaS_[celli]=gamma_[celli]*P1_[celli]*(1.0-P1_[celli])/P3_[celli];
	if(term2 < term1)
	{
	  gammaS_[celli]=gamma_[celli]*P2_[celli]*(1.0-P2_[celli])/P3_[celli];
	}
      }
      if(gammaS_[celli] < 0.0) {gammaS_[celli] = 0.0;}
    }
    else {gammaS_[celli] = 0.0;}  
  }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFModel::PDFModel
(
    const volVectorField& U,
    const volScalarField& rho,
    const surfaceScalarField& rhoPhi,
    const volScalarField& k,
    const volScalarField& epsilon
)
:
    IOdictionary
    (
        IOobject
        (
            "PDFdict",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    U_(U),
    rho_(rho),
    rhoPhi_(rhoPhi),
    k_(k),
    epsilon_(epsilon),

    PDFdict_(*this),

    C_phi_("C_phi",dimensionSet (0, 0, 0, 0, 0, 0, 0), PDFdict_.lookup("C_phi")),
    Sc_t_("Sc_t", dimensionSet (0, 0, 0, 0, 0, 0, 0), PDFdict_.lookup("Sc_t")),
    P_tol_("P_tol", dimensionSet (0, 0, 0, 0, 0, 0, 0), PDFdict_.lookup("P_tol")),
    Cmu_("Cmu", dimensionSet (0, 0, 0, 0, 0, 0, 0), PDFdict_.lookup("Cmu")),
    lamDiff_("lamDiff", dimensionSet (0, 2, -1, 0, 0, 0, 0), PDFdict_.lookup("lamDiff")),
    S1_("S1", dimensionSet (0, 0, 0, 0, 0, 0, 0), scalarList(PDFdict_.lookup("phi1"))[0]),
    S2_("S2", dimensionSet (0, 0, 0, 0, 0, 0, 0), scalarList(PDFdict_.lookup("phi2"))[0]),    
     
   
    
    P1_
    (
        IOobject
        (
            "P1",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    
    P2_
    (
        IOobject
        (
            "P2",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    
    P3_
    (
        IOobject
        (
	    "P3",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        scalar(1.0)-P1_-P2_
    ),
    S3_
    (
        IOobject
        (
            "S3",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    gamma_
    (
        IOobject
        (
            "gamma",
            U_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("gamma",dimensionSet(0, 0, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),
    gammaS_
    (
        IOobject
        (
            "gammaS",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("gammaS",dimensionSet(0, 0, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),
    mixFrac_
    (
        IOobject
        (
	    "mixFrac",
            U_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        S3_
    ),
    mixVar_
    (
        IOobject
        (
	    "mixVar",
            U_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mixFrac_
    ),
    magSqrMixFrac_
    (
        IOobject
        (
	    "magSqrMixFrac",
            U_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        magSqr(fvc::grad(mixFrac_))
    ),
    effDiff_
    (
        IOobject
        (
            "effDiff",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("effDiff",dimensionSet(1, -1, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ) 
{
    //update function values
  calcParameters(); 
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::PDFModel::specieSourceRx(dimensionedScalar& SE1, dimensionedScalar& SE2, volScalarField& sSource) const
{
    return
    (
       rho_*(P3_*sSource-gammaS_*P3_*(SE1+SE2)+gamma_*(P1_*(scalar(1.0)-P1_)*SE1+P2_*(scalar(1.0)-P2_)*SE2))
    );
}
Foam::tmp<Foam::volScalarField>
Foam::PDFModel::specieSourceEnv1(volScalarField& SE1) const
{
    return
    (
       rho_*(-gammaS_*P3_*SE1+gamma_*P1_*(scalar(1.0)-P1_)*SE1)
    );
}
Foam::tmp<Foam::volScalarField>
Foam::PDFModel::specieSourceEnv2(volScalarField& SE2) const
{
    return
    (
       rho_*(-gammaS_*P3_*SE2+gamma_*P2_*(scalar(1.0)-P2_)*SE2)
    );
}

Foam::tmp<Foam::volScalarField>
Foam::PDFModel::specieVolSourceEnv12Rx(volScalarField& SE1, volScalarField& SE2, volScalarField& dHrx) const
{
    return
    (
       (-gammaS_*P3_*(SE1+SE2)+gamma_*(P1_*(scalar(1.0)-P1_)*SE1+P2_*(scalar(1.0)-P2_)*SE2))*dHrx
    );
}
Foam::tmp<Foam::volScalarField>
Foam::PDFModel::meanSpecie(dimensionedScalar& SE1, dimensionedScalar& SE2, volScalarField& SE3) const
{
    return (P1_*SE1+P2_*SE2+P3_*SE3);
}

/*void Foam::PDFModel::specieEnv3(volScalarField wS3,volScalarField& SE3) const
{
  scalarField& SCells = wS3.internalField();  
  forAll(SCells, celli)
  {
    if(P3_[celli] >= P_tol_.value())
    {
      SE3[celli] = wS3[celli]/P3[celli];
    }
    else
    {
      SE3[celli] = 0.0;
    }
  }
}*/

void Foam::PDFModel::correct()
{
  calcParameters();  
  fvScalarMatrix P1Eqn
  (
    fvm::ddt(rho_,P1_)
    +fvm::div(rhoPhi_, P1_)
    -fvm::laplacian(effDiff_,P1_)
    ==
    rho_*(gammaS_*P3_-gamma_*P1_*(scalar(1.0)-P1_))
  );
  P1Eqn.relax();
  P1Eqn.solve();
  fvScalarMatrix P2Eqn
  (
    fvm::ddt(rho_, P2_)
    +fvm::div(rhoPhi_, P2_)
    -fvm::laplacian(effDiff_,P2_)
    ==
    rho_*(gammaS_*P3_-gamma_*P2_*(scalar(1.0)-P2_))
  );
  P2Eqn.relax();
  P2Eqn.solve();

  P1_ = min(max(P1_,scalar(0.0)),scalar(1.0));
  P2_ = min(max(P2_,scalar(0.0)),scalar(1.0));

  P1_.correctBoundaryConditions();
  P2_.correctBoundaryConditions();

  P3_ = scalar(1.0)-P1_-P2_;

  forAll(P1_.boundaryFieldRef(), patchi)
  {
    P3_[patchi]= scalar(1.0)-P1_[patchi]-P2_[patchi]; 
  }

  fvScalarMatrix S3Eqn
  (
    fvm::ddt(rho_, S3_)
    + fvm::div(rhoPhi_, S3_)
    - fvm::laplacian(effDiff_, S3_)
    == 
      rho_*(gamma_*(P1_*(scalar(1.0)-P1_)*S1_+P2_*(scalar(1.0)-P2_)*S2_)-gammaS_*P3_*(S1_+S2_))
  );
  S3Eqn.relax();
  S3Eqn.solve();
  S3_=min(max(S3_,scalar(0.0)),scalar(1.0));
  S3_.correctBoundaryConditions();
}

bool Foam::PDFModel::read()
{
    if (regIOobject::read())
    {	  
	  PDFdict_.lookup("C_phi") >> C_phi_;
	  PDFdict_.lookup("Sc_t") >> Sc_t_;
	  PDFdict_.lookup("P_tol") >> P_tol_;
	  PDFdict_.lookup("Cmu") >> Cmu_;
	  PDFdict_.lookup("lamDiff") >> lamDiff_;
	  S1_ = List<dimensionedScalar>(PDFdict_.lookup("phi1"))[0];
	  S2_ = List<dimensionedScalar>(PDFdict_.lookup("phi2"))[0];

          return true;
    }
    else
    {
        return false;
    }
}
// ************************************************************************* //
