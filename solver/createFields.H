//Create main variables of the model
IOdictionary transportProperties
(
  IOobject
  (
    "transportProperties",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
  )
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
  IOobject
  (
    "p_rgh",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
  IOobject
  (
    "U",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);
#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
twoPhaseMixtureThermoTransportModel thermoMixture(U, phi);

volScalarField& alpha1(thermoMixture.alpha1());
volScalarField& alpha2(thermoMixture.alpha2());


dimensionedScalar Dab(thermoMixture.lookup("Dab"));

// Read the reciprocal of the turbulent Schmidt number
dimensionedScalar alphatab(thermoMixture.lookup("alphatab"));

//thermoModel
volScalarField& rho1(thermoMixture.rho1());
volScalarField& rho2(thermoMixture.rho2());
volScalarField& rho(thermoMixture.rho());
//volScalarField& Cp(thermoMixture.Cp());

volScalarField& K(thermoMixture.K());
volScalarField& he(thermoMixture.he());
volScalarField& T(thermoMixture.T());

volScalarField dHmix
(
  IOobject
  (
    "dHmix",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("dHmix", dimensionSet (0,2,-2,0,0,0,0), 0.0)
);


surfaceScalarField rhoPhi
(
  IOobject
  (
    "rhoPhi",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  fvc::interpolate(rho)*phi
);

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
  incompressible::turbulenceModel::New(U, phi, thermoMixture)
);
const volScalarField& k(turbulence->k());
const volScalarField& epsilon(turbulence->epsilon());
//declares and upadates the value of alphaEff
volScalarField& alphaEff(thermoMixture.alphaEff());

Info<< "Calculating field g.h\n" << endl;
volScalarField gh("gh", g & mesh.C());
surfaceScalarField ghf("ghf", g & mesh.Cf());

volScalarField p
(
  IOobject
  (
    "p",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
  p,
  p_rgh,
  mesh.solutionDict().subDict("PIMPLE"),
  pRefCell,
  pRefValue
);

if (p_rgh.needReference())
{
  p += dimensionedScalar
  (
    "p",
    p.dimensions(),
    pRefValue - getRefCellValue(p, pRefCell)
  );
  p_rgh = p - rho*gh;
}
