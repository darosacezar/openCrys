Info<< "Reading PDFdict\n" << endl;
//Create PDF dictionary
IOdictionary PDFdict
(
  IOobject
  (
    "PDFdict",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
  )
);
//Read data from PDF dictionary

List<dimensionedScalar> phi1(PDFdict.lookup("phi1"));
List<dimensionedScalar> phi2(PDFdict.lookup("phi2"));
//Create probabilities for the environments
forAll(phi1,i)
{
  if(i > 0)
  {
    phi1[i].dimensions().reset(dimMass/dimVolume);
    phi2[i].dimensions().reset(dimMass/dimVolume);
  }
}

PDFModel pdf(U,rho,rhoPhi, k, epsilon);
volScalarField& P1(pdf.P1());
volScalarField& P2(pdf.P2());
volScalarField& P3(pdf.P3());
volScalarField& effDiff(pdf.effDiff());

dimensionedScalar& P_tol(pdf.P_tol());

volScalarField S_solute
(
  IOobject
  (
    "S_solute",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);
//create the weighted concentration of the solvent in environment 3
volScalarField S_solvent
(
  IOobject
  (
    "S_solvent",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);
//create the weighted concentration of the antisolvent in environment 3
volScalarField S_antiSolvent
(
  IOobject
  (
    "S_antiSolvent",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);
volScalarField S_soluteSource
(
  IOobject
  (
    "S_soluteSource",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("S_soluteSource", dimensionSet (1,-3,-1,0,0,0,0), 0.0)
);
//Create turbulent effective diffusivity
volScalarField hmixSource
(
  IOobject
  (
    "hmixSource",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("hmixSource", dimensionSet (1,-1,-3,0,0,0,0), 0.0)
);


