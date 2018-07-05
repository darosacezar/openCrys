Info<< "Reading PBEdict\n" << endl;

IOdictionary PBEdict
(
  IOobject
  (
    "PBEdict",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
  )
);
//create and read bins for population balance
List<dimensionedScalar> rend(PBEdict.lookup("rend"));
List<dimensionedScalar> particleSize(rend.size()-1);
//read laminar diffusivity for population balance
dimensionedScalar diffLam(PBEdict.lookup("diffusivity"));

dimensionedScalar rhoCrystal(PBEdict.lookup("density"));
dimensionedScalar kv(PBEdict.lookup("volumeShapeFactor"));
dimensionedScalar dH_crys(PBEdict.lookup("dH_crystallization"));
dimensionedScalar theta(PBEdict.lookup("minmodTheta"));
List<word> inletPatchName(PBEdict.lookup("inletPatchName"));
List<word> outletPatchName(PBEdict.lookup("outletPatchName"));
List<dimensionedScalar> r4_diff(rend.size()-1);
PtrList<volScalarField> Y(rend.size()-1);
PtrList<volScalarField> Yw(rend.size()-1);
PtrList<volScalarField> fj(rend.size()+3);
PtrList<volScalarField> fR(rend.size()+1);
PtrList<volScalarField> S_crystalSource(rend.size()-1);
word name;
for(int i = 0; i < (rend.size()-1); i++)
{
  name = itoa(i+1);
  Y.set
  (
    i,
    new volScalarField
    (
      IOobject
      (
	"C"+name,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
      ),
      mesh
    )
  );
  Yw.set
  (
    i,
    new volScalarField
    (
      IOobject
      (
	"Cw"+name,
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
      ),
      Y[i] 
    )
  );
  S_crystalSource.set
  (
    i,
    new volScalarField
    (
      IOobject
      (
	"CS"+name,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("CS"+name, dimensionSet (1,-3,-1,0,0,0,0), 0.0)      
    )
  );
  particleSize[i] = (rend[i+1].value()+rend[i].value())/2.0;
  rend[i].dimensions().reset(dimLength);
  r4_diff[i] = std::pow(rend[i+1].value(),4.0)-std::pow(rend[i].value(),4.0);
  r4_diff[i].dimensions().reset(dimArea*dimArea);
  particleSize[i].dimensions().reset(dimLength);
}
forAll(fj,i)
{
  name = "fj"+itoa(i);
  fj.set
  (
    i,
    new volScalarField
    (
      IOobject
      (
	name,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar(name, dimensionSet (0,-4,0,0,0,0,0), 0.0)      
    )
  );
}
forAll(fR,i)
{
  name = "fR"+itoa(i);
  fR.set
  (
    i,
    new volScalarField
    (
      IOobject
      (
	name,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar(name, dimensionSet (0,-5,0,0,0,0,0), 0.0)      
    )
  );
}
dimensionedScalar dr("dr", rend[2]-rend[1]);
volScalarField csat
(
  IOobject
  (
    "csat",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("csat", dimensionSet (0,0,0,0,0,0,0), 0.0)
);
volScalarField supersat
(
  IOobject
  (
    "supersat",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  csat
);
volScalarField relSupersat
(
  IOobject
  (
    "relSupersat",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("relSupersat", dimensionSet (0,0,0,0,0,0,0), 0.0)
);
volScalarField GR
(
  IOobject
  (
    "GR",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("GR", dimensionSet (0,1,-1,0,0,0,0), 0.0)
);
volScalarField NR
(
  IOobject
  (
    "NR",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("NR", dimensionSet (0,-3,-1,0,0,0,0), 0.0)
);
volScalarField meanSauterPartSize
(
  IOobject
  (
    "meanSauterPartSize",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("meanSauterPartSize", dimensionSet (0,1,0,0,0,0,0), 0.0)
);
volScalarField solute3
(
  IOobject
  (
    "solute3",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("solute3", dimensionSet (1,-3,0,0,0,0,0), 0.0)
);
volScalarField solvent3
(
  IOobject
  (
    "solvent3",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("solvent3", dimensionSet (1,-3,0,0,0,0,0), 0.0)
);
volScalarField antiSolvent3
(
  IOobject
  (
    "antiSolvent3",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("antiSolvent3", dimensionSet (1,-3,0,0,0,0,0), 0.0)
);
volScalarField solute
(
  IOobject
  (
    "solute",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("solute", dimensionSet (1,-3,0,0,0,0,0), 0.0)
);
volScalarField solvent
(
  IOobject
  (
    "solvent",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("solvent", dimensionSet (1,-3,0,0,0,0,0), 0.0)
);
volScalarField antiSolvent
(
  IOobject
  (
    "antiSolvent",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("antiSolvent", dimensionSet (1,-3,0,0,0,0,0), 0.0)
);
volScalarField totalParticleConc
(
  IOobject
  (
    "totalParticleConc",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("totalParticleConc", dimensionSet (1,-3,0,0,0,0,0), 0.0)
);
