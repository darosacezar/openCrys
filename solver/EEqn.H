// Solve Energy Equation
thermoMixture.correct();

scalar f1;
forAll(mesh.C(), celli)
{
  if((P3[celli] > P_tol.value()) && (antiSolvent3[celli]+solvent3[celli] > P_tol.value()))
  {
    f1 = antiSolvent3[celli]/(antiSolvent3[celli]+solvent3[celli]);      
    dHmix[celli] = -5.90806E+05*std::pow(f1,6.0)+1.45022E+06*std::pow(f1,5.0)-9.73106E+05*std::pow(f1,4.0)+3.90128E+04*std::pow(f1,3.0)+2.27166E+05*std::pow(f1,2.0)-1.52557E+05*f1+6.62235E+01;
  }
  else
  {
    dHmix[celli] = 0.0;
  }
}
fvScalarMatrix EEqn
(
  fvm::ddt(rho,he) + fvm::div(rhoPhi, he)
  +fvc::ddt(rho,K) + fvc::div(rhoPhi, K)
  -fvm::laplacian(alphaEff,he)
 ==
  fvc::div(rho*turbulence->nuEff()*(fvc::grad(U)&U))
  +fvc::div(rho*turbulence->nuEff()*(dev2(fvc::grad(U)().T())&U))
  -S_soluteSource*dH_crys
  -pdf.specieVolSourceEnv12Rx(rho1,rho2,dHmix)
  //+fvOptions(rho,he)
);

EEqn.relax();

EEqn.solve();

he.correctBoundaryConditions();
