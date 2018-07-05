// Solve Species Transport Equations on Environment 3
fvScalarMatrix S_soluteEqn
(
  fvm::ddt(rho, S_solute)
  + fvm::div(rhoPhi, S_solute)
  - fvm::laplacian(effDiff, S_solute)
  == 
    pdf.specieSourceRx(phi1[1],phi2[1],S_soluteSource)
);
S_soluteEqn.relax();
S_soluteEqn.solve();
S_solute.correctBoundaryConditions();
// Solve Solvent Equation on Environment 3
fvScalarMatrix S_solventEqn
(
  fvm::ddt(rho, S_solvent)
  + fvm::div(rhoPhi, S_solvent)
  - fvm::laplacian(effDiff, S_solvent)
  == 
    pdf.specieSourceEnv1(rho1)
);
S_solventEqn.relax();
S_solventEqn.solve();
S_solvent.correctBoundaryConditions();
// Solve AntiSolvent Equation on Environment 3
fvScalarMatrix S_antiSolventEqn
(
  fvm::ddt(rho, S_antiSolvent)
  + fvm::div(rhoPhi, S_antiSolvent)
  - fvm::laplacian(effDiff, S_antiSolvent)
  == 
    pdf.specieSourceEnv2(rho2)
);
S_antiSolventEqn.relax();
S_antiSolventEqn.solve();
//S_antiSolvent = max(S_antiSolvent, phi1[6]);
S_antiSolvent.correctBoundaryConditions();

solute = P1*phi1[1]+P2*phi2[1]+P3*solute3;
solvent = P1*rho1+P3*solvent3;
antiSolvent = P2*rho2+P3*antiSolvent3;

alpha1 = min(max(solvent/rho1,scalar(0.0)),scalar(1.0));
alpha1.correctBoundaryConditions();
alpha2 = scalar(1.0)-alpha1;

solute.correctBoundaryConditions();
