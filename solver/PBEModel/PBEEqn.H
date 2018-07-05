// Solve PBE equations
forAll(mesh.C(), celli)
{
  S_solvent[celli] = std::min(std::max(S_solvent[celli], 0.0), rho1[celli]);
  S_antiSolvent[celli] = std::min(std::max(S_antiSolvent[celli], 0.0), rho2[celli]);
  S_solute[celli] = std::min(std::max(S_solute[celli],0.0),std::max(phi1[1].value(),phi2[1].value()));
  //test P3 value to avoid undterminate values
  if(P3[celli] > P_tol.value())
  {
    solute3[celli] = S_solute[celli]/P3[celli];
    solvent3[celli] = S_solvent[celli]/P3[celli];
    antiSolvent3[celli] = S_antiSolvent[celli]/P3[celli];
    forAll(Yw,i)
    {
      Yw[i][celli] = Y[i][celli]/P3[celli];
    }
  }
  else
  {
    solute3[celli] = 0.0;
    solvent3[celli] = 0.0;
    antiSolvent3[celli] = 0.0;
    forAll(Yw,i)
    {
      Yw[i][celli] = (phi1[i+2].value()+phi2[i+2].value())/2.0;
    }
  }
  if((antiSolvent3[celli]+solvent3[celli]) > P_tol.value() && (antiSolvent3[celli] > 0.0))
  {
    csat[celli] = solubility(antiSolvent3[celli], solvent3[celli], T[celli]); //saturation condition
    //nfo << "csat" << endl;
    if(csat[celli] < 0.0) {csat[celli] = 0.0;}    
    if((csat[celli] > 0.0) && (solute3[celli] > 0.0))
    {
      supersat[celli] = solute3[celli]/(antiSolvent3[celli]+solvent3[celli])-csat[celli]; //supersaturation
      //Info << "end supersat" << endl;
      relSupersat[celli] = solute3[celli]/(antiSolvent3[celli]+solvent3[celli])/csat[celli]; //relative supersaturation
      GR[celli] = growthRate(supersat[celli], relSupersat[celli]); //growth rate
      //Info << "GR" << endl;
      NR[celli] = nucRate(supersat[celli], relSupersat[celli]);//nucleation rate
      //Info << "NR" << endl;     
    }
    else
    {
      relSupersat[celli] = 0.0;
      supersat[celli] = 0.0;
      GR[celli] = 0.0;
      NR[celli] = 0.0;
    }

  }
  else
  {
    csat[celli] = 0.0;
    supersat[celli] = 0.0;
    relSupersat[celli] = 0.0;
    GR[celli] = 0.0;
    NR[celli] = 0.0;
  }
  totalParticleConc[celli] = 0.0;
  forAll(Yw,i)
  {
    fj[i+2][celli] = scalar(4.0)*Yw[i][celli]/rhoCrystal.value()/kv.value()/r4_diff[i].value(); 
    totalParticleConc[celli] += P3[celli]*Yw[i][celli];
  }
  fj[0][celli] = 0.0;
  fj[1][celli] = 0.0;
  fj[r4_diff.size()+2][celli] = fj[r4_diff.size()+1][celli];
  fj[r4_diff.size()+3][celli] = fj[r4_diff.size()+1][celli];
  
  //calculate the gradiant using minmod limiter
  for(int i = 0; i<= rend.size(); i++)
  {
    fR[i][celli] = minmod(fj[i][celli],fj[i+1][celli],fj[i+2][celli],dr.value(),theta.value());
  }
  
  S_soluteSource[celli] = 0.0;
  forAll(Yw,i)
  {
    if(supersat[celli] > 0.0)
    {
      if(i == 0)
      {
	S_crystalSource[i][celli] = (rhoCrystal.value()*kv.value()*0.25/dr.value())*r4_diff[i].value()*(-GR[celli]*(fj[i+2][celli]+dr.value()*fR[i+1][celli]/2.0)+GR[celli]*(fj[i+1][celli]+dr.value()*fR[i][celli]/2.0)+NR[celli]);                       
      }
      else
      {
	S_crystalSource[i][celli] = (rhoCrystal.value()*kv.value()*0.25/dr.value())*r4_diff[i].value()*(-GR[celli]*(fj[i+2][celli]+dr.value()*fR[i+1][celli]/2.0)+GR[celli]*(fj[i+1][celli]+dr.value()*fR[i][celli]/2.0));
      }
    }
    else
    {
      S_crystalSource[i][celli] = (rhoCrystal.value()*kv.value()*0.25/dr.value())*r4_diff[i].value()*(-GR[celli]*(fj[i+3][celli]-dr.value()*fR[i+2][celli]/2.0)+GR[celli]*(fj[i+2][celli]-dr.value()*fR[i+1][celli]/2.0));
    }
    S_soluteSource[celli] += -S_crystalSource[i][celli];
  }
}
dictionary pbeSolveDict = mesh.solutionDict().subDict("solvers").subDict("Y");
forAll(Y, i)
{
  fvScalarMatrix MassTransport
  (
    fvm::ddt(rho, Y[i])
    + fvm::div(rhoPhi, Y[i], "div(phi,Y)")
    - fvm::laplacian(effDiff, Y[i], "laplacian(D,Y)")
    ==
    pdf.specieSourceRx(phi1[i+2],phi2[i+2],S_crystalSource[i])
  );
  MassTransport.relax();
  MassTransport.solve(pbeSolveDict);
  Y[i].correctBoundaryConditions();
  Yw[i].correctBoundaryConditions();
}
//solverPerformance::debug = 1;
