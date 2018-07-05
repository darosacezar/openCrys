/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "PBESource.H"
#include <math.h>
#include <algorithm>
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

double minmod(double fL, double f, double fR, double dr, double theta)
{
  double grad1, grad2, grad3, grad;
  grad1 = theta*(f-fL)/dr;
  grad2 = 0.5*(fR-fL)/dr;
  grad3 = theta*(fR-f)/dr;
  if ((grad1>0.0) && (grad2>0.0) && (grad3>0.0)) 
  {
    grad = std::min(grad1,grad2);
    grad = std::min(grad,grad3);
  }
  else if ((grad1<0.0) && (grad2<0.0) && (grad3<0.0)) 
  {
    grad = std::max(grad1,grad2);
    grad = std::max(grad, grad3);
  }
  else {grad = 0.0;}
  return grad;  
}
double growthRate(double supersat, double relSupersat)
{
  double G;
  if(supersat > 0.0)
  {
    G = 8.333e-30*pow((2.4623e3*log(relSupersat)),6.7);
  }
  else if (supersat < 0.0)
  {
    G = 0.0;
  }
  else
  {
    G = 0.0;
  }
  return G;
}
double nucRate(double supersat, double relSupersat)
{
  double B_homo, B_hetero, B;
  if(supersat > 0.0)
  {
    B_homo = 6.9656e14 * exp(-15.7647/((log(relSupersat))*(log(relSupersat))));
    B_hetero = 2.19196e8 * exp(-0.994387/((log(relSupersat))*(log(relSupersat))));
    B = B_homo+B_hetero;
  }
  else if (supersat < 0.0)
  {
    B = 0.0;
  }
  else
  {
    B = 0.0;
  }
  return B;
}
double solubility(double antiSolvent, double solvent, double T)
{
  double csat, Was, theta;
  if((antiSolvent+solvent)>0.0)
  {
    Was = 100.0*antiSolvent/(antiSolvent+solvent);
  }
  else{Was = 0.0;}
  Was = std::min(std::max(Was, 0.0),99.0);
  theta = T/296.0;
  if(Was <= 45.67)
  {
    csat = 0.001*exp(15.45763*(1.0-1.0/theta))*(-2.7455e-4*pow(Was,3.0)+3.3716e-2*pow(Was,2)-1.6704*Was+33.089);
  }
  else 
  {
    csat = 0.001*exp(15.45763*(1.0-1.0/theta))*(-1.7884e-2*Was+1.7888);
  }
  return csat;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //