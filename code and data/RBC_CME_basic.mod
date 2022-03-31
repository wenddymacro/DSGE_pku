/*
 * This file accompoanies Johannes Pfeifer (2013): "A Guide to Specifying Observation 
 * Equations for the Estimation of DSGE Models". It is "Listing 1: Basic RBC Classical Monetary Economy Model"
 *
 * This file was written by Johannes Pfeifer. In
 * case you spot mistakes, email me at jpfeifer@gmx.de
 *
 * The model is written in Dynare's end of period stock notation.
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model
 */

/*
 * Copyright (C) 2013 Johannes Pfeifer
 *
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This file is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You can receive a copy of the GNU General Public License
 * at <http://www.gnu.org/licenses/>.
 */

var y c k A z R Pi;
varexo eps_z;

parameters beta delta alpha rhoz phi_pi Pibar;

alpha   = 0.33;   // capital share
delta   = 0.025;  //deprecation rate
beta    = 0.99;   //discount factor
Pibar   = 1;
phi_pi  = 1.5;
rhoz    = 0.97;   //TFP autocorr. from linearly detrended Solow residual

model;
#Rbar    = 1/beta;
1/c=beta*1/c(+1)*(alpha*A(+1)*k^(alpha-1)+(1-delta));
1/c=beta*1/c(+1)*(R/Pi(+1));
A*k(-1)^alpha=c+k-(1-delta)*k(-1);
y=A*k(-1)^alpha;
R/Rbar=(Pi/Pibar)^phi_pi;
A=exp(z);
z=rhoz*z(-1)+eps_z;
end;

steady_state_model;
k=((1/beta-(1-delta))/alpha)^(1/(alpha-1));
y=k^alpha;
c=y-delta*k;
R=1/beta;
Pi=Pibar;
A=1;
z=0;
end;

shocks;
var eps_z=0.0068^2; //estimated value
end;

steady;
check;

stoch_simul(order=1,irf=20,periods=250);
rplot y;