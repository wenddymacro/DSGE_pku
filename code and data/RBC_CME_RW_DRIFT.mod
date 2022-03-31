/*
 * This file accompoanies Johannes Pfeifer (2013): "A Guide to Specifying Observation 
 * Equations for the Estimation of DSGE Models". It is "Listing 9: Nonlinear Model with Explicit Trend for Log-Linearization"
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

var y_tilde c_tilde k_tilde z R_tilde Pi_tilde mu_tilde Y;
varexo eps_z eps_x;

parameters beta delta alpha rhoz phi_pi Pibar Lambda_x;

alpha   = 0.33;     // capital share
delta   = 0.025;    //deprecation rate
beta    = 0.99;     //discount factor
Pibar   = 1;
phi_pi  = 1.5;
rhoz    = 0.97;     //TFP autocorr. from linearly detrended Solow residual
Lambda_x= 0.0055;   //2.2% output growth per year

model;
#Rbar    = exp(Lambda_x)/beta;
1/exp(c_tilde)=beta/exp(mu_tilde(+1))*1/exp(c_tilde(+1))*(alpha*exp(z(+1))*(exp(k_tilde)/exp(mu_tilde(+1)))^(alpha-1)+(1-delta));
1/exp(c_tilde)=beta/exp(mu_tilde(+1))*1/exp(c_tilde(+1))*(exp(R_tilde)/exp(Pi_tilde(+1)));
exp(z)*(exp(k_tilde(-1))/exp(mu_tilde))^alpha=exp(c_tilde)+exp(k_tilde)-(1-delta)*exp(k_tilde(-1))/exp(mu_tilde);
exp(y_tilde)=exp(z)*(exp(k_tilde(-1))/exp(mu_tilde))^alpha;
exp(R_tilde)/Rbar=(exp(Pi_tilde)/Pibar)^phi_pi;
z=rhoz*z(-1)+eps_z;
mu_tilde=Lambda_x+eps_x;
Y=y_tilde*exp(mu_tilde);
end;

steady_state_model;
mu_tilde=Lambda_x;
k_tilde=log(exp(mu_tilde)*((exp(mu_tilde)/beta-(1-delta))/alpha)^(1/(alpha-1)));
y_tilde=log((exp(k_tilde)/exp(mu_tilde))^alpha);
c_tilde=log(exp(y_tilde)-(1-(1-delta)/exp(mu_tilde))*exp(k_tilde));
R_tilde=log(exp(mu_tilde)/beta);
Pi_tilde=log(Pibar);
z=0;
Y=y_tilde*exp(mu_tilde);
end;


shocks;
var eps_z=0.0068^2;
var eps_x=0.005^2; //just some number
end;

steady;
check;

stoch_simul(order=1,irf=20,periods=200);

rplot Y;