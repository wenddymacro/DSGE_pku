var c g i l;

varexo ei;

parameters beta piss rho phi1 phi2 lss chi omega psi;

beta=0.9;
piss=1;
rho=0.9;
phi1=1.5;
phi2=1;
lss=1;
chi=0.1;
omega=0.1;
psi=1;

model;

1+i=(piss/beta)*piss^phi1*(l/lss)^phi2+ei;

c=g*c(1)*piss/(beta*(piss/beta)*piss^phi1*(l/lss)^phi2);

g(1)=beta*(c/c(1))*(chi*omega*l(1)+1);

psi*l=c+(g(1)-1)/chi;

end;

steady_state_model;

i=piss/beta-1;
l=lss;
g=beta*(chi*omega*l+1);
c=psi*l-(g-1)/chi;
end;

shocks;
var ei=0.01^2;
end;

steady;
check;
model_diagnostics;

stoch_simul;