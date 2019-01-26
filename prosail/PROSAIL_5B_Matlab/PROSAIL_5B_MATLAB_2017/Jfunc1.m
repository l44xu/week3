function Jout=Jfunc1(k,l,t)
%	J1 function with avoidance of singularity problem
del=(k-l)*t;
Jout(abs(del)>1e-3)=(exp(-l(abs(del)>1e-3)*t)-exp(-k*t))./(k-l(abs(del)>1e-3));
Jout(abs(del)<=1e-3)=0.5*t*(exp(-k*t)+exp(-l(abs(del)<=1e-3)*t)).*(1-del(abs(del)<=1e-3).*del(abs(del)<=1e-3)/12);
Jout=transpose(Jout);