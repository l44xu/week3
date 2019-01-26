function Jout=Jfunc2(k,l,t)
%	J2 function
Jout=(1.-exp(-(k+l)*t))./(k+l);
