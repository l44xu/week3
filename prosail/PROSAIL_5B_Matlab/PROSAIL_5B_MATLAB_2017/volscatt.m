function [chi_s,chi_o,frho,ftau]=volscatt(tts,tto,psi,ttl)

%********************************************************************************
%*	tts		= solar zenith
%*	tto		= viewing zenith
%*	psi		= azimuth
%*	ttl		= leaf inclination angle
%*	chi_s	= interception functions
%*	chi_o	= interception functions
%*	frho	= function to be multiplied by leaf reflectance rho
%*	ftau	= functions to be multiplied by leaf transmittance tau
%********************************************************************************

%	Compute volume scattering functions and interception coefficients
%	for given solar zenith, viewing zenith, azimuth and leaf inclination angle.

%	chi_s and chi_o are the interception functions.
%	frho and ftau are the functions to be multiplied by leaf reflectance rho and
%	leaf transmittance tau, respectively, in order to obtain the volume scattering
%	function.

%	Wout Verhoef, april 2001, for CROMA

% REAL(KIND=8),INTENT(in) :: tts
% REAL(KIND=8),INTENT(in) :: tto
% REAL(KIND=8),INTENT(in) :: psi
% REAL(KIND=8),INTENT(in) :: ttl
% REAL(KIND=8),INTENT(inout) :: chi_s
% REAL(KIND=8),INTENT(inout) :: chi_o
% REAL(KIND=8),INTENT(inout) :: frho
% REAL(KIND=8),INTENT(inout) :: ftau
% 
% REAL(KIND=8) costs,costo,sints,sinto,cospsi
% REAL(KIND=8) psir
% REAL(KIND=8) costl,sintl,cs,co,ss,so,ds
% REAL(KIND=8) cosbts,cosbto,bts,bto
% REAL(KIND=8) btran1,btran2,bt1,bt2,bt3,t1,t2
% REAL(KIND=8) doo
% REAL(KIND=8) denom
rd=pi/180;
costs=cos(rd*tts);
costo=cos(rd*tto);
sints=sin(rd*tts);
sinto=sin(rd*tto);
cospsi=cos(rd*psi);
psir=rd*psi;
costl=cos(rd*ttl);
sintl=sin(rd*ttl);
cs=costl*costs;
co=costl*costo;
ss=sintl*sints;
so=sintl*sinto;

%c ..............................................................................
%c     betas -bts- and betao -bto- computation
%c     Transition angles (beta) for solar (betas) and view (betao) directions
%c     if thetav+thetal>pi/2, bottom side of the leaves is observed for leaf azimut 
%c     interval betao+phi<leaf azimut<2pi-betao+phi.
%c     if thetav+thetal<pi/2, top side of the leaves is always observed, betao=pi
%c     same consideration for solar direction to compute betas
%c ..............................................................................


cosbts=5;
if (abs(ss)>1e-6)
	cosbts=-cs/ss;
end
cosbto=5;
if (abs(so)>1e-6)
	cosbto=-co/so;
end

if (abs(cosbts)<1)
	bts=acos(cosbts);
	ds=ss;
else
	bts=pi;
	ds=cs;
end

chi_s=2./pi*((bts-pi*.5)*cs+sin(bts)*ss);

if (abs(cosbto)<1)
	bto=acos(cosbto);
	doo=so;
elseif(tto<90)
	bto=pi;
	doo=co;
else
	bto=0;
	doo=-co;
end
chi_o=2./pi*((bto-pi*.5)*co+sin(bto)*so);

%c ..............................................................................
%c   Computation of auxiliary azimut angles bt1, bt2, bt3 used          
%c   for the computation of the bidirectional scattering coefficient w              
%c .............................................................................

btran1=abs(bts-bto);
btran2=pi-abs(bts+bto-pi);

if (psir<=btran1)
	bt1=psir;
	bt2=btran1;
	bt3=btran2;
else
	bt1=btran1;
	if (psir<=btran2)
		bt2=psir;
		bt3=btran2;
    else
		bt2=btran2;
		bt3=psir;
    end
end

t1=2.*cs*co+ss*so*cospsi;
t2=0;
if (bt2>0)
	t2=sin(bt2)*(2.*ds*doo+ss*so*cos(bt1)*cos(bt3));
end

denom=2.*pi*pi;
frho=((pi-bt2)*t1+t2)/denom;
ftau=    (-bt2*t1+t2)/denom;

if (frho<0)
	frho=0;
end

if (ftau<0)
	ftau=0;
end
