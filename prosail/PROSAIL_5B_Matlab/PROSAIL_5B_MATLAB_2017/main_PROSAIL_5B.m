% %*************************************************************************
% %*                                                                       *
% 	main_PROSAIL
% 
% 2017-12-28
% This program allows modeling reflectance data from canopy
% - modeling leaf optical properties with PROSPECT-5 (feret et al. 2008)
% - modeling leaf inclination distribution function with the subroutine campbell
% (Ellipsoidal distribution function caracterised by the average leaf 
% inclination angle in degree), or dladgen (2 parameters LIDF)
% - modeling canopy reflectance with 4SAIL (Verhoef et al., 2007)

% This version has been implemented by Jean-Baptiste Feret
% Jean-Baptiste Feret takes the entire responsibility for this version 
% All comments, changes or questions should be sent to:
% Jean-Baptiste FERET
% UMR-TETIS, IRSTEA Montpellier
% Maison de la Télédétection
% 500 rue Jean-Fracois Breton
% 34093 Montpellier cedex 5
% E-mail: jb.feret@teledetection.fr

% References:
% 	Verhoef et al. (2007) Unified Optical-Thermal Four-Stream Radiative
% 	Transfer Theory for Homogeneous Vegetation Canopies, IEEE TRANSACTIONS 
% 	ON GEOSCIENCE AND REMOTE SENSING, VOL. 45, NO. 6, JUNE 2007
% 	Féret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
% 	Properties Model Separating Photosynthetic Pigments, REMOTE SENSING OF 
% 	ENVIRONMENT
% The specific absorption coefficient corresponding to brown pigment is
% provided by Frederic Baret (EMMAH, INRA Avignon, baret@avignon.inra.fr)
% and used with his autorization.
% the model PRO4SAIL is based on a version provided by
%	Wout Verhoef
%	NLR	
%	April/May 2003

% The original 2-parameter LIDF model is developed by and described in:
% 	W. Verhoef, 1998, "Theory of radiative transfer models applied in 
%	optical remote sensing of vegetation canopies", Wageningen Agricultural
%	University,	The Netherlands, 310 pp. (Ph. D. thesis)
% the Ellipsoidal LIDF is taken from:
%   Campbell (1990), Derivtion of an angle density function for canopies 
%   with ellipsoidal leaf angle distribution, Agricultural and Forest 
%   Meteorology, 49 173-176
%*                                                                       *
%*************************************************************************
%** Bug correction
%** - lidf values harmonized among datasets: modification in campbell.m
%** this problem was not present in the fortran version 
%*************************************************************************

clc
TypeLidf=2;
% if 2-parameters LIDF: TypeLidf=1
if (TypeLidf==1)
    % LIDFa LIDF parameter a, which controls the average leaf slope
    % LIDFb LIDF parameter b, which controls the distribution's bimodality
    %	LIDF type 		a 		 b
    %	Planophile 		1		 0
    %	Erectophile    -1	 	 0
    %	Plagiophile 	0		-1
    %	Extremophile 	0		 1
    %	Spherical 	   -0.35 	-0.15
    %	Uniform 0 0
    % 	requirement: |LIDFa| + |LIDFb| < 1	
    LIDFa	=	-0.35;
    LIDFb	=	-0.15;
% if ellipsoidal LIDF: TypeLidf=2
elseif (TypeLidf==2)
    % 	LIDFa	= average leaf angle (degrees) 0 = planophile	/	90 = erectophile
    % 	LIDFb = 0
    LIDFa	=	30;
    LIDFb	=	0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEAF CHEM & STR PROPERTIES	%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cab		= 40;		% chlorophyll content (µg.cm-2) 
Car		= 8;		% carotenoid content (µg.cm-2)
Cbrown	= 0.0;	% brown pigment content (arbitrary units)
Cw		= 0.01;	% EWT (cm)
Cm		= 0.009;	% LMA (g.cm-2)
N		= 2.5;	% structure coefficient
data=dataSpec_P5B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Soil Reflectance Properties	%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rsoil1 = dry soil
% rsoil2 = wet soil
Rsoil1  = data(:,10);Rsoil2=data(:,11);
psoil	= 1;		% soil factor (psoil=0: wet soil / psoil=1: dry soil)
rsoil0  = psoil*Rsoil1+(1-psoil)*Rsoil2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	4SAIL canopy structure parm	%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAI		=	5.;     	% leaf area index (m^2/m^2)
hspot	=	0.01;       % hot spot
tts		=	30.;		% solar zenith angle (°)
tto		=	10.;		% observer zenith angle (°)
psi		=	90.;         % azimuth (°)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          CALL PRO4SAIL       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rdot: hemispherical-directional reflectance factor in viewing direction    
% rsot: bi-directional reflectance factor
% rsdt: directional-hemispherical reflectance factor for solar incident flux
% rddt: bi-hemispherical reflectance factor
[rdot,rsot,rddt,rsdt]=PRO4SAIL(N,Cab,Car,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,LAI,hspot,tts,tto,psi,rsoil0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	direct / diffuse light	%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the direct and diffuse light are taken into account as proposed by:
% Francois et al. (2002) Conversion of 400–1100 nm vegetation albedo 
% measurements into total shortwave broadband albedo using a canopy 
% radiative transfer model, Agronomie
% Es = direct
% Ed = diffuse
Es  = data(:,8);
Ed  = data(:,9);
rd  = pi/180;
skyl	=	0.847- 1.61*sin((90-tts)*rd)+ 1.04*sin((90-tts)*rd)*sin((90-tts)*rd); % % diffuse radiation
PARdiro	=	(1-skyl)*Es;
PARdifo	=	(skyl)*Ed;

% resv : directional reflectance
resv	= (rdot.*PARdifo+rsot.*PARdiro)./(PARdiro+PARdifo);
dlmwrite('Refl_CAN.txt',[data(:,1),resv],'delimiter','\t','precision',5)
figure,plot(data(:,1),resv)
axis([400 2500 0 0.6])
