%----- This code implements the Transfer Matrix Method  used to calculate propagation in stratified media (a stack of layers) in the Formulation
%of Azzam and Bashara. 1. Azzam, R. M. & Bashara, N. M. Ellipsometry and Polarized Light. vol. 1 (North-Holland, 1977).

% Juan P. Martinez, reference as Martinez, J. P. Light propagation in multilayered nanostructures. (2024) doi:10.13140/RG.2.2.30332.96640.

%--- Note that this code requires the function fresnel_interface.m

%---- INPUTS---------------

% --> n, is a vector with the refractive index (possibly
% complex) of each layer in order of approach of light rays. Notice that the first and last elements of the vector are the incoming and outgoing media If light is
% incident in a stack of 4 layers, which have air (n=1) on one side and water on another (n=1.3) n=[1,n1,n2,n3,n4,1.3]. 

% CONVENTION: NEGATIVE extinction coefficient k for absorption. complex refractive
% index is n_c=n-ik

% --> phi0 is the angle of incidence in degrees

% --> e is a vector with the thickness of each layer in nanometers (nm).
% Notice that the incoming and outgoing medium are semi-infinite, i.e. they
% do not have a thickness, and thus if length (n) = y then length (e) = y-2

%--> Lambda is the wavelength of light in nanometers

%----- OUTPUTS ------

% capital leters are quotient of intensities (R,T), lower case are of
% fields (r,t)

%--> Rs,Ts,Rp,Tp Reflectance (R) and Transmittance (T) of the stratified
%media, for the s (perpendicular to plane of incidence, also called TE) polarization and the p (parallel to
%the plane of incidence, also called TM) polarization.

%--> Ms and Mp : global transfer matrix of the stratified media for p and s
%polarizations

%--> rs,ts,rp,tp reflection (r) and transmission (t) coefficients of the stratified
%media, for the s (perpendicular to plane of incidence, also called TE) polarization and the p (parallel to
%the plane of incidence, also called TM) polarization.

function [Rs,Ts,Rp,Tp,Ms,Mp,rs,ts,rp,tp]=TMM_fresnel(n,e,phi0,lambda)%e y lambda en mismas unidades
%ko=2*pi/lambda;%calculo ko
%phi=deg2rad(phi0);%defino el primer phi (angulo de incidencia)
Ms=eye(2);
Mp=eye(2);
%-----Fresnel Coeficientes-----------------
[rs,ts,rp,tp,phi]=fresnel_interface(n,phi0);
for k=1:length(e)
    %------For s-polarized-----
    Ds=(1/ts(k))*[1,rs(k);rs(k),1];
    %path=(2*pi/lambda)*n(k+1)*(e(k)/cos(phi(k+1)));
    %path=(2*pi/lambda)*e(k)*n(k+1);
    path=(2*pi/lambda)*(n(k+1))*(e(k)*cos(phi(k+1)));
    phasedif=-1i*path;
    Ps=[1/(exp(phasedif)),0;0,exp(phasedif)];
    Ms=Ms*Ds*Ps;
    
    %-------For p-polarized------
    Dp=(1/tp(k))*[1,rp(k);rp(k),1];
    %----Aca agregaria otro path si tiene otro indice de refraccion
    Mp=Mp*Dp*Ps;
    
end
%----The last interface (with outgoing media)
%-----For s-polarized--------
Dsf=(1/ts(end))*[1,rs(end);rs(end),1];
Ms=Ms*Dsf;

%------For p-polarized-------
Dpf=(1/tp(end))*[1,rp(end);rp(end),1];
Mp=Mp*Dpf;

%------Now obtain r and t----------
% t=1/m11 y r=m21/m11

%---- s-polarized -----
ts=1/Ms(1,1);
rs=Ms(2,1)/Ms(1,1);

Rs=(abs(rs))^2;
%Ts=real(((n(end)*cos(phi(end)))/(n(1)*cos(deg2rad(phi0)))))*(abs(ts))^2;
Ts=(real(cos(phi(end)))/real(cos(deg2rad(phi0))))*(abs(ts))^2;

%---- p-polarized -----
tp=1/Mp(1,1);
rp=Mp(2,1)/Mp(1,1);

Rp=(abs(rp))^2;
%Tp=real(((n(end)*cos(phi(end)))/(n(1)*cos(deg2rad(phi0)))))*(abs(tp))^2;
Tp=(real(cos(phi(end)))/real(cos(deg2rad(phi0))))*(abs(tp))^2;
    
    
    