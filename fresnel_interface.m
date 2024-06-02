%---- Juan P. Martinez, reference as 1. Martinez, J. P. Light propagation in multilayered nanostructures. (2024) doi:10.13140/RG.2.2.30332.96640.
% This code calculates the Fresnel coefficients for s (perpendicular to the
% plane of incidence and p (parallel to plane of incidence) polarizations
% at each hinterface of a stack of layers (also called stratified media).

% ------The INPUTS are--------------:
% 
% --> n, which is a vector with the refractive index (possibly
% complex) of each layer in order of approach of light rays. Notice that the first and last elements of the vector are the incoming and outgoing media If light is
% incident in a stack of 4 layers, which have air (n=1) on one side and water on another (n=1.3) n=[1,n1,n2,n3,n4,1.3]. 

% CONVENTION: NEGATIVE extinction coefficient k for absorption. complex refractive
% index is n_c=n-ik

% --> phi0 is the angle of incidence in degrees

%------ OUTPUTS --------------
% --> rs,ts,rp,tp are vectors of each fresnel reflection (r) and
% transmission (t) coefficients at the interfaces. The order is the same as
% the n vector. Notice that if n has y elements there are y-1 elements, and
% thus the length of the output vectors is y-1

% --> phi is a vector with the angle in each media in RADIANS!. The first element is
% the phi0 converted to radians, and the values in each media are
% calculated using Snell's law. Note that the angles can be complex if the
% refractive index are complex.

function [rs,ts,rp,tp,phi]=fresnel_interface(n,phi0)
rs=zeros([1,length(n)-1]);
ts=zeros([1,length(n)-1]);
rp=zeros([1,length(n)-1]);
tp=zeros([1,length(n)-1]);
phi=zeros([1,length(n)]);
phi(1)=deg2rad(phi0);
for k=1:length(n)-1
  
    a=n(k)*cos(phi(k));
    phi(k+1)=asin((n(k)/n(k+1))*sin(phi(k)));
    b=n(k+1)*cos(phi(k+1));
    rs(k)=(a-b)/(a+b);
    ts(k)=(2*a)/(a+b);
    %------
    c=n(k)*cos(phi(k+1));%this is the refractive index in the incident media times cosine of the angle transmitted
    d=n(k+1)*cos(phi(k));% this is the refractive index in the transmitted media times cosine of the angle of incidence
    
    rp(k)=(d-c)/(d+c);
    tp(k)=(2*a)/(d+c);
   
end
end

    
  
    
    