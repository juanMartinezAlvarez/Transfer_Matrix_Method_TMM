%---- Juan P. Martinez, reference as 1. Martinez, J. P. Light propagation in multilayered nanostructures. (2024) doi:10.13140/RG.2.2.30332.96640.

% This function calculates the Effective Medium Approximation (EMA) using
% Bruggeman formulation (the effective medium is itself the host). It works for a mixture of TWO MATERIALS The work
% by Rouseel et.al is used to obtain the correct root of the second degree
% polynomial. If more than two materials are present, the mixture can be
% performed recursively (i.e mix 2 first and the resulting EMA with the
% third, etc). See Nazarov et.al.

%-- 2. Rouseel, Ph. J., Vanhellemont, J. & Maes, H. E. Numerical aspects of the implementation of effective-medium approximation models in spectroscopic ellipsometry regression software. Thin Solid Films 234, 423–427 (1993).
%-- 3. Nazarov, R., Zhang, T. & Khodzitsky, M. Effective medium theory for multi-component materials based on iterative method. Photonics 7, 1–8 (2020).

%---- INPUTS ------

%--> eps1 is the Dielectric function (epsilon) of consituent 1
%--> eps2 is the Dielectric function (epsilon) of consituent 2
%--> c is the VOLUME FRACTION of the SECOND constituent. It is a value
%between 0 and 1.
%--> Mode is a vector which depends on the type of inclusions. In Bruggeman
%formula it determines the shape factor L (the depolarization factor). In
%Rouseel formulation only two types of inclusions are possible

% 1) Homogeneous mixture of spheres --> input 'HME'
% 2) Spherical inclusions in a matrix --> input 'SIM'

%--- OUTPUT ---

% ema --> The effective dielectric function of the mixture

function [ema]=Bruggeman_Rous(eps1,eps2,c,mode)

% Output and inputs are in term of dielectric functions

n1=sqrt(eps1);
n2=sqrt(eps2);

%-- If the mixture is a Homogeneous Mixture of Spheres (HME)
if contains(mode,'HME')
p=n1./n2;
b=(1/4)*((3*c-1)*((1./p)-p)+p);
z=b+sqrt(b.^2+0.5);
ema=z.*n1.*n2;

%--- If the mixture is Spherical Inclusions in a MAtrix (SIM)
elseif contains(mode,'SIM')
    e1=n1^2;
    e2=n2^2;
    p=e1/e2;
    A=(1-c)*(1-p)/(p^(1/3));
    B=((1/2)*(1+(1+(4/27)*A^3)^(1/2)))^(1/3);
    z=A*(B-A/(3*B));
    ema=(1-z)*e2;
end
