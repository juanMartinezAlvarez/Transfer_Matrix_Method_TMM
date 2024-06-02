%---- Juan P. Martinez, reference as 1. Martinez, J. P. Light propagation in multilayered nanostructures. (2024) doi:10.13140/RG.2.2.30332.96640.

% This function calculates the Effective Medium Approximation (EMA) using
% Bruggeman formulation (the effective medium is itself the host). It works for a mixture of TWO MATERIALS The work
% by Rouseel et.al is used to obtain the correct root of the second degree
% polynomial. If more than two materials are present, the mixture can be
% performed recursively (i.e mix 2 first and the resulting EMA with the
% third, etc). See Nazarov et.al.

%-- 2. Rouseel, Ph. J., Vanhellemont, J. & Maes, H. E. Numerical aspects of the implementation of effective-medium approximation models in spectroscopic ellipsometry regression software. Thin Solid Films 234, 423–427 (1993).
%-- 3. Nazarov, R., Zhang, T. & Khodzitsky, M. Effective medium theory for multi-component materials based on iterative method. Photonics 7, 1–8 (2020).

%--> this code requires Bruggeman_Rous in the path.

% This code adds the possibility to calculate the spectral dependence of
% the EMA, by allowing for the dielectric function of the constituents to
% be vectors.

%---- INPUTS ------

%--> eps1 is the Dielectric function (epsilon) of consituent 1
%--> eps2 is the Dielectric function (epsilon) of consituent 2

% ## IMPORTANT ## eps1 and eps2 can be vectors, which allows to mix many
% materials at the same time (usually the spectral dependence of the
% dielectric function). The code mixes first value of each together, then
% second, etc. They must be in SAME ORDER and have SAME DIMENSIONS

%--> c is the VOLUME FRACTION of the SECOND constituent. It is a value
%between 0 and 1.
%--> Mode is a vector which depends on the type of inclusions. In Bruggeman
%formula it determines the shape factor L (the depolarization factor). In
%Rouseel formulation only two types of inclusions are possible

% 1) Homogeneous mixture of spheres --> input 'HME'
% 2) Spherical inclusions in a matrix --> input 'SIM'

%--- OUTPUT ---

% ema --> The effective dielectric function of the mixture. If eps1 and
% eps2 are vectors, ema is a vector of the same dimensions

function [ema]=Bruggeman_Rous_vec(eps1,eps2,c,mode)

for h=1:length(eps1)
    ema(h)=Bruggeman_Rous(eps1(h),eps2(h),c,mode);

end

%---Vectors eps1 and eps2 must be of same dimensions
