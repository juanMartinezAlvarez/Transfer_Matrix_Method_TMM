%---- Juan P. Martinez, reference as 1. Martinez, J. P. Light propagation in multilayered nanostructures. (2024) doi:10.13140/RG.2.2.30332.96640.

% This function calculates reflection and transmission coefficients allowing for an input vector of either angles or wavelengths. Thus, it
% calculates the angular or spectral dependence of the fields

%--- See the documentation in TMM_fresnel_inc_file_fields for details. The only difference in this function is that the inputs:

%--> lambda: wavelength of light in nm. Now this can be a vector of
%wavelengths. In this case the spectral dependence of n has to be given as
%a file, and all the wavelengths to be calculated must be on that file.

%--> phi0: angle of incidence in degrees. Now it can be a vector of angles of incidence, which results in angular dependence of the coefficients

%---> ### IMPORTANT ### ONLY ONE can be a vector!

%---- OUTPUTS -----

%--> rs,ts,rp,tp reflection (r) and transmission (t) of the stratified
%media, for the s (perpendicular to plane of incidence, also called TE) polarization and the p (parallel to
%the plane of incidence, also called TM) polarization.

%-- ### Each output is now a COLUMN VECTOR, which each row corresponding to
%the wavelength or the angle (in order). So if lambda is [600,700,800], rs
%will give a column vector [rs@600;rs@700:rs@800]


function [rs,ts,rp,tp]=TMM_fresnel_inc_file_fields_vec(n,e,phi0,lambda,coh)%n es cell array coh=1 si es incoherente

%----If lambda is a vector----
if length(phi0)==1
    rs=zeros(length(lambda),1);
    ts=zeros(length(lambda),1);
    rp=zeros(length(lambda),1);
    tp=zeros(length(lambda),1);
    
    for k=1:length(lambda)
        [rs(k),ts(k),rp(k),tp(k)]=TMM_fresnel_inc_file_fields(n,e,phi0,lambda(k),coh);
    end
    
    %----If phi0 is a vector vector
    
elseif length(lambda)==1
    rs=zeros(length(phi0),1);
    ts=zeros(length(phi0),1);
    rp=zeros(length(phi0),1);
    tp=zeros(length(phi0),1);
    
    for k=1:length(phi0)
        [rs(k,1),ts(k,1),rp(k,1),tp(k,1)]=TMM_fresnel_inc_file_fields(n,e,phi0(k),lambda,coh);
    end
end
end