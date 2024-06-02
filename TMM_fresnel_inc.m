%---- Juan P. Martinez, reference as 1. Martinez, J. P. Light propagation in multilayered nanostructures. (2024) doi:10.13140/RG.2.2.30332.96640.

%----- This code implements the Transfer Matrix Method  used to calculate propagation in stratified media (a stack of layers) in the Formulation
%of Azzam and Bashara. 1. Azzam, R. M. & Bashara, N. M. Ellipsometry and Polarized Light. vol. 1 (North-Holland, 1977).

%-- It requires fresnel_interface.m and TMM_fresnel.m, available in the repository

%-- It includes the posibility to work with both coherent and incoherent
%layers (where light propagated incoherently, i.e. it is a straight
%addition of intensities instead of fields, and there is no interference).

%-----The idea to integrate incoherent layers appears in:

% 1. Centurioni, E. Generalized matrix method for calculation of internal light energy flux in mixed coherent and incoherent multilayers. Applied Optics 44, 7532 (2005).
% 2. Katsidis, C. C. & Siapkas, D. I. General transfer-matrix method for optical multilayer systems with coherent, partially coherent, and incoherent interference. Applied Optics 41, 3978 (2002).

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

% --> Lambda is the wavelength of light in nanometers

% --> coh is a vector that indexes which layers are coherent (index is 0) and which ones are incoherent (index is 1). 
% So if 3 layers are present, the middle is incoherent coh=[0,1,0]. Notice that length of coh is the same as vector e 

%----- OUTPUTS ------

% capital leters are quotient of intensities (R,T), lower case are of
% fields (r,t)

%--> Rs,Ts,Rp,Tp Reflectance (R) and Transmittance (T) of the stratified
%media, for the s (perpendicular to plane of incidence, also called TE) polarization and the p (parallel to
%the plane of incidence, also called TM) polarization.

%--> Msi and Mpi : global intensity transfer matrix of the stratified media for p and s
%polarizations. Intensity matrix as output (elements square of the regular
%matrix, see papers above)


function [Rs,Ts,Rp,Tp,Msi,Mpi]=TMM_fresnel_inc(n,e,phi0,lambda,coh)%e y lambda en mismas unidades
%ko=2*pi/lambda;%calculo ko

%incoherent if coh=1

%----Find position of incoherent layers----------
indices=find(coh==1);

%-----Fresnel coefficients-----------------
[rs,ts,rp,tp,phi]=fresnel_interface(n,phi0);

%-- First case: all layers are coherent, and regular application of the TMM---

if all(indices==0) %
Ms=eye(2);
Mp=eye(2);
   
for k=1:length(e)
    %------For s-polarized-----
    Ds=(1/ts(k))*[1,rs(k);rs(k),1];
    path=(2*pi/lambda)*(n(k+1))*(e(k)*cos(phi(k+1)));
    phasedif=-1i*path;
    Ps=[1/(exp(phasedif)),0;0,exp(phasedif)];
    Ms=Ms*Ds*Ps;
    
    %-------For p-polarized------
    Dp=(1/tp(k))*[1,rp(k);rp(k),1];
    %----Aca agregaria otro path si tiene otro indice de refraccion
    Mp=Mp*Dp*Ps;
    
end
%----Last Interface
%-----For s-polarized--------
Dsf=(1/ts(end))*[1,rs(end);rs(end),1];
Ms=Ms*Dsf;

%------For p-polarized-------
Dpf=(1/tp(end))*[1,rp(end);rp(end),1];
Mp=Mp*Dpf;

%------now obtain r and t----------
% t=1/m11 y r=m21/m11

%---- s-polarized -----
ts=1/Ms(1,1);
rs=Ms(2,1)/Ms(1,1);

Rs=(abs(rs))^2;
Ts=(real(cos(phi(end)))/real(cos(deg2rad(phi0))))*(abs(ts))^2;

%---- p-polarized -----
tp=1/Mp(1,1);
rp=Mp(2,1)/Mp(1,1);

Rp=(abs(rp))^2;
Tp=(real(cos(phi(end)))/real(cos(deg2rad(phi0))))*(abs(tp))^2;

Msi=Ms;
Mpi=Mp;

%---- Second case: single incoherent layer. size(e)=1 incoh
elseif size(e)==1
    Msi=eye(2);
    Mpi=eye(2);
    
    %---- The matrices ------
    %----Extract coherent matrices----
    [~,~,~,~,Ms1,Mp1]=TMM_fresnel(n,e,phi0,lambda);

    %----Now calculate intensity matrix for incoherent layer from elements of coherent matrix. Calculations are in the papers of Centurioni or Katsidis (see above)

    %----- s polarized ---------
    rs1=Ms1(2,1)/Ms1(1,1);
    ts1=1/Ms1(1,1);
    rs1b=-Ms1(1,2)/Ms1(1,1);
    ts1b=det(Ms1)/Ms1(1,1);

    Msi=(1/(abs(ts1)^2))*[1,-(abs(rs1b))^2;abs(rs1)^2,abs(ts1*ts1b)^2-abs(rs1*rs1b)^2];

    %-------- p-polarized ---------------
    
    rp1=Mp1(2,1)/Mp1(1,1);
    tp1=1/Mp1(1,1);
    rp1b=-Mp1(1,2)/Mp1(1,1);
    tp1b=det(Mp1)/Mp1(1,1);
    
    Mpi=(1/(abs(tp1)^2))*[1,-(abs(rp1b))^2;abs(rp1)^2,abs(tp1*tp1b)^2-abs(rp1*rp1b)^2];
    %------

   %{
    Di=(1/abs(ts(1)))^2*[1,abs(rs(1))^2;abs(rs(1))^2,1];%s
    Dp=(1/abs(tp(1)))^2*[1,abs(rp(1))^2;abs(rp(1))^2,1];%p
    path=(2*pi/lambda)*(n(2))*(e)*(cos(phi(2)));
    phasedif=-1i*path;
    P=[abs(1/(exp(phasedif)))^2,0;0,abs(exp(phasedif))^2];
    Dse=(1/abs(ts(2)))^2*[1,abs(rs(2))^2;abs(rs(2))^2,1];%s end
    Dpe=(1/abs(tp(2)))^2*[1,abs(rp(2))^2;abs(rp(2))^2,1];%p end
    Msi=Msi*Di*P*Dse;
    Mpi=Mpi*Dp*P*Dpe;
   %}

    %----- Extract transmittance and reflectance coefficients ------
       
       %------ s-polarized -------
       rsi=Msi(2,1)/Msi(1,1);
       tsi=1/Msi(1,1);

       %----- Notice that transmittance and reflectance are given directly
       %by the elements of the matrix because for we converted the matrix
       %of the incoherent layer to an INTENSITY matrix

       Rs=rsi;
       Ts=tsi*real(n(end)*cos(phi(end))/(n(1)*cos(phi(1))));
       
       %----- p-polarized -----
       rpi=Mpi(2,1)/Mpi(1,1);
       tpi=1/Mpi(1,1);
       %-----
       Rp=rpi;
       Tp=tpi*real(n(end)*cos(phi(end))/(n(1)*cos(phi(1))));

%--- Case 3: other cases --------
else
    %------ First step: Calculate intensity global transfer matrix for the
    %stacks (groups) of coherent layers---
    
    %-------First stack of coherent layers-------
    [~,~,~,~,Ms1,Mp1]=TMM_fresnel(n(1:indices(1)+1),e(1:indices(1)-1),phi0,lambda);

   
    
    %------Extract r, t, r' y t', which are the Fresnel coefficients in both directions, see papers above -------
    %r=T21/T11 , t=1/T11  , r'= -T21/T11  t'= Det(T)/T11 % Expression for the coefficients
    
    %------- s-polarized ------------------
    rs1=Ms1(2,1)/Ms1(1,1);
    ts1=1/Ms1(1,1);
    rs1b=-Ms1(1,2)/Ms1(1,1);
    ts1b=det(Ms1)/Ms1(1,1);
    
    Tscoh1=(1/(abs(ts1)^2))*[1,-(abs(rs1b))^2;abs(rs1)^2,abs(ts1*ts1b)^2-abs(rs1*rs1b)^2];
    
    %--- The matrices are saved as a cell array. Each element of Ts and Tp
    %is a matrix
    
    Ts(1)={Tscoh1};
    %-------- p-polarized ---------------
    
    rp1=Mp1(2,1)/Mp1(1,1);
    tp1=1/Mp1(1,1);
    rp1b=-Mp1(1,2)/Mp1(1,1);
    tp1b=det(Mp1)/Mp1(1,1);
    
    Tpcoh1=(1/(abs(tp1)^2))*[1,-(abs(rp1b))^2;abs(rp1)^2,abs(tp1*tp1b)^2-abs(rp1*rp1b)^2];
    
   
    
    Tp(1)={Tpcoh1};

    %-------- Last Stack of Coherent layers--------
    [~,~,~,~,Msend,Mpend]=TMM_fresnel(n(indices(end)+1:end),e(indices(end)+1:end),rad2deg(phi(indices(end)+1)),lambda);
    
     %------- s-polarized ------------------
    rsend=Msend(2,1)/Msend(1,1);
    tsend=1/Msend(1,1);
    rsendb=-Msend(1,2)/Msend(1,1);
    tsendb=det(Msend)/Msend(1,1);
    
    Tscohend=(1/(abs(tsend)^2))*[1,-(abs(rsendb))^2;abs(rsend)^2,abs(tsend*tsendb)^2-abs(rsend*rsendb)^2];
     
    
    
    Ts(length(indices)+1)={Tscohend};
    
    %-------- p-polarized ---------------
    rpend=Mpend(2,1)/Mpend(1,1);
    tpend=1/Mpend(1,1);
    rpendb=-Mpend(1,2)/Mpend(1,1);
    tpendb=det(Mpend)/Mpend(1,1);
    
    Tpcohend=(1/(abs(tpend)^2))*[1,-(abs(rpendb))^2;abs(rpend)^2,abs(tpend*tpendb)^2-abs(rpend*rpendb)^2];
    
    
    Tp(length(indices)+1)={Tpcohend};
    
    %--- Now the rest of the stacks of incoherent layers : -----------

    for k=1:length(indices)-1
        [~,~,~,~,Ms,Mp]=TMM_fresnel(n(indices(k)+1:indices(k+1)+1),e(indices(k)+1:indices(k+1)-1),rad2deg(phi(indices(k)+1)),lambda);
        
        %------- s-polarized -----------
       rs=Ms(2,1)/Ms(1,1);
       ts=1/Ms(1,1);
       rsb=-Ms(1,2)/Ms(1,1);
       tsb=det(Ms)/Ms(1,1);
    
        Ts(k+1)={(1/(abs(ts)^2))*[1,-(abs(rsb))^2;abs(rs)^2,abs(ts*tsb)^2-abs(rs*rsb)^2]};
        
       
       %--------- p-polarized -----------
       
       rp=Mp(2,1)/Mp(1,1);
       tp=1/Mp(1,1);
       rpb=-Mp(1,2)/Mp(1,1);
       tpb=det(Mp)/Mp(1,1);
    
        Tp(k+1)={(1/(abs(tp)^2))*[1,-(abs(rpb))^2;abs(rp)^2,abs(tp*tpb)^2-abs(rp*rpb)^2]};
       
      
    end
    
    %------Now the propagation matrices in the incoherent layers-----
       for q=1:length(indices)
            path=(2*pi/lambda)*(n(indices(q)+1))*(e(indices(q))*cos(phi(indices(q)+1)));
            phasedif=-1i*path;
            P(q)={[(abs(1/(exp(phasedif))))^2,0;0,(abs(exp(phasedif)))^2]};
       end
       
       %-----### Multiplication of all the matrices of the whole stack ###-----
       
       %--- Inicialization with the first matrix (between the entrance and the first incoherent slab) 
         Msi=Ts{1};
         Mpi=Tp{1};
         
       %----- Now multiply each matrix in order --------------
       
       for k=1:length(indices)
           Msi=Msi*P{k}*Ts{k+1};
           Mpi=Mpi*P{k}*Tp{k+1};
       end
       
       %----- Extract values ------
       
       %------ s-polarized -------
       rsi=Msi(2,1)/Msi(1,1);
       tsi=1/Msi(1,1);
       %-----
       Rs=rsi;
       Ts=tsi*real(n(end)*cos(phi(end))/(n(1)*cos(phi(1))));
       
       %----- p-polarized -----
       rpi=Mpi(2,1)/Mpi(1,1);
       tpi=1/Mpi(1,1);
       %-----
       Rp=rpi;
       Tp=tpi*real(n(end)*cos(phi(end))/(n(1)*cos(phi(1))));
end
