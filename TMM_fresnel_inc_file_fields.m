%---- Juan P. Martinez, reference as 1. Martinez, J. P. Light propagation in multilayered nanostructures. (2024) doi:10.13140/RG.2.2.30332.96640.

%----- This code implements the Transfer Matrix Method  used to calculate propagation in stratified media (a stack of layers) in the Formulation
%of Azzam and Bashara. 1. Azzam, R. M. & Bashara, N. M. Ellipsometry and
%Polarized Light. vol. 1 (North-Holland, 1977). 
% It does the same as TMM_fresnel_inc_file but the output is in terms of
% fields instead of intensities.  It is useful to work with polarization
% measurements or calculate the Stokes parameters.

% --> ### IMPORTANT ### output of fields only works well if all layers are coherent, because to use incoherent layers fields are converted to intensity (field information is lost!)

%-- It requires fresnel_interface.m and TMM_fresnel.m, available in the repository

%-- It includes the posibility to work with both coherent and incoherent
%layers (where light propagated incoherently, i.e. it is a straight
%addition of intensities instead of fields, and there is no interference).

%-----The idea to integrate incoherent layers appears in:

% 1. Centurioni, E. Generalized matrix method for calculation of internal light energy flux in mixed coherent and incoherent multilayers. Applied Optics 44, 7532 (2005).
% 2. Katsidis, C. C. & Siapkas, D. I. General transfer-matrix method for optical multilayer systems with coherent, partially coherent, and incoherent interference. Applied Optics 41, 3978 (2002).

%--- It also includes the option to use a file for the refractive index of
%some of the materials, which has the spectral dependence of the complex
%refractive index

%---- INPUTS---------------

% --> n, is a cell array (in matlab {}) with the refractive index (possibly
% complex) of each layer in order of approach of light rays. Notice that the first and last elements of the vector are the incoming and outgoing media If light is
% incident in a stack of 4 layers, which have air (n=1) on one side and water on another (n=1.3) n={1,n1,n2,n3,n4,1.3}. If n2 is a file it would be n={1,n1,'n2.txt',n3,n4,1.3}

% The elements of the cell array can be either straight complex numbers or the name of a file containing the spectral dependence of the complex
% refractive index. The file must have the first column the wavelength, the second the real part of the refractive index and the third column the
% extinction coefficient k. Do not include headers in the file. The wavelength under study must be the exact wavelength in the file, this code does NOT interpolate. 
% The file must be accesible in the path. If there is no absorption (k=0), the value still has to be in the file (third column will be all 0)

% CONVENTION: NEGATIVE extinction coefficient k for absorption. complex refractive index is n_c=n-ik

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

%--> rs,ts,rp,tp reflection (r) and transmission (t) of the stratified
%media, for the s (perpendicular to plane of incidence, also called TE) polarization and the p (parallel to
%the plane of incidence, also called TM) polarization.


function [rs,ts,rp,tp]=TMM_fresnel_inc_file_fields(n,e,phi0,lambda,coh)%n cell cell array, coh=1 if incoherent

%---- 1st take the refractive index from the files, for the layers were a file is given
%file is included in the cell array
num_index=zeros(1,length(n));
%---- Check if all the elements are numbers or not
for k=1:length(n)
    num_index(k)=isnumeric(n{k});
end
%--- Next do it only if there is at least one file name in n (take values from file)
if any(num_index==0)
    index_string=find(num_index==0);% aca encuentro en que lugar de n estan los strings
    for k=1:length(index_string)

        %{
        filename=n{index_string(k)};
        fid=fopen(filename);
        data=textscan(fid,'%f%f%f','Headerlines',0);
        fclose(fid);
        n(index_string(k))={data{1,2}(data{1,1}==lambda)-1i*data{1,3}(data{1,1}==lambda)};
        
        %}

       data=readmatrix(n{index_string(k)});
       n(index_string(k))={data(data(:,1)==lambda,2)-1i*data(data(:,1)==lambda,3)};
         %-- Substitutes in n, in the places where there was a file, a
       %complex number taken from the file, at the wavelength of light
    end
end
    n=cell2mat(n); %Convert to a regular matrix to continue the calculations
    
    %----- Now starts TMM like in other codes without the file -------------------
    
    
      %----Find place of INCOHERENT layers----------
indices=find(coh==1);

%-----Fresnel Coefficients-----------------
[rs,ts,rp,tp,phi]=fresnel_interface(n,phi0);

%--- First case, all layers are coherent, regular TMM ----

if all(indices==0)
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
   
    Mp=Mp*Dp*Ps;
    
end
%----Now the last interface
%-----For s-polarized--------
Dsf=(1/ts(end))*[1,rs(end);rs(end),1];
Ms=Ms*Dsf;

%------For p-polarized-------
Dpf=(1/tp(end))*[1,rp(end);rp(end),1];
Mp=Mp*Dpf;

%------Now extract r and t----------
% t=1/m11 and r=m21/m11 % expression of coefficients, see references above

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
    
   %---- Set up the dynamic, D (interface) and Propagation (P) (also called layer) matrices------
   
    Di=(1/abs(ts(1)))^2*[1,abs(rs(1))^2;abs(rs(1))^2,1];%s-polarized
    Dp=(1/abs(tp(1)))^2*[1,abs(rp(1))^2;abs(rp(1))^2,1];%p-polarized
    path=(2*pi/lambda)*(n(2))*(e)*(cos(phi(2)));
    phasedif=-1i*path;
    P=[abs(1/(exp(phasedif)))^2,0;0,abs(exp(phasedif))^2];
    Dse=(1/abs(ts(2)))^2*[1,abs(rs(2))^2;abs(rs(2))^2,1];%s-polarized end matrix, last
    Dpe=(1/abs(tp(2)))^2*[1,abs(rp(2))^2;abs(rp(2))^2,1];%p-polarized end matrix, last
    Msi=Msi*Di*P*Dse;
    Mpi=Mpi*Dp*P*Dpe;
    
     %----- Extract values ------
       
         %----- Notice that transmittance and reflectance are given directly
       %by the elements of the matrix because for we converted the matrix
       %of the incoherent layer to an INTENSITY matrix

       
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
       
       %-----
       rs=rsi;
       ts=tsi;
       rp=rpi;
       tp=tpi;
     

else
    %---- Third case, not any of the above

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
    %-------- Last Stack of COHERENT layers --------
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

       rs=rsi;
       ts=tsi;
       rp=rpi;
       tp=tpi;
end
end

    
    
    