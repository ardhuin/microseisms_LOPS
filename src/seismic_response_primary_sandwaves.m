function [F_delta,Fx_delta,Fy_delta,pos,freqp,df,times]=seismic_response_primary_sandwaves(wave_spectrum, ...
    bottom_topography_spectrum, ...
    dpt, lat,lon, CgR,  Q, date1, date2, dt, lono, lato, dA, statname, ndcalc)


     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This functions computes synthetic seismic spectra due to primary mechanism
% WITH SOME ASSUMPTIONS ON THE DIRECTIONAL SPECTRUM... 
% at the point of longitude lono and latitude lato
%
% Last modified: March 1, 2017 
%
% History     : created March 1, 2017 by F. Ardhuin
%             : from seismic_response_primary_from_spec.m 
%
% Input parameters: 
% - pathNC    : file path where the *_p1.nc NetCDF files are stored
% - dpt       : array of water depths, should have the same size as the
%               wave data
% - Q         : quality factor (function of frequency)
% - date1     : start date in matlab form =datenum(iy,im,id,ih,imin,isec)
% - date2     : end date in matlab form   =datenum(iy,im,id,ih,imin,isec)
% - dt        : time step of computation, in fraction of day
% - lono      : longitude of receiving station, where spectrum is estimated
% - lato      : latitude of receiving station, where spectrum is estimated
% - deep      : reference depth
% - statname  : name of seismic station for which the computation is done
% - ndcalc    : number of directions used for subgrid islands 
%
% Output parameters: 
% - F_delta   : Spectral density of surface ground displacement in m^2/Hz
%               normalized so that sum(Ef*df) is the ground displacement
%               variance
% - pos       : centroid of sources 
% - freqp     : seismic frequencies
% - df        : frequency increment
% - times     : matlab dates at which the spectrum is computed
%
% Method:
% Uses Hasselmann (1963) as corrected by Ardhuin et al. (2015): 
%   noise source for a straight coastline with 
%   constant bottom slope (defined by Hp). 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1. Defines constants and coefficients
%
g    =9.81;     % gravity          in m/s^2
rhow =1026;     % density of water in kg/m^3
rhos =2600;     % density of rock  in kg/m^3
betas=2800;     % shear wave velocity of rock in m/s
R_E=4E7/(2*pi); % Earth radius     in m
pos=0.;
lg10=log(10);


fid=fopen(bottom_topography_spectrum,'r');
%         OPEN(183,FILE= 'bottomspectrum.inp', status='old')
%         READ(183,*) nkbx, nkby
%         READ(183,*) dkbx, dkby
%         WRITE(*,*) 'Bottom spec. dim.:', nkbx, nkby, dkbx, dkby
%         ALLOCATE(BOTSPEC(nkbx, nkby))
%         DO I=1, nkbx
%            READ(183,*) BOTSPEC(I,:)
%            END DO
%         CLOSE(183)
%         variance=0
%         DO i=1,nkbx
%            DO j=1,nkby
%               variance=variance+BOTSPEC(i,j)*dkbx*dkby
%               END DO
%            END DO
%         WRITE(*,*) 'Bottom variance:', variance

A=fscanf(fid,'%f');
fclose(fid);
nkbx=A(1);
nkby=A(2);
dkbx=A(3);
dkby=A(4);
botspec=reshape(A(5:end),nkby,nkbx)';
% Computes bottom wavenumbers   
   kbxmax=((nkbx-mod(nkbx,2))/2.).*dkbx;
   kbotx=linspace(-kbxmax,kbxmax,nkbx);
   kbymax=((nkby-mod(nkby,2))/2.).*dkby;
   kboty=linspace(-kbymax-(1-mod(nkby,2))*dkby,kbymax,nkby);

%figure(100)   
%pcolor(kbotx-dkbx/2,kboty-dkby/2,10.*log10(botspec'));shading flat;colorbar;


%
%    Get pressure to seismic coefficients (C^2 in Longuet-Higgins 1950)
%    Be careful that it may not be consistent with betas and other rheology
%    factors
[C2,nC2,c]=seismic_define_C2(0,CgR);
%
%
% NB: the F_p1(K,f) is computed from E(f,theta) in the NetCDF files 
%     spectral density of surface pressure (in meters)
%     such that (rhow*g)^2*F_p1(K,f)*df*dKx*dKy is
%     the variance of the pressure (in Pa^2) with df the 
%     ocean wave frequency interval.
%
%
% factor1 is the constant part of the integrand, 
% factor1=W^2*2*pi^2/(c^2*sigma)*sin(alpha)
%
factor1=2.*pi*(1/rhos)^2/(betas^5*R_E); 
%
%
% This source is attenuated during propagation 
% 

%
% 2. Reads netcdf file to get dimensions ... .
%
nt=floor((date2-date1)/dt)+1;     % number of time steps
times=linspace(date1,date2,nt);


%
% Reads NetCDF *spec.nc files or *p1.nc file  
%
  fid=netcdf.open(wave_spectrum,'NC_NOWRITE');
%
% Reads dimensions from wave spectrum NetCDF file
%
dimp = netcdf.inqDimID(fid,'station');
dimtime = netcdf.inqDimID(fid,'time');
dimf = netcdf.inqDimID(fid,'frequency');
[d0,np]=netcdf.inqDim(fid,dimp);
[d1,nf]=netcdf.inqDim(fid,dimf);
[d3,nts]=netcdf.inqDim(fid,dimtime);
%
% Reads variables from NetCDF file
%
vartime= netcdf.inqVarID(fid,'time');
varlon = netcdf.inqVarID(fid,'longitude');
varlat = netcdf.inqVarID(fid,'latitude');
varf   = netcdf.inqVarID(fid,'frequency');
time0=datenum(1990,1,1);
time=netcdf.getVar(fid,vartime)+time0;
freqp=netcdf.getVar(fid,varf);

xfr=exp(log(freqp(nf)/freqp(1))/(nf-1));  % determines the xfr geometric progression factor
df=freqp*0.5*(xfr-1/xfr);% frequency intervals in wave model 
   dimd = netcdf.inqDimID(fid,'direction');
   [d2,nd]=netcdf.inqDim(fid,dimd);
   varid   = netcdf.inqVarID(fid,'efth');
   vard   = netcdf.inqVarID(fid,'direction');
   vardpt   = netcdf.inqVarID(fid,'dpt');
   theta=netcdf.getVar(fid,vard);
   lonp2=netcdf.getVar(fid,varlon);
   latp2=netcdf.getVar(fid,varlat);
   dptpall=netcdf.getVar(fid,vardpt);
   lonp=lonp2(:,1);
   latp=latp2(:,1);


nx=length(lon);
ny=length(lat);
la2=repmat(lat',nx,1);
lo2=repmat(lon,1,ny);


nf=length(freqp);
omega=(2*pi).*freqp;     % seismic radian frequency omega=2*pi/T




%
% 4. Computes distance from station
%
   alpha=dist_sphere(lono,lonp,lato,latp).*(pi/180);  % spherical distance along shortest arc

%
% 5. Computes Efth spectrum to pressure transformation  
%

   Eftop1f=zeros(nd,nf);
   botspeci=zeros(nd,nf);

   dtor = pi/180;
   
   k=dispNewtonTH(freqp,dpt);
    
   for i=1:nf
       for j=1:nd
          kbx=k(i).*sin(theta(j).*dtor);
          kby=k(i).*cos(theta(j).*dtor);
          kbotxi=1+(kbxmax+kbx)/dkbx; % wavenumber in grid units 
          kbotyi=1+(kbymax+kby)/dkby;
          ibk=max(min(floor(kbotxi),nkbx-1),1);
          xbk=mod(kbotxi,1.0);
          jbk=max(min(floor(kbotyi),nkby-1),1);        
          ybk=mod(kbotyi,1.0);
%  bilinear interpolation of bottom spectrum onto wave wavenumber 
           botspeci(j,i)=(botspec(ibk,jbk)*(1-ybk)+botspec(ibk,jbk+1)*ybk)*(1-xbk) ...
                +(botspec(ibk+1,jbk)*(1-ybk)+botspec(ibk+1,jbk+1)*ybk)*xbk;   
          %fprintf('botspec: %d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f, \n',i,j,ibk,jbk,freqp(i),k(i),dkbx,dkby,kbx,kby,botspeci(j,i))
             
       end
   end



%
% 6. Computation of seismic sources 
%


%
% Loop on time steps from 1 to nt
%
efthall =netcdf.getVar(fid,varid);
dth=2*pi./nd;
       
for it=1:nt
   itf= 1+round((times(it)-times(1))/dt);
   efth_file=squeeze(efthall(:,:,1,itf));
   efth=efth_file;
   alphas=zeros(nf,1);
   for i=1:nf
       %if (i < nf-4)
       %    efth(:,i)=efth_file(:,i+4).*2;
       %end
       om=omega(i);k0=k(i);h0=dpt; %pall(1,itf);
        
  
       kh=k(i)*dpt; % should update later with tidal change of water depth
       if (kh < 7)
           dkdD=-2*k0^2./(2*kh+ sinh(2*kh));
           dkDdD=k0.*sinh(2*kh)./(2*kh+ sinh(2*kh));
           Cg=om./k0.*0.5*(1+(2*kh)/sinh(2*kh));
           dCgdD=om/k(i)*(dkDdD/sinh(2*kh)-2*kh*dkDdD*cosh(2*kh)/(sinh(2*kh)^2))-dkdD*om/k0.^2*0.5*(1+(2*kh)/sinh(2*kh));
           alpha1=-(k0+h0*dkdD).*tanh(k0.*h0)./k0;
           alpha2=-0.5.*dCgdD./(Cg.*k0);
           alpha3=-dkdD./(k0.^2);
           alphas(i) = alpha1+alpha2+alpha3;
       else
           alphas(i) = -1;
       end
       Eftop1f(:,i)=botspeci(:,i).*(rhow.*g.*k(i).*alphas(i)/cosh(kh)).^2.*dth;
       %fprintf('Eftop1fg: %d, %f, %f, %f, %g, %g, %g, %g, %g,\n',i,freqp(i),om.^2,9.81*k0*tanh(kh),(rhow.*g.*k(i).*alphas(i)/cosh(kh)).^2.*dth,k(i),alphas(i),sum(botspeci(:,i)),sum( Eftop1f(:,i)))
       if it==1 & i==10
           botspeci(:,i)
           botspeci(:,i).*(k(i)./cosh(kh)).^2.*dth
       end
   end
%
% Computes amplification and attenuation for the first time step
%
    if it==1
        nf=length(freqp);
        coeff=zeros(nf,1);
        attenuation=zeros(nf,1);
        coeff_love=zeros(nf,1);
        attenuation_love=zeros(nf,1);
        
        for j=1:nf
            omehoverbeta=0.*omega(j).*dpt./betas;     % set to zero because interaction in shallow water
            i=1+int16(100.*omehoverbeta);                % gets the nearest discretized omega*h/beta
            I=find(i>nC2);
            i(I)=nC2;
            coeff(j)=(factor1.*omega(j).*C2(i))./sin(alpha);
  
            coeff_love(j)=coeff(j); % to be corrected ... 
%
% correction factor for multiple orbits
%
            b=exp(-1.*omega(j).*(2.*pi).*(R_E./(abs(CgR)*Q(j))));  % attenuation for one orbit
%
%
% terms for shorter arc
%
            attenuation(j)=exp(-1.*omega(j).*alpha.*(R_E./(abs(CgR)*Q(j))))./(1-b);
%
% terms for longer arc
%
            attenuation(j)=attenuation(j)+exp(-1.*omega(j).*(2.*pi-alpha).*(R_E./(abs(CgR)*Q(j))))./(1-b);
            attenuation_love(j)=attenuation(j); % to be corrected ... 
            %fprintf('coeff: %d, %d, %f, %f, %f, %f, C2:%f, %f \n',it,j,alpha,alpha*R_E,coeff(j)*1E30,attenuation(j),alphas(j),C2(i))
               
            I2=find(isfinite(c)==0 | C2 == 0);
            c(I2)=0;
            coeff_all(j)=coeff(j).*attenuation(j).*dA;
            coeff_all_love(j)=coeff_love(j).*attenuation_love(j).*dA;
            %fprintf('coeff: %d, %g, %f, %f, %g,\n',j,coeff(j),attenuation(j),dA,coeff_all(j))
        end
    end 
    
   if (mod(it,50) ==1 ) 
       msg=sprintf(': %d out of %d ', it, nt);
       fprintf('computing for time step %s\n',msg); 
  
   end
   for j=1:nf
       p1f(j)= sum( Eftop1f(:,j).*efth(:,j));  % this is the spectrum at K =0
       %fprintf('I... %d, %d, %g,  %g, %g, %g\n',it,j,p1f(j),sum( Eftop1f(:,j)),sum(efth(:,j)),coeff_all(j))
       t1yf(j)= sum( Eftop1f(:,j).*efth(:,j).*(cos(theta.*dtor).^2)./(alphas(j).^2));  % this is the spectrum at K =0
       t1xf(j)= sum( Eftop1f(:,j).*efth(:,j).*(sin(theta.*dtor).^2)./(alphas(j).^2));  % this is the spectrum at K =0
       t1xyf(j)= sum( Eftop1f(:,j).*efth(:,j).*(sin(theta.*dtor).*cos(theta.*dtor))./(alphas(j).^2));  % this is the spectrum at K =0
   end;
   
   source=coeff_all(:).*p1f(:);             
   sourcex=coeff_all_love(:).*t1xf(:);             
   sourcey=coeff_all_love(:).*t1yf(:);             

%
% Only takes points that are in the ocean (source is finite) and not too
% close (singularity for alpha=0) 
%
    F_delta(:,it)= source(:);
    Fx_delta(:,it)= sourcex(:);
    Fy_delta(:,it)= sourcey(:);

    if ( it == nt  )
%
% Save table at end of month
%
        dateendfile=datevec(now);
        force=p1f.*0;
        for j=1:nf
          force(:,j)=(2*pi).*sqrt(p1f(:,j).*dA.*df(j));
        end
       %figure(101)
       %j=11;
       %set(gca,'FontSize',16);
       %set(gca,'XTick',[-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180],'LineWidth',1);
       %title('p1f from spec');
       % scatter(lo2(indb),la2(indb),20,p1f(:,j),'LineWidth',2);
        
        
        filenames=[ 'Pf_sandwave_'  sprintf('%04d',dateendfile(1))  sprintf('%02d',dateendfile(2)) ];
        filenamef=[ 'force_sandwave_'  sprintf('%04d',dateendfile(1))  sprintf('%02d',dateendfile(2)) ];
        %save(filenames, 'freqp','time', 'los', 'las', 'indb' , 'Pfo', ...
        %                'statname' ,  'coeff_all');
        
        %save('tempo','fullnam2','Pfo','freqp','las','los','time','dptp1p','dAp');
        %seismic_write%_Pfo_NC(fullnam2,Pfo,freqp,las,los,time,dptp1p,dAp);
        % 
        %force_small=single(force(:,1:25));            
        %freq_small=freqp(1:25);            
        %save(filenamef, 'times', 'freq_small',  'force_small');
    end
end
nf=length(freqp);

end

