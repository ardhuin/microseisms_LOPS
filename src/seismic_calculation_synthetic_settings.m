%
% Defines start and and date for seismic response calculation
%
iy=2023;
date1=datenum(iy,01,01,00,00,00);
date2=datenum(iy,12,31,21,00,00);
%date2=datenum(2008,02,29,21,00,00);
dt=1/8;                     % time step in days: 1/8 is the usual time resolution 
                            % of the wave model output

channel='LHZ';                            
% Choice of seismic station and Q factors 
statname='EDA'; QQ=[400];   % The array QQ can contain several values of Q
%statname='ESK'; QQ=[400];   % The array QQ can contain several values of Q
%statname='SSB'; QQ=[240];   % The array QQ can contain several values of Q
%
%statname='SNZO'; QQ=[240];   % The array QQ can contain several values of Q
%statname='IGLA'; QQ=[400];   % The array QQ can contain several values of Q
%statname='ELSH'; channel='BHZ';   QQ=[240];   % The array QQ can contain several values of Q
%statname='PPTF'; QQ=[500];   % The array QQ can contain several values of Q
%statname='MBWA'; QQ=[400 500 600];   % The array QQ can contain several values of Q
%statname='BFO'; QQ=[300];   % The array QQ can contain several values of Q
%statname='KIP'; QQ=[800];   % The array QQ can contain several values of Q
%statname='MBO'; QQ=[250];   % The array QQ can contain several values of Q
%statname='GRFO'; QQ=[300];   % The array QQ can contain several values of Q
%statname='KONO'; QQ=[300];   % The array QQ can contain several values of Q
%statname='PAB'; QQ=[300];   % The array QQ can contain several values of Q
%statname='PVAQ'; QQ=[300];   % The array QQ can contain several values of Q
%statname='BBSR'; QQ=[300];   % The array QQ can contain several values of Q
%statname='KIEV'; QQ=[300];   % The array QQ can contain several values of Q
%statname='PEL'; QQ=[400];   % The array QQ can contain several values of Q
%statname='KIP'; QQ=[800];   % The array QQ can contain several values of Q
%statname='TAM'; QQ=[200];   % The array QQ can contain several values of Q
%statname='CAN'; QQ=[400];   % The array QQ can contain several values of Q
%statname='RSSD'; QQ=[400];   % The array QQ can contain several values of Q
%statname='RER'; QQ=[400];   % The array QQ can contain several values of Q
%statname='NWAO'; QQ=[400];   % The array QQ can contain several values of Q
%statname='SCZ'; QQ=[300];   % The array QQ can contain several values of Q
%statname='BKS'; QQ=[300];   % The array QQ can contain several values of Q
%statname='KBS'; QQ=[400];   % The array QQ can contain several values of Q
%statname='COR'; QQ=[300];   % The array QQ can contain several values of Q
%statname='PAS'; QQ=[300];   % The array QQ can contain several values of Q
%statname='KIP'; QQ=[800];   % The array QQ can contain several values of Q
%statname='NNA'; QQ=[400];   % The array QQ can contain several values of Q
%statname='PPT'; QQ=[800];   % The array QQ can contain several values of Q
%statname='KDAK'; QQ=[400];   % The array QQ can contain several values of Q
CgR  =1800;                 % group speed of seismic waves, assumed constant

%bottom_topography_spectrum='spectrum_smoothPDC1_sq2.bsp';
bottom_topography_spectrum='spectrum_PDC1.bsp';
%bottom_topography_spectrum='spectrum_StMalo.bsp';
%bottom_topography_spectrum='spectrum_Ireland_rock_lr.bsp';
bottom_topography_spectrum='spectrum_Ireland_shallow_rocks.bsp';
%bottom_topography_spectrum='spectrum_UK_lr.bsp';


fid=netcdf.open(wave_spectrum,'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid]=netcdf.inq(fid); 
dimf = netcdf.inqDimID(fid,'frequency'); 
[dimd] = netcdf.inqDimID(fid,'direction'); 
[d2,ndir]=netcdf.inqDim(fid,dimd);
varf   = netcdf.inqVarID(fid,'frequency');
freqp=netcdf.getVar(fid,varf);
freq=freqp;

%
%  Choice of model grid 
%  1) global or regional zoom: will use only one grid if modgrid1b=''
%  2) Both global AND regional zoom if modgrid1b~=''
%
modgrid1b='';modgrid2b='';


%modgrid1='GLOBAL05'; reftag = '_RMAP';modgrid1b='';modgrid2b='';vtag='_TEST451_RMAP_Q';   

nmodes=4;     % Number of Rayleigh modes used in calculation 
nmodesout=0;  % Number of modes for separate output 

date=datevec(date1);
iy=date(1);
date0=datenum(iy,01,01,00,00,00);
yname=sprintf('%.4d',date(1));
im=sprintf('%.2d',date(2));

% Sets the frequency indices for comparison with data      
   ifmin=30; %36; %30; %25; %1;
   ifmax=34; %18; %35; %40; %16;
%



path1='/media/ardhuin/ondadelmar/';      %path1=['/media/VERBATIM4/SISMO/' modgrid1 '_'];
onepath='/media/ardhuin/FabLinux/GLOBAL_IG_2013/';modgrid1='GLOBAL05'; reftag = '_REF10_20_40';vtag='_TEST471_REF102040Q';  
vtag='_CCI_WW3_Ekoundou_Q'; 
wave_spectrum='CCI_WW3-GLOB-30M-28606_2023_spec.nc';
wave_spectrum='CCI_WW3-GLOB-30M-Akpo_2023_spec.nc';
wave_spectrum='CCI_WW3-GLOB-30M-Ekoundou_2023_spec.nc';

pathref1=[path1  modgrid1 '_' yname reftag '/'];
pathnoref1=[path1  modgrid1 '_' yname '_NOREF/'];
%pathnoref2=[path1  yname '_NOREF/' modgrid2b  '/'];


%onepath='36D_r1178/';
%pathnoref1=onepath;

pathref1=onepath;

%
% Reads frequency array from first file
%
varp2l='p2l';  % For "old" files this is 'Fp3D' 
filename=sprintf('ww3.%.4d%.2d_p2l.nc',date(1),mod(date(2),100))
%varef='ef'; filename=sprintf('ww3.%.4d%.2d_ef.nc',date(1),mod(date(2),100))
fullfile=[pathref1 filename]

if exist(fullfile,'file') 
   fid=netcdf.open(fullfile,'NC_NOWRITE');
   [ndims,nvars,ngatts,unlimdimid]=netcdf.inq(fid); 
   dimf = netcdf.inqDimID(fid,'f'); 
   varf   = netcdf.inqVarID(fid,'f');
   freqp=netcdf.getVar(fid,varf);
   freq=2.*freqp;
end

allspec=0

if allspec==1
%
% Test if we have the mat file for primary microseisms 
%
%filenamema=sprintf('Pf_%.4d%.2d.mat',date(1),mod(date(2),100));
filenamema=sprintf('ww3.%.4d%.2d_p1.nc',date(1),mod(date(2),100));
fullfima=[pathref1 filenamema]
resmat=exist(fullfima,'file')
resmat2=exist(fullfile,'file')

if resmat2 == 2 
[lat,lon,time,freq,P2f,unit,scale,fillv,MAPSTA] ...
         =read_WWNCf_noscale(fullfile,varp2l,date1,1);
end
resmat
if resmat < 2 
   %filename=sprintf('glob_30m.%.4d%.2d_spec.nc',date(1),mod(date(2),100))
   filename=sprintf('ww3.%.4d%.2d_spec.nc',date(1),mod(date(2),100))
   resmat3=exist(filename,'file');
   %if resmat3 > 0
   fullname=[pathref1 filename]
   fid=netcdf.open(fullname,'NC_NOWRITE');
   [ndims,nvars,ngatts,unlimdimid]=netcdf.inq(fid); 
   dimf = netcdf.inqDimID(fid,'frequency'); 
   [dimd] = netcdf.inqDimID(fid,'direction'); 
   [d2,ndir]=netcdf.inqDim(fid,dimd);
   varf   = netcdf.inqVarID(fid,'frequency');
   freqp=netcdf.getVar(fid,varf);
   freq=freqp;
else
   fid=netcdf.open(fullfima,'NC_NOWRITE');
   [ndims,nvars,ngatts,unlimdimid]=netcdf.inq(fid); 
   dimf = netcdf.inqDimID(fid,'frequency'); 
   varf   = netcdf.inqVarID(fid,'frequency');
   freqp=netcdf.getVar(fid,varf);
   ndir=36; % change this ... 
   freq=freqp;

%         load(fullfima,'freqp');
end

         
filename=sprintf('ww3.%.4d%.2d_p2l.nc',date(1),mod(date(2),100))
fullfile=[pathref1 filename]
resmat2=exist(fullfile,'file')
end

[lato,lono,thresh1,thresh2,thresh3,maxd,maxdp,r] = seismic_stations(freq,statname);
readdata=0
if readdata==1     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gets seismic data for validation
%
% here the data are supposed to be in the folder mydatapath/YEAR/STAT
% where year is '2000', ... 
% stat is a string identifying the seismic station ... 
% if, instead, you have the data in a folder 'DATAFOLDER', replace the 
% line 
%
iy1=iy;
date2v=datevec(date2);
iy2=date2v(1);

mydatapath='/home/ardhuin/DATA/';
%
date0_old=[];
for iydata=iy1:iy2;
   ynamed=sprintf('%.4d',iydata);
   if length( mydatapath) == 0
     seismicdatafile=[statname '.' channel '.' ynamed '_03h_depl.mat' ];
   else
     seismicdatafile=[mydatapath ynamed '/' statname '/' statname  '.' channel '.' ynamed '_03h_depl.mat' ];
   end
   load(seismicdatafile);
   spectre_year=spectre0;
   date0_year=date0;
   frq0_year=frq0;
   if length(date0_old) > 0 
       if (length(frq0_old)==length(frq0) & frq0_old(1)==frq0(1) )
           nold=length(date0_old);
           nnew=length(date0);
           spectre_temp=spectre0_old;
           date0_temp=date0_old;
           spectre_old=zeros(nold+nnew,length(frq0));
           date0_old=zeros(nold+nnew);
           spectre_old(1:nold,:)=spectre_temp;
           spectre_old(1:nnew,:)=spectre0;
           date0_old(1:nold)=date0_temp;
           date0_old(1:nnew)=date0;
       end
   else
      spectre_old=spectre0;
      date0_old=date0;
      frq0_old=frq0;
   end
end
spectre0=spectre_old;
date0=date0_old;
%
% filters out earthquakes
%
   [ spectre1, Er, Er2, Er3, Esum ] = sismo_filterseisms( date0,frq0,spectre0,thresh1,thresh2,thresh3 );
   spectrea=spectre1;
   [ spectre1, Er, Er2, Er3, Esum ] = sismo_filterseisms( date0,frq0,spectrea,thresh1,thresh2,thresh3 );
%
% this assumes thate date0 is a matlab date ... not true for seismic
% spectra provided by IPGP, in that case, use  datesi=date0+datenum(iy,1,1)-1;
%
   datesi=date0;
   nts=length(date0);
   frq2=repmat(frq0,nts,1);
   dfsi=frq0(2)-frq0(1);
   Indexf=find(frq0 > freqp(ifmin)*0.95 & frq0 < freqp(ifmax)*1.05);
   Esi=sum(spectre0(:,Indexf),2).*dfsi;
   Esi1=sum(spectre1(:,Indexf),2).*dfsi;
   dsi1=sqrt(Esi1)*1E6;
%
% End of seismic data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end

%
% Frequency dependance for Q. 
Qf=freq.*0+1;
load('UQPREM');

Qf=0.4+0.4*(1-tanh(15*(2*freq-0.14)));  
Qfp=0.4+0.4*(1-tanh(15*(freqp-0.14)));  

%%% WARNING !!!!
Qf=Qf.*0+1;
Qfp=Qfp.*0+1;
% Adjustment of reflection coefficient... 
%r=max(0, (1.4-8*freq)).*0+2.;  % this is for RMAP runs
%r=max(0, (1.4-8*freq)); %.*0+1;  % 


% Reads in the depths used in the seismic source calculation
% Here is an exemple with an ASCII file format
file_bathy1='ww3.07121700.dpt'; 
% Here is an exemple with an NetCDF
file_bathy2='ww3.PACE_dpt.nc';
[lat,lon,dpt1,mat2,norm,var,nodata]=readWW31(file_bathy1);dpt1=dpt1'; 
% Here is an exemple with an NetCDF
%file_bathy2='ww3.PACE_dpt.nc';
%[lat2,lon2,time,dpt,mat2,mat3,var1,var2,var3,unit1,unit2,unit3]=read_WWNC(file_bathy2);
%dpt2=double(dpt);
