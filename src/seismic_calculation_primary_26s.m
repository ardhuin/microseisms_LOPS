if (exist('freqp') == 0) 
    seismic_calculation_synthetic_settings
    ndir
end

nQ=length(QQ);
for iQ=1:nQ
   Q=QQ(iQ).*Qfp;
   OK=exist('UQ');
   if (OK == 1) 
       fUQ=1./TUQ;
       QP=interp1(fUQ,UQ,freqp);
       I=find(freqp > 0.075);
       QP(I)=QP(I(1)-1);
       Q=QQ(iQ).*Qfp;
       Q=QP.*(1-tanh((freqp-0.05)./0.02))./2+Q.*(1+tanh((freqp-0.05)./0.02))./2;
   end

%
% Gets the modelled seismic responses 
%
calcspec=0
if calcspec > 0     
% Method 2: uses actual frequency-direction spectra
   [Ef_pri2 ,pos,freqp,dfps,times ]= ...
         seismic_response_primary_from_spec_v2 (pathref1,  dpt1, lat,lon, CgR, ...
         Q, date1,date2,dt,lono,lato, 4000., statname, ndir);

     
%   seismic_response_primary_from_spec

% Method 1: assumes isotropic wave spectrum
%   [Ef_pri1 ,pos,freqp,dfps,times ]= ...
%         seismic_response_primary_list (pathref1, ...
%            dpt1, lat,lon, CgR, Q, date1,date2,dt,lono,lato, varef, 4000.);


   prim=10-6*tanh((freqp-0.01)*20);
   prim=prim.*0+7.7; % adjustment of bottom slope here 
   prim2=repmat(prim,1,length(times));
   
   
   nf=length(freqp);
   nt=length(times);
   df2=repmat(dfps,1,nt);
   freq=freqp;


   I=find(isnan(Ef_pri2));
   Ef_pri2(I)=0;
%
% Saves synthetic spectra to file for later use
%
   Qstr=sprintf('%.3d',QQ(iQ));
   oname=[statname '_' modgrid1 '_' yname im vtag Qstr '_primary.mat' ];
%   save(oname, 'Ef_pri1','Ef_pri2','freqp','dfps','times', 'Q',  'date1', 'date2', 'dt', 'lono', 'lato');
   save(oname, 'Ef_pri2','freqp','dfps','times', 'Q',  'date1', 'date2', 'dt', 'lono', 'lato','-v7.3');
   
   Ef_pric=prim2.*Ef_pri2;
end
% Testing contribution of sandwaves .. .

deps=25; dA=2E9; % this area is 10 by 10 km  
deps=15; dA=1E8; % this area is 10 by 10 km  
deps=150; dA=1E8; % this area is 10 by 10 km  
deps=30; dA=1E8; % this area is 10 by 10 km  
%dA=0;
ifmins=4;ifmaxs=8;
[Ef_pris , Efx_pris , Efy_pris ,pos,freqps,dfps,times ]= ...
         seismic_response_primary_sandwaves (wave_spectrum, bottom_topography_spectrum,deps, lat,lon, CgR, ...
         Q, date1,date2,dt,lono,lato, dA, statname, ndir);
Qstr=sprintf('%.3d',QQ(iQ));
dstr=sprintf('%.3d',deps);
oname=[statname '_IRL_' modgrid1 '_' yname im vtag Qstr '_' dstr '_primary.mat' ];
save(oname, 'Ef_pris','freqp','dfps','times', 'Q',  'date1', 'date2', 'dt', 'lono', 'lato','-v7.3');

[nf,nt]=size(Ef_pris)


df2s=repmat(dfps,1,nt);
      
%
% Computes standard deviation of the vertical displacement in microns
%
   %Ef_pri2=Ef_pric;
   %Ef_pri2(ifmin:ifmax,:) = Ef_pri2(ifmin:ifmax,:)+Ef_pris(ifmins:ifmaxs,:);
   %delta_pri2=1.E6*sqrt(sum(Ef_pri2(ifmin:ifmax,:).*df2(ifmin:ifmax,:),1));
   delta_pris=1.E6*sqrt(sum(Ef_pris(ifmins:ifmaxs,:).*df2s(ifmins:ifmaxs,:),1));
% Interpolates measurements on the model time steps  
   dsi1i=interp1(datesi,dsi1,times);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plots time series
%

d1=date1;
d2=date2;
i1s=min(find(date0 > d1));
i2s=max(find(date0 < d2));

i1h=min(find(times >= d1));
i2h=min(find(times >= d2));



   figure(1)
   clf
   hold on;
   
   I=find(isfinite(dsi1i) & isfinite(delta_pri2));
   X=dsi1i(I);
   Y=delta_pris(I); 

   P=mean(X)./mean(Y);
   
   hp=plot(times,delta_pris,'r-','LineWidth',2);
   set(gca,'LineWidth',2,'FontSize',16,'Position',[0.05 0.165312 0.93 0.756978]);
   datetick('x',19);
   grid on;
   xlabel(['time (mm/dd of 2013), UTC'],'FontSize',20);
%
% Overlays seismic data
%
   I1=find(datesi > times(1));
   I2=find(datesi < times(end));
   I1=I1(1);
   I2=I2(end);
      hp=plot(times,delta_pri2,'k-','LineWidth',2);
   plot(datesi(I1:I2),dsi1(I1:I2),'b-','LineWidth',2);
   datetick('x',7)


  dskip=8;
  set(gca,'XTick',(times(i1h):dskip:times(i2h)), ...
      'XMinorTick','on','TickDir','out','XTickLabel',datestr(times(i1h):dskip:times(i2h),'mm/dd'), ...
      'FontSize',16,'Box','on');
   
   ylabel('<{\delta^2}>^{0.5} ({\mu}m)','FontSize',18)
   hl=legend('sandwave patch','7.7% slope', ...
             'measured');  %'observed raw data'); %,
   set(hl,'FontSize',14);
   hold off;
%  axis([times(i1h),times(i2h),0 12]);
%
% Saves figure to file
%
   oname=[statname '_' modgrid1 '_' yname im vtag Qstr 'pri_timeseries.png' ];
   saveas(gcf, oname, 'png')
   oname=[statname  '_' modgrid1 '_' yname im vtag Qstr 'pri_timeseries.fig' ];
   saveas(gcf, oname, 'fig')
%
% Compute frequency-dependent ratios
%
   ratio=zeros(nf,5);
   finf=(1+1/1.1)/2;
   fsup=(1+1.1)/2;
   df0=dfsi*0.5;
   for i=1:nf
    f1=freq(i)*finf;
    f2=freq(i)*fsup;
    Indexf=find(frq0+df0 > f1 & frq0-df0 < f2);
    nind=length(Indexf);
    Esi1f=Esi1.*0;
    for j=1:nind
        wf=min([f2,frq0(Indexf(j))+df0])-max([f1,frq0(Indexf(j))-df0]);
        Esi1f=Esi1f+spectre1(:,Indexf(j))*wf;
    end
    dsi1f=sqrt(Esi1f)*1E6;
    dsi1if=interp1(datesi,dsi1f,times);
    I=find(isfinite(dsi1if));
    delf=1.E6*sqrt(Ef_pri2(i,:).*df2(i));
    ratio(i,1)=sum(delf(I))/sum(dsi1if(I));
    ratio(i,3)=median(dsi1if(I)).^2/(df2(i));
    ratio(i,4)=median(delf(I)).^2/(df2(i));
    delf=1.E6*sqrt(Ef_pri2(i,:).*df2(i));
    ratio(i,5)=median(delf(I)).^2/(df2(i));
   end
%
%  Plots average spectra
%
   figure(22);
   semilogy(freq,ratio(:,3),'b-',freq,ratio(:,4),'k-', ...
       freq,ratio(:,5),'r-','LineWidth',2)    
   xlabel('seismic frequency f_s (Hz)','FontSize',16)
   ylabel('median spectral density ({\mu}m^2 / Hz)','FontSize',16)
   legend('observations','model from 1D','model from 2D')
   set(gca,'XLim',[0.003 0.1])
   %set(gca,'XTick',linspace(0.08,0.32,7))
   set(gca,'FontSize',16)
%
% Saves figure to file
%
   oname=[statname '_' modgrid1 '_' yname im vtag Qstr 'pri_spectra.png' ];
   saveas(gcf, oname, 'png')
   oname=[statname '_' modgrid1 '_' yname im vtag Qstr 'pri_spectra.fig' ];
   saveas(gcf, oname, 'fig')
%
%  Computes mismatch statistics
%
    J=[];
   dsi1i(J)=NaN;
  
   
   I=find(isfinite(dsi1i) & isfinite(delta_pri2));
   X=dsi1i(I);
   Y=delta_pri2(I); %.*P; 
   obs_rms=sum(X.^2);
   obs_mean=mean(X);
   obs_scat=sqrt(sum((X-obs_mean).^2));
   mod_mean=mean(Y);
   mod_scat=sqrt(sum((Y-mod_mean).^2));
   nrmse=sqrt(sum((X-Y).^2)./obs_rms);
   corr=sum((X-obs_mean).*(Y-mod_mean))./(obs_scat*mod_scat);
   ymax=max(Y);
   xmax= max(X);
%
%  Makes a scatter plot
%
   figure(23)
   plot(X,Y,'k+','Linewidth',2);
   axis equal;
   maxd=maxdp;
   axis([0 maxd 0 maxd])
   ypos1=maxd*0.9;
   ypos2=maxd*0.8;
   ypos3=maxd*0.7;
   xmax=maxd;
   xpos1=0.05.*xmax;
   str1=sprintf('Corr (r): %6.4f',corr);
   str2=sprintf('NRMSE(%%): %5.2f',100.*nrmse);
   str3=sprintf('P : %5.3f',P);
   ht=text(xpos1,ypos1,str1);
   set(ht,'Fontsize',14)
   ht=text(xpos1,ypos2,str2);
   set(ht,'Fontsize',14)
   ht=text(xpos1,ypos3,str3);
   set(ht,'Fontsize',14)
   set(gca,'LineWidth',1 ,'FontSize',14);
   xlabel('observed <{\delta^2}>^{0.5} ({\mu}m)','FontSize',16)
   ylabel('modelled <{\delta^2}>^{0.5} ({\mu}m)','FontSize',16)
   grid on;
%
% Saves figure to file
%
   oname=[statname '_' modgrid1 '_' yname im vtag Qstr 'pri_scatter.png' ];
   saveas(gcf, oname, 'png')
   oname=[statname '_' modgrid1 '_' yname im vtag Qstr 'pri_scatter.fig' ];
   saveas(gcf, oname, 'fig')
end

%
% plots seismic spectra
%
figure(5)
set(gcf, 'Renderer', 'painters');
d1=date1;
%d1=datenum(2002,05,30);
d2=date2;
%d2=datenum(2002,06,02);
i1s=min(find(date0 > d1));
i2s=max(find(date0 < d2));

pcolor(date0(i1s:i2s),frq0,10*log10(spectre0(i1s:i2s,:)'));
shading flat;colorbar;caxis([-200,-50]);
axis([date0(i1s)+datenum(iy,1,1)-1,date0(i2s)+datenum(iy,1,1)-1,0.005,0.1]);
set(gca,'FontSize',16)
datetick('x',7)
dskip=3;
i1h=min(find(times >= d1));
i2h=min(find(times >= d2));

set(gca,'XTick',(times(i1h):dskip:times(i2h)),'XMinorTick','on','TickDir','out','XTickLabel',datestr(times(i1h):dskip:times(i2h),'dd'),'FontSize',16,'Box','on');
title(['Recorded noise spectrum (BHZ) at ' statname ])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
clf
set(gcf, 'Renderer', 'painters');
pcolor(times(i1h:i2h),freqp,10*log10(double(Ef_pri2(:,i1h:i2h))));
shading flat;colorbar;caxis([-200,-50]);datetick('x',19);
set(gca,'FontSize',16)
axis([times(i1h),times(i2h),0.005,0.1]);
datetick('x',7)

dskip=3;
set(gca,'XTick',(times(i1h):dskip:times(i2h)),'XMinorTick','on','TickDir','out','XTickLabel',datestr(times(i1h):dskip:times(i2h),'dd'),'FontSize',16,'Box','on');
title(['Modeled primary noise spectrum (BHZ) at ' statname ])
