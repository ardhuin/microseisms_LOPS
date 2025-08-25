
function [lato,lono,thresh1,thresh2,thresh3,maxd,maxdp,r,network,startdate] = seismic_stations(freq,statname)
network='undefined';
startdate='undefined';
switch statname
    case 'ADK'
        lono = -176.6844
        lato = 51.8837
        network='IU';
        startdate='2009/07/19';
        thresh1 = 1E-14;
        thresh2 = 10;
        thresh3 = 1E-8         
        maxd=10; 
        r=max(0, (0.8-5*freq));
    case 'AIS'
        lono = 	77.57;
        lato = 	-37.80;
        network='G';
        startdate='1993/12/25';
        thresh1 = 1E-14;
        thresh2 = 10;
        thresh3 = 1E-8         
        maxd=10; 
        r=max(0, (0.8-5*freq));
    case 'ALOHA'
        lato=22+45/60.;
        lono=-158.0;
        network='UH';
        startdate='2007/02/14';
        thresh1=1E-13; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=6; 
        r=max(0, (0.8-5*freq));
    case 'ANMO'
        lato=34.946;
        lono=-106.457;
        network='IU';
        startdate='2007/02/14';
        thresh1=1E-13; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=1.5; 
        r=max(0, (0.8-5*freq));
    case 'ATD'
        lono = 42.85;
        lato = 11.53; % Arta Cave - Arta, Republic of Djibouti 
        network='G';
        startdate='2010/01/01';
        thresh1 = 1E-14;
        thresh2 = 10;
        thresh3 = 1E-8         
        maxd=10; 
        r=max(0, (0.8-5*freq));
    case 'BBSR'
        lato = 32.3712;
        lono = -64.6962;
        network='IU';
        startdate='2002/07/19';
        thresh1 = 1E-14; % SAME AS CMLA
        thresh2 = 10;    % SAME AS CMLA
        thresh3 = 1E-9;  % SAME AS CMLA
        maxd=8;         % SAME AS CMLA
        r=max(0, (0.8-5*freq));
    case 'BFO'
        station_full_name='Black Forest Observatory, Schiltach, Germany';  
        lato=48.33; 
        lono=8.33;         
        network='II';
        startdate='1986/01/01';
        thresh1=1E-15; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=3; 
        r=max(0, (0.8-5*freq));        
    case 'BKS'
        lato    = 37.8762;
        lono    = -122.2356;
        network='BK';
        startdate='1991/05/01';
        thresh1 = 1E-14; 
        thresh2 = 10; 
        thresh3 = 5.E-8;
        maxd    = 5;
        r       = max(0, (1.4-8*freq));
    case 'BORG'
        lato=64.7474; 
        lono=-21.3268;   
        network='IU';
        startdate='1994/07/30';
        thresh1=3E-12;
        thresh2=8;
        thresh3=1.E-8;
        maxd=48;
        r=max(0, (0.8-5*freq));
        picture='http://www.iris.edu/hq/gallery/photo/443';
   case 'CAN'
        lono = 	149.00;
        lato = 	-35.32;
        network='G';
        startdate='1987/11/27';
        thresh1 = 1E-14; 
        thresh2 = 10;
        thresh3 = 1E-9;
        maxd=8; 
        r=max(0, (0.8-5*freq)); % to be adjusted ... 
   case 'CMLA'
        lono = -25.5243;
        lato = 37.7637;
        network='II';
        startdate='1996/03/10';
        thresh1 = 1E-14; 
        thresh2 = 10;
        thresh3 = 1E-9;
        maxd=8; 
        r=max(0, (0.8-5*freq));
   case 'CCM'
        lono = -91.245;
        lato = 38.056;
        network='IU';
        startdate='1989/01/01';
        thresh1 = 1E-14; 
        thresh2 = 10;
        thresh3 = 1E-9;
        maxd=2; 
        r=max(0, (0.8-5*freq));
    case 'COR'
        lato=44.586; 
        lono=-123.305;
        network='IU';
        startdate='1992/01/01';
        thresh1=1E-14;
        thresh2=8;
        thresh3=1.E-8;
        maxd=10;
        r=max(0, (0.8-5*freq))./6;
    case 'COYC'
        lato=-45.57; 
        lono=-72.08;
        network='G';
        startdate='2004/12/17';
        thresh1=1E-14;
        thresh2=8;
        thresh3=1.E-8;
        maxd=10;
        r=max(0, (0.8-5*freq))./6;
    case 'CRZF'
        lato=-46.43; 
        lono=51.86;
        network='G';
        startdate='1986/02/01';
        thresh1=1E-13;
        thresh2=4;
        thresh3=1.E-8;
        maxd=10;
        r=max(0, (0.8-5*freq))./6;
   case 'DBG'   
        station_full_name='Daneborg, Greenland';  
        lato=74.31; 
        lono= -20.22;      
        startdate='2010/08/11';
        thresh1=1E-14;
        thresh2=10;
        thresh3=5.E-8;
        network='II';
        maxd=5;
        r=max(0, (0.8-5*freq)); 
    case 'DRV'
        lato=-66.66; 
        lono=140.00;
        network='G';
        startdate='1986/02/01';
        thresh1=1E-14;
        thresh2=8;
        thresh3=1.E-8;
        maxd=10;
        r=max(0, (0.8-5*freq))./6;
    case 'DWPF'
        lato=	28.11; 
        lono=	-81.43;
        network='IU';
        startdate='2010/09/24';
        thresh1=1E-14;
        thresh2=8;
        thresh3=1.E-8;
        maxd=2;
        r=max(0, (0.8-5*freq))./6;
    case 'EDA'
        lato=	3.778868; 	
        lono=	10.153427;
        network='G';
        startdate='2019/05/04';
        thresh1=1E-14;
        thresh2=8;
        thresh3=1.E-8;
        maxd=2;
        maxdp=1;           
        r=max(0, (0.8-5*freq))./6;
    case 'ELSH'
        lono = 1.136600;
        lato = 51.147600;
        network='GB';
        startdate='2008/06/21';
        thresh1=1E-15; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=10; 
        maxdp=0.3; 
        r=max(0, (0.8-5*freq));
    case 'ERM'
        lono = 143.1571
        lato = 42.0150
        network='IU';
        startdate='1990/05/21';
        thresh1 = NaN;
        thresh2 = NaN;
        thresh3 = NaN;         
        maxd=10; 
        r=max(0, (0.8-5*freq));
    case 'ESK'
        lato    = 55.317;
        lono    = -3.205;
        network='II';
        startdate='1987/11/13';
        thresh1=1E-14; 
        thresh2=30; 
        thresh3=5.E-8;
        maxd=5;
        r=max(0, (0.8-5*freq));
    case 'FDF'
        lato    = 14.7349711;
        lono    = -61.146311;
        network='G';
        startdate='2008/09/01';
        thresh1=1E-14; 
        thresh2=10; 
        thresh3=1.E-9;
        maxd=10;
        r=max(0, (0.8-5*freq));
    case 'FOMA'
        lato    = 	-24.98;
        lono    = 	46.98;
        network='G';
        startdate='2008/09/01';
        thresh1=1E-14; 
        thresh2=30; 
        thresh3=5.E-8;
        maxd=5;
        r=max(0, (0.8-5*freq));
    case 'GRFO'    
        station_full_name='Grafenberg, Germany';
        lono =49.69;
        lato =11.22;
        network='IU';
        startdate='1994/01/26';
        thresh1=1E-13; % same as KIP
        thresh2=8; % same as KIP
        thresh3=1.E-11; % same as KIP         
         maxd=5; 
        r=max(0, (0.8-5*freq));
    case 'H2O'
        lono = -141.991736
        lato = 27.881910
        network='H2';
        startdate='1999/10/04';
        thresh1 = 1E-14;
        thresh2 = 10;
        thresh3 = 1E-8         
        maxd=10; 
        r=max(0, (0.8-5*freq));
    case 'HRV' % FAIT MO 09/06/11
        lato = 42.5064;
        lono = -71.5583;
        network='IU';
        startdate='1988/01/01';
        thresh1=0.5E-14; 
        thresh2=3;     
        thresh3=1.E-8; 
        maxd=10;       
        r=max(0, (0.8-5*freq));
    case 'ICC';
        network='G';
        startdate='1996/04/01';
        lato= -20.28 ; 
        lono= -70.03;  % Iquique, Chile
        thresh1=1E-14;  % TO BE ADJUSTED
        thresh2=8; thresh3=1.E-9;maxd=4;
        r=max(0, (0.8-5*freq));
    case 'IDGL'
        lono = 	-7.51
        lato = 	55.07 % INCH ISLAND, CO DONEGAL, IRELAND
        network='EI';
        startdate='2011/01/01';
        thresh1 = 1e-13;    
        thresh2 = 4;        
        thresh3 = 1E-11;    
        maxd=10;            
        r=max(0, (0.8-5*freq));
    case 'IGLA'
        lono = 	-9.375000
        lato = 	53.419500 % Station Glengowla, Ireland
        network='EI';
        startdate='2011/01/01';
        thresh1 = 3e-15;    
        thresh2 = 4;        
        thresh3 = 1E-11;    
        maxd=10;  
        maxdp=0.6;           
        r=max(0, (0.8-5*freq));
    case 'INU';
        network='G';
        lato= 35.350 ; 
        lono= 137.029 ;  % Inuyama, Japon
        thresh1=1E-14; 
        thresh2=8; thresh3=1.E-9;maxd=4;
        r=max(0, (0.8-5*freq));
    case 'JOHN'
        lono = -169.5292
        lato = 16.7329
        network='IU';
        startdate='1998/08/04';
        thresh1=1E-13; % same as KIP
        thresh2=8; % same as KIP
        thresh3=1.E-11; % same as KIP         
         maxd=10; 
        r=max(0, (0.8-5*freq));
    case 'JTS'
        lono = -84.953;
        lato =  10.291;
        network='IU';
        thresh1=1E-15; 
        thresh2=10; 
        thresh3=5.E-8;
        maxd=3; 
        r=max(0, (0.8-5*freq));
    case 'KBS'    
        station_full_name='Ny-Alesund, Spitzbergen, Norway';
        lono =11.94;
        lato =78.92;
        network='IU';
        startdate='2010/06/08';
        thresh1=1E-13; % same as KIP
        thresh2=8; % same as KIP
        thresh3=1.E-11; % same as KIP         
         maxd=5; 
        r=max(0, (0.8-5*freq));
    case 'KEKH'
        lato=21.98;
        lono=-159.71;
        thresh1=1E-13; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=6; 
        r=max(0, (0.8-5*freq));
        station_full_name='Kekaha, Kauai, Hawaii';
        network='PT';
    case 'KDAK'
        lato=57.7828;
        lono= -152.5835;
        network='II';
        startdate='1997/06/09';
        thresh1=1E-14; 
        thresh2=10; 
        thresh3=5.E-8;
        maxd=10;
        r=max(0, (0.8-5*freq));
    case 'KIP'
        lato=21.4233;
        lono=-158.015;
        network='G';
        startdate='1986/04/17';
        thresh1=1E-13; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=10; 
        r=max(0, (0.8-5*freq));
    case 'KIEV'
        lato=50.70;
        lono=29.22;
        network='IU';
        startdate='1995/01/30';
        thresh1=1E-13; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=10; 
        r=max(0, (0.8-5*freq));
    case 'KNTN'
        lato=-2.77;
        lono=-171.72;
        network='IU';
        startdate='2007/12/04';
        thresh1=1E-13; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=6; 
        r=max(0, (0.8-5*freq));
    case 'KONO' % FAIT MO 10/06/11
        lato    = 59.6491;
        lono    = 9.5982 ;
        thresh1=1E-14;  
        thresh2=10;     
        thresh3=5.E-8;  
        maxd=10;        
        r=max(0, (0.8-5*freq));
% for the station below, we need to define thresholds
    case 'KWAJ'
        lono = 167.6130
        lato = 8.8019
        thresh1 = 1e-13;    %SAME AS WAKE
        thresh2 = 4;        %SAME AS WAKE
        thresh3 = 1E-11;    %SAME AS WAKE
        maxd=10;            %SAME AS WAKE
        r=max(0, (0.8-5*freq));
    case 'LCO'
        lato=-29.01;
        lono=-70.70; % Las Campanas Astronomical Observatory, Chile
        network='IU';
        startdate='2009/07/30';
        thresh1=1E-13; % 'to be defined';
        thresh2=8; 
        thresh3=1.E-11;
        maxd=6; 
        r=max(0, (0.8-5*freq));
    case 'LVC'
        lato=-22.61;
        lono=-68.91;
        network='IU';
        startdate='2009/04/10';
        thresh1=1E-13; % 'to be defined';
        thresh2=8; 
        thresh3=1.E-11;
        maxd=2; 
        r=max(0, (0.8-5*freq));
    case 'MAJO'
        lono = 138.2070
        lato = 36.5427
        thresh1 = NaN;
        thresh2 = NaN;
        thresh3 = NaN;         
        maxd=3; 
        r=max(0, (0.8-5*freq));
    case 'MBO'
        station_full_name='M Bour, Senegal'; 
        lato=14.39;
        lono=-16.96; % 
        network='G';
        startdate='1985/09/01';
        thresh1=1E-13; % 'to be defined';
        thresh2=8; 
        thresh3=1.E-11;
        maxd=0.5; 
        r=max(0, (0.8-5*freq));
    case 'MBWA'
        station_full_name='Marble Bar, Western Australia'; 
        lato=-21.159;
        lono=119.731; % 
        network='IU';
        startdate='2001/09/01';
        thresh1=1E-13; % 'to be defined';
        thresh2=8; 
        thresh3=1.E-11;
        maxd=4; 
        r=max(0, (0.8-5*freq));
    case 'MIDW'
        lono = -177.3697
        lato = 28.2157
        thresh1=1E-13; % same as KIP
        thresh2=8; % same as KIP
        thresh3=1.E-11; % same as KIP         
        maxd=6; 
        r=max(0, (0.8-5*freq));
    case 'NACB'
        lato=-24.1738;
        lono=121.5947; 
        network='TW';
        thresh1=1e-14;
        thresh2=5;
        thresh3=1;
        maxd=3;
        r=max(0, (0.8-5*freq)); 
    case 'NE71'
        lato=31.68973; 
        lono=-115.90526;      
        thresh1=1e-14;
        thresh2=5;
        thresh3=1;
        maxd=3;
        r=max(0, (0.8-5*freq)); 
    case 'NE75'     % comme NE71
        lato=27.29334; 
        lono=-112.85649;      
        thresh1=1E-14; 
        thresh2=5;
        thresh3=1;
        maxd=3;
        r=max(0, (0.8-5*freq)); 
    case 'NE79'     % comme NE71
        lato=23.11937; 
        lono=-109.75611;      
        thresh1=1E-14;
        thresh2=5;
        thresh3=1;
        maxd=3;
        r=max(0, (0.8-5*freq)); 
    case 'NNA'     
        lato=-11.988; 
        lono= -76.842;      
        thresh1=1E-14;
        thresh2=10;
        thresh3=5.E-8;
        network='II';
        maxd=1;
        r=max(0, (0.8-5*freq)); 
    case 'NOR'   
        station_full_name='Station Nord, Greenland';  
        lato=81.60; 
        lono= -16.66;      
        startdate='2010/07/06';
        thresh1=1E-14;
        thresh2=10;
        thresh3=5.E-8;
        network='II';
        maxd=2;
        r=max(0, (0.8-5*freq)); 
  case 'NWAO';
        station_full_name='Narrogin, Australia';  
        lato=-32.93;
        lono=117.24;  
        network='IU';
        startdate='1991/11/25';
        thresh1=1E-14; thresh2=10; thresh3=5.E-8;maxd=4;
        r=max(0, (0.8-5*freq));
        %QQ=[40 50 60 70 80 200 400 600 ];Qbest=80;

    case 'OTAV'
        lato=0.24;
        lono=-78.45; % Otavalo, Ecuador
        network='IU';
        startdate='2009/04/01';
        thresh1=1E-13; % 'to be defined';
        thresh2=8; 
        thresh3=1.E-9;
        maxd=0.5; 
        r=max(0, (0.8-5*freq));
    case 'PFO'
        lato=33.6092; 
        lono=-116.4553; % same as PAS      
        thresh1=1E-13;
        thresh2=4;
        thresh3=1.E-8;
        maxd=3;
        r=max(0, (0.8-5*freq));
    case 'POHA'
        lato=19.7575;
        lono=-155.5325;
        thresh1=1E-13;  %SAME AS KIP
        thresh2=8;      %SAME AS KIP
        thresh3=1.E-11; %SAME AS KIP
        maxd=9;        %SAME AS KIP
        r=max(0, (0.8-5*freq));
        
    case  'PAB' % FAIT MO 09/06/11
        lato = 39.5458;
        lono = -4.3483;
        thresh1=0.5E-14;
        thresh2=10;    
        thresh3=1.E-8;
        maxd=10;
        r=max(0, (0.8-5*freq));
    case 'PAF'
        lato=	-49.35; 
        lono=	70.21;   
        network='G';
        thresh1=1E-13;
        thresh2=4;
        thresh3=1.E-8;
        maxd=3;
        r=max(0, (0.8-5*freq));
    case 'PAS'
        lato=34.15; 
        lono=-118.15;      
        thresh1=1E-13;
        thresh2=4;
        thresh3=1.E-8;
        maxd=3;
        r=max(0, (0.8-5*freq));
    case 'PAYG';
        network='IU';
        startdate='2010/05/05';
        lato= -0.67 ; 
        lono= -90.29 ;  %Puerto Ayora, Galapagos Islands 
        thresh1=1E-14;  % TO BE ADJUSTED
        thresh2=8; thresh3=1.E-9;maxd=4;
        r=max(0, (0.8-5*freq));

    case 'PEL';
        network='G';
        startdate='1995/10/04';
        lato= -33.14 ; 
        lono= -70.67 ;  % Peldehue, Chile 
        thresh1=1E-14;  % TO BE ADJUSTED
        thresh2=8; thresh3=1.E-9;maxd=2;
        r=max(0, (0.8-5*freq));
    case 'PET'
        lono = 158.6531
        lato = 53.0239
        thresh1 = 1E-14;
        thresh2 = 10;
        thresh3 = 1E-8;         
        maxd=4; 
        r=max(0, (0.8-5*freq));     
    case 'PTCN';
        network='IU';
        startdate='1996/12/29';
        lato= -25.07 ; 
        lono= -130.10 ;  % Pitcairn Island, South Pacific
        thresh1=1E-14;  % TO BE ADJUSTED
        thresh2=8; thresh3=1.E-9;maxd=4;
        r=max(0, (0.8-5*freq));
    case  'PPT' % done by FA in 2010 
        lato=-17.569;
        lono= -149.576 ;
        network='G';
        startdate='1986/05/31';
        thresh1=1E-14; 
        thresh2=10; 
        thresh3=5.E-8;
        maxd=10;
        r=max(0, (0.8-5*freq)).*0.2;
    case  'PPTF' % done by FA in 2010 
        lato=-17.569;
        lono= -149.576 ;
        network='G';
        startdate='2009/05/27';
        thresh1=1E-14; 
        thresh2=10; 
        thresh3=5.E-8;
        maxd=10;
        r=max(0, (1.4-8*freq));
    case 'NOR'   
        station_full_name='Station Nord, Greenland';  
        lato=81.60; 
        lono= -16.66;      
        startdate='2010/07/06';
        thresh1=1E-14;
        thresh2=10;
        thresh3=5.E-8;
        network='II';
        maxd=2;
        r=max(0, (0.8-5*freq)); 
     case 'PVAQ'
        lato=37.40;
        lono=-7.72;
        network='PM';
        startdate='	2006/12/23';
         thresh1=1E-14; 
        thresh2=10; 
        thresh3=1.E-9;
        maxd=10; 
        r=max(0, (0.8-5*freq));
   case 'RAR'
        lato=-21.21;
        lono=-159.77;
        network='IU';
        startdate='1992/03/07';
         thresh1=1E-14; 
        thresh2=10; 
        thresh3=1.E-9;
        maxd=30; 
        r=max(0, (0.8-5*freq));
     case 'RER'
        station_full_name='Riviere de l Est - Sainte Rose - La Reunion island, France';  
        lato=-21.17;
        lono=55.74;
        network='G';
        startdate='1986/02/10';
         thresh1=1E-14; 
        thresh2=10; 
        thresh3=1.E-9;
        maxd=6; 
        r=max(0, (0.8-5*freq));
    case 'ROCAM'
        lato=	-19.76;
        lono=63.37;
        network='G';
        startdate='2012/12/15';
         thresh1=1E-14; 
        thresh2=10; 
        thresh3=1.E-9;
        maxd=50; 
        r=max(0, (0.8-5*freq));
    case 'RODM'
        lato=		-19.70; % Rodrigues Island, Republic of Mauritius
        lono=63.44;
        network='G';
        startdate='2012/12/15';
         thresh1=1E-14; 
        thresh2=10; 
        thresh3=1.E-9;
        maxd=50; 
        r=max(0, (0.8-5*freq));
    case 'ROSA'
        lato=38.72; 
        lono=-28.25;         
        network='PM';
        startdate='2008/03/06';
        thresh1=1E-15; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=3; 
        r=max(0, (0.8-5*freq));
    case 'RPN'
        lato=-27.13;
        lono=-109.33; %Rapanui, Easter Island, Chile 
        network='II';
        thresh1=1E-13; % TO BE ADJUSTED...
        thresh2=8; 
        thresh3=1.E-11;
        maxd=1.5; 
        r=max(0, (0.8-5*freq));
    case 'RSSD'
        lato=44.121;
        lono=-104.036;
        network='IU';
        startdate='1993/01/01';
        thresh1=1E-13; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=1.5; 
        r=max(0, (0.8-5*freq));
    case 'SACV'
        station_full_name='Santiago Island, Cape Verde';
        lato=14.97; 
        lono=-23.61;         
        network='II';
        startdate='2000/05/29';
        thresh1=1E-15; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=3; 
        r=max(0, (0.8-5*freq));
    case 'SCO' % FAIT FA 08/12/2013
        station_full_name='Ittoqqortoormiit, Greenland';  
        lato=70.49; 
        lono= -21.95;      
        startdate='2010/08/05';
        thresh1=1E-12;
        thresh2=10;
        thresh3=1.E-6;
        network='II';
        maxd=5;
        r=max(0, (0.8-5*freq)); 
    case 'SCZ'
        lato    = 36.59798; 
        lono    = -121.40481;
        thresh1 = 1E-14;  
        thresh2 = 8; 
        thresh3 = 1.E-8;
        maxd    = 5;
        network='G';
        r       = max(0, (0.8-5*freq));
    case 'SFJD' % FAIT MO 09/06/11
        lono = -50.6215;
        lato = 66.9960;
        thresh1=5E-15; 
        thresh2=10;     
        thresh3=1.E-8; 
        maxd=10;       
        r=max(0, (0.8-5*freq));
     case 'SNZO' 
        lono = 174.7046;
        lato = -41.3101;
        thresh1=5E-15; 
        thresh2=10;     
        thresh3=1.E-8; 
        maxd=10;       
        maxdp=0.4; 
        r=max(0, (0.8-5*freq));
        network='IU';
  case 'SOEG'   
        station_full_name='Sodalen, Greenland';  
        lato=68.20; 
        lono= -31.38;      
        startdate='2011/07/27';
        thresh1=1E-14;
        thresh2=10;
        thresh3=5.E-8;
        network='II';
        maxd=5;
        r=max(0, (0.8-5*freq)); 
    case 'SPB'
        station_full_name='Sao Paulo, Brazil';  
        lato=-23.59; 
        lono=-47.43;         
        network='G';
        startdate='1996/06/17';
        thresh1=1E-15; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=5; 
        r=max(0, (0.8-5*freq));
        % QQ=[120 180 220 280 300 320 260 ];Qbest=260;
    case 'SSB'
        lato=45.279; 
        lono=4.542;         
        network='G';
        startdate='1982/05/02';
        thresh1=1E-15; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=3; 
        maxdp=0.2; 
        r=max(0, (0.8-5*freq));
        % QQ=[120 180 220 280 300 320 260 ];Qbest=260;
    case 'SSPA' % FAIT MO 09/06/11
        lato = 40.6358;
        lono = -77.8880;
        thresh1=0.5E-14; 
        thresh2=10;     
        thresh3=1.E-8; 
        maxd=10;       
        r=max(0, (0.8-5*freq));
    case 'SUR' % 
        lato =-32.38;
        lono =20.81;
        network='II';
        startdate='1990/10/30';
        thresh1=1.E-13; 
        thresh2=8;     
        thresh3=1.E-8; 
        maxd=3;       
        r=max(0, (0.8-5*freq));
    case 'PET'
        lono = 158.6531
        lato = 53.0239
        thresh1 = 1E-14;
        thresh2 = 10;
        thresh3 = 1E-8;         
        maxd=10; 
        r=max(0, (0.8-5*freq));     
    case 'TAM'
        lono = 5.53; % Tamanrasset, Algeria 
        network='G';
        lato = 	22.79;
        thresh1 = 1E-14;
        thresh2 = 10;
        thresh3 = 1E-8;         
        maxd=10; 
        r=max(0, (0.8-5*freq));    
        startdate='1983/11/16';
        thresh1=1E-15; thresh2=10; thresh3=5.E-9;maxd=0.5;
        r=max(0, (0.8-5*freq));
        %QQ=[40 50 60 70 80 200 400 600 ];Qbest=80;
     case 'TAOE';
        lato=-8.85;
        lono=-140.15; 
        station_full_name='Taiohae - Marquesas islands, France';
        network='G';
        startdate='2004/11/01';
        thresh1=1E-14; 
        thresh2=10; 
        thresh3=5.E-8;
        r=max(0, (0.8-5*freq));
        maxd=10;
        %QQ=[40 50 60 70 80 200 400 600 ];Qbest=80;
     case 'TEIG';
        lato=	20.23;
        lono=	-88.28; % Tamanrasset, Algeria 
        network='IU';
        startdate='2009/09/18';
        thresh1=1E-14; thresh2=10; thresh3=5.E-8;maxd=1;
        r=max(0, (0.8-5*freq));
        %QQ=[40 50 60 70 80 200 400 600 ];Qbest=80;
    case 'TRIS';
        lato=-37.07;
        lono= -12.32 ;  
        network='G';
        startdate='2004/03/03';
        thresh1=1E-14; thresh2=10; thresh3=5.E-8;maxd=4;
        r=max(0, (0.8-5*freq));
        %QQ=[40 50 60 70 80 200 400 600 ];Qbest=80;
    case 'TRQA';
        station_full_name='Tornquist, Argentina';  
        lato=-38.06;
        lono=-61.98;  
        network='IU';
        startdate='2000/10/28';
        thresh1=1E-14; thresh2=10; thresh3=5.E-8;maxd=4;
        r=max(0, (0.8-5*freq));
        %QQ=[40 50 60 70 80 200 400 600 ];Qbest=80;
    case 'UCC'
        lato=50.279; 
        lono=1.5;         
        network='G';
        startdate='1982/05/02';
        thresh1=1E-15; 
        thresh2=8; 
        thresh3=1.E-11;
        maxd=5; 
        r=max(0, (0.8-5*freq));
        % QQ=[120 180 220 280 300 320 260 ];Qbest=260;
    case 'UNM';
        lato=19.329662;lono= -99.178065 ;
        network='G';
        startdate='1990/06/06';
        thresh1=1E-14; thresh2=10; thresh3=5.E-8;maxd=4;
        r=max(0, (0.8-5*freq));
        %QQ=[40 50 60 70 80 200 400 600 ];Qbest=80;
    case 'WAKE'
        lono = 166.6536;
        lato = 19.2833;
        thresh1 = 1e-13;
        thresh2 = 4;
        thresh3 = 1E-11;         
        maxd=10; 
        r=max(0, (0.8-5*freq)); 
    case 'WKER1'
        thresh1 = 1e-13;
        thresh2 = 4;
        thresh3 = 1E-11;         
        maxd=10; 
        network='IUEM';
        r=max(0, (0.8-5*freq)); 
        lato=-46.634166;
        lono=60.1303;
    case 'XMAS'
        thresh1 = 1e-13;
        thresh2 = 4;
        thresh3 = 1E-11;         
        maxd=10; 
        r=max(0, (0.8-5*freq)); 
        lono=-153.910;    
        lato=0.000;
    case 'YSS'
        lono = 142.7550
        lato = 46.9539
        thresh1 = NaN;
        thresh2 = NaN;
        thresh3 = NaN;         
        maxd=10; 
        r=max(0, (0.8-5*freq));
    otherwise
        disp('STATION NOT FOUND')
        lato    = NaN; 
        lono    = NaN;
        thresh1 = NaN;
        thresh2 = NaN;
        thresh3 = NaN;
        maxd    = NaN;
        r       = NaN(1,length(freq));
end
        
        


% 
% 
% station='KIP';lato=21.4233; lono=-158.015;         
% thresh1=1E-13; thresh2=8; thresh3=1.E-11;maxd=27; %Q=600; 
% r=max(0, (0.8-5*freq));
% 
% 
% station='UNM';lato=19.329662;lono= -99.178065 ;
% thresh1=1E-14; thresh2=10; thresh3=5.E-8;maxd=4;
% r=max(0, (0.8-5*freq));
% 
% station='RPN';lato=-27.127;lono= -109.334 ;
% thresh1=1E-14; thresh2=10; thresh3=5.E-8;maxd=30;
% r=max(0, (0.8-5*freq))./1.5^2;
% 
% station='JTS';lato=10.291;lono= -84.953 ;
% thresh1=1E-15; thresh2=10; thresh3=5.E-8;maxd=3;
% r=max(0, (0.8-5*freq));
% 

% 
% station='TRI';lato=-37.07 ;lono=	-12.32 ;  % Tristan da Cunha: TRIS
% thresh1=1E-14; thresh2=10; thresh3=5.E-8;maxd=30;
% r=max(0, (0.8-5*freq));
% 
% station='CML';lato=37.7665 ;lono=-25.5225 ;  % % Chã de Macela, Azores
% thresh1=1E-14; thresh2=10; thresh3=1.E-9;maxd=30;
% r=max(0, (0.8-5*freq));
% 
% station='NOU';lato= -22.101;lono=166.303 ;  % % Noumea, NC
% thresh1=1E-14; thresh2=10; thresh3=1.E-9;maxd=30;
% r=max(0, (0.8-5*freq));

% station='TAT';lato= 24.974 ; lono= 121.497 ;  % Inuyama, Japon
% thresh1=1E-14; thresh2=8; thresh3=1.E-9;maxd=10;
% r=max(0, (0.8-5*freq));
% 

