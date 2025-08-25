bottom_topography_spectrum='spectrum_Ireland_shallow_rocks.bsp';D=17;
bottom_topography_spectrum='spectrum_PDC1.bsp';D=29;
%bottom_topography_spectrum='spectrum_Ireland_rock_lr.bsp';
%bottom_topography_spectrum='spectrum_UK_lr.bsp';
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

figure(104)   
hold on
set(gcf, 'Renderer', 'painters');
pcolor(kbotx-dkbx/2,kboty-dkby/2,10.*log10(botspec'));shading flat;colorbar;
axis equal;
caxis([-0 40]);
axis equal
axis([-0.1,0.1,-0.1,0.1]);
set(gca,'FontSize',16);
D=200
k1=dispNewtonTH(0.1,D);
theta=linspace(0,2*pi,37)
plot(k1.*cos(theta),k1.*sin(theta),'k-','LineWidth',2)
 k1=dispNewtonTH(1/16,D);
plot(k1.*cos(theta),k1.*sin(theta),'k-','LineWidth',2)
k1=dispNewtonTH(1/26,D);
plot(k1.*cos(theta),k1.*sin(theta),'k-','LineWidth',2)
