function dispNewtonTH=dispNewtonTH(f,dep)
% inverts the linear dispersion relation (2*pi*f)^2=g*k*tanh(k*dep) to get 
% k from f and dep. 2 Arguments: f and dep. 
eps=0.000001;	
g=9.81;
   sig=2.*pi.*f;
   Y=dep.*sig.^2./g ;   %a is the squared adimensional frequency
	X=sqrt(Y);
	I=1;
   F=1.;
   while abs(max(F)) > eps
		H=tanh(X);
		F=Y-X.*H;
		FD=-H-X./cosh(X).^2;
		X=X-F./FD;
   end
   dispNewtonTH=X./dep;
   

