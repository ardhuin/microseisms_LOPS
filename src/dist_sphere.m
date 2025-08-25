function d=dist_sphere(lo1,lo2,la1,la2);
dtor=pi/180.;
d=acos(sin(la2.*dtor).*sin(la1.*dtor)+ ...
    cos(la2.*dtor).*cos(la1.*dtor).*cos((lo2-lo1).*dtor))./dtor;
end