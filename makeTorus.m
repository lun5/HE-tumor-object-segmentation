     a=5;
     c=10;
     [u,v]=meshgrid(0:10:360);
     x=(c+a*cosd(v)).*cosd(u);
     y=(c+a*cosd(v)).*sind(u);
     z=a*sind(v);
     surfl(x,y,z)
     axis equal;