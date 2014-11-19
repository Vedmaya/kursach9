function [ n ] = nc_bolc777(t,j,x,o)
plank=6.6260695729*10^-34;
light_speed=3*10^10;
bolcman=1.380648813*10^-23;
system_koef=(light_speed*plank)/bolcman;
y=zeros(1,j+1);
z=0;
e=zeros(1,j+1);
Z_el=1+6*exp(-48686.7*system_koef/t);
n0_1=1/Z_el*exp(-o*system_koef/t);
for i=1:1:j+1
    e(i)=(x(1)*(i-1/2)-x(2)*(i-1/2)^2)*plank*light_speed;
    z=z+exp(-e(i)/(bolcman*t));
end;
for i=1:1:j+1
    y(i)=exp(-e(i)/(bolcman*t))/z;
end;
n=y*n0_1;
end

