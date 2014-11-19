function [z] =energy_el(c)
global bolcman;
plank=6.6260695729*10^-34;
light_speed=3*10^10;
system_koef=(light_speed*plank)/bolcman; 
if(c==1)
    x=[2169.81 , 13.288];
    e=zeros(1,69);
    for(j=0:1:68)
        e(j+1)=(x(1)*(j+1/2)-x(2)*(j+1/2)^2)*system_koef*bolcman;
    end;
else
    x=[1743.41 , 14.36];
    Te=48686.7;
    e=zeros(1,35);
    for(j=0:1:34)
        e(j+1)=(Te+x(1)*(j+1/2)-x(2)*(j+1/2)^2)*system_koef*bolcman;
    end;
end;
z=e;
end

