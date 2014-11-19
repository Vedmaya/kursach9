function [z] = R_2(j,T,n_1,n,en,en1,tau0,n0)
global bolcman;
plank=6.6260695729*10^-34;      %ƒж*с    
light_speed=3*10^10;            %см/с
system_koef=(light_speed*plank)/bolcman;
x=[1743.41 , 14.36];            %спектроскопические посто€нные дл€ a3Pi
%--------------------------------------------------------------------------R_VV
L=[0:1:33];
k_1=VV2_invers(j+1,L,T,x)*tau0*n0;
if(j==0)
    k_2=zeros(1,length(L));
else
    k_2=VV2_invers(j,L,T,x)*tau0*n0;
end;
k_3=k_1;
k_4=k_2;
for i=1:1:length(k_1)
    k_3(i)=k_3(i)*exp(2*x(2)*(j-i+1)*system_koef/T);
    k_4(i)=k_4(i)*exp(2*x(2)*(j-i)*system_koef/T);
end;
n1=[n(1:34)];
r_1=0;
if(j+2<36)
    r_1=(k_1*n1)*n(j+2);
end;
r_2=0;
if(j+1<36)
    r_2=(k_2*n1)*n(j+1);
end;
n2=[n(2:35)];
r_3=0;
if(j+1<36)
    r_3=(k_3*n2)*n(j+1);
end;
if(j==0)
    r_VV=r_1-r_3;
else
    r_4=(k_4*n2)*n(j);
    r_VV=r_1-r_2-r_3+r_4;
end;
R_VV=r_VV;
%--------------------------------------------------------------------------
R=R_VV;
z=R*10^(-6);
end

