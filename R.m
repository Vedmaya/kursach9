function [z] = R(j,T,n,n_1,en,en1,tau0,n0)
plank=6.6260695729*10^-34;      %Дж*с    
bolcman=1.380648813*10^-23;     %Дж/К
light_speed=3*10^10;            %см/с
system_koef=(light_speed*plank)/bolcman;
x=[2169.81 , 13.288];           %спектроскопические постоянные для X1Sigma+
%-------------------------------------------------------------------------- R_VV
L=[0:1:67]; 
k_1=VV2_invers(j+1,L,T,x)*tau0*n0;
if(j==0)
    k_2=zeros(1,length(L));
else
    k_2=VV2_invers(j,L,T,x)*tau0*n0;
end;
k_3=k_1;
k_4=k_2;
for i=1:1:length(k_1)
    k_3(i)=k_3(i)*exp((2*x(2)*(j-i+1))*system_koef/T);
    k_4(i)=k_4(i)*exp((2*x(2)*(j-i))*system_koef/T);
end;
n1=[n(1:68)]; 
r_1=0;
if(j+2<70)  
      r_1=(k_1*n1)*n(j+2);
end;
r_2=0;
if(j+1<70)  
    r_2=(k_2*n1)*n(j+1);
end;
n2=[n(2:69)]; 
r_3=0;
if(j+1<70)
    r_3=(k_3*n2)*n(j+1);
end;
if(j==0)
    r_VV=(r_1-r_3);
else
    r_4=(k_4*n2)*n(j);
    r_VV=(r_1+r_4-r_2-r_3);
end;
R_VV=r_VV;
%--------------------------------------------------------------------------
R=R_VV;
z=R*10^(-6);
end

