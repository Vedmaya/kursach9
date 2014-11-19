function [ z ] = VE3_invers2_all( T, n,en,en1)       %зависимость k_VE от j без учёта выбора наилучшей возможность VE перехода
s0=5*10^(-13);
c=13.58*10^(-2);
x=[1743.41 , 14.36];
x2=[2169.81 , 13.288];
o1=48686.7;
o2=0;
k=1.3807*10^(-23);
j=[0:1:68];
plank=6.6260695729*10^-34;
light_speed=3*10^10;
system_koef=(light_speed*plank);
bolcman=1.380648813*10^-23;     %Дж/К
E2=(x(1)*(n+1/2)-x(2)*(n+1/2)^2)+o1;
E=zeros(1,length(j));
deltaE=zeros(1,length(j));
deltaE1=zeros(1,length(j));
w=zeros(1,length(j));
for l=1:1:length(j)
    E(l)=(x2(1)*(j(l)+1/2)-x2(2)*(j(l)+1/2)^2)+o2;
    deltaE(l)=-E(l)+E2;
    deltaE1(l)=deltaE(l)*system_koef;
    w(l)=s0*exp(-(deltaE(l)/(c*x(1)))^2)*exp(-deltaE1(l)/(2*k*T));
end;
u=w;
% Be1=1.9313*6;
% Be2=1.6912;
% Be=Be2/Be1;
for i=1:69
    u(i)=u(i)*exp((en1(n+1)-en(i))/(T*bolcman));
end;
z=[w,u];
end


