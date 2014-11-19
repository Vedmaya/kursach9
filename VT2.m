function [z] =VT2(T,j,x)
plank=6.6260695729*10^-34;
bolcman=1.380648813*10^-23;
light_speed=3*10^10;
system_koef=(light_speed*plank)/bolcman; 

c=0.456;
A=-15.23;
B=280.5;
C=-549.6;

deltaE0=(x(1)-x(2)*2)*system_koef;
lambda0=2^(-3/2)*sqrt(c/T)*abs(deltaE0);
F0=(3-exp(-2*lambda0/3))*exp(-2*lambda0/3)/2;
b=exp(A+B*T^(-1/3)+C*T^(-2/3));
p=10^7*T*bolcman/(b*F0*(1-exp(-(x(1)*system_koef)/T)));

deltaE=zeros(1,length(j));
lambda=zeros(1,length(j));
F=zeros(1,length(j));
P=zeros(1,length(j));

for i=1:1:length(j)
    deltaE(i)=(x(1)-x(2)*2*j(i))*system_koef;
    lambda(i)=2^(-3/2)*sqrt(c/T)*abs(deltaE(i));
    F(i)=(3-exp(-2*lambda(i)/3))*exp(-2*lambda(i)/3)/2;
    P(i)=p*j(i)/(1-(x(2)/x(1))*j(i))*F(i);
end;
%display(P);
z=P;
end

