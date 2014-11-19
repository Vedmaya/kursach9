function [z] =VV2_invers(j,k,T,x)
plank=6.6260695729*10^-34;
bolcman=1.380648813*10^-23;
light_speed=3*10^10;
system_koef=(light_speed*plank)/bolcman; 

c=0.456;
b=40.36;
a=x(1)/x(2);
D2=((a+1)/(a+3-2*j))^2*j*(a+2-2*j)*(a+4-2*j)/(a*(a+3-j));

deltaE=zeros(1,length(k));
D1=zeros(1,length(k));
lambda=zeros(1,length(k));
F=zeros(1,length(k));
s=zeros(1,length(k));
l=zeros(1,length(k));
Q=zeros(1,length(k));

for i=1:1:length(k)
    deltaE(i)=-2*x(2)*(j-k(i))*system_koef;
    D1(i)=((a+1)/(a+3-2*(k(i))))^2*(k(i))*(a+2-2*(k(i)))*(a+4-2*(k(i)))/(a*(a+3-(k(i))));
    lambda(i)=2^(-3/2)*sqrt(c/T)*abs(deltaE(i));
    F(i)=(3-exp(-2*lambda(i)/3))*exp(-2*lambda(i)/3)/2;
    s(i)=4.93*10^(-4)*(T/300)*k(i)*j/((1-(x(2)/x(1))*k(i))*(1-x(2)/x(1)*j))*F(i);  
    l(i)=5.37*10^(-3)*(300/T)*D1(i)*D2*exp(-deltaE(i)^2/(b*T));   
    Q(i)=3*10^(-10)*(T/300)^(1/2)*(s(i)+l(i))*exp(-deltaE(i)/(2*T)); 
end;
%display(Q);
z=Q;
end


