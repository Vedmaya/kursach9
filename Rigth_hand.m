function [z] = Rigth_hand(y,en,en1,tau0,n0,T0)   
bolcman=1.380648813*10^-23;     %Дж/К

w=zeros(105,1);
n=[y(1:69)];
n2=[y(70:104)];
T=y(105)*T0;
R0=zeros(1,104);
for i=1:1:69
    R0(i)=R(i-1,T,n,n2,en,en1,tau0,n0);
    w(i)=R0(i);
end;
for i=1:1:35
    R0(69+i)=R_2(i-1,T,n,n2,en,en1,tau0,n0);
    w(69+i)=R0(69+i);
end;
S_R=(sum(R0(1:69).*en)+sum(R0(70:104).*en1))*(n0/tau0);
N=sum(y(1:104))*n0;
w(105)=-tau0/T0*(2/(5*bolcman*N))*S_R;
z=w;    
end

