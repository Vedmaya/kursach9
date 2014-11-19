function [ z ] = VE3_invers( T, n,en,en1)       %зависимость k_VE от j без учёта выбора наилучшей возможность VE перехода
s0=5*10^(-13);
c=13.58*10^(-2);
x=[1743.41 , 14.36 , -0.045 ];
x2=[2169.81 , 13.288 , 0.0105 ];
o1=48686.7;
o2=0;
k=1.3807*10^(-23);
j=[0:1:34];
plank=6.6260695729*10^-34;
light_speed=3*10^10;
system_koef=(light_speed*plank);
bolcman=1.380648813*10^-23;     %Дж/К
E2=(x2(1)*(n+1/2)-x2(2)*(n+1/2)^2)+o2;
E=zeros(1,length(j));
deltaE=zeros(1,length(j));
deltaE1=zeros(1,length(j));
w=zeros(1,length(j));
for l=1:1:length(j)
    E(l)=(x(1)*(j(l)+1/2)-x(2)*(j(l)+1/2)^2)+o1;
    deltaE(l)=E(l)-E2;
    deltaE1(l)=deltaE(l)*system_koef;
    w(l)=s0*exp(-(deltaE(l)/(c*x2(1)))^2)*exp(-deltaE1(l)/(2*k*T));
end;
[x,y]=max(w);
Be1=1.9313*6;
Be2=1.6912;
Be=Be2/Be1;
if(y==1)
    figny1=w(1)*Be*exp((en1(1)-en(n+1))/(bolcman*T));
    figny2=w(2)*Be*exp((en1(2)-en(n+1))/(bolcman*T));
    g=[0,w(1),w(2),0,figny1,figny2,y];
    
else if(y==length(j))
        figny1=w(y-1)*Be*exp((en1(y-1)-en(n+1))/(bolcman*T));
        figny2=w(y)*Be*exp((en1(y)-en(n+1))/(bolcman*T));
        g=[w(y-1),w(y),0,figny1,figny2,0,y];
    else
        figny1=w(y-1)*Be*exp((en1(y-1)-en(n+1))/(bolcman*T));
        figny2=w(y)*Be*exp((en1(y)-en(n+1))/(bolcman*T));
        figny3=w(y+1)*Be*exp((en1(y+1)-en(n+1))/(T*bolcman));
        g=[w(y-1),w(y),w(y+1),figny1,figny2,figny3,y];
    end;
end;
z=g;
end

