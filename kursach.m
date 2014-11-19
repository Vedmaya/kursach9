% clear all;
format long e;

global e;
global en;
global I;
global index;
index=0;
delete('backup/time','backup/N','backup/E','backup/T');
I=0;
%--������������������ ����������-------------------------------------------
x1=[2169.81 , 13.288 , 0.0105 ];                            %��� ���� �: �=[we, wexe, weye] X1sigma+
x2=[1743.41 , 14.36 , -0.045 ];                             %a3Pi
x3=[1518.24 , 19.4, 0.766];                                 %A1Pi
%--���-�� ������������� ������� � ������ ����������� ���������-------------
% j1=68;
% j2=34;
% j3=18;
%--------------------------------------------------------------------------
global T0;
global T1;
global bolcman;
global system_koef;
global myu;
global n0;
global tau0;
%--------------------------------------------------------------------------
T0=500;                         %K
T1=5000;                        %K
O=[0,48686.7,65075.77];         %������� ����������� ������� X1Sigma+,a3Pi,A1Pi cm^(-1)
bolcman=1.380648813*10^-23;     %��/�
plank=6.6260695729*10^-34;      %��*�    
light_speed=3*10^10;            %��/�
system_koef=(light_speed*plank)/bolcman;
diam=3.5*10^(-10);              %m
R_gas=8.31;                     %��/�/����
myu=0.028;                      %��/����
p0=10^5;                        %Pa
n0=p0/(bolcman*T0);             %m^(-3)
tau0=bolcman*sqrt(myu*T0/(6*R_gas))/(pi*diam^2*p0);%�
%--------------------------------------------------------------------------
e=energy_el(1);
en=energy_el(2);
n_EL=nc_bolc_El(T0);
n_EL(1)=n_EL(1)-9.97911636*10^(-1);
n_EL(75)=n_EL(75)+9.97911636*10^(-1);
n=n_EL';
y333=[n;1];

% % n40=(nc_2(T0, T1, O(1), 68, x1))';
% % a=sum(n40);
% % n50=(nc_2(T0, T1, O(2), 34, x2))';
% % b=sum(n50);
% y333=[n4;n5;1];

buratino=Rigth_hand(y333,e,en,tau0,n0,T0);
options=odeset('InitialStep',10e-30,'MaxStep',10e-2);
[t, Y]=ode15s(@Solve,[0,1e+6],y333,options);


