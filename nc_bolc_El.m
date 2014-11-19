function [ n ] = nc_bolc_El(t)  %��� ������������ ������ ������� ��� ��������������� �������������, 
                               %� ������ ������������ ����������,������ � � �������, n-������������
global e;
global en;
global bolcman;
Be=[1.9313,1.6912];
Z=(sum(exp(-e/(bolcman*t)))/Be(1)+6*sum(exp(-en/(bolcman*t)))/Be(2));
n=[exp(-e/(bolcman*t))/Be(1), 6*exp(-en/(bolcman*t))/Be(2)]/Z;
end

