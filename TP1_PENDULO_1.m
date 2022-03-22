clc, clear all, close all;

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas

-------------------------------------------------------------------------
%}

m=0.1; F=0.1; l=0.6; g=9.8; M=0.5; 
X01=[0 0 -0.01 0];
X02=[0 0 3.01 0];
deltat=10^-4;
ts=10;
u=1;

t=0:deltat:(ts-deltat);
valores=zeros(4,ts/deltat);
valores(1,1)=X01(1);
valores(2,1)=X01(2);
valores(3,1)=X01(3);
valores(4,1)=X01(4);

for i=2:1:(ts/deltat-1)
    valores(1,i)=valores(1,i-1)+deltat*valores(2,i-1);
    valores(2,i)=valores(2,i-1)+deltat*(u-F*valores(2,i-1)+m*l*valores(4,i-1)*valores(4,i-1)*sin(valores(3,i-1))-m*g*cos(valores(3,i-1))*sin(valores(3,i-1)))/(M+m-m*cos(valores(3,i-1))*cos(valores(3,i-1)));
    valores(3,i)=valores(3,i-1)+deltat*valores(4,i-1);
    valores(4,i)=valores(4,i-1)+deltat*(g*sin(valores(3,i-1))/l-valores(2,i)*cos(valores(3,i-1))/l);
end

figure
plot(t,valores(3,:))

A=[0 1 0 0; 0 -F/M -m*g/M 0; 0 0 0 1; 0 F/(l*M) g*(1+m/M)/l 0];
B=[0; 1/M; 0; -1/(M*l)];
C=[0 0 1 0];
D=0;

sys=ss(A,B,C,D);