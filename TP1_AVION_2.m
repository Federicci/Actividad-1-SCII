clc, clear all, close all;

%Motor con carga
%Variables de estado: x1=ia, x2=wr, x3=titat

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas
Agregar en el informe la obtenci�n del sistema linealizado en el equilibrio
X=[0 0 0 0]t
-------------------------------------------------------------------------
%}

w=2;
a=0.05;
b=5;
c=50;
deltat=10^-3;
ts=20;
u=0;
h0=1000; %Altura inicial

alpha_0=0; tita_0=0; titadif_0=0; h_0=0; %Cond. iniciales

t=0:deltat:(ts-deltat);
valores=zeros(4,ts/deltat);
valores(1,1)=alpha_0;
valores(2,1)=tita_0;
valores(3,1)=titadif_0;
valores(4,1)=h_0;

for i=2:1:(ts/deltat-1)
    valores(1,i)=valores(1,i-1)+deltat*(a*valores(2,i-1)-a*valores(1,i-1));
    valores(2,i)=valores(2,i-1)+deltat*(valores(3,i-1));
    valores(3,i)=valores(3,i-1)+deltat*(-w^2*valores(2,i-1)+w^2*valores(1,i-1)+b*w^2*u);
    valores(4,i)=valores(4,i-1)+deltat*(c*valores(1,i-1));
end

figure
plot(t,valores(4,:)+h0);