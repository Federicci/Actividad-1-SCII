clc, clear all, close all;

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas

-------------------------------------------------------------------------
%}

m=0.1; F=0.1; l=0.6; g=9.8; M=0.5; 

X01=[0 0 -0.01 0];
% X02=[0 0 3.01 0];
deltat=10^-4;
ts=10;
u=0;

t=0:deltat:(ts-deltat);
valores=zeros(4,ts/deltat);
valores(1,1)=X01(1);
valores(2,1)=X01(2);
valores(3,1)=X01(3);
valores(4,1)=X01(4);
delta_dd=0;
phi_dd=0;

for i=1:1:(ts/deltat-1)
     delta_dd=(u-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*l*phi_dd*cos(valores(3,i)))/(M+m);
     phi_dd=(g*sin(valores(3,i))-delta_dd*cos(valores(3,i)))/l;
     valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
     valores(2,i+1)=valores(2,i)+deltat*delta_dd;
     valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
     valores(4,i+1)=valores(4,i)+deltat*phi_dd;
%      valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
%      valores(2,i+1)=valores(2,i)+deltat*(u-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*g*cos(valores(3,i))*sin(valores(3,i)))/(M+m-m*cos(valores(3,i))*cos(valores(3,i)));
%      valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
%      valores(4,i+1)=valores(4,i)+deltat*(g*sin(valores(3,i))/l-valores(2,i+1)*cos(valores(3,i))/l);
%esta version de Euler no funciona, no se por qué, del sistema no lineal
end

figure
plot(t,valores(3,:))
hold on;

A=[0 1 0 0; 0 -F/M -m*g/M 0; 0 0 0 1; 0 F/(l*M) g*(1+m/M)/l 0];
B=[0; 1/M; 0; -1/(M*l)];
C=[0 0 1 0];
D=0;

valoreslin(1,1)=X01(1);
valoreslin(2,1)=X01(2);
valoreslin(3,1)=X01(3);
valoreslin(4,1)=X01(4);

for i=2:1:(ts/deltat-1)
    valoreslin(1,i)=valoreslin(1,i-1)+deltat*valoreslin(2,i-1);
    valoreslin(2,i)=valoreslin(2,i-1)+deltat*(-F*valoreslin(2,i-1)/M-m*g*valoreslin(3,i-1)/M);
    valoreslin(3,i)=valoreslin(3,i-1)+deltat*valoreslin(4,i-1);
    valoreslin(4,i)=valoreslin(4,i-1)+deltat*(F*valoreslin(2,i-1)/(l*M)+(1+m/M)*g*valoreslin(3,i-1)/l);
%esta forma de Euler no anda, no se por qué, del sistema linealizado
end

