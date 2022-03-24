clc, clear all, close all;

%Motor con carga
%Variables de estado: x1=ia, x2=wr, x3=titat

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas

-------------------------------------------------------------------------
%}

mediciones=xlsread('Curvas_Medidas_Motor.xls');
%Columnas: Tiempo wr ia
ti=0.02;
x1=mediciones(3664,2);
t1=0.02022;
x2=mediciones(3671,2);
t2=0.02043;
x3=mediciones(3679,2);
t3=0.02067;

%ALGORITMO DE CHEN para aproximación de FT de la forma
%G(s)=K*(T3*s+1)/((T1*s+1)*(T2*s+1))
ganancia=198.248802;
k1=x1/ganancia-1;
k2=x2/ganancia-1; 
k3=x3/ganancia-1; 
b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta=(2*k1^3+3*k1*k2+k3-sqrt(b))/(sqrt(b));
T1=-0.001/log(alfa1);
T2=-0.001/log(alfa2);
T3=beta*(T1-T2)+T1; %No hay cero, no se usa

s=tf('s');
G=ganancia/((T1*s+1)*(T2*s+1));
[num,den]=tfdata(G,'v');
num=num/12;

t=0:0.001:0.6;
u=zeros(1,0.6/0.001+1);
for j=20:1:600
    u(1,j)=1;
end

figure
lsim(G,u,t);
hold on;
plot(mediciones(:,1),mediciones(:,2),'r');

%Prueba con valores arbitrarios (?)
Ki=num(3);
Laa=1.2;
J=den(1)/Laa;
Bm=0.0001;
Ra=(den(2)-Laa*Bm)/J;
Km=(den(3)-Ra*Bm)/Ki;

ia_0=0; wr_0=0; titat_0=0; %Cond. iniciales
ts=0.6;
deltat=10^-7;

Va=zeros(1,ts/deltat);
for i=200000:1:6000000
    Va(1,i)=12;
end

T=zeros(1,ts/deltat);
for i=1000000:1:6000000
    T(1,i)=0.075;
end    

t=0:deltat:(ts-deltat);
variables=zeros(3,ts/deltat);
variables(1,1)=ia_0;
variables(2,1)=wr_0;
variables(3,1)=titat_0;

for i=2:1:(ts/deltat)
    variables(1,i)=variables(1,i-1)+deltat*(-Ra*variables(1,i-1)/Laa-Km*variables(2,i-1)/Laa+Va(1,i-1)/Laa);
    variables(2,i)=variables(2,i-1)+deltat*(Ki*variables(1,i-1)/J-Bm*variables(2,i-1)/J-T(1,i-1)/J);
    variables(3,i)=variables(3,i-1)+deltat*variables(2,i-1);
end

plot(t,variables(2,:),'g');

