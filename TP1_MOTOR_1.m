clc, clear all, close all;

%Motor con carga
%Variables de estado: x1=ia, x2=wr, x3=titat

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas

-------------------------------------------------------------------------
%}

Laa=366e-6;
J=5e-9;
Ra=55.6;
Bm=0;
Ki=6.49e-3;
Km=6.53e-3;

%Simulación por integración de Euler
deltat=10^-7;
Va=12;
Tl=0; %Se itera hasta que ia tenga un valor coherente ~500mA
ia_0=0; wr_0=0; titat_0=0; %Cond. iniciales


iai_1=ia_0+deltat*(-Ra*ia_0/Laa-Km*wr_0/Laa+Va/Laa);
wri_1=wr_0+deltat*(Ki*ia_0/J-Bm*wr_0/J-Tl/J);
titati_1=wr_0+deltat*(wr_0);

t=0:deltat:(1-deltat);
variables=zeros(3,10000000,10);
for j=1:1:10
    Tl=(j-1)*0.001;
    iai_1=ia_0+deltat*(-Ra*ia_0/Laa-Km*wr_0/Laa+Va/Laa);
    wri_1=wr_0+deltat*(Ki*ia_0/J-Bm*wr_0/J-Tl/J);
    titati_1=wr_0+deltat*(wr_0);
    variables(1,1,j)=ia_0;
    variables(2,1,j)=wr_0;
    variables(3,1,j)=titat_0;
    variables(1,2,j)=iai_1;
    variables(2,2,j)=wri_1;
    variables(3,2,j)=titati_1;
end

for j=1:1:10
    Tl=(j-1)*0.001;
    for i=3:1:(10000000-3)
       variables(1,i,j)=variables(1,i-1,j)+deltat*(-Ra*variables(1,i-1,j)/Laa-Km*variables(2,i-1,j)/Laa+Va/Laa);
       variables(2,i,j)=variables(2,i-1,j)+deltat*(Ki*variables(1,i-1,j)/J-Bm*variables(2,i-1,j)/J-Tl/J);
       variables(3,i,j)=variables(3,i-1,j)+deltat*variables(2,i-1,j);
    end
end

%gráfico de corrientes para distintos torques
figure
hold on;
for k=1:1:10
    plot(t, variables(1,:,k))
end

figure
plot(t, variables(1,:,1))

%Para un torque de 0.003Nm el motor tiene un pico de ~450mA, se considera
%aceptable