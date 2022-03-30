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
deltat=10^-6;
ts=5;

Va=zeros(1,ts/deltat);
for i=0.5/deltat:1:ts/deltat
    Va(i)=12;
end

cant_torques=8;
Tl=zeros(ts/deltat,cant_torques);
for j=1:1:cant_torques
    for i=2/deltat:1:ts/deltat
        Tl(i,j)=0+j*0.0002;
    end
end %Genero una matriz de distintos torques

t=0:deltat:(ts-deltat);
variables=zeros(3,ts/deltat,cant_torques);
ia_0=0; wr_0=0; titat_0=0; %Cond. iniciales

for j=1:1:cant_torques
    for i=2:1:(ts/deltat-1)
        variables(1,i,j)=variables(1,i-1,j)+deltat*(-Ra*variables(1,i-1,j)/Laa-Km*variables(2,i-1,j)/Laa+Va(1,i-1)/Laa);
        variables(2,i,j)=variables(2,i-1,j)+deltat*(Ki*variables(1,i-1,j)/J-Bm*variables(2,i-1,j)/J-Tl(i-1,j)/J);
        variables(3,i,j)=variables(3,i-1,j)+deltat*variables(2,i-1,j);
    end
end

figure
plot(t,variables(2,:,1));
hold on;
xlabel('Tiempo');
ylabel('wr');
title('Respuesta del motor con entrada de 12 v ante distintas perturbaciones');
a='\downarrow Torque: ';
b=num2str(1*0.0002);
c=[a b];
text(3, variables(2,3000000,1)+100, c);
for i=2:1:cant_torques
    plot(t,variables(2,:,i));
    b=num2str(i*0.0002);
    c=[a b];
    text(3, variables(2,3000000,i)+100, c);
end


figure
plot(t,variables(1,:,1));
hold on;
xlabel('Tiempo');
ylabel('ia');
title('Respuesta del motor con entrada de 12 v ante distintas perturbaciones');
a='\uparrow Torque: ';
b=num2str(1*0.0002);
c=[a b];
text(3, variables(1,3000000,1)-0.01, c);
for i=2:1:cant_torques
    plot(t,variables(1,:,i));
    b=num2str(i*0.0002);
    c=[a b];
    text(3, variables(1,3000000,i)-0.01, c);
end  

% iai_1=ia_0+deltat*(-Ra*ia_0/Laa-Km*wr_0/Laa+Va/Laa);
% wri_1=wr_0+deltat*(Ki*ia_0/J-Bm*wr_0/J-Tl/J);
% titati_1=wr_0+deltat*(wr_0);


% for j=1:1:10
%     Tl=(j-1)*0.001;
%     iai_1=ia_0+deltat*(-Ra*ia_0/Laa-Km*wr_0/Laa+Va/Laa);
%     wri_1=wr_0+deltat*(Ki*ia_0/J-Bm*wr_0/J-Tl/J);
%     titati_1=wr_0+deltat*(wr_0);
%     variables(1,1,j)=ia_0;
%     variables(2,1,j)=wr_0;
%     variables(3,1,j)=titat_0;
%     variables(1,2,j)=iai_1;
%     variables(2,2,j)=wri_1;
%     variables(3,2,j)=titati_1;
% end
% 
% for j=1:1:10
%     Tl=(j-1)*0.001;
%     for i=3:1:(10000000-3)
%        variables(1,i,j)=variables(1,i-1,j)+deltat*(-Ra*variables(1,i-1,j)/Laa-Km*variables(2,i-1,j)/Laa+Va/Laa);
%        variables(2,i,j)=variables(2,i-1,j)+deltat*(Ki*variables(1,i-1,j)/J-Bm*variables(2,i-1,j)/J-Tl/J);
%        variables(3,i,j)=variables(3,i-1,j)+deltat*variables(2,i-1,j);
%     end
% end
% 
% %gráfico de corrientes para distintos torques
% figure
% hold on;
% for k=1:1:10
%     plot(t, variables(1,:,k))
% end
% 
% figure
% plot(t, variables(1,:,1))
% 
% %Para un torque de 0.003Nm el motor tiene un pico de ~450mA, se considera
% %aceptable

