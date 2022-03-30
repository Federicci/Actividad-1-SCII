clc, clear all, close all;

%Definicion del sistema RLC serie en variables de estado
%Variables de estado: x1=i, x2=vc

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas
Mismo ejercicio pero cambiando los valores de las constantes
Los comentarios de esta sección son los realizados en el archivo
TP1_1_RLC_1
-------------------------------------------------------------------------
%}

%Constantes
L=0.00001; %[H]
C=0.0000001; %[F]
R=4700; %[ohm]

%Matrices
A=[-R/L -1/L; 1/C 0];
B=[1/L; 0];
C=[0; 1]';
D=0;

%Definicion de la ecuación de estado y de salida
sys=ss(A,B,C,D)

%Definicion de la entrada
u=zeros(1,1000);
paso=0.01/1000;
t=0:paso:(0.01-paso);

signo=false;
for i=100:1:1000
    if mod(i,100)==0
       signo=not(signo);
    end
    if signo==1
        u(1,i)=12;
    end
    if signo==0
        u(1,i)=-12;
    end
end
%plot(t,u)

%Simulación
figure
lsim(sys,u,t);
title('Voltaje en el capacitor para 12v fluctuantes');
legend({'vc(t)'},'Location','southeast')
%Se observa la carga del capacitor que no llega a completarse

%Para ver la corriente:
C=[1; 0]';
sys=ss(A,B,C,D);
figure
[ysim tsim]=lsim(sys,u,t); %el orden de magnitud de la corriente es mucho
%menor al voltaje de entrada, grafico por separado
plot(tsim,ysim);
legend({'i(t)'},'Location','southeast')
title('Corriente para 12v fluctuantes');
xlabel('Tiempo');
%se observa la caida exponencial de la corriente debido a la carga del
%capacitor 

