clc, clear all, close all;

%Definicion del sistema RLC serie en variables de estado
%Variables de estado: x1=i, x2=vc

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas

-------------------------------------------------------------------------
%}

%Importo gráficos de mediciones y grafico
mediciones=xlsread('Curvas_Medidas_RLC.xls');
%Columnas: Tiempo Corriente Vcap

%Definicion de la entrada
u=zeros(1,1000);
paso=0.1/1000;
t=0:paso:(0.1-paso);

signo=true;
for i=100:1:1000
    if mod(i,500)==0
       signo=not(signo);
    end
    if signo==1
        u(1,i)=12;
    end
    if signo==0
        u(1,i)=-12;
    end
end
figure
plot(t,u)

%Obtención de R L y C

%ALGORITMO DE CHEN para aproximación de FT de la forma
%G(s)=K*(T3*s+1)/((T1*s+1)*(T2*s+1))

k1=1.429467279/12-1; %t1=0.011 y(t1)=1,429467279 TENER EN CUENTA EL DELAY!
k2=3.467338816/12-1; %2t1=0.012 y(2t1)=3,467338816
k3=5.20835075/12-1; %3t1=0.013 y(3t1)=5,20835075
b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta=(2*k1^3+3*k1*k2+k3-sqrt(b))/(sqrt(b));
T1=-0.001/log(alfa1);
T2=-0.001/log(alfa2);
T3=beta*(T1-T2)+T1; %No hay cero, no se usa

s=tf('s');
G=12/((T1*s+1)*(T2*s+1));
[num,den]=tfdata(G,'v');

%Se sabe que Vc/vin (s) = (1/(L*C)) / (s^2+R*s/L+1/(L*C))
den_norm=den/(den(1))
num_norm=num/(12*den(1))
L=0.1;
R=L*den_norm(2);
Cap=1/(L*den_norm(3));

figure
[yaprox,taprox]=lsim(G,u/12,t);
plot(0.011,1.429467279,'x');
hold on;
plot(0.012,3.467338816,'x');
plot(0.013,5.20835075,'x');
plot(taprox,yaprox);
title('Respuesta del voltaje del capacitor aproximada');
xlabel('Tiempo');
legend({'vc(t1)','vc(2t1)','vc(3t1','vc(t)'},'Location','southeast')

%verifico que el sistema está bien aproximado
figure
plot(mediciones(:,1),mediciones(:,3));
hold on;
plot(taprox,yaprox); %Son muy parecidos, la aproximación es válida
legend({'vc(t) real','vc(t) aprxoimada'},'Location','southeast')
title('Contraste de respuestas');
xlabel('Tiempo');

%se verifica la corriente?

%Matrices
A=[-R/L -1/L; 1/Cap 0];
B=[1/L; 0];
C=[1; 0]';
D=0;

%Definicion de la ecuación de estado y de salida (salida de corriente)
sys1=ss(A,B,C,D)
figure
[yout,yt]=lsim(sys1,u,t);
plot(yt,yout);
hold on;
plot(mediciones(:,1),mediciones(:,2)); %se verifica con los componentes obtenidos
legend({'i(t) aproximada','i(t) real'},'Location','southeast')
title('Contraste de respuestas');
xlabel('Tiempo');

