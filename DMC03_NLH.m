%% NLH-DMC - Nonlinear Hammerstein Dynamic Matrix Control


clc
clear all
close all
%% Modelo  do  Sistema
T=1;                          %Tempo da Amostragem

%Processo nao-linear 1 (KATAYAMA; YAMAMOTO, 2004):
d=0;
num=[0 1.2 -0.1];
B=[zeros(1,d) num];
A=[1 -0.6 0.1];
ftz=tf(B,A,T);
ftz.iodelay=d;
gp=ftz;
ftz2=tf([0 0.6 -0.1],A,T);
ftz2.iodelay=d;

%% Define os ajustes do Controle Preditivo

gpi=stepinfo(gp);       %Armazeno informacao da planta
P=35;                   %Horizonte de predição (2 constantes de tiempo puede ser)
N=5;                    %Horizonte de controle (5 - 10)
delta=1*eye(P);         %Matriz de ponderacion del error
lamda=10;                %(Grande control é pequeno; Pequeno o control é livre) 

Ql=eye(N)*lamda;
Qd=eye(P)*delta;

%% Obter o vetor Gi
%Armazeno a respota ao degrau.
gi=step(ftz,P+d);  
Nm=length(gi)-1;        %Numero de amostras que vou pegar do degrau (Pontos finitos)

%% Calculo da Matriz G
G=zeros(P,N);
G(:,1)=gi(1+d:P+d);  %encho a primeira coluna com os valores
                 %da resposta ao degrau até horizonte de
                 %predição, quitando o retrasso
 
 for i=2:N
     for j=2:P
        G(j,i)=G(j-1,i-1); 
     end
 end

 %% calcula matriz Mn Matriz Ganho
% Mn=inv(G'*delta*G + lamda*eye(N))*G'*delta;
Mn=inv(G'*Qd*G+Ql)*G'*Qd; %Calculo de la Funcion de Costo sin Restriccion


%Calculo do controldor K1 (Primeira linha de Mn)
K1=Mn(1,:);

%% Loop de Controle

% inicializa parametros
nit=140;
inc_x=0;
u_ant(1:10) = 0;
u(1:nit) = 0; ym(1:nit) = 0; r(1:nit)=0; x(1:nit)=0; uf(1:nit)=0;
xt(1:nit)=0;
% Referência
r(10:70) = 1; r(71:nit)=1.2;
do(1:179)= 0;do(180:nit) = 0.1;

dx=zeros(1,Nm); %Acao de controle livre (Delta U Free)

w=0; %Termino para colocar referencias futuras     
for k=2:nit-w
     % calcula salida proceso 
      %x=1.5.*[0 u(1:k-1)]-1.5.*[0 u(1:k-1)].^2+0.5.*[0 u(1:k-1)].^3;
      x(k)=1.5*u(k-1)-1.5*u(k-1)^2+0.5*u(k-1)^3;
      t = 0:T:(k-1)*T;
      ym=lsim(ftz,x(1:k),t,'zoh')';     
      

% perturbacion deterministica a la salida
    if k>120
     ymr(k)=ym(k)+0.1;
    else
     ymr(k)=ym(k);
    end
     
     %% CALCULO DA RESPOSTA LIVRE
     f=zeros(P,1); % Vetor f (free) Resposta livre
     for kk=1:P
    % monta un vector con las gkk+i - gkk
         for i=1:Nm-P
             vect_g(i)=gi(kk+i)-gi(i);
         end
         for i=Nm-P+1:Nm
             vect_g(i)=gi(Nm)-gi(i);
         end
         f(kk)=ymr(k)+vect_g*dx';  %Calculo da resposta livre
         %f= Vector de respuesta libre con tamaño P
         %duf= (du libre) es la u correspondiente a la respuesta libre
               %ese vector siempre esta en el pasado. Es cero en el futuro
               %y vale unicamente en el pasado
         %ym= Salida de la planta
     
    end   

     %Calculo do Controle
     %Projeto onde nao tenho as referencias futuras
     inc_x=K1*(r(k+w)*ones(1,P)'-f);
     %actualiza dx
     aux_2=dx(1:Nm-1);
     dx=[inc_x aux_2];
    
     %Encuentra la Pseudo-Salida x(k)
     if k==1
        xt(k)=inc_x;
     else
        xt(k)=xt(k-1)+ inc_x;
     end
     
     %Crea el Polinomio de tercer orden x(k)
     xp=[0.5 -1.5 1.5 -xt(k)];
     %Busca las Raices para encontrar la u(k)
     ur=roots(xp);
     %Busca la Raiz que sea REAL
     ww=find(imag(ur)==0);
     %Aplica esa Raiz al proceso
     u(k)=real(ur(ww));
     
     %Satura el Controlador
     if u(k)>2.5
         u(k)=2.5;
     end
     if u(k)<-2.5
         u(k)=-2.5;
     end
    
end

nm=nit;
t = 0:T:(nm-1)*T;
figure
subplot(2,1,1),plot(t(1:nit-w),r(1:nit-w),'--k',t(1:nit-w),ymr,'-r','Linewidth',2)
xlabel('Tempo (s)');
ylabel('Saida');
legend('y_r','y','Location','SouthEast')
grid on;
hold
subplot(2,1,2),plot(t(1:nit-w),u,'b','Linewidth',3)
xlabel('Tempo (s)');
ylabel('Controle');
legend('u')
grid on;

