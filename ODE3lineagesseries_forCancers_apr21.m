%%%%%%%Algoritm_Cellular_Automaton_Tumor%%%%%%%%%
%                                               %
%%%               %3D - MODEL%                %%%
%                                               %
%            vers. for Cancers (MDPI)           %
%                with lineages                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Authors:       Luca Meacci and Mario Primicerio
%University:    ICMC Sao Carlos - Universidade de Sao Paulo
%Date:          Apr. 26, 2021

%Legend
% 2   -> empty
% 0.7 -> CC newborn (1)
% 1.0 -> CC juvenile (2)
% 1.3 -> CC adult (3)
% 1.6 -> CC senescent (4)
% 0   -> CSC

%setup
n=50; N=n^3;
d=0.8;
rho=1;
rho1=[0.6,0.5,0.4]; 
rho2=[0.8,0.7,0.6]; 
rho3=[0.9,0.8,0.7]; 
rho4=[0.4,0.3,0.2];
gamma1=[1,1,1]; gamma2=[1,1,1]; gamma3=[1,1,1];
tSpan=[0 200];
y0 =[9/N,0/N,0/N,0/N,0/N,0/N,0/N,0/N,0/N,0/N,0/N,0/N,0/N,9/N,9/N]; %[CSC-L1, CC1-L1, CC2, CC3, CC4, ..., CSC-L2, CSC-L3]

xi=[0.6,1.5,0.6,0.8,1.0,1.5,2.0,2.5,3.0,3.5];
indexM=size(xi,2);
for index=1:indexM
mu1=xi(index)*[0.09,0.1,0.11]; mu2=xi(index)*[0.09,0.1,0.11]; mu3=xi(index)*[0.09,0.1,0.11]; mu4=xi(index)*[0.14,0.15,0.16];  
 
options=odeset('RelTol',1E-8,'AbsTol',1E-8,'InitialStep',1E-8,'MaxStep',10000000000);
[t,y]=ode45(@(t,y)Sisdif3lineages_forCancers(t,y,d,rho,rho1,rho2,rho3,rho4,mu1,mu2,mu3,mu4,gamma1,gamma2,gamma3),tSpan,y0,options);
%p=sum(y,2);

csc1=y(:,1);
cc11=y(:,2);
cc21=y(:,3);
cc31=y(:,4);
cc41=y(:,5);
cc12=y(:,6);
cc22=y(:,7);
cc32=y(:,8);
cc42=y(:,9);
cc13=y(:,10);
cc23=y(:,11);
cc33=y(:,12);
cc43=y(:,13);
csc2=y(:,14);
csc3=y(:,15);

p=csc1+cc11+cc21+cc31+cc41+cc12+cc22+cc32+cc42+cc13+cc23+cc33+cc43+csc2+csc3;

%Plotting
figure(1)
hold on
plot(t,p,'LineWidth',1.5)
set(gca,'FontSize',10)
xlabel('time')
title('p')
ylim([0 1])

EX=0;

if EX ==1
figure(2)
hold on
plot(t,csc1,'k','LineWidth',1.5)
plot(t,cc11,'b','LineWidth',1.5)
plot(t,cc21,'g','LineWidth',1.5)
plot(t,cc31,'y','LineWidth',1.5)
plot(t,cc41,'r','LineWidth',1.5)
xlabel('time')
ylabel('densities')
ylim([0 0.35])
legend('CSC-L1','CC1-L1','CC2-L1','CC3-L1','CC4-L1')

figure(3)
hold on
plot(t,csc2,'k','LineWidth',1.5)
plot(t,cc12,'b','LineWidth',1.5)
plot(t,cc22,'g','LineWidth',1.5)
plot(t,cc32,'y','LineWidth',1.5)
plot(t,cc42,'r','LineWidth',1.5)
xlabel('time')
ylabel('densities')
legend('CSC-L2','CC1-L2','CC2-L2','CC3-L2','CC4-L2')

figure(4)
hold on
plot(t,csc3,'k','LineWidth',1.5)
plot(t,cc13,'b','LineWidth',1.5)
plot(t,cc23,'g','LineWidth',1.5)
plot(t,cc33,'y','LineWidth',1.5)
plot(t,cc43,'r','LineWidth',1.5)
xlabel('time')
ylabel('densities')
legend('CSC-L3','CC1-L3','CC2-L3','CC3-L3','CC4-L3')
pause

end

end
legend('ξ=0.2','ξ=0.4','ξ=0.6','ξ=0.8','ξ=1.0','ξ=1.5','ξ=2.0','ξ=2.5','ξ=3.0','ξ=3.5')
% %Plotting
% figure(1)
% plot(t,p,'m','LineWidth',1.5)
% set(gca,'FontSize',12)
% xlabel('time')
% title('p')
% 
% figure(2)
% hold on
% plot(t,csc1,'k','LineWidth',1.5)
% plot(t,cc11,'b','LineWidth',1.5)
% plot(t,cc21,'g','LineWidth',1.5)
% plot(t,cc31,'y','LineWidth',1.5)
% plot(t,cc41,'r','LineWidth',1.5)
% xlabel('time')
% ylabel('densities')
% legend('CSC','CC1','CC2','CC3','CC4')
% hold off




%%%%



