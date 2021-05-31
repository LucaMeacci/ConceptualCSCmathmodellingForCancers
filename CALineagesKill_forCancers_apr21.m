%%%%%%%Algoritm_Cellular_Automaton_Tumor%%%%%%%%%
%                                               %
%%%               %3D - MODEL%                %%%
%                                               %
%            vers. for Cancers (MDPI)           %
%                with lineages                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Authors:       Luca Meacci and Mario Primicerio
%University:    ICMC Sao Carlos - Universidade de Sao Paulo
%Date:          May 31, 2021

%Legend
% 2   -> empty
% 0.7 -> CC newborn (1)
% 1.0 -> CC juvenile (2)
% 1.3 -> CC adult (3)
% 1.6 -> CC senescent (4)
% 0   -> CSC

clear all

%Parameters setup
n=50;
N=n*n*n;
tMAX=250;

ALLdens_t=[]; CSCdens_t=[]; 
CC1dens_t=[]; CC2dens_t=[]; CC3dens_t=[]; CC4dens_t=[];

NL1dens_t=[]; NL2dens_t=[]; NL3dens_t=[];

mh=10;

delta=0; %dormancy
delta1=0; delta2=0; delta3=0; delta4=0;

rho=[mh*0.1, mh*0.1, mh*0.1, mh*0.1]; %replication
rho1=[mh*0.08, mh*0.05, mh*0.05]; 
rho2=[mh*0.10, mh*0.07, mh*0.07]; 
rho3=[mh*0.11, mh*0.08, mh*0.08]; 
rho4=[mh*0.06, mh*0.03, mh*0.03]; 

d=0.8; %CSC -> CC (1)

fm=3.5;
mu1=[fm*mh*0.009, fm*mh*0.010, fm*mh*0.010, 2]; %mortality (1)
mu2=[fm*mh*0.009, fm*mh*0.010, fm*mh*0.010, 2]; %mortality (2)
mu3=[fm*mh*0.009, fm*mh*0.010, fm*mh*0.010, 2]; %mortality (3) 
mu4=[fm*mh*0.014, fm*mh*0.015, fm*mh*0.015, 2]; %mortality (4)

p1=mh*0.1; 
p2=mh*0.1; 
p3=mh*0.1; %probability age progression

%KillTime
tK=300;

%Lineage Matrix
L=zeros(n,n,n);

%Initialization of Matrix
casei=1;
MCC1=[]; MCC2=[]; MCC3=[]; MCC4=[];
MCSC=[];
if casei == 0
G=ones(n,n,n).*2; %empty
G(2:4,2:4,2:4)=0.7; %CC (1)
for i=2:4
    for j=2:4
        for l=2:4
MCC1=[MCC1;i,j,l];
        end
    end
end
G(2:4,2:4,2:4)=1.0; %CC (2)
for i=40:42
    for j=40:42
        for l=40:42
MCC2=[MCC2;i,j,l];
        end
    end
end
G(2:4,42:44,42:44)=1.3; %CC (3)
for i=2:4
    for j=42:44
        for l=42:44
MCC3=[MCC3;i,j,l];
        end
    end
end
G(32:34,2:4,32:34)=0.7; %CC 4)
for i=32:34
    for j=2:4
        for l=32:34
MCC4=[MCC4;i,j,l];
        end
    end
end
G(24:26,24:26,24:26)=0; %CCS
for i=24:26
    for j=24:26
        for l=24:26
MCSC=[MCSC;i,j,l];
        end
    end
end
else
G=ones(n,n,n).*2; %empty
G(24:26,24:26,24:26)=0; %CSC
L(24:26,24:26,24)=1; %Lineage 1
L(24:26,24:26,25)=2; %Lineage 2
L(24:26,24:26,26)=3; %Lineage 3
%MCC1=[25,25,25];
end

%Plotting intial condition
plot0=0; %flag
if plot0 == 1
t=0;
plot3([1,1,1,1,n,n,n,n],[1,1,n,n,1,1,n,n],[1,n,1,n,1,n,1,n],'Owhite'); %box
title(['t = 0']);
hold on
if size(MCC1,1)>0
plot3(MCC1(:,1),MCC1(:,2),MCC1(:,3),'.b','MarkerSize',30);
end
if size(MCC2,1)>0
plot3(MCC2(:,1),MCC2(:,2),MCC2(:,3),'.g','MarkerSize',30);
end
if size(MCC3,1)>0
plot3(MCC3(:,1),MCC3(:,2),MCC3(:,3),'.y','MarkerSize',30);
end
if size(MCC4,1)>0
plot3(MCC4(:,1),MCC4(:,2),MCC4(:,3),'.r','MarkerSize',30);
end
if size(MCSC,1)>0
plot3(MCSC(:,1),MCSC(:,2),MCSC(:,3),'.k','MarkerSize',30);
end
saveas(gcf, ['CA_srnsh_time' ,num2str(t), '.png']);
pause
hold off
end

%Plotting initial lineages
plotLi=0;
if plotLi==1
t=0;
figure(4)
plot3([1,1,1,1,n,n,n,n],[1,1,n,n,1,1,n,n],[1,n,1,n,1,n,1,n],'Owhite'); %box
title(['t = ' num2str(t) ' ' ])
hold on
[LL1x,LL1yz]=find((~(L-1)));
LL1y=mod((LL1yz-1),n)+1;
LL1z=(LL1yz-LL1y)/n+1;
plot3(LL1x,LL1y,LL1z,'.m','MarkerSize',20)
%
[LL2x,LL2yz]=find((~(L-2)));
LL2y=mod((LL2yz-1),n)+1;
LL2z=(LL2yz-LL2y)/n+1;
plot3(LL2x,LL2y,LL2z,'.c','MarkerSize',20)
%
[LL3x,LL3yz]=find((~(L-3)));
LL3y=mod((LL3yz-1),n)+1;
LL3z=(LL3yz-LL3y)/n+1;
plot3(LL3x,LL3y,LL3z,'.','color',[0.9100    0.4100    0.1700],'MarkerSize',20);
hold off
pause
end

G_prec=G; %G prev. state
L_prec=L;
UV=[];
for t=1:tMAX
    t
    %reset MCC and MCS
    MCC1=[]; MCC2=[]; MCC3=[]; MCC4=[];
    MCSC=[];
    L1=[]; L2=[]; L3=[];
    for i=1:n    
        for j=1:n
            for l=1:n
            if G_prec(i,j,l)==0 %CSC case
                if rand<(1-delta) % no dormancy condition
                    n_vuote=0; %reset counters
                    M_vuote=[];
                    for k=i-1:i+1
                        for v=j-1:j+1
                            for m=l-1:l+1
                            if k==0
                                k=n;
                            elseif k==n+1
                                k=1;
                            end
                            if v==0
                                v=n;
                            elseif v==n+1
                                v=1;
                            end
                            if m==0
                                m=n;
                            elseif m==n+1
                                m=1;
                            end
                            if G(k,v,m)==2
                                n_vuote=n_vuote+1;
                                M_vuote=[M_vuote;k,v,m];
                            end
                        end
                        end
                    end
                    Nlib=n_vuote/26;                    
                    if rand<(Nlib*rho(L_prec(i,j,l))) %replication process
                        xyz_r=M_vuote(floor(rand*n_vuote)+1,:);
                        x_r=xyz_r(1,1);
                        y_r=xyz_r(1,2);
                        z_r=xyz_r(1,3);
                        if rand<(1-d) 
                            G(x_r,y_r,z_r)=0;
                            L(x_r,y_r,z_r)=L_prec(i,j,l);
                            MCSC=[MCSC;x_r,y_r,z_r]; %write new MCSC
                            
                            G(i,j,l)=0;
                            L(i,j,l)=L_prec(i,j,l);
                            MCSC=[MCSC;i,j,l]; %write MCSC
                        else
                            G(x_r,y_r,z_r)=0.7;
                            L(x_r,y_r,z_r)=L_prec(i,j,l);
                            MCC1=[MCC1;x_r,y_r,z_r]; %write new MCC1
                            G(i,j,l)=0;
                            L(i,j,l)=L_prec(i,j,l);
                            MCSC=[MCSC;i,j,l]; %write MCSC 
                        end
                    end
                else
                            G(i,j,l)=0;
                            MCSC=[MCSC;i,j,l]; %write MCSC
                            L(i,j,l)=L_prec(i,j,l);
                end
            elseif G_prec(i,j,l)==0.7 %CC1 case
                if rand<mu1(L_prec(i,j,l)) %apoptosis condition
                    G(i,j,l)=2; %CC1 dies
                    L(i,j,l)=0;
                else
                    if rand<(1-delta1) % no dormancy condition
                    n_vuote=0; 
                    M_vuote=[];
                    for k=i-1:i+1
                        for v=j-1:j+1
                            for m=l-1:l+1
                            if k==0
                                k=n;
                            elseif k==n+1
                                k=1;
                            end
                            if v==0
                                v=n;
                            elseif v==n+1
                                v=1;
                            end
                            if m==0
                                m=n;
                            elseif m==n+1
                                m=1;
                            end
                            if G(k,v,m)==2
                                n_vuote=n_vuote+1;
                                M_vuote=[M_vuote;k,v,m];
                            end
                        end
                        end
                    end
                    Nlib=n_vuote/26;
                    if rand<(Nlib*rho1(L_prec(i,j,l))) %replication process CC1
                        xyz_r=M_vuote(floor(rand*n_vuote)+1,:);
                        x_r=xyz_r(1,1);
                        y_r=xyz_r(1,2);
                        z_r=xyz_r(1,3);
                        G(x_r,y_r,z_r)=0.7;
                        L(x_r,y_r,z_r)=L_prec(i,j,l);
                        MCC1=[MCC1;x_r,y_r,z_r]; %write new MCC1
                        G(i,j,l)=0.7;
                        L(i,j,l)=L_prec(i,j,l);
                        MCC1=[MCC1;i,j,l]; %write MCC1
                    else % progression age stage
                        if rand<p1
                        G(i,j,l)=1;
                        L(i,j,l)=L_prec(i,j,l);
                        MCC2=[MCC2;i,j,l]; %write MCC2
                        else
                        G(i,j,l)=0.7;
                        L(i,j,l)=L_prec(i,j,l);
                        MCC1=[MCC1;i,j,l]; %write MCC1
                        end
                    end
                    else
                        G(i,j,l)=0.7;
                        MCC1=[MCC1;i,j,l]; %write MCC1
                        L(i,j,l)=L_prec(i,j,l);
                    end
                end
                elseif G_prec(i,j,l)==1 %CC2 case
                if rand<mu2(L_prec(i,j,l)) %apoptosis condition
                    G(i,j,l)=2; %CC1 dies
                    L(i,j,l)=0;
                else
                    if rand<(1-delta2) % no dormancy condition
                    n_vuote=0; 
                    M_vuote=[];
                    for k=i-1:i+1
                        for v=j-1:j+1
                            for m=l-1:l+1
                            if k==0
                                k=n;
                            elseif k==n+1
                                k=1;
                            end
                            if v==0
                                v=n;
                            elseif v==n+1
                                v=1;
                            end
                            if m==0
                                m=n;
                            elseif m==n+1
                                m=1;
                            end
                            if G(k,v,m)==2
                                n_vuote=n_vuote+1;
                                M_vuote=[M_vuote;k,v,m];
                            end
                        end
                        end
                    end
                    Nlib=n_vuote/26;
                    if rand<(Nlib*rho2(L_prec(i,j,l))) %replication process CC2
                        xyz_r=M_vuote(floor(rand*n_vuote)+1,:);
                        x_r=xyz_r(1,1);
                        y_r=xyz_r(1,2);
                        z_r=xyz_r(1,3);
                        G(x_r,y_r,z_r)=0.7;
                        L(x_r,y_r,z_r)=L_prec(i,j,l);
                        MCC1=[MCC1;x_r,y_r,z_r]; %write new MCC1
                        G(i,j,l)=0.7;
                        L(i,j,l)=L_prec(i,j,l);
                        MCC1=[MCC1;i,j,l]; %write MCC1
                    else % progression age stage
                        if rand<p2
                        G(i,j,l)=1.3;
                        MCC3=[MCC3;i,j,l]; %write MCC3
                        L(i,j,l)=L_prec(i,j,l);
                        else
                        G(i,j,l)=1;
                        MCC2=[MCC2;i,j,l]; %write MCC2
                        L(i,j,l)=L_prec(i,j,l);
                        end
                    end
                    else
                        G(i,j,l)=1;
                        MCC2=[MCC2;i,j,l]; %write MCC2
                        L(i,j,l)=L_prec(i,j,l);
                    end
                end
            elseif G_prec(i,j,l)==1.3 %CC3 case
                if rand<mu3(L_prec(i,j,l)) %apoptosis condition
                    G(i,j,l)=2; %CC3 dies
                    L(i,j,l)=0;
                else
                    if rand<(1-delta3) % no dormancy condition
                    n_vuote=0; 
                    M_vuote=[];
                    for k=i-1:i+1
                        for v=j-1:j+1
                            for m=l-1:l+1
                            if k==0
                                k=n;
                            elseif k==n+1
                                k=1;
                            end
                            if v==0
                                v=n;
                            elseif v==n+1
                                v=1;
                            end
                            if m==0
                                m=n;
                            elseif m==n+1
                                m=1;
                            end
                            if G(k,v,m)==2
                                n_vuote=n_vuote+1;
                                M_vuote=[M_vuote;k,v,m];
                            end
                        end
                        end
                    end
                    Nlib=n_vuote/26;
                    if rand<(Nlib*rho3(L_prec(i,j,l))) %replication process CC3
                        xyz_r=M_vuote(floor(rand*n_vuote)+1,:);
                        x_r=xyz_r(1,1);
                        y_r=xyz_r(1,2);
                        z_r=xyz_r(1,3);
                        G(x_r,y_r,z_r)=0.7;
                        L(x_r,y_r,z_r)=L_prec(i,j,l);
                        MCC1=[MCC1;x_r,y_r,z_r]; %write new MCC1
                        G(i,j,l)=0.7;
                        L(i,j,l)=L_prec(i,j,l);
                        MCC1=[MCC1;i,j,l]; %write MCC1
                    else % progression age stage
                        if rand<p3
                        G(i,j,l)=1.6;
                        L(i,j,l)=L_prec(i,j,l);
                        MCC4=[MCC4;i,j,l]; %write MCC4
                        else
                        G(i,j,l)=1.3;
                        L(i,j,l)=L_prec(i,j,l);
                        MCC3=[MCC3;i,j,l]; %write MCC3
                        end
                    end
                    else
                        G(i,j,l)=1.3;
                        L(i,j,l)=L_prec(i,j,l);
                        MCC3=[MCC3;i,j,l]; %write MCC3
                    end
                end
            elseif G_prec(i,j,l)==1.6 %CC4 case
                if rand<mu4(L_prec(i,j,l)) %apoptosis condition
                    G(i,j,l)=2; %CC1 dies
                    L(i,j,l)=0;
                else
                    if rand<(1-delta4) % no dormancy condition
                    n_vuote=0; 
                    M_vuote=[];
                    for k=i-1:i+1
                        for v=j-1:j+1
                            for m=l-1:l+1
                            if k==0
                                k=n;
                            elseif k==n+1
                                k=1;
                            end
                            if v==0
                                v=n;
                            elseif v==n+1
                                v=1;
                            end
                            if m==0
                                m=n;
                            elseif m==n+1
                                m=1;
                            end
                            if G(k,v,m)==2
                                n_vuote=n_vuote+1;
                                M_vuote=[M_vuote;k,v,m];
                            end
                        end
                        end
                    end
                    Nlib=n_vuote/26;
                    if rand<(Nlib*rho4(L_prec(i,j,l))) %replication process CC4
                        xyz_r=M_vuote(floor(rand*n_vuote)+1,:);
                        x_r=xyz_r(1,1);
                        y_r=xyz_r(1,2);
                        z_r=xyz_r(1,3);
                        G(x_r,y_r,z_r)=0.7;
                        L(x_r,y_r,z_r)=L_prec(i,j,l);
                        MCC1=[MCC1;x_r,y_r,z_r]; %write new MCC1
                        G(i,j,l)=0.7;
                        L(i,j,l)=L_prec(i,j,l);
                        MCC1=[MCC1;i,j,l]; %write MCC1
                    else % no progression age stage
                        G(i,j,l)=1.6;
                        L(i,j,l)=L_prec(i,j,l);
                        MCC4=[MCC4;i,j,l]; %write MCC4
                    end 
                    else %dormancy
                    G(i,j,l)=1.6;
                    L(i,j,l)=L_prec(i,j,l);
                    MCC4=[MCC4;i,j,l]; %write MCC4
                    end
                end
            end
            end
            end
    end

%Plotting screenshots
plotS=0; %flag
if plotS == 1
figure(1)
plot3([1,1,1,1,n,n,n,n],[1,1,n,n,1,1,n,n],[1,n,1,n,1,n,1,n],'Owhite'); %box
title(['t = ' num2str(t) ' ' ])
hold on
if size(MCC1,1)>0
plot3(MCC1(:,1),MCC1(:,2),MCC1(:,3),'.b','MarkerSize',20);
end
if size(MCC2,1)>0
plot3(MCC2(:,1),MCC2(:,2),MCC2(:,3),'.g','MarkerSize',20);
end
if size(MCC3,1)>0
plot3(MCC3(:,1),MCC3(:,2),MCC3(:,3),'.y','MarkerSize',20);
end
if size(MCC4,1)>0
plot3(MCC4(:,1),MCC4(:,2),MCC4(:,3),'.r','MarkerSize',20);
end
if size(MCSC,1)>0
plot3(MCSC(:,1),MCSC(:,2),MCSC(:,3),'.k','MarkerSize',20);
end
%saveas(gcf, ['CA_srnsh_time' ,num2str(t), '.png']);
%pause(0.1)
hold off
end

%Plotting lineages
plotL=0;
if plotL==1
figure(4)
plot3([1,1,1,1,n,n,n,n],[1,1,n,n,1,1,n,n],[1,n,1,n,1,n,1,n],'Owhite'); %box
title(['t = ' num2str(t) ' ' ])
hold on
[LL1x,LL1yz]=find((~(L-1)));
LL1y=mod((LL1yz-1),n)+1;
LL1z=(LL1yz-LL1y)/n+1;
plot3(LL1x,LL1y,LL1z,'.m','MarkerSize',20);
%
[LL2x,LL2yz]=find((~(L-2)));
LL2y=mod((LL2yz-1),n)+1;
LL2z=(LL2yz-LL2y)/n+1;
plot3(LL2x,LL2y,LL2z,'.c','MarkerSize',20);
%
[LL3x,LL3yz]=find((~(L-3)));
LL3y=mod((LL3yz-1),n)+1;
LL3z=(LL3yz-LL3y)/n+1;
plot3(LL3x,LL3y,LL3z,'.','color',[0.9100    0.4100    0.1700],'MarkerSize',20);
pause(0.1)
hold off
end


%saving over time
ALLdens=size(find(G-2),1)/N; 
NCSC=size(find(G),1);
NCC1=size(find(G-0.7),1);
NCC2=size(find(G-1),1);
NCC3=size(find(G-1.3),1);
NCC4=size(find(G-1.6),1);
CSCdens=(N-(NCSC))/N;
CC1dens=(N-(NCC1))/N;
CC2dens=(N-(NCC2))/N;
CC3dens=(N-(NCC3))/N;
CC4dens=(N-(NCC4))/N;

NL1dens=size(find(~(L-1)),1)/N;
NL2dens=size(find(~(L-2)),1)/N;
NL3dens=size(find(~(L-3)),1)/N;

ALLdens_t=[ALLdens_t; ALLdens];
CSCdens_t=[CSCdens_t; CSCdens];
CC1dens_t=[CC1dens_t; CC1dens];
CC2dens_t=[CC2dens_t; CC2dens];
CC3dens_t=[CC3dens_t; CC3dens];
CC4dens_t=[CC4dens_t; CC4dens];

NL1dens_t=[NL1dens_t; NL1dens];
NL2dens_t=[NL2dens_t; NL2dens];
NL3dens_t=[NL3dens_t; NL3dens];

%save prev. time step
G_prec=G;
L_prec=L;

if t==(tK-1)
    lK=1; %lineage to kill
   L_prec(find(~(L_prec-lK)))=4; 
end
if t==(tK)
    lK=1; %lineage to kill
   L_prec(find(~(L_prec-4)))=lK; 
end

end

%Post processing plotting

figure(2)
plot(ALLdens_t,'m','LineWidth',1.5)
xlabel('time')
ylabel('density p')

figure(3)
hold on
plot(CSCdens_t,'k','LineWidth',1.5)
plot(CC1dens_t,'b','LineWidth',1.5)
plot(CC2dens_t,'g','LineWidth',1.5)
plot(CC3dens_t,'y','LineWidth',1.5)
plot(CC4dens_t,'r','LineWidth',1.5)
xlabel('time')
ylabel('densities')
legend('CSC','CC1','CC2','CC3','CC4')
hold off

figure(5)
hold on
plot(NL1dens_t,'m','LineWidth',1.5)
plot(NL2dens_t,'c','LineWidth',1.5)
plot(NL3dens_t,'color',[0.9100    0.4100    0.1700],'LineWidth',1.5)
xlabel('time')
ylabel('lineages')
legend('L1','L2','L3')
hold off



