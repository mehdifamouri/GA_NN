clc
clear all
format long eng
%========== GENETIC  PARAMETER =====================
NS=40;
NG=100;
NP=2;
zigma=0;
q=0;

if (NP == 6);NN=30;end;
if (NP == 4);NN=15;end;
if (NP == 2);NN=1;end;
%NN=1
if (NP == 6);Rc=.95;end;
if (NP == 4);Rc=.93;end;
if (NP == 2);Rc=.9;end;
PR=.05;
PM=.12;
PC=1.;

if (NP==2)
    AU(1)=100;AL(1)=0;
    AU(2)=1000;AL(2)=0;
 end
if (NP==4)
    AU(1)=100;AL(1)=0;
    AU(2)=1;AL(2)=-1;
    AU(3)=1000;AL(3)=0;
    AU(4)=10;AL(4)=-10;
end
if (NP==6)
    AU(1)=100;AL(1)=0;
    AU(2)=1;AL(2)=-1;        %.0166
    AU(3)=0.01;AL(3)=-0.01;  %-0.000021
    AU(4)=1000;AL(4)=0;
    AU(5)=10;AL(5)=-10;      %.291
    AU(6)=0.1;AL(6)=-0.1;    %-0.0005597
end
NSPR=round(NS*PR);
%=====================================
LP1=1;LP2=1;
if (q == 0)
    DATACONDUC(1)=14.1;         % K1
    DATACONDUC(2)=0;      % K2
    DATACONDUC(3)=0;           % K3
    DATACONDUC(4)=448.;         % CP1
    DATACONDUC(5)=0;       % CP2
    DATACONDUC(6)=0;           % CP3
end
if (q == 1)
    DATACONDUC(1)=14.1;         % K1
    DATACONDUC(2)=.0166;      % K2
    DATACONDUC(3)=0;           % K3
    DATACONDUC(4)=448.;         % CP1
    DATACONDUC(5)=.291;       % CP2
    DATACONDUC(6)=0;           % CP3
end
if (q == 2)
    DATACONDUC(1)=14.1;         % K1
    DATACONDUC(2)=0.0200;      % K2
    DATACONDUC(3)=-0.000021;           % K3
    DATACONDUC(4)=448.;         % CP1
    DATACONDUC(5)=0.5455;       % CP2
    DATACONDUC(6)=-0.0005597;           % CP3
end

% DATACONDUC(1)=14.1;         % K1
% DATACONDUC(2)=.0166*LP1;      % K2
% DATACONDUC(3)=-0.000021*LP2;           % K3
% DATACONDUC(4)=448.;         % CP1
% DATACONDUC(5)=.291*LP1;       % CP2
% DATACONDUC(6)=-0.0005597*LP2;           % CP3
DATACONDUC(7)=6288.25;        % Ro
DATACONDUC(8)=23.;         % T Initial
DATACONDUC(9)=1000;         % NT
DATACONDUC(10)=41;          % NX
DATACONDUC(11)=41;          % Ny
DATACONDUC(12)=2;        % DT
DATACONDUC(13)=100*10^3;    %q
DATACONDUC(14)=.3;           % th
DATACONDUC(15)=.1;           % Lx
DATACONDUC(16)=.1;           % Ly
DATACONDUC(17)=6;           % Error
DATACONDUC(18)=DATACONDUC(13)*DATACONDUC(15)/DATACONDUC(1);
NT=DATACONDUC(9);
NX=DATACONDUC(10);
NY=DATACONDUC(11);
TIME=DATACONDUC(12)*(NT-1.);
TIMEPLUS=TIME*(DATACONDUC(1)/(DATACONDUC(7)*DATACONDUC(4)))/(DATACONDUC(15)^2);
tic
T=DIRECTSOLUTION(DATACONDUC);
TE=REFRENCESOLUTION(T,DATACONDUC,zigma);
toc
DATAPRINT(1,NP+5+1)=toc;
Tmax=max(max(max(TE)));
ij=0;

%GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
A=INTITALPOPULATION(TE,DATACONDUC,NS,NP,NN,AU,AL);
toc
DATAPRINT(2,NP+5+1)=toc;
for K=1:NG
    A=REALTIVEFITTNESS(A,NS,NP);
    A=TOURNOMENTSELECTION(NS,NP,A);
    ANEW=CROOSOVER(A,NS,NP,PC,AU,AL);
    ANEW=MOTTATION(ANEW,AU,AL,PM,NS,NP);
    ANEW=FITTNESSCACULATION(TE,DATACONDUC,ANEW,NS,NSPR,NP,AU,AL);
    ANEW=FINALSELECTION(A,ANEW,NS,NP,NSPR);
    %==========================================================================
    q=0;q(1:NP)=A(1,1:NP);MAX=A(1,NP+3);AVE=A(2,NP+3); q=[K,q,AVE,MAX,(1/MAX-.001)^2, toc];q=q'
    A=0;A=ANEW;
    q=q';
    DATAPRINT(K,1:NP+5)=q(1:NP+5);
    %======================================================================
    if (K > 2 )
        AU(1:NP)=(1-Rc)*A(1,1:NP)+Rc*AU(1:NP);
        AL(1:NP)=(1-Rc)*A(1,1:NP)+Rc*AL(1:NP);
        qq=[AU;AL]
    end
    save
    toc
    %======================================================================
end
%GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
%AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA