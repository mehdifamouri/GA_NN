clc
load
AK=K+1
format long eng
for K=AK:NG
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