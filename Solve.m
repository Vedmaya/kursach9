function [ F ] = Solve( X, Y )
    global e;
    global en;
    global I;
    global index;
    global myu;
    global bolcman;
    global T0;
    global n0; 
    global tau0;
    %----------------------------------------------------------------------
    Navagadr= 6.022*10^23;
    m=myu/Navagadr;
    %----------------------------------------------------------------------
    if(I>=1000)
        index=index+1;        
        T=Y(105)*T0;
        N=sum(Y(1:104))*n0;
        E=5/2*bolcman*T/m+n0/(m*N)*(sum(e*Y(1:69))+sum(en*Y(70:104)));
        display(T);
        display(X*tau0);
        display(N);
        display(E);
        display(min(Y(1:104)));
        saveIndex(index);
        saveData('T',T);
        saveData('time',X*tau0);
        saveData('N',N);
        saveData('E',E);
        saveData('Y',Y);
        if(min(Y(1:104))<0)
            error('kapec');
        end;
        I=0; 
    end
    I=I+1;
    F=Rigth_hand(Y);
end

