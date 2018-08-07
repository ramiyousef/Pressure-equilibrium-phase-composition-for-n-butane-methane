function [x,y,check]=flash(P,z,T)
z=z(:);
R=8.314;%pa.m^3/(mol.k)
Tc=[ 190.6; 425.1]; %Tc(1) is n-butae Tc(2) is methane
Pc=[ 4.599e6;3.796e6 ];
w=[0.012;0.2]; %butane and methane respectively
ac=[0.42747*R^2*Tc(1)^2/Pc(1); 0.42747*R^2*Tc(2)^2/Pc(2)]; %equation 5 in paper
b=[0.08664*R*Tc(1)/Pc(1); 0.08664*R*Tc(2)/Pc(2)]; %equation 6 in paper
Tr=[T/Tc(1) ; T/Tc(2)];%finding the reduced temperature at 100 fehrenheit
mbutane=0.480+1.574*w(1)-0.176*w(1)^2; %calculating the m parameter to find alpha that will be used to find ai(T)
mmethane=0.480+1.574*w(2)-0.176*w(2)^2;
alphabutane=(1+mbutane*(1-sqrt(Tr(1))))^2;%calculating alpha
alphamethane=(1+mmethane*(1-sqrt(Tr(2))))^2;
alpha=[alphabutane  ;  alphamethane];
ai =[alpha(1)*ac(1) ;  alpha(2)*ac(2)];
kappa=[0 0;0 0];
%% initial guess of K equilibrium constants using the wilson relation
K=[exp(log(Pc(1)/P)+5.37*(1+w(1))*(1-Tc(1)/T)); exp(log(Pc(2)/P)+5.37*(1+w(2))*(1-Tc(2)/T))];
kmax=100;
k=1;
while k<kmax %flash calculations are performed inside the loop
    check=0;
    %% calculating P at alpha=1 & alha=0
    
    Pzero=z(1)*(K(1)-1)+z(2)*(K(2)-1); % @alpha=0
    
    Pone=1-(z(1)/K(1)+z(2)/K(2)); % @alpha=1
    
    %% check if there are any roots between zero and 1
    if Pzero*Pone<=0
        NPhase=2;
        check=1;
        [alphastar,nrfailnewton] = newtonraphsontermproject(0,1,z,K,2);
        if nrfailnewton==1
            %              display('newton raphson failed, bisection used instead')
            [alphastar,nrfailbisection]=bisection_termproject(@Palpha,z,K);
            if nrfailbisection==1
                %                 display('bisection  failed');
                break
            end
        end
    else
        %if p(0) and P(1) are negative, alpha* is negative and feed is all liquid
        if Pzero<=0
            alphastar=0;
            NPhase=1;
        else
            if Pone>=0
                alphastar=1;
                NPhase=1;
            else
                error('No results can be obtained! check your function')
            end
        end
        
    end
    %% calculating x and y
    x=z./(1+alphastar.*(K-1));
    y=K.*x;
    %% normalizing x and y
    x=x./(x(1)+x(2));
    y=y./(y(1)+y(2));
    %% calculating bv and bl and av and al
    bvap=y(1)*b(1)+y(2)*b(2);
    bliq=x(1)*b(1)+x(2)*b(2);
    avap=0; aliq=0;
    %% calculating avapor and aliquid
    for i=1:2
        for j=1:2
            avap=avap+y(i)*y(j)*sqrt(ai(i))*sqrt(ai(j))*(1-kappa(i,j));
            aliq=aliq+x(i)*x(j)*sqrt(ai(i))*sqrt(ai(j))*(1-kappa(i,j));
        end
    end
    %% calculating capital A and B for each phase to find the Z-factor
    Al=0.42747*P/T^2*(x(1)*Tc(1)*sqrt(alpha(1)/Pc(1)) +x(2)*Tc(2)*sqrt(alpha(2)/Pc(2)))^2;
    Bl=0.08664*P/T*(x(1)*Tc(1)/Pc(1)+x(2)*Tc(2)/Pc(2));
    Av=0.42747*P/T^2*(y(1)*Tc(1)*sqrt(alpha(1)/Pc(1)) +y(2)*Tc(2)*sqrt(alpha(2)/Pc(2)))^2;
    Bv=0.08664*P/T*(y(1)*Tc(1)/Pc(1)+y(2)*Tc(2)/Pc(2));
    cvap=[1 -1 (Av-Bv-Bv^2) (-Av*Bv)];
    cliq=[1 -1 (Al-Bl-Bl^2) (-Al*Bl)];
    %% finding the roots of the two Z- equations and taking the real part of them
    Zvap=roots(cvap);
    index=find(imag(Zvap)==0);
    if length(index)>1
        Zv=Zvap(index);
    else
        Zv=real(Zvap(index));%taking the real part of the roots
    end
    Zliq=roots(cliq);
    index=find(imag(Zliq)==0);
    if length(index)>1
        Zl=Zliq(index);
    else
        Zl=real(Zliq(index));%taking the real part of the roots
    end
    vliquid=min(Zl)*R*T/P; vvapor=max(Zv)*R*T/P; % calculating the molar volumes for the obtained Z factors
    %% calculation of bi/b to be used in fugacity calculations
    bioverbliquid=[(Tc(1)/Pc(1))/(x(1)*Tc(1)/Pc(1)+x(2)*Tc(2)/Pc(2)) ; (Tc(2)/Pc(2))/(x(1)*Tc(1)/Pc(1)+x(2)*Tc(2)/Pc(2))];
    bioverbvapor=[(Tc(1)/Pc(1))/(y(1)*Tc(1)/Pc(1)+y(2)*Tc(2)/Pc(2)) ; (Tc(2)/Pc(2))/(y(1)*Tc(1)/Pc(1)+y(2)*Tc(2)/Pc(2))];
    sqrtaioveraliquid=[(sqrt(alpha(1))*Tc(1)/sqrt(Pc(1)))/(x(1)*sqrt(alpha(1))*Tc(1)/sqrt(Pc(1))+x(2)*sqrt(alpha(2))*Tc(2)/sqrt(Pc(2)))...
        ; (sqrt(alpha(2))*Tc(2)/sqrt(Pc(2)))/(x(1)*sqrt(alpha(1))*Tc(1)/sqrt(Pc(1))+x(2)*sqrt(alpha(2))*Tc(2)/sqrt(Pc(2)))];
    sqrtaioveravopor=[(sqrt(alpha(1))*Tc(1)/sqrt(Pc(1)))/(y(1)*sqrt(alpha(1))*Tc(1)/sqrt(Pc(1))+y(2)*sqrt(alpha(2))*Tc(2)/sqrt(Pc(2)))...
        ; (sqrt(alpha(2))*Tc(2)/sqrt(Pc(2)))/(y(1)*sqrt(alpha(1))*Tc(1)/sqrt(Pc(1))+y(2)*sqrt(alpha(2))*Tc(2)/sqrt(Pc(2)))];
    
    %% fugacity calculations f1v f2v f1l f2l using the previously found values
    
    fugacityvapor=P*y.*exp((bioverbvapor.*(max(Zv)-1))-log(max(Zv)-Bv)-Av/Bv.*(2.* sqrtaioveravopor-bioverbvapor).*log(1+Bv/max(Zv)));
    fugacityliquid=P*x.*exp((bioverbliquid.*(min(Zl)-1))-log(min(Zl)-Bl)-Al/Bl.*(2.* sqrtaioveraliquid-bioverbliquid).*log(1+Bl/min(Zl)));
    sum=0;
    for q=1:2
        sum =sum+y(q)*log(fugacityvapor(q)/fugacityliquid(q));
    end
    thetanew=sum;
    
    if abs(thetanew)<=1e-7
        %         display('flash converged')
        break
    else
        if k==1
            K=K.*(fugacityliquid./fugacityvapor);%updating the K factor
        elseif k~=1
            if abs(thetanew-thetaold)<=1e-8&&Nphase==1
                if thetanew>0
                    display('Feed is liquid like check if alpha=0')
                    break
                else
                    display('Feed is vapor like check if alpha=1')
                    break
                end
            else
                K=K.*(fugacityliquid./fugacityvapor);
            end
            
            if NPhase==2
                if k>=kmax
                    break
                else
                    k=k+1;
                end % return to begining of loop
            else
                K=K.*exp(thetanew);%updating K if 2 phases are present
                if k>=kmax
                    %                 display('100 iterations no convergence')
                    break
                else
                    k=k+1;
                end % return to begining of loop
            end
        end
    end
    thetaold=thetanew;
    
end %end of while loop
