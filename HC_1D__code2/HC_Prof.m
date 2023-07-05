clear variables, close all, clc
% Hydro-Chemical code for hydration of zero-porosity rock
% Stefan Schmalholz, 31.03.2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look up related stuff; approxmation of LU variables with tanh, and Pf as
% "reverse" ln-function of total density
% Look up table for atg + bru = fo + H2O reaction
load LOOK_UP_atg
Pf_lu       = P_vec*1e8/1e9;    % Pf in GPa
Pf_re       = 1.265;            % Reaction Pf in GPa
rhos_lp     = Rho_s(1);
rhos_hp     = Rho_s(end);
rhof        = mean(Rho_f);
xs_lp       = X_s_vec(1);
xs_hp       = X_s_vec(end);

% Fitting rhos and xs from LU with tanh; rhof is assumed constant
rw              = 0.1; % Width of step in tanh
rhos            = rhos_lp + (rhos_hp-rhos_lp) *(0.5*(tanh((Pf_lu-Pf_re)*2*pi/rw)+1));
xs              = xs_lp   + (xs_hp-xs_lp)     *(0.5*(tanh((Pf_lu-Pf_re)*2*pi/rw)+1));
rhostot         = min(rhos)*(1-0)*min(xs);
Phi_lu          = 1-(rhostot./(rhos.*xs));
rhoTOT_lu       = (1-Phi_lu).*rhos + Phi_lu.*rhof;
rT_lp           = rhoTOT_lu(1);
rT_hp           = rhoTOT_lu(end);
rhoTOT          = rT_lp   + (rT_hp-rT_lp)     *(0.5*(tanh((Pf_lu-Pf_re)*2*pi/rw)+1));
Pf_reverse      = real( (rw.*log((rhoTOT-rT_lp)./(-rhoTOT+rT_hp))+4*pi.*Pf_re)./4/pi );
figure(1)
subplot(211)
hold on
% plot(Pf_lu,lup(:,[3]),'oc')
% plot(Pf_lu,lup(:,[2]),'or')
plot(Pf_lu,rhos,'-b')
% plot(Pf_lu,lup(:,[4]),'xc')
plot(Pf_lu,rhof*ones(size(Pf_lu)),'-r')
% plot(Pf_lu,rhoTOT_lu,'-xk')
plot(Pf_lu,rhoTOT,'-k')
plot(Pf_reverse,rhoTOT,'ok')
legend('\rho_s','\rho_f','\rho_T','location','southeast'), xlabel('P_f'), ylabel('\rho')
title('Look up and fitting: densities'),
% axis([1.7 2.2 1000 3400]), grid on
subplot(212)
plot(Pf_lu,xs,'-b'), xlabel('P_f'), ylabel('X_s'),
% axis([1.7 2.2 0.99 1.001]), grid on
title('Look up and fitting: solid mass fraction, X_s')
set(gcf,'Position',[10.6000 193.8000 539.2000 580])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of HC model, configuration and intialization
% Physics
Lx              = 1;
betaf           = 1e11/1e9;
k0_etaf         = 1e-19/(1e-3/1e9);
n_poro          = 2;
% Numerics
nx              = 100;
dx              = Lx/(nx-1);
X               = 0:dx:Lx;
nt              = 1e5;
CFL             = 2e-2/betaf/1;
rel             = 0.25;
% Initial values, in detail. First define ambient zero-porosity values and
% then add pertubation at boundary assuming total solid mass is constant.
np                      = 1;
Phi_amb                 = zeros(size(X));
Delta_Pf                = 0.1;
Pf_amb                  = Pf_re*ones(size(X)) + Delta_Pf/2;
rhos_amb                = rhos_lp + (rhos_hp-rhos_lp)*(0.5*(tanh( (Pf_amb-Pf_re)*2*pi/rw ) + 1));
xs_amb                  = xs_lp + (xs_hp-xs_lp)*(0.5*(tanh( (Pf_amb-Pf_re)*2*pi/rw ) + 1));
rhostot                 = max( rhos_amb.*(1-Phi_amb).*xs_amb );
Pf_amb(1:np)            = Pf_amb(1) - Delta_Pf;
Pf_ini                  = Pf_amb;
rhos_ini                = rhos_lp + (rhos_hp-rhos_lp)*(0.5*(tanh( (Pf_ini-Pf_re)*2*pi/rw ) + 1));
xs_ini                  = xs_lp + (xs_hp-xs_lp)*(0.5*(tanh( (Pf_ini-Pf_re)*2*pi/rw ) + 1));
Phi_ini                 = 1-(rhostot./(rhos_ini.*xs_ini));
phi0                    = max(Phi_ini);
rhoTOT                  = (1-Phi_ini).*rhos_ini + Phi_ini.*rhof;
rhoTOT_ini              = rhoTOT;
Phi                     = Phi_ini;
rhos                    = rhos_ini;
Pf_reverse              = real( (rw.*log((rhoTOT-rT_lp)./(-rhoTOT+rT_hp))+4*pi.*Pf_re)./4/pi );
Pf                      = Pf_reverse;
Pf_ini                  = Pf_reverse;
% figure(2)
% subplot(211)
% plot(X,Pf_ini,'-k'), hold on
% plot(X,Pf_reverse,'o')
% plot([X(1) X(end)],[1 1]*Pf_re,'--xr')
% subplot(212)
% plot(X,rhoTOT,'-k')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time loop and solver.
Time                    = [0];
Dehy_front              = [ ];
Max_rhoTOT              = [ ];
for it = 1:nt
    Dehy_front(it)      = X(max(find(Pf<Pf_re)));
    Max_rhoTOT(it)      = max(real(rhoTOT));
    rhoTOT_old          = rhoTOT;
    Diffmax             = max(k0_etaf .* (Phi/phi0).^n_poro) ./ betaf;
    dt                  = CFL*dx^2*min(rhof)/max(real(rhoTOT))/Diffmax;
    if it==1;dt_ini = dt; end
    error_it            = 1;
    iter                = 0;
    while error_it > 1e-6
        iter            = iter+1;
        rhoTOT_it       = rhoTOT;
        Perm            = rhof*k0_etaf .* (Phi/phi0).^n_poro;
        Flux            = (Perm(1:end-1) + Perm(2:end))/2 .* diff(Pf)/dx;
        rhoTOT(2:end-1) = real( rhoTOT_old(2:end-1) + dt*diff(Flux)/dx );
        error_rho       = max(abs(rhoTOT_it-rhoTOT));
        Pfr             = real( (rw.*log((rhoTOT-rT_lp)./(-rhoTOT+rT_hp))+4*pi.*Pf_re)./4/pi );
        Pf              = (1-rel)*Pf + rel*Pfr;
        Pf(1:np)        = Pf_ini(1:np);
        error_Pf        = max(abs(Pf-Pfr));
        error_it        = max([error_rho error_Pf]);
        rhos            = rhos_lp + (rhos_hp-rhos_lp)*(0.5*(tanh( (Pf-Pf_re)*2*pi/rw ) + 1));
        xs              = xs_lp   + (xs_hp-xs_lp)    *(0.5*(tanh( (Pf-Pf_re)*2*pi/rw ) + 1));
        Phi             = 1-rhostot./(rhos.*xs);
    end
    if it>1, Time = [Time Time(end)+dt]; end
    if it==1, Flux_ini = Flux; end
    if mod(it,10000)== 0 || it==1
        figure(3)
        subplot(411), hold off
        plot(X,Pf,'-ok',X,Pf_ini,'--k'),hold on
        plot([X(1) X(end)],[1 1]*Pf_re,'--r')
        legend('P_f','P_{f0}','P_{f reaction}'), set(gca,'FontSize',12)
        title(['Fluid pressure evolution. Time step: ',num2str(it),'. Iterations: ',num2str(iter)])
        ylabel('P_f'),xlabel('X'),axis([-0.01 Lx min(Pf_ini)*0.98 max(Pf_ini)*1.02]), grid on
        subplot(412)
        plot(X,Phi,'-ok',X,Phi_ini,'--k')
        legend('\phi','\phi_0'), set(gca,'FontSize',12)
        ylabel('\phi'),xlabel('X'),axis([-0.01 Lx -0.01 max(Phi_ini)*1.2]), grid on
        title('Porosity evolution')
        subplot(413)
        plot(X,rhoTOT,'-ok',X,rhoTOT_ini,'--k',X,rhos,'-ob',X,rhos_ini,'--b')
        legend('\rho_T','\rho_{T 0}','\rho_s','\rho_{s 0}'), set(gca,'FontSize',12)
        ylabel('\rho'),xlabel('X'),axis([-0.01 Lx 2200 3300]), grid on
        title(['Density evolution. Maximal \rho_T: ',num2str(max(rhoTOT))])
        subplot(414)
        Phi_plot                = zeros(length(X),3);
        Phi_plot(:,3)           = Phi;
        Phi_plot(Pf>Pf_re,2)    = 1-Phi(Pf>Pf_re);
        Phi_plot(Pf<=Pf_re,1)   = 1-Phi(Pf<=Pf_re);
        ar = area(X,Phi_plot);ar(1).FaceColor=[140 56 28]/255;ar(2).FaceColor=[183 167 148]/255;ar(3).FaceColor=[30 143 195]/255;
        xlabel('X'),legend('Olivine', 'Serpentinite', 'H2O'), set(gca,'FontSize',12)
        title('Relative spatial distribution of phases')
        %         plot((X(1:end-1)+X(2:end))/2,-Flux,'-k',(X(1:end-1)+X(2:end))/2,-Flux_ini,'--k')
        %         legend('q','q_{0}'), set(gca,'FontSize',12)
        %         ylabel('Total mass flux'),xlabel('X'),axis([-0.01 Lx -50 500]), grid on
        %         title('Total mass flux evolution')
        set(gcf,'Position',[553.8000 50.6000 968.0000 722.4000])
        drawnow
    end
end

figure(4)
subplot(211)
plot(sqrt(Time),Dehy_front,'-xk')
xlabel('Time^{1/2}'), ylabel('Dehydration front')
subplot(212)
plot(Time,Max_rhoTOT,'-xk')
xlabel('Time'), ylabel('Maximal \rho_T')

