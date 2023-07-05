clear variables, close all, clf
%% INPUT %%
%1)import of 10 simulations with different porosity exponent: 
load('vein_width_npow1','vein_width')
vein_width1=vein_width;
load('vein_width_npow2','vein_width')
vein_width2=vein_width;
load('vein_width_npow3','vein_width')
vein_width3=vein_width;
load('vein_width_npow4','vein_width')
vein_width4=vein_width;
load('vein_width_npow5','vein_width')
vein_width5=vein_width;
load('vein_width_npow6','vein_width')
vein_width6=vein_width;
load('vein_width_npow7','vein_width')
vein_width7=vein_width;
load('vein_width_npow8','vein_width')
vein_width8=vein_width;
load('vein_width_npow9','vein_width')
vein_width9=vein_width;
load('vein_width_npow10','vein_width')
vein_width10=vein_width;

load('t_vein1','t_vein')
t_vein1=t_vein;
load('t_vein2','t_vein')
t_vein2=t_vein;
load('t_vein3','t_vein')
t_vein3=t_vein;
load('t_vein4','t_vein')
t_vein4=t_vein;
load('t_vein5','t_vein')
t_vein5=t_vein;
load('t_vein6','t_vein')
t_vein6=t_vein;
load('t_vein7','t_vein')
t_vein7=t_vein;
load('t_vein8','t_vein')
t_vein8=t_vein;
load('t_vein9','t_vein')
t_vein9=t_vein;
load('t_vein10','t_vein')
t_vein10=t_vein;

%2)import of 10 simulations with different porosity exponent: 
load('vein_width_npow1_2','vein_width')
vein_width1_2=vein_width;
load('vein_width_npow2_2','vein_width')
vein_width2_2=vein_width;
load('vein_width_npow3_2','vein_width')
vein_width3_2=vein_width;
load('vein_width_npow4_2','vein_width')
vein_width4_2=vein_width;
load('vein_width_npow5_2','vein_width')
vein_width5_2=vein_width;
load('vein_width_npow6_2','vein_width')
vein_width6_2=vein_width;
load('vein_width_npow7_2','vein_width')
vein_width7_2=vein_width;
load('vein_width_npow8_2','vein_width')
vein_width8_2=vein_width;
load('vein_width_npow9_2','vein_width')
vein_width9_2=vein_width;
load('vein_width_npow10_2','vein_width')
vein_width10_2=vein_width;

load('t_vein1_2','t_vein')
t_vein1_2=t_vein;
load('t_vein2_2','t_vein')
t_vein2_2=t_vein;
load('t_vein3_2','t_vein')
t_vein3_2=t_vein;
load('t_vein4_2','t_vein')
t_vein4_2=t_vein;
load('t_vein5_2','t_vein')
t_vein5_2=t_vein;
load('t_vein6_2','t_vein')
t_vein6_2=t_vein;
load('t_vein7_2','t_vein')
t_vein7_2=t_vein;
load('t_vein8_2','t_vein')
t_vein8_2=t_vein;
load('t_vein9_2','t_vein')
t_vein9_2=t_vein;
load('t_vein10_2','t_vein')
t_vein10_2=t_vein;

Pmetam     = 12.65e8; % Reaction pressure
Pfamb      = 12.75e8; % Ambient pressure
Pamp       = 5e8; % Pressure perturbation amplitude
wini       = 0.01;   % Initial perturbation width
k_etaf0    = 1e-19/1e-3;   % m^2/s/Pa - k = 1e-19 m2, etaf = 1e3 Pa.s
Lx         = 1;          % model length, m
npow       = 3;    % Porosity nonlinearity in permeability
% each of these simulations are done separetely by replacing npow  
% npow1      = 1;    % Porosity nonlinearity in permeability
% npow2      = 2;    % Porosity nonlinearity in permeability
% npow3      = 3;    % Porosity nonlinearity in permeability
% npow4      = 4;    % Porosity nonlinearity in permeability
% npow5      = 5;    % Porosity nonlinearity in permeability
% npow6      = 6;    % Porosity nonlinearity in permeability
% npow7      = 7;   % Porosity nonlinearity in permeability
% npow8      = 8;   % Porosity nonlinearity in permeability
% npow9      = 9;   % Porosity nonlinearity in permeability
% npow10     = 10;   % Porosity nonlinearity in permeability
%% Data from look-up table %%
load LOOK_UP_atg
Pf_lt       = P_vec*1e8;    % Corresponding fluid pressure array; converted to Pa
rhos_lt     = Rho_s;
rhof_lt     = Rho_f;
dP          = Pf_lt(2) - Pf_lt(1);
Pf_lt_min   = Pf_lt(1);
Pf_ltc      = ( Pf_lt(1:end-1) + Pf_lt(2:end))/2;
dPc         = Pf_ltc(2) - Pf_ltc(1);
Pf_ltc_min  = Pf_ltc(1);
ieclo_in    = find(rhos_lt < (max(rhos_lt)+min(rhos_lt))/2, 1, 'first');
Pmetam      = Pf_lt(ieclo_in);
Pf_eclo_in  = Pf_lt(ieclo_in);
Pf_gr_out   = Pf_lt(ieclo_in-1);
% Densities and beta
betaf_lt    = diff(log(rhof_lt))./diff(Pf_lt);
% Total solid mass
Xs_lt       = X_s_vec;
rhos_tot    = min(Xs_lt)*min(rhos_lt); % Ambient density of atg with zero porosity
phi_lt      = 1 - rhos_tot ./ rhos_lt ./ Xs_lt;
rhot_lt     = (1 - phi_lt) .* rhos_lt + phi_lt .* rhof_lt;
% Gridded interpolant function to replace reverse interpolation %
rhotpfGItp  = griddedInterpolant(rhot_lt,Pf_lt);

%% Model numerics %%
% numerics
CFL0        = 0.1;%2e-2/betaf/1;
niter       = 1e5;%3e5
rel0        = 0.1;
eps_err     = 1e-6;
nx          = 100;  % number of grid points

% preprocessing
dx       = Lx/(nx-1);     % grid spacing
x        = 0:dx:Lx;   % grid points cooarrdinates

%% Initialisation %%
% Fluid pressure %
Pf_bg    = Pfamb * ones(nx,1); % Background pressure without fluid pressure perturbation
Pf       = Pf_bg;
Pf(x<=wini) = Pf_bg(x<=wini) - Pamp;
% Pf    = Pf_bg   + Pamp*exp( -(x/wg).^2 )';
Pf0      = Pf;
% Density look up - Initial perturbation is fully eclogitised from the start
rhos    = Itp1D(Pf_lt, rhos_lt, Pf, zeros(nx,1), dP, Pf_lt_min);
rhof    = Itp1D(Pf_lt, rhof_lt, Pf, zeros(nx,1), dP, Pf_lt_min);
Xs      = Itp1D(Pf_lt, Xs_lt  , Pf, zeros(nx,1), dP, Pf_lt_min);
% Porosity
phi     = 1 - rhos_tot ./ rhos ./ Xs;
phi0    = max(phi);  
% Preallocation
beta       = zeros(nx, 1);
Pfr        = zeros(nx ,1);
rhog       = zeros(nx, 1);
rhoe       = zeros(nx, 1);
niters     = [];
err        = [];
vein_width = [];
time_it    = [];

%% Time and saving stuff %%
beta       = Itp1D(Pf_ltc, betaf_lt, Pf, beta, dPc, Pf_ltc_min);
time_sc    = Lx^2/k_etaf0*max(beta);
t_tot      = 2e2*time_sc;   % total time, s
t_pert     = t_tot;
time       = 0;
nout       = 2500;
nitp       = 10;
dp         = [];
%% Action %%
it = 0;
while time < t_tot,  it = it + 1;
    Pf_old = Pf;
    CFL      = CFL0;
    rel      = rel0;
    for iter = 1:niter
        if iter > 400, rel = 0.1; end
        errorPf = Pf;
        % Reverse lut for fluid pressure
        if it+iter > 2
            Pfr = rhotpfGItp(rhot);
            Pf = (1-rel)*Pf + rel*Pfr;
        end
        % Density look up and porosity
        if mod(iter, nitp) == 0 || iter == 1
            rhos = Itp1D(Pf_lt, rhos_lt, Pf, rhos, dP, Pf_lt_min);
            rhof = Itp1D(Pf_lt, rhof_lt, Pf, rhof, dP, Pf_lt_min);
        end
        Xs      = Itp1D(Pf_lt, Xs_lt  , Pf, zeros(nx,1), dP, Pf_lt_min);
        phi     = 1 - rhos_tot ./ rhos ./ Xs;
        beta    = Itp1D(Pf_ltc, betaf_lt, Pf, beta, dPc, Pf_ltc_min);
        if it   == 1, rhot     = (1 - phi) .* rhos + phi .* rhof; end          % forward
        if iter == 1, rhot_old = rhot; end

        % boundary conditions and averaging
        phiex        = [phi(1); phi; phi(end)];
        rhofex       = [rhof(1); rhof; rhof(end)];
        if time < t_pert
            Pfex     = [2*Pf0(1)-Pf(1); Pf; Pf(end)];
        else
            Pfex     = [Pf(1); Pf; Pf(end)];
        end
        phic         = (phiex(1:end-1) + phiex(2:end))/2;
        rhofc        = (rhofex(1:end-1) + rhofex(2:end))/2;

        % total mass conservation with darcy flux
        k_etaf       = k_etaf0*(phic/phi0).^npow;%
        drhot_dt     = diff(rhofc.*k_etaf.*diff(Pfex)/dx)/dx;
        % update
        Dcmax        = max(max(k_etaf(1:end-1),k_etaf(2:end)) ./ beta);
        dt           = CFL*dx*dx*min(rhof)/max(rhot)/Dcmax/rel;
        rhot         = rhot_old + drhot_dt * dt;

        errorPf      = max(abs(Pf-Pfr));
        max_err      = max(abs(errorPf))/max(abs(Pf));
        if max_err < eps_err && iter>1, break, end
    end
    time             = time + dt;
    % postprocessing
    niters(it)       = iter;
    err(it)          = max_err;
    vein_width(it)   = 1-sum(Pf>Pf_gr_out) * dx;
    time_it(it)      = time;
    t_vein           = sqrt(time_it);
    dp(it)           = max(Pf - Pmetam);
    if mod(it, 100000) == 0, fprintf('Completed iteration %i - time = %i seconds\n', it, time), end
    if max(isnan(Pf))==1,fprintf("nan at time step %i, iteration %i\n", it, iter),break,end
    if mod(it,nout) == 1 || it == 1

        phi_plot = zeros(nx,4);
        phi_plot(:,4) = phi;
        phi_plot(Pf < Pf_gr_out, 1) = 1 - phi(Pf < Pf_gr_out);
        phi_plot(Pf >= Pf_eclo_in, 2) = 1 - phi(Pf >= Pf_eclo_in);
        phi_plot(Pf > Pf_gr_out & Pf < Pf_eclo_in, 3) = 1 - phi(Pf > Pf_gr_out & Pf < Pf_eclo_in);
        figure(1)
%         subplot(322),plot(time_it, vein_width), title('Front width = ' + string(vein_width(end)))         
%         subplot(324),plot(niters), ylabel('niters'),title(max(niters))
        subplot(411),plot(x,Pf0,'b',x,Pf,'--r', [x(1) x(end)], [Pf_eclo_in Pf_eclo_in], ':k'),ylabel('Pf'), legend('Initial','Current'), xlabel('x'),title('Fluid pressure evolution')
        subplot(412),plot(x,phi), ylabel('\phi'), title('Porosity evolution')
        if time<2e5
        subplot(413),ar = area(x,phi_plot);ylabel('Proportions'),legend('Olivine', 'Serpentine', 'Transition', 'Water'), title('Volumetric proportions')
        ar(1).FaceColor=[10 141 10]/255;ar(2).FaceColor=[186 252 228]/255;ar(3).FaceColor = [0.5 1 0];ar(4).FaceColor=[14 176 246]/255;
        end
        subplot(414),ar = area(x,phi_plot);ylabel('Proportions'),legend('Olivine', 'Serpentine', 'Transition', 'Water'), title('Volumetric proportions')
        ar(1).FaceColor=[10 141 10]/255;ar(2).FaceColor=[186 252 228]/255;ar(3).FaceColor = [0.5 1 0];ar(4).FaceColor=[14 176 246]/255;
        %         subplot(326),plot(log10(err)),xlabel('log10(iter)'), ylabel('Error')
        drawnow
        
                % search of slope and constant  
        line_equation1     = polyfit(t_vein, vein_width,1);
        line_equation2    = polyfit(t_vein2, vein_width2,1);
        line_equation3    = polyfit(t_vein3, vein_width3,1);
        line_equation4    = polyfit(t_vein4, vein_width4,1);
        line_equation5    = polyfit(t_vein5, vein_width5,1);
        line_equation6    = polyfit(t_vein6, vein_width6,1);
        line_equation7    = polyfit(t_vein7, vein_width7,1);
        line_equation8    = polyfit(t_vein8, vein_width8,1);
        line_equation9    = polyfit(t_vein9, vein_width9,1);
        line_equation10    = polyfit(t_vein10, vein_width10,1);

        % creation of an equation for the fitted line with the slope multiplied by x and
        % constant found before
        y_line_eq1   =  line_equation1(1)*t_vein + line_equation1(2);    
        y_line_eq2  =  line_equation2(1)*t_vein2 + line_equation2(2);    
        y_line_eq3  =  line_equation3(1)*t_vein3 + line_equation3(2);    
        y_line_eq4  =  line_equation4(1)*t_vein4 + line_equation4(2);    
        y_line_eq5  =  line_equation5(1)*t_vein5 + line_equation5(2);    
        y_line_eq6  =  line_equation6(1)*t_vein6 + line_equation6(2);    
        y_line_eq7  =  line_equation7(1)*t_vein7 + line_equation7(2);    
        y_line_eq8  =  line_equation8(1)*t_vein8 + line_equation8(2);    
        y_line_eq9  =  line_equation9(1)*t_vein9 + line_equation9(2);    
        y_line_eq10  =  line_equation10(1)*t_vein10 + line_equation10(2);    
        
        % calculation of effective diffusivity (m^2/s) 
        K_eff1       =   (line_equation1(1))^2; 
        K_eff2      =   (line_equation2(1))^2; 
        K_eff3      =   (line_equation3(1))^2; 
        K_eff4      =   (line_equation4(1))^2; 
        K_eff5      =   (line_equation5(1))^2; 
        K_eff6      =   (line_equation6(1))^2; 
        K_eff7      =   (line_equation7(1))^2; 
        K_eff8      =   (line_equation8(1))^2; 
        K_eff9      =   (line_equation9(1))^2; 
        K_eff10      =   (line_equation10(1))^2;
        K_eff_tot   = [K_eff1,K_eff2,K_eff3,K_eff4,K_eff5,K_eff6,K_eff7,K_eff8,K_eff9,K_eff10];
        n_por       = (1:10);
        
%         figure(2)
%         subplot(311),ar = area(x,phi_plot);xlabel('x[m]'),ylabel('Volumetric proportions'),legend('olivine', 'serpentine', 'transition', 'water')
%         ar(1).FaceColor = [0.7 0 0]; ar(2).FaceColor = [0 0.6 0.2]; ar(3).FaceColor = [0.5 1 0]; ar(4).FaceColor = [0.1 0.2 0.4];
%         subplot(312),plot(x,phi), ylabel('Proportion of water'), title(it)
%         subplot(313),plot(x,Pf0,'b',x,Pf,'--r', [x(1) x(end)], [Pf_eclo_in Pf_eclo_in], ':k'),ylabel('Fluid pressure'), legend('initial','current'), xlabel('x[m]')
        t_vein=sqrt(time_it); 
        if time>2e6
        figure (3)
        subplot(521), plot(t_vein, vein_width,"k",t_vein, y_line_eq1,"r"), title('Front reaction width with porosity exponent = 1'); xlabel('Time(s)^{1/2}');ylabel('Vein width (m)'),equ = sprintf('y = %.2ex %.4f', line_equation1(1), line_equation1(2));text(100,0.6,equ,'FontSize',11,'Color','r') 
        subplot(522), plot(t_vein2, vein_width2,"k",t_vein2, y_line_eq2,"r"), title('With porosity exponent = 2'); xlabel('Time(s)^{1/2}');ylabel('Vein width (m)'),equ2 = sprintf('y = %.2ex %.4f', line_equation2(1), line_equation2(2));text(100,0.6,equ2,'FontSize',11,'Color','r') 
        subplot(523), plot(t_vein3, vein_width3,"k",t_vein3, y_line_eq3,"r"), title('With porosity exponent = 3'); xlabel('Time(s)^{1/2}');ylabel('Vein width (m)'),equ3 = sprintf('y = %.2ex %.4f', line_equation3(1), line_equation3(2));text(250,0.3,equ3,'FontSize',11,'Color','r') 
        subplot(524), plot(t_vein4, vein_width4,"k",t_vein4, y_line_eq4,"r"), title('With porosity exponent = 4'); xlabel('Time(s)^{1/2}');ylabel('Vein width (m)'),equ4 = sprintf('y = %.2ex %.4f', line_equation4(1), line_equation4(2));text(250,0.2,equ4,'FontSize',11,'Color','r') 
        subplot(525), plot(t_vein5, vein_width5,"k",t_vein5, y_line_eq5,"r"), title('With porosity exponent = 5'); xlabel('Time(s)^{1/2}');ylabel('Vein width (m)'),equ5 = sprintf('y = %.2ex %.4f', line_equation5(1), line_equation5(2));text(250,0.2,equ5,'FontSize',11,'Color','r')
        subplot(526), plot(t_vein6, vein_width6,"k",t_vein6, y_line_eq6,"r"), title('With porosity exponent = 6'); xlabel('Time(s)^{1/2}');ylabel('Vein width (m)'),equ6 = sprintf('y = %.2ex %.4f', line_equation6(1), line_equation6(2));text(500,0.2,equ6,'FontSize',11,'Color','r')
        subplot(527), plot(t_vein7, vein_width7,"k",t_vein7, y_line_eq7,"r"), title('With porosity exponent = 7'); xlabel('Time(s)^{1/2}');ylabel('Vein width (m)'),equ7 = sprintf('y = %.2ex %.4f', line_equation7(1), line_equation7(2));text(500,0.2,equ7,'FontSize',11,'Color','r')
        subplot(528), plot(t_vein8, vein_width8,"k",t_vein8, y_line_eq8,"r"), title('With porosity exponent = 8'); xlabel('Time(s)^{1/2}');ylabel('Vein width (m)'),equ8 = sprintf('y = %.2ex %.4f', line_equation8(1), line_equation8(2));text(1000,0.2,equ8,'FontSize',11,'Color','r')
        subplot(529), plot(t_vein9, vein_width9,"k",t_vein9, y_line_eq9,"r"), title('With porosity exponent = 9'); xlabel('Time(s)^{1/2}');ylabel('Vein width (m)'),equ9 = sprintf('y = %.2ex %.4f', line_equation9(1), line_equation9(2));text(1000,0.1,equ9,'FontSize',11,'Color','r')
        subplot(5,2,10), plot(t_vein10, vein_width10,"k",t_vein10, y_line_eq10,"r"), title('With porosity exponent = 10'); xlabel('Time(s)^{1/2}');ylabel('Vein width (m)'),equ10 = sprintf('y = %.2ex %.4f', line_equation10(1), line_equation10(2));text(1000,0.04,equ10,'FontSize',11,'Color','r')
        drawnow 
        
        %porosity exponent vs diffusivity
        figure(4)
        plot(n_por,log(K_eff_tot)),semilogy(n_por,log(K_eff_tot)),title("Diffusivities in function of porositiy exponents"),xlabel("Porosity exponents"),ylabel("Log(diffusivity [m^{2}/s])")
        end
%         
        if time>7e5; break 
        end

    end
end
        
%         save('t_vein10_2', 't_vein')
%         save('vein_width_npow10_2','vein_width')
%%

%% ============================= Functions ============================= %%
function var = Itp1D(xlt, varlt, xdata, var, dx, xmin)
for i = 1:length(xdata)
    iW = floor((xdata(i) - xmin) / dx) + 1;
    wW = 1 - (xdata(i) - xlt(iW)) / dx;
    var(i) = wW * varlt(iW) + (1 - wW) * varlt(iW + 1);
end
end