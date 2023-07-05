clear variables, close all, clf
%% INPUT %%
Pmetam     = 12.65e8; % Reaction pressure
Pfamb      = 12.75e8; % Ambient pressure
Pamp       = 5e8; % Pressure perturbation amplitude
wini       = 0.01;   % Initial perturbation width
k_etaf0    = 1e-19/1e-3;   % m^2/s/Pa - k = 1e-19 m2, etaf = 1e3 Pa.s
Lx         = 1;          % model length, m
npow       = 8;    % Porosity nonlinearity in permeability

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
CFL0        = 0.1;
niter       = 100000;
rel0        = 1;
eps_err     = 1e-6;
nx          = 400;  % number of grid points
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
% Preallocation
beta       = zeros(nx, 1);
Pfr        = zeros(nx ,1);
rhog       = zeros(nx, 1);
rhoe       = zeros(nx, 1);
niters     = [];
err        = [];
eclo_width = [];
time_it    = [];

%% Time and saving stuff %%
beta       = Itp1D(Pf_ltc, betaf_lt, Pf, beta, dPc, Pf_ltc_min);
time_sc    = Lx^2/k_etaf0*max(beta);
t_tot      = 2e2*time_sc;   % total time, s
t_pert     = 1e1*time_sc;
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
        k_etaf       = k_etaf0*phic.^npow;
        drhot_dt     = diff(rhofc.*k_etaf.*diff(Pfex)/dx)/dx;
        % update
        Dcmax        = max(max(k_etaf(1:end-1),k_etaf(2:end)) ./ beta);
        dt           = CFL*dx*dx*min(rhof)/max(rhot)/Dcmax/rel;
        rhot         = rhot_old + drhot_dt * dt;

        errorPf      = Pf-errorPf;
        max_err      = max(abs(errorPf))/max(abs(Pf));
        if max_err < eps_err && iter>1, break, end
    end
    time             = time + dt;
    % postprocessing
    niters(it)       = iter;
    err(it)          = max_err;
    eclo_width(it)   = 1-sum(Pf>Pf_gr_out) * dx;
    time_it(it)      = time;
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
        subplot(321),plot(x,phi), ylabel('phi H2O'), title(it)
        subplot(322),plot(time_it, eclo_width), title('Front width = ' + string(eclo_width(end)))
        subplot(323),ar = area(x,phi_plot);ylabel('phi'),legend('ol', 'serp', 'tr', 'H2O')
        ar(1).FaceColor = [0.7 0 0]; ar(2).FaceColor = [0 0.6 0.2]; ar(3).FaceColor = [0.5 1 0]; ar(4).FaceColor = [0.1 0.2 0.4];
        subplot(324),plot(niters), ylabel('niters'),title(max(niters))
        subplot(325),plot(x,Pf0,'g',x,Pf,'--r', [x(1) x(end)], [Pf_eclo_in Pf_eclo_in], ':k'),ylabel('Pf'), legend('initial','current'),grid, xlabel('x')
        subplot(326),plot(log10(err)),xlabel('log10(iter)'), ylabel('Error')
        drawnow
        t_vein=sqrt(time_it);
        figure(2)
        plot(t_vein, eclo_width), title('Front reaction width with porosity exponent = 1'); xlabel("Time(s)^{1/2}");ylabel("vein width (m)")
    end
end
%%

%% ============================= Functions ============================= %%
function var = Itp1D(xlt, varlt, xdata, var, dx, xmin)
for i = 1:length(xdata)
    iW = floor((xdata(i) - xmin) / dx) + 1;
    wW = 1 - (xdata(i) - xlt(iW)) / dx;
    var(i) = wW * varlt(iW) + (1 - wW) * varlt(iW + 1);
end
end