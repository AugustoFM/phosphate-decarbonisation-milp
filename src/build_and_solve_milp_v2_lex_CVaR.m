function sol = build_and_solve_milp_v2_lex_CVaR(TimeUTC, G_POA, DNI_kWhm2, ...
                                           Site_Load_kWh, H2_need_kg, ...
                                           unit, econ, p_s)
% Time-series inputs are [T × S]; p_s is S×1 probability; sum(p_s)=1
% MILP for PV + CSP (field+TES+turbine) + BESS + Electrolyser
% Two-pass, lexicographic:
%   Pass 1: Minimize sum(uSite_t)
%   Pass 2: Fix sum(uSite_t)=U* and minimize GA NetCost (CAPEX+OPEX+water - revenues - credits + penalties)
% No binaries; integers only for sizes -> robust and convergent.
%
% Electrolyser only consumes energy that is explicitly routed (pvE+cspE+b2E).
% Site is prioritized by the 2-pass scheme (global site-first).

%% ----------- dimensions & constants -------------------------------------
[T,S] = size(G_POA);                               % ✱ CVaR ADD
assert(all(size(DNI_kWhm2)==[T,S]));


if nargin < 8 || isempty(p_s)
    % if caller did not supply probabilities, assume equal weight
    p_s = ones(S,1) / S;
end
assert(isvector(p_s) && numel(p_s)==S, ...
      'p_s must be an S×1 vector of scenario probabilities');
p_s   = p_s(:);        % enforce column shape


%T = numel(TimeUTC);
%assert(isvector(G_POA)         && numel(G_POA)==T);
%assert(isvector(DNI_kWhm2)     && numel(DNI_kWhm2)==T);
%assert(isvector(Site_Load_kWh) && numel(Site_Load_kWh)==T);
%assert(isvector(H2_need_kg)    && numel(H2_need_kg)==T);

% Electrolyser conversion
kWh2kgH2    = unit.elec_eta / 33.3;                 % kg-H2 per kWh (LHV)
H2_need_kWh = H2_need_kg ./ max(kWh2kgH2, eps);     % kWh demand each hour

% Battery efficiencies (fixed to keep linear)
eta_rt = unit.batt_eta;
eta_c  = sqrt(eta_rt);
eta_d  = sqrt(eta_rt);

% CSP "electric-equivalent" availability per m^2
csp_e_per_m2 = unit.CSP_opt_eff * unit.CSP_therm_eff * unit.CSP_cycle_eff;

%% ----------- variables & indexing --------------------------------------
% Hourly continuous variables (lb=0):
%  1 pvL    PV -> Site
%  2 pvB    PV -> BESS
%  3 pvX    PV -> Export
%  4 pvE    PV -> Electrolyser
%  5 f2T    CSP field -> Turbine
%  6 f2S    CSP field -> TES (charge)
%  7 s2T    TES -> Turbine (discharge)
%  8 yTurb  Turbine electric output
%  9 cspL   Turbine -> Site
% 10 cspB   Turbine -> BESS
% 11 cspX   Turbine -> Export
% 12 cspE   Turbine -> Electrolyser
% 13 Eb     Battery SOC (AC-equivalent kWh)
% 14 b2L    BESS -> Site
% 15 b2E    BESS -> Electrolyser
% 16 ETES   TES SOC (electric-equivalent kWh)
% 17 eEL    Electrolyser total kWh (= pvE + cspE + b2E)
% 18 uSite  Site unmet (kWh)
% 19 uH2    H2 unmet (kWh eq.)

%% ----------- variables & indexing --------------------------------------
nh  = 19;

% helper with explicit scenario index
idx = @(t,s,k) ((s-1)*T + (t-1))*nh + k;   % t = 1…T  ,  s = 1…S


NHV = nh*T*S;                                   % hourly variables
NV  = NHV;                                      % next free slot pointer

% ---------- CVaR auxiliaries ----------
z      = NV + 1;                 % VaR level
eps_s  = (NV+1 : NV+S).';        % S excess variables
NV     = NV + 1 + S;             % advance pointer past them

% ---------- design-size integers ----------
nPV    = NV + 1;
nBatt  = NV + 2;
nElec  = NV + 3;
nTurb  = NV + 4;
nStore = NV + 5;
nField = NV + 6;
NV     = NV + 6;                 % advance pointer again

intcon = [nPV nBatt nElec nTurb nStore nField];


%% ----------- defaults for missing econ fields ---------------------------
defaults = { ...
 'cost_of_CAPEX', 0.02; ...
 'Water_req_L_kgH2', 9; ...
 'Cost_Water_USD_m3', 4; ...
 'kgHFO_per_kgH2', 3.0; ...
 'HFO_Price_USD_kg', 0.40; ...
 'HFO_Truck_capacity_kg', 30000; ...
 'HFO_Transport_distance_km', 2400; ...
 'HFO_Trucking_cost_USD_km', 2.0; ...
 'Road_CO2_EF_kg_tkm', 2.68; ...
 'HFO_CO2_EF_kgkg', 3.2327; ...
 'HFO_CH4_EF_kgkg', 0.0030; ...
 'HFO_N2O_EF_kgkg', 0.0072; ...
 'price_CO2_emission', 0.010; ...
 'price_CH4_emission', 0.2; ...
 'price_N2O_emission', 0.1; ...
 'Price_export_energy_kWh_USD', 0.0; ...
 'unserved_cost_USD_per_kWh', 0; ...
 'Price_Excess_Elec_USD_kWh', 0.10; ...
 'CO2_red_factor_PVvsDiesel', 0.78 };
for i=1:size(defaults,1)
    if ~isfield(econ, defaults{i,1}), econ.(defaults{i,1}) = defaults{i,2}; end
end

%% ----------- bounds -----------------------------------------------------
lb = zeros(NV,1);               % default: all ≥ 0
ub =  inf(NV,1);

% VaR (z) is un-bounded below;  ε_s are ≥ 0
lb(z)      = -inf;
lb(eps_s)  = 0;                 % already 0, kept for clarity
ub(z)      =  inf;              % already inf
ub(eps_s)  =  inf;              % already inf

% Integer design bounds (align with GA, allow zeros)
lb([nPV nBatt nElec nTurb nStore nField]) = [0;0;0;0;0;0];
ub([nPV nBatt nElec nTurb nStore nField]) = [3e5; 2000; 500; 50; 2000; 2e6];


%% ----------- constraint builders ---------------------------------------
Ai=[]; Aj=[]; Av=[];  b=[];     ineq_row=0;
AEi=[];AEj=[];AEv=[]; beq=[];   eq_row=0;

    function addI(cols, vals, rhs)
        ineq_row = ineq_row + 1;
        Ai = [Ai; repmat(ineq_row, numel(cols), 1)];
        Aj = [Aj; cols(:)];
        Av = [Av; vals(:)];
        b(ineq_row,1) = rhs;
    end
    function addE(cols, vals, rhs)
        eq_row = eq_row + 1;
        AEi = [AEi; repmat(eq_row, numel(cols), 1)];
        AEj = [AEj; cols(:)];
        AEv = [AEv; vals(:)];
        beq(eq_row,1) = rhs;
    end

%% ----------- network & physics constraints -----------------------------
for s = 1:S
    for t=1:T
        % PV availability: pvL + pvB + pvX + pvE ≤ G_POA*inv_eff*(nPV*PV_kWp)
        addI([idx(t,s,1) idx(t,s,2) idx(t,s,3) idx(t,s,4) nPV], ...
             [1 1 1 1, -G_POA(t,s)*unit.inv_eff*unit.PV_kWp], 0);% [1 1 1 1, -G_POA(t)*unit.inv_eff*unit.PV_kWp], 0)

        % CSP field availability: f2T + f2S ≤ DNI * nField * m2 * eff
        addI([idx(t,s,5) idx(t,s,6) nField], ...
            [1 1, -DNI_kWhm2(t,s)*unit.CSP_field_unit_m2*csp_e_per_m2], 0);

        % Turbine composition & cap
        addE([idx(t,s,8) idx(t,s,5) idx(t,s,7)], [1 -1 -1], 0);                 % yTurb = f2T + s2T
        addI([idx(t,s,8) nTurb], [1, -unit.CSP_turb_unit_kW_e], 0);          % yTurb ≤ cap

        % Split turbine outputs: cspL + cspB + cspX + cspE = yTurb
        addE([idx(t,s,9) idx(t,s,10) idx(t,s,11) idx(t,s,12) idx(t,s,8)], [1 1 1 1 -1], 0);

        % TES dynamics (cyclic) + capacity + charge cap
        tnext = (t < T) * (t+1) + (t==T) * 1;
        addE([idx(tnext,s,16) idx(t,s,16) idx(t,s,6) idx(t,s,7)], [1 -1 -1 1], 0); % ETES_{t+1} = ETES_t + f2S - s2T
        addI([idx(t,s,16) nStore], [1, -unit.CSP_store_unit_kWh_e], 0);
        addI([idx(t,s,6)  nStore], [1, -unit.CSP_store_charge_kW_per_unit], 0);

        % Battery dynamics (cyclic), SOC bounds, power caps
        tnext = (t < T) * (t+1) + (t==T) * 1;
        addE([idx(tnext,s,13) idx(t,s,13) idx(t,s,2) idx(t,s,10) idx(t,s,14) idx(t,s,15)], ...
            [1, -1, -eta_c, -eta_c, 1/eta_d, 1/eta_d], 0);
        addI([idx(t,s,13) nBatt], [1, -unit.Batt_kWh*unit.Batt_SOC_max], 0);
        addI([idx(t,s,13) nBatt], [-1,  unit.Batt_kWh*unit.Batt_SOC_min], 0);
        addI([idx(t,s,14) idx(t,s,15) nBatt], [1 1, -unit.Batt_maxDis], 0);   % b2L+b2E ≤ dis cap
        addI([idx(t,s,2)  idx(t,s,10) nBatt], [1 1, -unit.Batt_maxChg], 0);   % pvB+cspB ≤ ch cap

        % Site energy balance (NO eEL term here)
        % pvL + cspL + b2L + uSite = Site_Load
        addE([idx(t,s,1) idx(t,s,9) idx(t,s,14) idx(t,s,18)], [1 1 1 1], Site_Load_kWh(t));

        % Electrolyser routing & demand/capacity
        % eEL = pvE + cspE + b2E
        addE([idx(t,s,17) idx(t,s,4) idx(t,s,12) idx(t,s,15)], [1 -1 -1 -1], 0);
        % eEL + uH2 = H2_need
        addE([idx(t,s,17) idx(t,s,19)], [1 1], H2_need_kWh(t));
        % eEL ≤ nameplate
        addI([idx(t,s,17) nElec], [1, -unit.Elec_kW], 0);
    end
end

A = sparse(Ai, Aj, Av, ineq_row, NV);
Aeq = sparse(AEi,AEj,AEv,eq_row,NV);
% optional sanity-check
assert(size(A,1)==ineq_row && size(Aeq,1)==eq_row, 'Row count mismatch');

%% ----------- PASS 1: minimize sum(uSite) -------------------------------
f1 = zeros(NV,1);
for s = 1:S
    for t = 1:T
        f1(idx(t,s,18)) = 1;          % pick uSite of every scenario
    end
end

opts = optimoptions('intlinprog','Display','iter','MaxTime',3600);
[xx1,fval1,exitflag1,output1] = intlinprog(f1, intcon, A, b, Aeq, beq, lb, ub, opts);

sol.pass1.exitflag = exitflag1;
sol.pass1.output   = output1;
sol.pass1.uSite_sum_min = fval1;

if exitflag1 <= 0
    warning('Pass-1 (min uSite) did not find an optimal solution. Exitflag=%d', exitflag1);
    sol.exitflag = exitflag1;
    return;
end

%% ----------- PASS 2: add sum(uSite) = U* and minimize GA NetCost -------
% --- Add equality:  sum_{t,s} uSite(t,s)  =  U* --------------------------
uSite_idx = [];                % collect all uSite indices (T*S of them)
for s = 1:S
    uSite_idx = [uSite_idx; idx(1:T, s, 18)];
end
row   = size(Aeq,1) + 1;       % new row in Aeq
AEi2  = [AEi ; repmat(row, numel(uSite_idx), 1)];
AEj2  = [AEj ; uSite_idx(:)];
AEv2  = [AEv ; ones(numel(uSite_idx), 1)];
beq2  = [beq ; fval1];
Aeq2  = sparse(AEi2, AEj2, AEv2, row, NV);


% eEL coefficient (water cost - revenues/credits)
usd_per_kWh_water   = kWh2kgH2 * (econ.Water_req_L_kgH2/1000) * econ.Cost_Water_USD_m3;
usd_per_kWh_fuelsav = kWh2kgH2 * econ.kgHFO_per_kgH2 * econ.HFO_Price_USD_kg;
usd_per_kWh_truck   = kWh2kgH2 * econ.kgHFO_per_kgH2 * (1/max(econ.HFO_Truck_capacity_kg,eps)) ...
                      * econ.HFO_Transport_distance_km * econ.HFO_Trucking_cost_USD_km;
usd_per_kWh_CO2comb = kWh2kgH2 * econ.kgHFO_per_kgH2 * econ.HFO_CO2_EF_kgkg * econ.price_CO2_emission;
usd_per_kWh_CH4comb = kWh2kgH2 * econ.kgHFO_per_kgH2 * econ.HFO_CH4_EF_kgkg * econ.price_CH4_emission;
usd_per_kWh_N2Ocomb = kWh2kgH2 * econ.kgHFO_per_kgH2 * econ.HFO_N2O_EF_kgkg * econ.price_N2O_emission;
usd_per_kWh_CO2road = kWh2kgH2 * econ.kgHFO_per_kgH2 * (econ.HFO_Transport_distance_km/1000) ...
                      * econ.Road_CO2_EF_kg_tkm * econ.price_CO2_emission;

coef_eEL =  + usd_per_kWh_water ...
            - usd_per_kWh_fuelsav ...
            - usd_per_kWh_truck ...
            - (usd_per_kWh_CO2comb + usd_per_kWh_CH4comb + usd_per_kWh_N2Ocomb + usd_per_kWh_CO2road);


% ------------------------------------------------------------------
% GA-style Net-Cost vector  ➜  call it f_cost for later CVaR use
% ------------------------------------------------------------------
f_cost = zeros(NV,1);          % start with an all-zero vector

% 1) Annualised CAPEX  (+)          (same 6-element slice)
capex_per_n = [ ...
    unit.PV_kWp   * unit.PV_Capex; ...
    unit.Batt_kWh * unit.Batt_Capex; ...
    unit.Elec_kW  * unit.Elec_Capex; ...
    unit.CSP_turb_unit_kW_e   * unit.CSP_Capex_turb_USD_per_kWe; ...
    unit.CSP_store_unit_kWh_e * unit.CSP_Capex_store_USD_per_kWh_e; ...
    unit.CSP_field_unit_m2    * unit.CSP_Capex_field_USD_per_m2];

f_cost([nPV nBatt nElec nTurb nStore nField]) = ...
        f_cost([nPV nBatt nElec nTurb nStore nField]) ...
      + econ.cost_of_CAPEX * capex_per_n;

% 2) Fixed O&M  (+)
opex_per_n = [ ...
    unit.PV_kWp   * unit.PV_Opex; ...
    unit.Batt_kWh * unit.Batt_Opex; ...
    unit.Elec_kW  * unit.Elec_Opex; ...
    unit.CSP_turb_unit_kW_e * unit.CSP_Opex_USD_per_kWe_yr; ...
    0; 0];                                % no O&M for TES tanks & field

f_cost([nPV nBatt nElec nTurb nStore nField]) = ...
        f_cost([nPV nBatt nElec nTurb nStore nField]) + opex_per_n;

% 3) Hourly terms  (water cost, export revenue, unserved penalties …)
for s = 1:S
    for t = 1:T
        f_cost(idx(t,s,17)) = f_cost(idx(t,s,17)) + coef_eEL;                 % eEL
        f_cost(idx(t,s, 3)) = f_cost(idx(t,s, 3)) - econ.Price_export_energy_kWh_USD;  % PV export
        f_cost(idx(t,s,11)) = f_cost(idx(t,s,11)) - econ.Price_export_energy_kWh_USD;  % CSP export
        f_cost(idx(t,s,18)) = f_cost(idx(t,s,18)) + econ.unserved_cost_USD_per_kWh;    % uSite
        f_cost(idx(t,s,19)) = f_cost(idx(t,s,19)) + econ.unserved_cost_USD_per_kWh;    % uH2
    end
end

%% ---------- CVaR objective vector ---------------------------------------
lambda = 0.5;      % risk-aversion weight 1
alpha  = 0.80;

% ----- ❶ Expected-value part (length = NV) -------------------------------
f_exp = zeros(NV,1);

for s = 1:S
    base = (s-1)*T*nh;                % first index of scenario-s hourly vars
    rng  = base + (1:T*nh);           % that whole block
    f_exp(rng) = p_s(s) * f_cost(rng);
end
f_exp(NHV+1:NV) = f_cost(NHV+1:NV);   % z, ε_s, and 6 size integers (no prob-weighting)

% ----- ❷ CVaR part (length = NV) -----------------------------------------
f_cvar = lambda * [ ...
        zeros(NHV,1) ;               % hourly vars
        0            ;               %      z
        (1/(1-alpha))*p_s ;          %  ε_s  (S entries)
        zeros(6,1) ];                %  6 size-integer vars

% ----- ❸ Full objective & loss copy --------------------------------------
f2     = f_exp + f_cvar;
f_loss = f_exp;                       % use in loss_s-z ≤ ε_s constraints



for s = 1:S
    for t=1:T
        f2(idx(t,s,17)) = f2(idx(t,s,17)) + coef_eEL;                       % eEL
        f2(idx(t,s,3))  = f2(idx(t,s,3))  - econ.Price_export_energy_kWh_USD; % pvX revenue
        f2(idx(t,s,11)) = f2(idx(t,s,11)) - econ.Price_export_energy_kWh_USD; % cspX revenue
        f2(idx(t,s,18)) = f2(idx(t,s,18)) + econ.unserved_cost_USD_per_kWh;   % uSite penalty (kept, but U fixed)
        f2(idx(t,s,19)) = f2(idx(t,s,19)) + econ.unserved_cost_USD_per_kWh;   % uH2 penalty
    end
end

%% ----------- CVaR "loss – z ≤ ε" constraints -----------------------
for s = 1:S
    ineq_row  = ineq_row + 1;

    % ------- include these two lines -------------------------------------
    base      = (s-1)*T*nh;          % first element of scenario-s block
    loss_vals = f_loss(base + (1:T*nh));   % slice of f_loss for scenario s
    % -------------------------------------------------------------------

    loss_cols = base + (1:T*nh);     % same indices as above

    Ai = [Ai; repmat(ineq_row,numel(loss_cols)+2,1)];
    Aj = [Aj; loss_cols(:); z; eps_s(s)];
    Av = [Av; loss_vals(:); -1; -1]; % −z  −εₛ
    b (ineq_row,1) = 0;
end

%% ----------- *NOW* build the sparse matrices ----------------------
A   = sparse(Ai,Aj,Av,ineq_row,NV);   % rebuilt with the new row count
Aeq = sparse(AEi,AEj,AEv,eq_row,NV);  % unchanged

[xx,fval,exitflag,output] = intlinprog(f2, intcon, A, b, Aeq2, beq2, lb, ub, opts);

%% ----------- unpack -----------------------------------------------------
sol.exitflag = exitflag;
sol.output   = output;
sol.fval_USD = fval;          % GA-like objective value excluding constants

% Sizes
sol.NumPV        = round(xx(nPV));
sol.NumBatt      = round(xx(nBatt));
sol.NumElec      = round(xx(nElec));
sol.NumCSP_Turb  = round(xx(nTurb));
sol.NumCSP_Store = round(xx(nStore));
sol.NumCSP_Field = round(xx(nField));

% Hourly series getter
get = @(k) reshape(xx(k:nh:nh*T*S), T, S);   % returns T×S matrix

sol.PV_to_Load_kWh    = get(1);
sol.PV_to_Batt_kWh    = get(2);
sol.PV_Export_kWh     = get(3);
sol.PV_to_Elec_kWh    = get(4);

sol.Field_to_Turb_kWh = get(5);
sol.Field_to_TES_kWh  = get(6);
sol.TES_to_Turb_kWh   = get(7);
sol.Turb_kWh          = get(8);

sol.CSP_to_Load_kWh   = get(9);
sol.CSP_to_Batt_kWh   = get(10);
sol.CSP_Export_kWh    = get(11);
sol.CSP_to_Elec_kWh   = get(12);

sol.Ebatt_kWh         = get(13);
sol.Batt_to_Load_kWh  = get(14);
sol.Batt_to_Elec_kWh  = get(15);

sol.ETES_kWh          = get(16);
sol.Elec_kWh          = get(17);     % = pvE + cspE + b2E
sol.Unserved_Site_kWh = get(18);
sol.Unserved_H2_kWh   = get(19);

% Sums
sol.Export_kWh          = sol.PV_Export_kWh + sol.CSP_Export_kWh;
sol.Unserved_total_kWh  = sum(sol.Unserved_Site_kWh) + sum(sol.Unserved_H2_kWh);
sol.annualH2_kg         = kWh2kgH2 * sum(sol.Elec_kWh);
sol.pass2_fval          = fval;

%% ----------- GA-style economics (same as your driver) -------------------
P_Turb_kW   = sol.NumCSP_Turb  * unit.CSP_turb_unit_kW_e;
E_TES_kWh   = sol.NumCSP_Store * unit.CSP_store_unit_kWh_e;
A_Field_m2  = sol.NumCSP_Field * unit.CSP_field_unit_m2;

Capex_USD = ...
  sol.NumPV       * unit.PV_kWp   * unit.PV_Capex + ...
  sol.NumBatt     * unit.Batt_kWh * unit.Batt_Capex + ...
  sol.NumElec     * unit.Elec_kW  * unit.Elec_Capex + ...
  P_Turb_kW       * unit.CSP_Capex_turb_USD_per_kWe + ...
  E_TES_kWh       * unit.CSP_Capex_store_USD_per_kWh_e + ...
  A_Field_m2      * unit.CSP_Capex_field_USD_per_m2;

AnnualisedCAPEX_USD = econ.cost_of_CAPEX * Capex_USD;

Fixed_Opex_USD = ...
  sol.NumPV     * unit.PV_kWp   * unit.PV_Opex + ...
  sol.NumBatt   * unit.Batt_kWh * unit.Batt_Opex + ...
  sol.NumElec   * unit.Elec_kW  * unit.Elec_Opex + ...
  P_Turb_kW     * unit.CSP_Opex_USD_per_kWe_yr;

H2_prod_kg_y    = sol.annualH2_kg;
HFO_displaced_kg= H2_prod_kg_y * econ.kgHFO_per_kgH2;
Revenue_FuelSaved   = HFO_displaced_kg * econ.HFO_Price_USD_kg;
nTrips              = HFO_displaced_kg / max(econ.HFO_Truck_capacity_kg,eps);
TruckingCost_USD    = nTrips * econ.HFO_Transport_distance_km * econ.HFO_Trucking_cost_USD_km;
ton_km              = (HFO_displaced_kg/1000) * econ.HFO_Transport_distance_km;
Road_CO2_kg         = ton_km * econ.Road_CO2_EF_kg_tkm;
CO2_avoided_kg      = HFO_displaced_kg * econ.HFO_CO2_EF_kgkg;
CH4_avoided_kg      = HFO_displaced_kg * econ.HFO_CH4_EF_kgkg;
N2O_avoided_kg      = HFO_displaced_kg * econ.HFO_N2O_EF_kgkg;

Water_m3            = H2_prod_kg_y * (econ.Water_req_L_kgH2/1000);
Cost_Water_USD      = Water_m3 * econ.Cost_Water_USD_m3;

ExportRevenue_USD   = sum(sol.Export_kWh) * econ.Price_export_energy_kWh_USD;

Unserved_cost_USD   = econ.unserved_cost_USD_per_kWh * sol.Unserved_total_kWh;

Total_y_Site_Load_kWh     = sum(Site_Load_kWh);
Total_y_Site_USD          = Total_y_Site_Load_kWh * econ.Price_Excess_Elec_USD_kWh;
CO2_Total_y_Site_Load_kWh = Total_y_Site_Load_kWh * econ.CO2_red_factor_PVvsDiesel;

Cash_Positive = TruckingCost_USD + Revenue_FuelSaved + Total_y_Site_USD + ExportRevenue_USD + ...
    (Road_CO2_kg + CO2_avoided_kg + CO2_Total_y_Site_Load_kWh) * econ.price_CO2_emission + ...
     CH4_avoided_kg * econ.price_CH4_emission + N2O_avoided_kg * econ.price_N2O_emission;

Cash_Negative = Cost_Water_USD + Fixed_Opex_USD + AnnualisedCAPEX_USD + Unserved_cost_USD;

sol.NetCost_USD_GA = Cash_Negative - Cash_Positive;

% fval excludes the two constant "site" terms below; compute offset if needed
Const_CO2_site_USD = CO2_Total_y_Site_Load_kWh * econ.price_CO2_emission;
Const_site_USD     = Total_y_Site_USD;
const_offset       = - (Const_site_USD + Const_CO2_site_USD);

sol.NetCost_USD_GAequiv = sol.pass2_fval + const_offset;

% Bookkeeping
sol.obj_components = struct( ...
    'AnnualisedCAPEX_USD',AnnualisedCAPEX_USD, ...
    'Fixed_Opex_USD',Fixed_Opex_USD, ...
    'Cost_Water_USD',Cost_Water_USD, ...
    'Unserved_cost_USD',Unserved_cost_USD, ...
    'Revenue_FuelSaved_USD',Revenue_FuelSaved, ...
    'TruckingCost_USD',TruckingCost_USD, ...
    'ExportRevenue_USD',ExportRevenue_USD, ...
    'Road_CO2_kg',Road_CO2_kg, ...
    'CO2_Total_y_Site_Load_kWh', CO2_Total_y_Site_Load_kWh, ...
    'CO2_avoided_kg',CO2_avoided_kg, ...
    'CH4_avoided_kg',CH4_avoided_kg, ...
    'N2O_avoided_kg',N2O_avoided_kg, ...
    'Const_site_USD',Const_site_USD, ...
    'Const_CO2_site_USD',Const_CO2_site_USD);

end
