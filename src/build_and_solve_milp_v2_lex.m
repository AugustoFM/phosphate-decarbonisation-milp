function sol = build_and_solve_milp_v2_lex(TimeUTC, G_POA, DNI_kWhm2, Site_Load_kWh, H2_need_kg, unit, econ)
% MILP for PV + CSP (field+TES+turbine) + BESS + Electrolyser
% Two-pass, lexicographic:
%   Pass 1: Minimize sum(uSite_t)
%   Pass 2: Fix sum(uSite_t)=U* and minimize GA NetCost (CAPEX+OPEX+water - revenues - credits + penalties)
% No binaries; integers only for sizes -> robust and convergent.
%
% Electrolyser only consumes energy that is explicitly routed (pvE+cspE+b2E).
% Site is prioritized by the 2-pass scheme (global site-first).

%% ----------- dimensions & constants -------------------------------------
T = numel(TimeUTC);
assert(isvector(G_POA)         && numel(G_POA)==T);
assert(isvector(DNI_kWhm2)     && numel(DNI_kWhm2)==T);
assert(isvector(Site_Load_kWh) && numel(Site_Load_kWh)==T);
assert(isvector(H2_need_kg)    && numel(H2_need_kg)==T);

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
nh  = 19;
idx = @(t,k) (t-1)*nh + k;
NHV = nh*T;

% Global integer sizing vars:
% nPV, nBatt, nElec, nTurb, nStore, nField
NV     = NHV + 6;
nPV    = NHV + 1;
nBatt  = NHV + 2;
nElec  = NHV + 3;
nTurb  = NHV + 4;
nStore = NHV + 5;
nField = NHV + 6;

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
lb = zeros(NV,1);
ub = inf(NV,1);

% Integer design bounds (align with your GA, but allow zeros except PV≥5e3, EL≥5)
lb([nPV nBatt nElec nTurb nStore nField]) = [ 0; 0; 0; 0; 0; 0];
ub([nPV nBatt nElec nTurb nStore nField]) = [3.0e5; 2000; 500; 50; 2000; 2.0e6];

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
for t=1:T
    % PV availability: pvL + pvB + pvX + pvE ≤ G_POA*inv_eff*(nPV*PV_kWp)
    addI([idx(t,1) idx(t,2) idx(t,3) idx(t,4) nPV], ...
         [1 1 1 1, -G_POA(t)*unit.inv_eff*unit.PV_kWp], 0);

    % CSP field availability: f2T + f2S ≤ DNI * nField * m2 * eff
    addI([idx(t,5) idx(t,6) nField], ...
         [1 1, -DNI_kWhm2(t)*unit.CSP_field_unit_m2*csp_e_per_m2], 0);

    % Turbine composition & cap
    addE([idx(t,8) idx(t,5) idx(t,7)], [1 -1 -1], 0);                 % yTurb = f2T + s2T
    addI([idx(t,8) nTurb], [1, -unit.CSP_turb_unit_kW_e], 0);          % yTurb ≤ cap

    % Split turbine outputs: cspL + cspB + cspX + cspE = yTurb
    addE([idx(t,9) idx(t,10) idx(t,11) idx(t,12) idx(t,8)], [1 1 1 1 -1], 0);

    % TES dynamics (cyclic) + capacity + charge cap
    tnext = (t < T) * (t+1) + (t==T) * 1;
    addE([idx(tnext,16) idx(t,16) idx(t,6) idx(t,7)], [1 -1 -1 1], 0); % ETES_{t+1} = ETES_t + f2S - s2T
    addI([idx(t,16) nStore], [1, -unit.CSP_store_unit_kWh_e], 0);
    addI([idx(t,6)  nStore], [1, -unit.CSP_store_charge_kW_per_unit], 0);

    % Battery dynamics (cyclic), SOC bounds, power caps
    tnext = (t < T) * (t+1) + (t==T) * 1;
    addE([idx(tnext,13) idx(t,13) idx(t,2) idx(t,10) idx(t,14) idx(t,15)], ...
         [1, -1, -eta_c, -eta_c, 1/eta_d, 1/eta_d], 0);
    addI([idx(t,13) nBatt], [1, -unit.Batt_kWh*unit.Batt_SOC_max], 0);
    addI([idx(t,13) nBatt], [-1,  unit.Batt_kWh*unit.Batt_SOC_min], 0);
    addI([idx(t,14) idx(t,15) nBatt], [1 1, -unit.Batt_maxDis], 0);   % b2L+b2E ≤ dis cap
    addI([idx(t,2)  idx(t,10) nBatt], [1 1, -unit.Batt_maxChg], 0);   % pvB+cspB ≤ ch cap

    % Site energy balance (NO eEL term here)
    % pvL + cspL + b2L + uSite = Site_Load
    addE([idx(t,1) idx(t,9) idx(t,14) idx(t,18)], [1 1 1 1], Site_Load_kWh(t));

    % Electrolyser routing & demand/capacity
    % eEL = pvE + cspE + b2E
    addE([idx(t,17) idx(t,4) idx(t,12) idx(t,15)], [1 -1 -1 -1], 0);
    % eEL + uH2 = H2_need
    addE([idx(t,17) idx(t,19)], [1 1], H2_need_kWh(t));
    % eEL ≤ nameplate
    addI([idx(t,17) nElec], [1, -unit.Elec_kW], 0);
end

A   = sparse(Ai,Aj,Av,ineq_row,NV);
Aeq = sparse(AEi,AEj,AEv,eq_row,NV);

%% ----------- PASS 1: minimize sum(uSite) -------------------------------
f1 = zeros(NV,1);
for t=1:T, f1(idx(t,18)) = 1; end

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
% Add equality: sum_t uSite_t = Ustar
row = size(Aeq,1) + 1;
AEi2 = [AEi; repmat(row, T, 1)];
uSite_idx = (idx(1,18):nh:idx(T,18)).';
AEj2 = [AEj; uSite_idx];
AEv2 = [AEv; ones(T,1)];
beq2 = [beq; fval1];
Aeq2 = sparse(AEi2,AEj2,AEv2,row,NV);

% Build GA NetCost objective vector f2
f2 = zeros(NV,1);

% Annualised CAPEX (+) on sizes
Capex_per_n = [ ...
  unit.PV_kWp   * unit.PV_Capex; ...
  unit.Batt_kWh * unit.Batt_Capex; ...
  unit.Elec_kW  * unit.Elec_Capex; ...
  unit.CSP_turb_unit_kW_e   * unit.CSP_Capex_turb_USD_per_kWe; ...
  unit.CSP_store_unit_kWh_e * unit.CSP_Capex_store_USD_per_kWh_e; ...
  unit.CSP_field_unit_m2    * unit.CSP_Capex_field_USD_per_m2];
f2([nPV nBatt nElec nTurb nStore nField]) = f2([nPV nBatt nElec nTurb nStore nField]) ...
    + econ.cost_of_CAPEX * Capex_per_n;

% Fixed O&M (+)
Opex_per_n = [ ...
  unit.PV_kWp   * unit.PV_Opex; ...
  unit.Batt_kWh * unit.Batt_Opex; ...
  unit.Elec_kW  * unit.Elec_Opex; ...
  unit.CSP_turb_unit_kW_e * unit.CSP_Opex_USD_per_kWe_yr; ...
  0; 0];
f2([nPV nBatt nElec nTurb nStore nField]) = f2([nPV nBatt nElec nTurb nStore nField]) + Opex_per_n;

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

for t=1:T
    f2(idx(t,17)) = f2(idx(t,17)) + coef_eEL;                       % eEL
    f2(idx(t,3))  = f2(idx(t,3))  - econ.Price_export_energy_kWh_USD; % pvX revenue
    f2(idx(t,11)) = f2(idx(t,11)) - econ.Price_export_energy_kWh_USD; % cspX revenue
    f2(idx(t,18)) = f2(idx(t,18)) + econ.unserved_cost_USD_per_kWh;   % uSite penalty (kept, but U fixed)
    f2(idx(t,19)) = f2(idx(t,19)) + econ.unserved_cost_USD_per_kWh;   % uH2 penalty
end

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
get = @(k) xx(k:nh:nh*T);

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
