% run_milp_milp_v2_lex.m — GA-compatible driver for the MILP (intlinprog)
clear; clc;

%% A. Load & align inputs (same behavior as your GA driver)
irr = load('irradiance_hourly.mat');  % expects TimeUTC + (G_POA or GHI_kWhm2)
assert(isfield(irr,'TimeUTC'), 'irradiance_hourly.mat must contain TimeUTC');
TimeUTC = irr.TimeUTC(:);
assert(isdatetime(TimeUTC) && iscolumn(TimeUTC), 'TimeUTC must be a datetime column vector');
nH = numel(TimeUTC);

if     isfield(irr,'G_POA'),       G_POA = irr.G_POA(:);
elseif isfield(irr,'GHI_kWhm2'),   G_POA = irr.GHI_kWhm2(:);
else,  error('irradiance_hourly.mat must have G_POA or GHI_kWhm2'); 
end
assert(numel(G_POA)==nH, 'Irradiance vector length must match TimeUTC');

% Hydrogen demand (hourly)
dem = load('hourly_H2_demand_100_level.mat');  % expects H2_need_hour
assert(isfield(dem,'H2_need_hour'), 'hourly_H2_demand_100.mat must contain H2_need_hour');
H2_need_hour = dem.H2_need_hour(:);
assert(numel(H2_need_hour)==nH || ~isfield(dem,'TimeUTC'), 'H2_need_hour length mismatch vs TimeUTC');

% Diesel hourly energy (site load proxy), prefer .mat, fallback .xlsx
Diesel_hourly_energy_kWh = [];
if exist('diesel_hourly_energy.mat','file')
    d = load('diesel_hourly_energy.mat');  % expects Energy_produced_Diesel_kWh (+ TimeUTC)
    assert(isfield(d,'Energy_produced_Diesel_kWh'), 'diesel_hourly_energy.mat must contain Energy_produced_Diesel_kWh');
    if isfield(d,'TimeUTC') && isequal(d.TimeUTC(:), TimeUTC)
        Diesel_hourly_energy_kWh = d.Energy_produced_Diesel_kWh(:);
    else
        if isfield(d,'TimeUTC')
            ttMaster = timetable(TimeUTC, 'VariableNames', {'dummy'});
            ttDiesel = timetable(d.TimeUTC(:), d.Energy_produced_Diesel_kWh(:), 'VariableNames', {'diesel'});
            tt = synchronize(ttMaster, ttDiesel, 'first', 'nearest');
            Diesel_hourly_energy_kWh = tt.diesel;
        else
            assert(numel(d.Energy_produced_Diesel_kWh)==nH, 'Diesel hourly length must match TimeUTC length');
            Diesel_hourly_energy_kWh = d.Energy_produced_Diesel_kWh(:);
        end
    end
elseif exist('diesel_hourly_energy.xlsx','file')
    Tdies = readtable('diesel_hourly_energy.xlsx','VariableNamingRule','preserve');
    if ismember('Energy_produced_Diesel_kWh', Tdies.Properties.VariableNames) && ...
       ismember('TimeUTC', Tdies.Properties.VariableNames)
        Diesel_hourly_energy_kWh = Tdies.Energy_produced_Diesel_kWh(:);
        if ~isequal(Tdies.TimeUTC(:), TimeUTC)
            ttMaster = timetable(TimeUTC, 'VariableNames', {'dummy'});
            ttDiesel = timetable(Tdies.TimeUTC(:), Diesel_hourly_energy_kWh, 'VariableNames', {'diesel'});
            tt = synchronize(ttMaster, ttDiesel, 'first', 'nearest');
            Diesel_hourly_energy_kWh = tt.diesel;
        end
    else
        % Fallback: first numeric column assumed to be hourly kWh already
        numMask = varfun(@isnumeric, Tdies, 'OutputFormat', 'uniform');
        assert(any(numMask), 'No numeric column found in diesel_hourly_energy.xlsx');
        c = find(numMask,1,'first');
        Diesel_hourly_energy_kWh = Tdies{:,c}(:);
        assert(numel(Diesel_hourly_energy_kWh)==nH, 'Hourly diesel length in XLSX must match TimeUTC length');
    end
else
    error('Cannot find diesel_hourly_energy.mat or diesel_hourly_energy.xlsx in the working folder.');
end
assert(numel(Diesel_hourly_energy_kWh)==nH, 'Diesel hourly vector must match TimeUTC length');

dni = load('dni_hourly.mat');   % must contain DNI_kWhm2 and optionally TimeUTC
assert(isfield(dni,'DNI_kWhm2'),'dni_hourly.mat must contain DNI_kWhm2');
if isfield(dni,'TimeUTC')
    assert(isequal(dni.TimeUTC(:), TimeUTC), 'DNI TimeUTC must match irradiance TimeUTC.');
end
DNI_kWhm2 = dni.DNI_kWhm2(:);

Site_Load_kWh = Diesel_hourly_energy_kWh(:);
H2_need_kg    = H2_need_hour(:);

%% B. Techno-economic parameters (copied from GA)
% PV
unit.PV_kWp     = 1;       unit.PV_Capex   =  550;   unit.PV_Opex   = 12;
unit.inv_eff    = 0.97;
unit.PV_area_m2_per_kWp = 6;      % for land area reporting

% BESS
unit.Batt_kWh   = 250;     unit.Batt_Capex =  250;   unit.Batt_Opex = 7.5;
unit.batt_eta   = 0.92;    unit.Batt_maxChg = 250;   unit.Batt_maxDis = 250; % kWh/h (AC)
unit.Batt_SOC_min = 0.10;  unit.Batt_SOC_max = 0.95;

% Electrolyser
unit.Elec_kW    = 500;     unit.Elec_Capex = 1950;   unit.Elec_Opex = 75;
unit.elec_eta   = 0.70;

% CSP
unit.CSP_turb_unit_kW_e             = 10000;
unit.CSP_store_unit_kWh_e           = 50000;
unit.CSP_field_unit_m2              = 1000;
unit.CSP_store_charge_kW_per_unit   = unit.CSP_store_unit_kWh_e / 8;
unit.CSP_opt_eff   = 0.80;
unit.CSP_therm_eff = 0.55;
unit.CSP_cycle_eff = 0.38;
unit.CSP_Capex_turb_USD_per_kWe     = 1500;
unit.CSP_Capex_store_USD_per_kWh_e  = 30;
unit.CSP_Capex_field_USD_per_m2     = 200;
unit.CSP_Opex_USD_per_kWe_yr        = 65;

%% C. Macro-economic / environmental (copied from GA)
econ.kgHFO_per_kgH2            = 3.0;
econ.HFO_Price_USD_kg          = 0.40;
econ.HFO_Truck_capacity_kg     = 30000;
econ.HFO_Transport_distance_km = 2400;
econ.HFO_Trucking_cost_USD_km  = 2.0;

econ.HFO_CO2_EF_kgkg = 3.2327; econ.HFO_CH4_EF_kgkg = 0.0030; econ.HFO_N2O_EF_kgkg = 0.0072;
econ.Road_CO2_EF_kg_tkm = 2.68;

econ.Cost_Water_USD_m3 = 4.0;  econ.Water_req_L_kgH2 = 9.0;

econ.Price_Excess_Elec_USD_kWh = 0.20; % Upper bound of cost of Energy generated by diesel
econ.CO2_red_factor_PVvsDiesel = 0.78;

econ.Price_export_energy_kWh_USD = 0.00;

econ.price_CO2_emission = 0.0; econ.price_CH4_emission = 0.0; econ.price_N2O_emission = 0; %Recommended prices for KAPSARC 2030 0.05 CO2, 0.54 CH4, 5.46 N20

econ.cost_of_CAPEX = 0.02;
econ.unserved_cost_USD_per_kWh = 0; %50USD suggested

%% D. (Optional) short horizon for a smoke test
use_short_horizon = false; hours_to_keep = 24*7;
if use_short_horizon
    idxKeep = 1:min(hours_to_keep, nH);
    TimeUTC = TimeUTC(idxKeep);
    G_POA   = G_POA(idxKeep);
    DNI_kWhm2 = DNI_kWhm2(idxKeep);
    Site_Load_kWh = Site_Load_kWh(idxKeep);
    H2_need_kg    = H2_need_kg(idxKeep);
    fprintf('Using short horizon of %d hours.\n', numel(idxKeep));
end

%% E. Solve MILP (objective matches GA NetCost)
fprintf('Solving MILP with intlinprog... (8760 h can take a while)\n');
tic;
sol = build_and_solve_milp_v2_lex(TimeUTC, G_POA, DNI_kWhm2, Site_Load_kWh, H2_need_kg, unit, econ);
if ~isfield(sol,'exitflag') || sol.exitflag <= 0
    warning('MILP infeasible or no solution returned (exitflag=%d).', ...
            getfield(sol,'exitflag',-Inf)); %#ok<GFLD>
    % Optional: write a minimal diagnostics file or bail out early:
    return;
end
t_solve = toc;

%% F. Derived sizes, capacities, areas (GA-compatible names)
NumPV        = sol.NumPV;
NumBatt      = sol.NumBatt;
NumElec      = sol.NumElec;
NumCSP_Turb  = sol.NumCSP_Turb;
NumCSP_Store = sol.NumCSP_Store;
NumCSP_Field = sol.NumCSP_Field;

Total_PV_kWp   = NumPV        * unit.PV_kWp;            Total_PV_MWp = Total_PV_kWp/1e3;
Total_Batt_kWh = NumBatt      * unit.Batt_kWh;          Total_Batt_MWh = Total_Batt_kWh/1e3;
Total_Elec_kW  = NumElec      * unit.Elec_kW;           Total_Elec_MW  = Total_Elec_kW/1e3;
Cap_turb_kW_tot= NumCSP_Turb  * unit.CSP_turb_unit_kW_e;
Cap_store_kWh  = NumCSP_Store * unit.CSP_store_unit_kWh_e;
Field_area_m2  = NumCSP_Field * unit.CSP_field_unit_m2;

% PV area (for reporting)
PV_Area_m2  = Total_PV_kWp * unit.PV_area_m2_per_kWp;
PV_Area_ha  = PV_Area_m2 / 1e4;
PV_Area_km2 = PV_Area_m2 / 1e6;

% Converter constants
kWh2kgH2    = unit.elec_eta / 33.3;
H2_need_kWh = H2_need_kg ./ max(kWh2kgH2, eps);

%% G. Hourly rollups for GA-like fields
% Potential PV & CSP field energy (as in GA)
PV_potential_kWh       = Total_PV_kWp .* G_POA .* unit.inv_eff;
CSP_Field_kWh_e_avail  = DNI_kWhm2 .* Field_area_m2 .* unit.CSP_opt_eff .* unit.CSP_therm_eff .* unit.CSP_cycle_eff;

% Served/unserved & attributions
Site_served_kWh  = Site_Load_kWh - sol.Unserved_Site_kWh;
Elec_served_kWh  = sol.Elec_kWh;                               % H2 actually produced
Served_Total_kWh = Site_served_kWh + Elec_served_kWh;
Total_Load_kWh   = Site_Load_kWh + H2_need_kWh;
Unserved_Total_kWh = sol.Unserved_Site_kWh + sol.Unserved_H2_kWh;

% Grid export
GridExport_kWh = sol.PV_Export_kWh + sol.CSP_Export_kWh;

% Electrolyser production in kg/h
H2_prod_kg = sol.Elec_kWh .* kWh2kgH2;

% Exact source attribution (replace the proportional version)
Elec_from_PV_kWh   = sol.PV_to_Elec_kWh;
Elec_from_CSP_kWh  = sol.CSP_to_Elec_kWh;
Elec_from_Batt_kWh = sol.Batt_to_Elec_kWh;

% "Turbines_On" estimate (continuous)
Turbines_On_est = sol.Turb_kWh ./ max(unit.CSP_turb_unit_kW_e, eps);

%% H. Print concise summary (GA-style)
fprintf('\n===== MILP SOLUTION (GA-compatible) =====\n');
fprintf('Status: exitflag=%d | Solve time: %.1f s\n', sol.exitflag, t_solve);
fprintf('Sizes  : PV=%d kWp | BESS=%d x%dkWh=%.1f MWh | EL=%d x%dkW=%.1f MW\n', ...
        NumPV, NumBatt, unit.Batt_kWh, Total_Batt_MWh, NumElec, unit.Elec_kW, Total_Elec_MW);
fprintf('CSP    : Turb=%d (%.1f MW_e) | TES=%d (%.1f MWh_e) | Field=%d (%.0f m^2)\n', ...
        NumCSP_Turb, Cap_turb_kW_tot/1e3, NumCSP_Store, Cap_store_kWh/1e3, NumCSP_Field, Field_area_m2);
fprintf('H2 prod: %.2f kt/yr | Export: %.2f GWh | Unserved: %.2f GWh\n', ...
        sum(H2_prod_kg)/1e6, sum(GridExport_kWh)/1e6, sum(Unserved_Total_kWh)/1e6);
fprintf('GA NetCost (reconstructed): $%.2f M | GA-equivalent objective (solver+const): $%.2f M\n', ...
        sol.NetCost_USD_GA/1e6, sol.NetCost_USD_GAequiv/1e6);

%% I. Write EVERYTHING to ONE Excel file (same sheets & fields as GA)
stamp   = char(datetime('now','Format','yyyyMMddHHmm'));
oneFile = sprintf('%s_H2_CSP_System_AllResults_MILP_lex_v2.xlsx', stamp);

% 1) Hourly sheet (column names mirror GA)
Hourly = table( ...
    TimeUTC, G_POA, DNI_kWhm2, ...
    repmat(NumPV,nH,1), PV_potential_kWh, ...
    repmat(NumCSP_Field,nH,1), CSP_Field_kWh_e_avail, ...
    repmat(NumCSP_Turb,nH,1),  sol.Turb_kWh, Turbines_On_est, ...
    sol.Field_to_TES_kWh, sol.TES_to_Turb_kWh, sol.ETES_kWh, ...
    repmat(NumElec,nH,1), sol.Elec_kWh, H2_need_kg, H2_prod_kg, ...
    Site_Load_kWh, Total_Load_kWh, Served_Total_kWh, ...
    sol.Unserved_Site_kWh, sol.Unserved_H2_kWh, Unserved_Total_kWh, ...
    Elec_from_PV_kWh, Elec_from_CSP_kWh, Elec_from_Batt_kWh, ...
    sol.PV_to_Batt_kWh, sol.Ebatt_kWh, GridExport_kWh, ...
    repmat(NumBatt,nH,1), ...
    'VariableNames', {'Datetime','Irrad_kWh_kWp','DNI_kWhm2', ...
        'NumPV','PV_kWh', ...
        'NumCSP_Field','CSP_Field_kWh_e', ...
        'NumCSP_Turb','CSP_Turbine_kWh','Turbines_On', ...
        'CSP_Charge_kWh','CSP_Discharge_kWh','CSP_SOC_kWh', ...
        'NumElec','Elec_E_kWh','H2_demand_kg','H2_prod_kg', ...
        'Site_Load_kWh','Total_Load_kWh','Served_Total_kWh', ...
        'Unserved_Site_kWh','Unserved_Elec_kWh','Unserved_Total_kWh', ...
        'Elec_from_PV_kWh','Elec_from_CSP_kWh','Elec_from_Batt_kWh', ...
        'PV_to_Batt_kWh','Batt_SOC_kWh','GridExport_kWh', ...
        'NumBatt'} );
writetable(Hourly, oneFile, 'Sheet', 'Hourly');

% 2) Summary sheet (same fields/order as GA)
% Economics recomputation to fill the GA metric fields
% (Most of these are also in sol.obj_components)
P_Turb_kW   = Cap_turb_kW_tot;
E_TES_kWh   = Cap_store_kWh;
A_Field_m2  = Field_area_m2;

Capex_USD = ...
  NumPV       * unit.PV_kWp   * unit.PV_Capex + ...
  NumBatt     * unit.Batt_kWh * unit.Batt_Capex + ...
  NumElec     * unit.Elec_kW  * unit.Elec_Capex + ...
  P_Turb_kW   * unit.CSP_Capex_turb_USD_per_kWe + ...
  E_TES_kWh   * unit.CSP_Capex_store_USD_per_kWh_e + ...
  A_Field_m2  * unit.CSP_Capex_field_USD_per_m2;

AnnualisedCAPEX_USD = econ.cost_of_CAPEX * Capex_USD;

Fixed_Opex_USD = ...
  NumPV     * unit.PV_kWp   * unit.PV_Opex + ...
  NumBatt   * unit.Batt_kWh * unit.Batt_Opex + ...
  NumElec   * unit.Elec_kW  * unit.Elec_Opex + ...
  P_Turb_kW * unit.CSP_Opex_USD_per_kWe_yr;

Total_Elec_E_kWh     = sum(sol.Elec_kWh);
Total_Load_kWh_yr    = sum(Total_Load_kWh);
Unserved_Total_kWh_y = sum(Unserved_Total_kWh);
Served_Total_kWh_y   = sum(Served_Total_kWh);
PV_export_kWh_y      = sum(GridExport_kWh);     % total export (PV+CSP)

H2_prod_kg_y         = sum(H2_prod_kg);
HFO_displaced_kg     = H2_prod_kg_y * econ.kgHFO_per_kgH2;
Revenue_FuelSaved    = HFO_displaced_kg * econ.HFO_Price_USD_kg;
nTrips               = HFO_displaced_kg / max(econ.HFO_Truck_capacity_kg,eps);
TruckingCost_USD     = nTrips * econ.HFO_Transport_distance_km * econ.HFO_Trucking_cost_USD_km;
ton_km               = (HFO_displaced_kg/1000) * econ.HFO_Transport_distance_km;
Road_CO2_kg          = ton_km * econ.Road_CO2_EF_kg_tkm;

CO2_avoided_kg       = HFO_displaced_kg * econ.HFO_CO2_EF_kgkg;
CH4_avoided_kg       = HFO_displaced_kg * econ.HFO_CH4_EF_kgkg;
N2O_avoided_kg       = HFO_displaced_kg * econ.HFO_N2O_EF_kgkg;

Water_m3             = H2_prod_kg_y * (econ.Water_req_L_kgH2/1000);
WaterCost_USD        = Water_m3 * econ.Cost_Water_USD_m3;
ExportRevenue_USD    = PV_export_kWh_y * econ.Price_export_energy_kWh_USD;

Total_y_Site_Load_kWh = sum(Site_Load_kWh);
Total_y_Site_USD      = Total_y_Site_Load_kWh * econ.Price_Excess_Elec_USD_kWh;
CO2_Total_y_Site_Load_kWh = Total_y_Site_Load_kWh * econ.CO2_red_factor_PVvsDiesel;

% NetCost in GA semantics (for the Summary line)
Unserved_Penalty_USD = econ.unserved_cost_USD_per_kWh * Unserved_Total_kWh_y;
NetCost_USD = (WaterCost_USD + Fixed_Opex_USD + AnnualisedCAPEX_USD + Unserved_Penalty_USD) ...
              - (TruckingCost_USD + Revenue_FuelSaved + Total_y_Site_USD + ExportRevenue_USD ...
                 + (Road_CO2_kg + CO2_avoided_kg + CO2_Total_y_Site_Load_kWh) * econ.price_CO2_emission ...
                 + CH4_avoided_kg * econ.price_CH4_emission + N2O_avoided_kg * econ.price_N2O_emission);

metricNames = { ...
  'H2_prod_kg'
  'HFO_displaced_kg'
  'Total_Elec_E_kWh'
  'Total_Load_kWh'
  'Total_y_Site_Load_kWh'
  'Total_y_Site_USD'
  'Unserved_Total_kWh'
  'Served_Total_kWh'
  'PV_export_kWh'
  'PV_area_m2'
  'PV_area_ha'
  'PV_area_km2'
  'Total_CAPEX_USD'
  'Annualised_CAPEX_USD'
  'Fixed_Opex_USD'
  'TruckingCost_USD'
  'WaterCost_USD'
  'FuelSaved_Revenue_USD'
  'PV_Export_Revenue_USD'    % NOTE: this is total export revenue (PV+CSP)
  'NetCost_USD'
  'CO2_avoided_kg'
  'Road_CO2_kg'
  'CO2_Total_y_Site_Load_kWh'
  'CH4_avoided_kg'
  'N2O_avoided_kg'};

metricValues = [ ...
  H2_prod_kg_y; ...
  HFO_displaced_kg; ...
  Total_Elec_E_kWh; ...
  Total_Load_kWh_yr; ...
  Total_y_Site_Load_kWh; ...
  Total_y_Site_USD; ...
  Unserved_Total_kWh_y; ...
  Served_Total_kWh_y; ...
  PV_export_kWh_y; ...
  PV_Area_m2; ...
  PV_Area_ha; ...
  PV_Area_km2; ...
  Capex_USD; ...
  AnnualisedCAPEX_USD; ...
  Fixed_Opex_USD; ...
  TruckingCost_USD; ...
  WaterCost_USD; ...
  Revenue_FuelSaved; ...
  ExportRevenue_USD; ...
  NetCost_USD; ...
  CO2_avoided_kg; ...
  Road_CO2_kg; ...
  CO2_Total_y_Site_Load_kWh; ...
  CH4_avoided_kg; ...
  N2O_avoided_kg ];

SummaryTbl = table(metricNames, metricValues, 'VariableNames', {'Metric','Value_yearly'});
writetable(SummaryTbl, oneFile, 'Sheet', 'Summary');

% 3) Cost breakdown (same items; PV export revenue is total export)
CAPEX_PV   = NumPV   * unit.PV_kWp   * unit.PV_Capex;
CAPEX_Elec = NumElec * unit.Elec_kW  * unit.Elec_Capex;
CAPEX_Batt = NumBatt * unit.Batt_kWh * unit.Batt_Capex;
CAPEX_CT   = Cap_turb_kW_tot * unit.CSP_Capex_turb_USD_per_kWe;
CAPEX_CS   = Cap_store_kWh   * unit.CSP_Capex_store_USD_per_kWh_e;
CAPEX_CF   = Field_area_m2   * unit.CSP_Capex_field_USD_per_m2;

OPEX_PV    = NumPV   * unit.PV_kWp   * unit.PV_Opex;
OPEX_Elec  = NumElec * unit.Elec_kW  * unit.Elec_Opex;
OPEX_Batt  = NumBatt * unit.Batt_kWh * unit.Batt_Opex;
OPEX_CT    = Cap_turb_kW_tot * unit.CSP_Opex_USD_per_kWe_yr;

Total_Water_m3 = H2_prod_kg_y * (econ.Water_req_L_kgH2/1000);

CostName  = { ...
  'Num_PV'; 'Total_PV_MWp'; ...
  'PV_Area_m2'; 'PV_Area_ha'; 'PV_Area_km2'; ...
  'Num_Elec'; 'Total_Elec_MW'; ...
  'Num_Batt'; 'Total_Batt_MWh'; ...
  'Num_CSP_Turb'; 'CSP_Turbine_MW_e'; ...
  'Num_CSP_Store'; 'CSP_Store_MWh_e'; ...
  'Num_CSP_Field'; 'CSP_Field_m2'; ...
  'CAPEX_PV_USD'; 'CAPEX_Elec_USD'; 'CAPEX_Batt_USD'; ...
  'CAPEX_CSP_Turb_USD'; 'CAPEX_CSP_Store_USD'; 'CAPEX_CSP_Field_USD'; ...
  'OPEX_PV_USD_per_year'; 'OPEX_Elec_USD_per_year'; 'OPEX_Batt_USD_per_year'; 'OPEX_CSP_Turb_USD_per_year'; ...
  'Total_Water_m3_per_year'; ...
  'CAPEX_total_USD'; 'Annualised_CAPEX_USD'; 'Fixed_Opex_USD'; ...
  'Water_cost_USD'; 'Trucking_cost_USD'; ...
  'Revenue_HFO_saved_USD'; 'Revenue_PV_export_USD'; ...
  'NetCost_USD'};

CostValue = [ ...
  NumPV; Total_PV_MWp; ...
  PV_Area_m2; PV_Area_ha; PV_Area_km2; ...
  NumElec; Total_Elec_MW; ...
  NumBatt; Total_Batt_MWh; ...
  NumCSP_Turb; Cap_turb_kW_tot/1e3; ...
  NumCSP_Store; Cap_store_kWh/1e3; ...
  NumCSP_Field; Field_area_m2; ...
  CAPEX_PV; CAPEX_Elec; CAPEX_Batt; ...
  CAPEX_CT; CAPEX_CS; CAPEX_CF; ...
  OPEX_PV; OPEX_Elec; OPEX_Batt; OPEX_CT; ...
  Total_Water_m3; ...
  Capex_USD; AnnualisedCAPEX_USD; Fixed_Opex_USD; ...
  WaterCost_USD; TruckingCost_USD; ...
  Revenue_FuelSaved; ExportRevenue_USD; ...
  NetCost_USD ];

CostTbl = table(CostName, CostValue, 'VariableNames', {'Item','USD_or_units_per_year_or_total'});
writetable(CostTbl, oneFile, 'Sheet', 'CostBreakdown');

% 4) Emissions sheet (same as GA)
EmissName  = {'CO2 avoided (combustion)'; 'Road transport CO2'; 'CH4 avoided'; ...
              'N2O avoided'; 'CO2_Total_y_Site_Load_kWh'};
EmissValue = [CO2_avoided_kg; Road_CO2_kg; CH4_avoided_kg; N2O_avoided_kg; CO2_Total_y_Site_Load_kWh];
EmissUnits = {'kg/yr'; 'kg/yr'; 'kg/yr'; 'kg/yr'; 'kg/yr'};
EmissTbl   = table(EmissName, EmissValue, EmissUnits, ...
                   'VariableNames', {'Pollutant_or_Item','Value','Units'});
writetable(EmissTbl, oneFile, 'Sheet', 'Emissions');

% 5) Assumptions: econ + unit (shows PV_area_m2_per_kWp too)
fnE  = fieldnames(econ);  valE = struct2cell(econ);  tagE = repmat({'econ'}, numel(fnE), 1);
fnU  = fieldnames(unit);  valU = struct2cell(unit);  tagU = repmat({'unit'}, numel(fnU), 1);
Parameter = [fnE; fnU];  Value = [valE; valU];  Source = [tagE; tagU];

try
    AssumpTbl = table(Parameter, cell2mat(Value), Source, ...
                      'VariableNames', {'Parameter','Value','Struct'});
catch
    toStr = @(x) string(x);
    AssumpTbl = table(Parameter, cellfun(toStr, Value, 'UniformOutput', false), Source, ...
                      'VariableNames', {'Parameter','Value','Struct'});
end
writetable(AssumpTbl, oneFile, 'Sheet', 'Assumptions');

%% 1b) Monthly sheet ----------------------------------------------------
% Monthly totals for: PV_kWh, CSP_Field_kWh_e, H2 demand & production, site load
monthlyVars = {'PV_kWh','CSP_Field_kWh_e', ...
               'H2_demand_kg','H2_prod_kg','Site_Load_kWh'};

% (i) Hours → timetable keyed by timestamp
ttHourly = table2timetable(Hourly(:, monthlyVars), ...
                           'RowTimes', Hourly.Datetime);

% (ii) Collapse to calendar-month sums
ttMonthly = retime(ttHourly, 'monthly', 'sum');

% ------- NEW: give clearer names for the two H2 fields
ttMonthly.Properties.VariableNames{'H2_demand_kg'} = 'H2_demand_total_kg';
ttMonthly.Properties.VariableNames{'H2_prod_kg'}   = 'H2_prod_total_kg';

% (iii) Back to a normal table
Monthly = timetable2table(ttMonthly,'ConvertRowTimes',true);
Monthly.Properties.VariableNames{'Time'} = 'Month';

% (iv) Write the sheet
writetable(Monthly, oneFile, 'Sheet', 'Monthly');

fprintf('✓ All results written to %s (Sheets: Hourly, Monthly, Summary, CostBreakdown, Emissions, Assumptions)\n', oneFile);


%fprintf('✓ All results written to %s (Sheets: Hourly, Summary, CostBreakdown, Emissions, Assumptions)\n', oneFile);
