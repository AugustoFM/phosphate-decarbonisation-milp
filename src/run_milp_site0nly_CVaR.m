% run_milp_siteonly_CVaR.m — GA-compatible driver for the MILP (intlinprog)
clear; clc;

% ======================================================================
% ✱ CVaR ADD – LOAD 10-YEAR SCENARIOS
% ======================================================================

% ======================================================================
% ✱ FIX – robust alignment of irradiance and DNI timelines
%      • keeps leap-days if present
%      • inserts NaN where either series is missing an hour
% ======================================================================
irr10 = load('irradiance_hourly_2013-2022.mat');   % TimeUTC + GHI_kWhm2 irradiance_hourly_2013-2022.mat
dni10 = load('dni_hourly_10Y.mat');                % TimeUTC + DNI_kWhm2

% --- build two single-column timetables
ttIrr = timetable(irr10.TimeUTC(:), irr10.GHI_kWhm2(:), ...
                  'VariableNames', {'GHI'});
ttDni = timetable(dni10.TimeUTC(:), dni10.DNI_kWhm2(:), ...
                  'VariableNames', {'DNI'});

% --- union-join (outer join) so every hour in either set is present
tt   = synchronize(ttIrr, ttDni, 'union', 'first');   % keeps time order

% --- report gaps (optional, but useful)
nGapIrr = sum(ismissing(tt.GHI));
nGapDni = sum(ismissing(tt.DNI));
fprintf('After union-join: %d GHI gaps, %d DNI gaps\n', nGapIrr, nGapDni);

% --- fill ≤3-hour holes by linear interpolation, ends by nearest
tt.GHI = fillmissing(tt.GHI,'linear','MaxGap',3,'EndValues','nearest');
tt.DNI = fillmissing(tt.DNI,'linear','MaxGap',3,'EndValues','nearest');

%% ✱ NEW ✱ — remove all leap-day (29-Feb) hours so every year = 8760 h
isLeapHour = (month(tt.Time)==2) & (day(tt.Time)==29);
if any(isLeapHour)
    fprintf('Removing %d leap-day hours (29-Feb) from the horizon…\n', ...
            sum(isLeapHour));
    tt(isLeapHour,:) = [];              % drop those rows
end


% --- replace the original vectors -------------------------------
TimeUTC   = tt.Time;
G_POA     = tt.GHI;          % or tt.GPOA if you stored POA irradiance
DNI_kWhm2 = tt.DNI;

% -------- ensure no residual NaNs (MILP cannot handle them) -------------
G_POA(isnan(G_POA))         = 0;
DNI_kWhm2(isnan(DNI_kWhm2)) = 0;

nH = numel(TimeUTC);         % <- RESTORED: total number of hours
% --- reshape as [T  x  S]  matrix so the MILP sees each year separately
years   = year(TimeUTC);                            % 2013 … 2022
yearList= unique(years,'stable');                   % 1×10
S       = numel(yearList);                          % 10 scenarios
TperYr  = sum(years==yearList(1));                  % 8 760 or 8 784
% ----------------------------------------------------------------


% ------------------------------------------------------------------
% >>> WHICH YEARS DO YOU WANT THE MILP TO SEE?
%      – leave empty  ( [] )  to keep **all** 10 scenarios
%      – list indices (1-based)  or the calendar years themselves
% ------------------------------------------------------------------

%fprintf('Available calendar years:\n');
%disp(yearList);
%userYearRequest = input('List the calendar years (e.g. [2013 2015 2019]) or [] for all: ');


userYearRequest = [2013 2014 2015 2016 2017 2018 2019 2020 2021 2022];                 % e.g.  [2015 2017 2021]  or  [1 3 8]
% ------------------------------------------------------------------

% Convert whatever the user typed into a clean index vector  selYears
if isempty(userYearRequest)
    selYears = 1:numel(yearList);                 % keep everything
elseif all(userYearRequest <= numel(yearList))    % treated as positions
    selYears = unique(userYearRequest(:)).';
else                                              % treated as calendar years
    [tf,pos] = ismember(userYearRequest, yearList);
    selYears = unique(pos(tf));                   % drop any year not present
    assert(~isempty(selYears), 'None of the requested years exist in yearList');
end


% Hydrogen demand (hourly)
dem = load('hourly_H2_demand_100.mat');  % expects H2_need_hour
assert(isfield(dem,'H2_need_hour'), 'hourly_H2_demand_100.mat must contain H2_need_hour');
H2_need_hour = dem.H2_need_hour(:);
% --- NEW: replicate to full 10-yr horizon if needed
if numel(H2_need_hour) == TperYr            % one calendar year only
    H2_need_hour = repmat(H2_need_hour, S, 1);    % S = 10
end

% now the usual length check
assert(numel(H2_need_hour) == nH, ...
      'H2_need_hour length mismatch vs TimeUTC');

% Diesel hourly energy (site load proxy), prefer .mat, fallback .xlsx
% ------------------------------------------------------------------
%  Diesel hourly energy (site load proxy) – prefer .mat, fallback .xlsx
% ------------------------------------------------------------------
Diesel_hourly_energy_kWh = [];

load_single_year = @(v) repmat(v(:), S, 1);   % helper

if exist('diesel_hourly_energy.mat','file')
    d = load('diesel_hourly_energy.mat');      % needs Energy_produced_Diesel_kWh
    diesel_vec = d.Energy_produced_Diesel_kWh(:);

    if numel(diesel_vec)==TperYr                     % *** one calendar year ***
        Diesel_hourly_energy_kWh = load_single_year(diesel_vec);

    elseif numel(diesel_vec)==nH                     % already multi-year
        Diesel_hourly_energy_kWh = diesel_vec;

    else
        error('Diesel profile length (%d) is neither 8 760 nor %d hours.', ...
              numel(diesel_vec), nH);
    end

elseif exist('diesel_hourly_energy.xlsx','file')
    Tdies = readtable('diesel_hourly_energy.xlsx','VariableNamingRule','preserve');
    numMask = varfun(@isnumeric,Tdies,'OutputFormat','uniform');
    c = find(numMask,1);
    diesel_vec = Tdies{:,c};

    if numel(diesel_vec)==TperYr
        Diesel_hourly_energy_kWh = load_single_year(diesel_vec);

    elseif numel(diesel_vec)==nH
        Diesel_hourly_energy_kWh = diesel_vec;

    else
        error('Diesel XLSX has %d rows – expecting 8 760 or %d.', ...
              numel(diesel_vec), nH);
    end

else
    error('Neither diesel_hourly_energy.mat nor .xlsx found.');
end


Site_Load_kWh = Diesel_hourly_energy_kWh(:);
H2_need_kg    = H2_need_hour(:);

% contiguous-year sanity check ------------------------------------------
dYear = diff(years).';                          % 1 × (nH-1)   row
mask  = mod((1:TperYr:S*TperYr)',TperYr)==1;    %  S × 1       column

assert(all(all(dYear==1 | mask)), ...
       'Non-contiguous yearly blocks – check leap-year trimming.');

G_POA   = reshape(G_POA  ,TperYr,S);
DNI_kWhm2 = reshape(DNI_kWhm2,TperYr,S);

Site_Load_kWh = reshape(Site_Load_kWh,TperYr,S);
H2_need_kg    = reshape(H2_need_kg   ,TperYr,S);

% ----------------------------------------------------------------
% >>> KEEP ONLY 3 SCENARIOS (years) instead of all 10 <<<
%     – here we keep the first three calendar years that appear in
%       yearList, but you can pick any three indices you like.
% ----------------------------------------------------------------
% ----------------------------------------------------------------
%  Generic trimming to the selected calendar years
% ----------------------------------------------------------------
% a) keep the  [T×S]  matrices
G_POA         = G_POA(:, selYears);
DNI_kWhm2     = DNI_kWhm2(:, selYears);
Site_Load_kWh = Site_Load_kWh(:, selYears);
H2_need_kg    = H2_need_kg(:, selYears);

% b) rebuild the *flat* vectors that drive the Hourly sheet later
idxSel        = ismember(years, yearList(selYears));  % logical mask
TimeUTC       = TimeUTC(idxSel);
G_POA_flat    = G_POA(:);      % already in correct order (column-wise)
DNI_flat      = DNI_kWhm2(:);
Site_flat     = Site_Load_kWh(:);
H2_flat       = H2_need_kg(:);

% c) update scenario count and probabilities
S             = numel(selYears);
p_s           = ones(S,1) / S;       % equal probability
T             = TperYr;              % name convenience − solver expects T×S
nH            = numel(TimeUTC);      % *** this is the new total hours ***

% ----------------------------------------------------------------


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

econ.price_CO2_emission = 0.00; econ.price_CH4_emission = 0.00; econ.price_N2O_emission = 0; %Recommended prices for KAPSARC 2030 0.05 CO2, 0.54 CH4, 5.46 N20

econ.cost_of_CAPEX = 0.02;
econ.unserved_cost_USD_per_kWh = 0; %50USD suggested

%% D. (Optional) short horizon for a smoke test the short horizon does not work because it has only one scenario, we might have to add at least one more scenario
use_short_horizon = false; hours_to_keep = 24*3;
if use_short_horizon
    idxKeep = 1:min(hours_to_keep, nH);
    TimeUTC = TimeUTC(idxKeep);
    G_POA   = G_POA(idxKeep);
    DNI_kWhm2 = DNI_kWhm2(idxKeep);
    Site_Load_kWh = Site_Load_kWh(idxKeep);
    H2_need_kg    = H2_need_kg(idxKeep);
    fprintf('Using short horizon of %d hours.\n', numel(idxKeep));
    nH = numel(TimeUTC);     % keep nH consistent after trimming
    S      = 1;                % one "scenario" only
    TperYr = nH;               % length of the toy horizon
    p_s    = 1;                % probability vector length = S

    fprintf('Using short horizon of %d hours (TperYr = %d, S = %d).\n', ...
            nH, TperYr, S);

end

%% E. Solve MILP (objective matches GA NetCost)
%fprintf('Solving MILP with intlinprog... (8760 h can take a while)\n');
%tic;
%sol = build_and_solve_milp_v2_lex_CVaR( ...
%    TimeUTC(1:TperYr), G_POA, DNI_kWhm2, ...
%     Site_Load_kWh, H2_need_kg, unit, econ,p_s); %econ, p_s 

%% E. Solve MILP (site-only CVaR)
fprintf('Solving MILP with intlinprog... (%d h × %d scenarios)\n', T, S);
tic;
sol = build_and_solve_milp_siteonly_CVaR( ...
        TimeUTC(1:T), ...      % single-year stamp axis
        G_POA, DNI_kWhm2, ...
        Site_Load_kWh, ...
        unit, econ, p_s);

if ~isfield(sol,'exitflag') || sol.exitflag <= 0
    warning('MILP infeasible or no solution returned (exitflag=%d).', ...
            getfield(sol,'exitflag',-Inf));
    return;
end
t_solve = toc;

%% F. Derived sizes, capacities, areas (site-only)
NumPV        = sol.NumPV;
NumBatt      = sol.NumBatt;
NumCSP_Turb  = sol.NumCSP_Turb;
NumCSP_Store = sol.NumCSP_Store;
NumCSP_Field = sol.NumCSP_Field;

Total_PV_kWp   = NumPV        * unit.PV_kWp;           Total_PV_MWp   = Total_PV_kWp/1e3;
Total_Batt_kWh = NumBatt      * unit.Batt_kWh;         Total_Batt_MWh = Total_Batt_kWh/1e3;
Cap_turb_kW_tot= NumCSP_Turb  * unit.CSP_turb_unit_kW_e;
Cap_store_kWh  = NumCSP_Store * unit.CSP_store_unit_kWh_e;
Field_area_m2  = NumCSP_Field * unit.CSP_field_unit_m2;

% PV area (for reporting)
PV_Area_m2  = Total_PV_kWp * unit.PV_area_m2_per_kWp;
PV_Area_ha  = PV_Area_m2 / 1e4;
PV_Area_km2 = PV_Area_m2 / 1e6;

% Served/unserved & export (T×S)
Site_served_kWh    = Site_Load_kWh - sol.Unserved_Site_kWh;
GridExport_kWh     = sol.Export_or_Dump_kWh;          % PV_Dump + CSP_Dump
Unserved_Total_kWh = sol.Unserved_Site_kWh;           % site-only

% Optional informational H2 from surplus (does NOT affect sizing)
kWh2kgH2        = unit.elec_eta / 33.3;
H2_from_surplus_kg = sol.H2_surplus_kg;               % from solver (kWh * η/33.3)

% "Turbines_On" estimate (continuous)
Turbines_On_est = sol.Turb_kWh ./ max(unit.CSP_turb_unit_kW_e, eps);

%% H. Print concise summary (site-only)
fprintf('\n===== MILP SOLUTION (Site-only, CVaR) =====\n');
fprintf('Status: exitflag=%d | Solve time: %.1f s\n', sol.exitflag, t_solve);
fprintf('Sizes  : PV=%d kWp | BESS=%d x%dkWh=%.1f MWh\n', ...
        NumPV, NumBatt, unit.Batt_kWh, Total_Batt_MWh);
fprintf('CSP    : Turb=%d (%.1f MW_e) | TES=%d (%.1f MWh_e) | Field=%d (%.0f m^2)\n', ...
        NumCSP_Turb, Cap_turb_kW_tot/1e3, NumCSP_Store, Cap_store_kWh/1e3, NumCSP_Field, Field_area_m2);
fprintf('Info   : Surplus-H2=%.2f kt/yr | Export=%.2f GWh | Unserved=%.2f GWh\n', ...
        sum(H2_from_surplus_kg,'all')/1e6, sum(GridExport_kWh,'all')/1e6, sum(Unserved_Total_kWh,'all')/1e6);



%% I. Write EVERYTHING to ONE Excel file (site-only shapes)
stamp   = char(datetime('now','Format','yyyyMMddHHmm'));
oneFile = sprintf('%s_H2_CSP_System_AllResults_SITEONLY_CVaR.xlsx', stamp);

%% G. Build Hourly sheet (site-only)
% Flatten all T×S matrices column-wise to match time stacking
Hourly = table( ...
    TimeUTC, ...
    G_POA_flat, ...
    DNI_flat, ...
    repmat(sol.NumPV       , nH,1), ...
    (sol.PV_to_Load_kWh(:) + sol.PV_to_Batt_kWh(:) + sol.PV_Dump_kWh(:)), ... % PV_kWh produced
    repmat(sol.NumCSP_Field, nH,1), ...
    (sol.Field_to_Turb_kWh(:) + sol.Field_to_TES_kWh(:)), ...                  % CSP field energy (e-eq)
    repmat(sol.NumCSP_Turb , nH,1), ...
    sol.Turb_kWh(:), ...
    sol.Turb_kWh(:)/unit.CSP_turb_unit_kW_e, ...
    sol.Field_to_TES_kWh(:), ...
    sol.TES_to_Turb_kWh(:), ...
    sol.ETES_kWh(:), ...
    repmat(sol.NumBatt     , nH,1), ...
    sol.Ebatt_kWh(:), ...
    sol.PV_to_Batt_kWh(:), ...
    sol.CSP_to_Batt_kWh(:), ...
    sol.Batt_to_Load_kWh(:), ...
    Site_flat, ...
    (Site_flat - sol.Unserved_Site_kWh(:)), ...
    sol.Unserved_Site_kWh(:), ...
    sol.Export_or_Dump_kWh(:), ...
    H2_from_surplus_kg(:), ...
    'VariableNames', { ...
      'Datetime','Irrad_kWh_kWp','DNI_kWhm2', ...
      'NumPV','PV_kWh', ...
      'NumCSP_Field','CSP_Field_kWh_e', ...
      'NumCSP_Turb','CSP_Turbine_kWh','Turbines_On', ...
      'CSP_Charge_kWh','CSP_Discharge_kWh','CSP_SOC_kWh', ...
      'NumBatt','Batt_SOC_kWh', ...
      'PV_to_Batt_kWh','CSP_to_Batt_kWh','Batt_to_Load_kWh', ...
      'Site_Load_kWh','Site_Served_kWh','Unserved_Site_kWh', ...
      'GridExport_kWh','H2_from_surplus_kg'} );

writetable(Hourly, oneFile, 'Sheet', 'Hourly');


%% Summary (site-only)
% Probability-weighted totals
Total_Load_kWh_yr    = sum((Site_Load_kWh   ) * p_s);             % scalar
Unserved_Total_kWh_y = sum((sol.Unserved_Site_kWh) * p_s);
Served_Total_kWh_y   = Total_Load_kWh_yr - Unserved_Total_kWh_y;
PV_export_kWh_y      = sum((GridExport_kWh) * p_s);
H2_surplus_kg_y      = sum((H2_from_surplus_kg) * p_s);

% Economics (no electrolyser CAPEX/OPEX/water; export optional)
P_Turb_kW   = Cap_turb_kW_tot;
E_TES_kWh   = Cap_store_kWh;
A_Field_m2  = Field_area_m2;

CAPEX_PV   = NumPV   * unit.PV_kWp   * unit.PV_Capex;
CAPEX_Batt = NumBatt * unit.Batt_kWh * unit.Batt_Capex;
CAPEX_CT   = P_Turb_kW * unit.CSP_Capex_turb_USD_per_kWe;
CAPEX_CS   = E_TES_kWh * unit.CSP_Capex_store_USD_per_kWh_e;
CAPEX_CF   = A_Field_m2 * unit.CSP_Capex_field_USD_per_m2;

Capex_USD          = CAPEX_PV + CAPEX_Batt + CAPEX_CT + CAPEX_CS + CAPEX_CF;
AnnualisedCAPEX_USD= econ.cost_of_CAPEX * Capex_USD;

OPEX_PV    = NumPV   * unit.PV_kWp   * unit.PV_Opex;
OPEX_Batt  = NumBatt * unit.Batt_kWh * unit.Batt_Opex;
OPEX_CT    = P_Turb_kW * unit.CSP_Opex_USD_per_kWe_yr;

Fixed_Opex_USD     = OPEX_PV + OPEX_Batt + OPEX_CT;

ExportRevenue_USD  = PV_export_kWh_y * econ.Price_export_energy_kWh_USD;
Unserved_Penalty_USD = econ.unserved_cost_USD_per_kWh * Unserved_Total_kWh_y;

% Optional informational CO2 credit vs diesel for site load (not in NetCost)
CO2_Total_y_Site_Load_kWh = Total_Load_kWh_yr * econ.CO2_red_factor_PVvsDiesel;

% Site-only NetCost (no H2 credits/costs)
NetCost_USD = AnnualisedCAPEX_USD + Fixed_Opex_USD + Unserved_Penalty_USD ...
              - ExportRevenue_USD;

metricNames = { ...
  'Total_Load_kWh'
  'Served_Total_kWh'
  'Unserved_Total_kWh'
  'PV_export_kWh'
  'H2_from_surplus_kg (info)'
  'PV_area_m2'
  'PV_area_ha'
  'PV_area_km2'
  'Total_CAPEX_USD'
  'Annualised_CAPEX_USD'
  'Fixed_Opex_USD'
  'Export_Revenue_USD'
  'Unserved_Penalty_USD'
  'NetCost_USD'
  'CO2_Total_y_Site_Load_kWh (info)'};  % informational, not in NetCost

metricValues = [ ...
  Total_Load_kWh_yr
  Served_Total_kWh_y
  Unserved_Total_kWh_y
  PV_export_kWh_y
  H2_surplus_kg_y
  PV_Area_m2
  PV_Area_ha
  PV_Area_km2
  Capex_USD
  AnnualisedCAPEX_USD
  Fixed_Opex_USD
  ExportRevenue_USD
  Unserved_Penalty_USD
  NetCost_USD
  CO2_Total_y_Site_Load_kWh ];

SummaryTbl = table(metricNames, metricValues, 'VariableNames', {'Metric','Value_yearly'});
writetable(SummaryTbl, oneFile, 'Sheet', 'Summary');


%% Per-scenario rollups (site-only)
slice  = @(M,s) M(:,s);
yearsS = yearList(selYears);

Scenario = {};  Metric = {};  Value = [];

for s = 1:S
    tag   = sprintf('Y%d', yearsS(s));
    loadS = slice(Site_Load_kWh, s);
    uS    = slice(sol.Unserved_Site_kWh, s);
    exS   = slice(GridExport_kWh, s);
    h2S   = slice(H2_from_surplus_kg, s);

    Scenario = [Scenario; repmat({tag}, 5,1)];
    Metric   = [Metric  ; {'Site_Load_kWh'; 'Site_Served_kWh'; ...
                           'Unserved_Site_kWh'; 'Export_or_Dump_kWh'; ...
                           'H2_from_surplus_kg'}];
    Value    = [Value   ; sum(loadS); sum(loadS - uS); sum(uS); sum(exS); sum(h2S)];
end

PerScenarioTbl = table(Scenario, Metric, Value);
writetable(PerScenarioTbl, oneFile, 'Sheet','PerScenario');


%% Cost breakdown (site-only)
CostName  = { ...
  'Num_PV'; 'Total_PV_MWp'; ...
  'PV_Area_m2'; 'PV_Area_ha'; 'PV_Area_km2'; ...
  'Num_Batt'; 'Total_Batt_MWh'; ...
  'Num_CSP_Turb'; 'CSP_Turbine_MW_e'; ...
  'Num_CSP_Store'; 'CSP_Store_MWh_e'; ...
  'Num_CSP_Field'; 'CSP_Field_m2'; ...
  'CAPEX_PV_USD'; 'CAPEX_Batt_USD'; ...
  'CAPEX_CSP_Turb_USD'; 'CAPEX_CSP_Store_USD'; 'CAPEX_CSP_Field_USD'; ...
  'OPEX_PV_USD_per_year'; 'OPEX_Batt_USD_per_year'; 'OPEX_CSP_Turb_USD_per_year'; ...
  'CAPEX_total_USD'; 'Annualised_CAPEX_USD'; 'Fixed_Opex_USD'; ...
  'Export_Revenue_USD'; 'Unserved_Penalty_USD'; ...
  'NetCost_USD'};

CostValue = [ ...
  NumPV; Total_PV_MWp; ...
  PV_Area_m2; PV_Area_ha; PV_Area_km2; ...
  NumBatt; Total_Batt_MWh; ...
  NumCSP_Turb; Cap_turb_kW_tot/1e3; ...
  NumCSP_Store; Cap_store_kWh/1e3; ...
  NumCSP_Field; Field_area_m2; ...
  CAPEX_PV; CAPEX_Batt; ...
  CAPEX_CT; CAPEX_CS; CAPEX_CF; ...
  OPEX_PV; OPEX_Batt; OPEX_CT; ...
  Capex_USD; AnnualisedCAPEX_USD; Fixed_Opex_USD; ...
  ExportRevenue_USD; Unserved_Penalty_USD; ...
  NetCost_USD ];

CostTbl = table(CostName, CostValue, 'VariableNames', {'Item','USD_or_units_per_year_or_total'});
writetable(CostTbl, oneFile, 'Sheet', 'CostBreakdown');


%% Emissions (simple, informational)
EmissName  = {'CO2_Total_y_Site_Load_kWh (info)'};
EmissValue = [CO2_Total_y_Site_Load_kWh];
EmissUnits = {'kg/yr'};
EmissTbl   = table(EmissName, EmissValue, EmissUnits, ...
                   'VariableNames', {'Pollutant_or_Item','Value','Units'});
writetable(EmissTbl, oneFile, 'Sheet', 'Emissions');

%% Monthly sheet (site-only)
monthlyVars = {'PV_kWh','CSP_Field_kWh_e', ...
               'Site_Load_kWh','Site_Served_kWh','GridExport_kWh','H2_from_surplus_kg'};

ttHourly  = table2timetable(Hourly(:, monthlyVars), 'RowTimes', Hourly.Datetime);
ttMonthly = retime(ttHourly, 'monthly', 'sum');

% friendlier name for H2 column
ttMonthly.Properties.VariableNames{'H2_from_surplus_kg'} = 'H2_from_surplus_total_kg';

Monthly = timetable2table(ttMonthly,'ConvertRowTimes',true);
Monthly.Properties.VariableNames{'Time'} = 'Month';
writetable(Monthly, oneFile, 'Sheet', 'Monthly');

fprintf('✓ All results written to %s (Sheets: Hourly, Monthly, Summary, PerScenario, CostBreakdown, Emissions, Assumptions)\n', oneFile);
