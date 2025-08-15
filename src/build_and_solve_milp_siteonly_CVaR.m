function sol = build_and_solve_milp_siteonly_CVaR(TimeUTC, G_POA, DNI_kWhm2, ...
                                                  Site_Load_kWh, unit, econ, p_s)
% Site-only, CVaR-aware MILP for PV + CSP (field+TES+turbine) + BESS.
% Sizes ONLY for Site_Load_kWh; hydrogen is NOT modelled in the optimization.
% Any H2 reported is post-hoc from surplus (export/dump) energy and cannot
% affect sizing.

% ---------------------- dimensions / probabilities -----------------------
[T,S] = size(G_POA);
assert(all(size(DNI_kWhm2)==[T,S]));
if nargin < 7 || isempty(p_s), p_s = ones(S,1)/S; end
p_s = p_s(:);  assert(numel(p_s)==S);

% ---------------------- default econ fields (used here) ------------------
defaults = { ...
 'cost_of_CAPEX', 0.02; ...
 'Price_export_energy_kWh_USD', 0.00; ...
 'unserved_cost_USD_per_kWh', 50; ...
 'CO2_red_factor_PVvsDiesel', 0.78 };
for i=1:size(defaults,1)
    if ~isfield(econ, defaults{i,1}), econ.(defaults{i,1}) = defaults{i,2}; end
end

% ---------------------- constants / efficiencies -------------------------
eta_rt = unit.batt_eta;   eta_c = sqrt(eta_rt);   eta_d = sqrt(eta_rt);
csp_e_per_m2 = unit.CSP_opt_eff * unit.CSP_therm_eff * unit.CSP_cycle_eff;

% ---------------------- variables & indexing -----------------------------
% Hourly variables (lb >= 0):
%  1 pvL      PV -> Site
%  2 pvB      PV -> BESS (charge)
%  3 pvDump   PV -> Dump/Export (surplus rail)
%  4 f2T      CSP field -> Turbine
%  5 f2S      CSP field -> TES (charge)
%  6 s2T      TES -> Turbine (discharge)
%  7 yTurb    Turbine electric output
%  8 cspL     Turbine -> Site
%  9 cspB     Turbine -> BESS (charge)
% 10 cspDump  Turbine -> Dump/Export
% 11 Eb       BESS SOC
% 12 b2L      BESS -> Site (discharge)
% 13 ETES     TES SOC (electric-equiv)
% 14 uSite    Site unmet (kWh)
nh = 14;

idx = @(t,s,k) ((s-1)*T + (t-1))*nh + k;     % helper

NHV = nh*T*S;
NV  = NHV;                                    % pointer to add more vars

% CVaR auxiliaries
z     = NV + 1;
eps_s = (NV+1 : NV+S).';
NV    = NV + 1 + S;

% Design-size integers (NO electrolyser here)
nPV    = NV + 1;
nBatt  = NV + 2;
nTurb  = NV + 3;
nStore = NV + 4;
nField = NV + 5;
NV     = NV + 5;

intcon = [nPV nBatt nTurb nStore nField];

% ---------------------- bounds ------------------------------------------
lb = zeros(NV,1);  ub = inf(NV,1);
lb(z)     = -inf;  % VaR level can be negative
lb(eps_s) = 0;

lb([nPV nBatt nTurb nStore nField]) = 0;
ub([nPV nBatt nTurb nStore nField]) = [3e5; 2000; 50; 2000; 2e6];

% ---------------------- sparse builders ----------------------------------
Ai=[]; Aj=[]; Av=[]; b=[]; ineq_row=0;
AEi=[];AEj=[];AEv=[]; beq=[]; eq_row=0;
    function addI(cols, vals, rhs)
        ineq_row = ineq_row + 1;
        Ai = [Ai; repmat(ineq_row,numel(cols),1)];
        Aj = [Aj; cols(:)];
        Av = [Av; vals(:)];
        b(ineq_row,1) = rhs;
    end
    function addE(cols, vals, rhs)
        eq_row = eq_row + 1;
        AEi = [AEi; repmat(eq_row,numel(cols),1)];
        AEj = [AEj; cols(:)];
        AEv = [AEv; vals(:)];
        beq(eq_row,1) = rhs;
    end

% ---------------------- physics & network --------------------------------
for s = 1:S
  for t = 1:T
    % PV availability
    addI([idx(t,s,1) idx(t,s,2) idx(t,s,3) nPV], ...
         [1 1 1, -G_POA(t,s)*unit.inv_eff*unit.PV_kWp], 0);

    % CSP field availability
    addI([idx(t,s,4) idx(t,s,5) nField], ...
         [1 1, -DNI_kWhm2(t,s)*unit.CSP_field_unit_m2*csp_e_per_m2], 0);

    % Turbine: power and split
    addE([idx(t,s,7) idx(t,s,4) idx(t,s,6)],[1 -1 -1],0);            % yTurb=f2T+s2T
    addI([idx(t,s,7) nTurb],[1, -unit.CSP_turb_unit_kW_e],0);        % yTurb ≤ cap
    addE([idx(t,s,8) idx(t,s,9) idx(t,s,10) idx(t,s,7)],[1 1 1 -1],0);% cspL+cspB+cspDump=yTurb

    % TES dyn/cap/charge-cap (cyclic)
    tnext = (t<T)*(t+1) + (t==T)*1;
    addE([idx(tnext,s,13) idx(t,s,13) idx(t,s,5) idx(t,s,6)],[1 -1 -1 1],0);
    addI([idx(t,s,13) nStore],[1, -unit.CSP_store_unit_kWh_e],0);
    addI([idx(t,s,5)  nStore],[1, -unit.CSP_store_charge_kW_per_unit],0);

    % Battery dyn/cap/power caps (cyclic)
    addE([idx(tnext,s,11) idx(t,s,11) idx(t,s,2) idx(t,s,9) idx(t,s,12)], ...
         [1 -1  eta_c      eta_c     -1/eta_d], 0);
    addI([idx(t,s,11) nBatt],[ 1, -unit.Batt_kWh*unit.Batt_SOC_max],0);
    addI([idx(t,s,11) nBatt],[-1,  unit.Batt_kWh*unit.Batt_SOC_min],0);
    addI([idx(t,s,12) nBatt],[ 1, -unit.Batt_maxDis],0);              % b2L ≤ dis cap
    addI([idx(t,s,2)  nBatt],[ 1, -unit.Batt_maxChg],0);              % pvB ≤ ch cap
    addI([idx(t,s,9)  nBatt],[ 1, -unit.Batt_maxChg],0);              % cspB ≤ ch cap

    % Site balance
    addE([idx(t,s,1) idx(t,s,8) idx(t,s,12) idx(t,s,14)], [1 1 1 1], Site_Load_kWh(t));
  end
end

A   = sparse(Ai,Aj,Av,ineq_row,NV);
Aeq = sparse(AEi,AEj,AEv,eq_row,NV);

% ---------------------- PASS 1: minimize total unmet site ----------------
f1 = zeros(NV,1);
for s=1:S, for t=1:T, f1(idx(t,s,14)) = 1; end, end    % sum uSite

opts = optimoptions('intlinprog','Display','off','MaxTime',3600);
[xx1,fval1,exitflag1,output1] = intlinprog(f1,intcon,A,b,Aeq,beq,lb,ub,opts);

sol.pass1.exitflag   = exitflag1;
sol.pass1.output     = output1;
sol.pass1.uSite_sum_min = fval1;
if exitflag1 <= 0
    warning('Pass-1 (min uSite) failed (exitflag=%d).', exitflag1);
    sol.exitflag = exitflag1;  return;
end

% ---------------------- PASS 2: cost + λ·CVaR(unserved) ------------------
% Fix total unmet site energy to U* (from pass 1)
uSite_all = [];
for s=1:S, uSite_all = [uSite_all; idx(1:T,s,14)]; end %#ok<AGROW>
row   = size(Aeq,1)+1;
AEi2  = [AEi ; repmat(row,numel(uSite_all),1)];
AEj2  = [AEj ; uSite_all(:)];
AEv2  = [AEv ; ones(numel(uSite_all),1)];
beq2  = [beq ; fval1];
Aeq2  = sparse(AEi2,AEj2,AEv2,row,NV);

% --- Expected cost vector (sizes + hourly revenue/penalties) -------------
f_cost = zeros(NV,1);

% Annualised CAPEX (+) and fixed O&M (+)
capex_per_n = [ ...
  unit.PV_kWp   * unit.PV_Capex; ...
  unit.Batt_kWh * unit.Batt_Capex; ...
  unit.CSP_turb_unit_kW_e   * unit.CSP_Capex_turb_USD_per_kWe; ...
  unit.CSP_store_unit_kWh_e * unit.CSP_Capex_store_USD_per_kWh_e; ...
  unit.CSP_field_unit_m2    * unit.CSP_Capex_field_USD_per_m2 ];
opex_per_n  = [ ...
  unit.PV_kWp   * unit.PV_Opex; ...
  unit.Batt_kWh * unit.Batt_Opex; ...
  unit.CSP_turb_unit_kW_e   * unit.CSP_Opex_USD_per_kWe_yr; ...
  0; 0];

f_cost([nPV nBatt nTurb nStore nField]) = ...
    econ.cost_of_CAPEX * capex_per_n + opex_per_n;

% Hourly export revenue (dump rails) and optional unmet penalty (USD)
for s=1:S
  for t=1:T
    f_cost(idx(t,s, 3)) = f_cost(idx(t,s, 3)) - econ.Price_export_energy_kWh_USD;  % PV dump
    f_cost(idx(t,s,10)) = f_cost(idx(t,s,10)) - econ.Price_export_energy_kWh_USD;  % CSP dump
    f_cost(idx(t,s,14)) = f_cost(idx(t,s,14)) + econ.unserved_cost_USD_per_kWh;    % unmet site
  end
end

% Expected value over scenarios
f_exp = zeros(NV,1);
for s=1:S
  base = (s-1)*T*nh;
  rng  = base + (1:T*nh);
  f_exp(rng) = p_s(s) * f_cost(rng);
end
f_exp(NHV+1:NV) = f_cost(NHV+1:NV);   % sizes, z, eps_s untouched by p_s

% --- CVaR on monetized unmet site energy per scenario --------------------
lambda = 1.0;     % risk aversion
alpha  = 0.95;    % CVaR level

% loss selector: only uSite (monetized) contributes to scenario "loss"
f_loss = zeros(NV,1);
for s=1:S, for t=1:T, f_loss(idx(t,s,14)) = econ.unserved_cost_USD_per_kWh; end, end

% CVaR linear term
f_cvar = lambda * [ ...
   zeros(NHV,1);       % hourly vars
   1;                  % z
   (1/(1-alpha))*p_s;  % eps_s
   zeros(5,1) ];       % 5 integer sizes

% Full objective (EV + CVaR)
f2 = f_exp + f_cvar;

% CVaR constraints:  (loss_s - z ≤ eps_s)
for s=1:S
  ineq_row = ineq_row + 1;
  base     = (s-1)*T*nh;
  loss_cols = base + (1:T*nh);
  loss_vals = f_loss(loss_cols);
  Ai = [Ai; repmat(ineq_row,numel(loss_cols)+2,1)];
  Aj = [Aj; loss_cols(:); z; eps_s(s)];
  Av = [Av; loss_vals(:); -1; -1];
  b (ineq_row,1) = 0;
end
A2 = sparse(Ai,Aj,Av,ineq_row,NV);

[xx,fval,exitflag,output] = intlinprog(f2,intcon,A2,b,Aeq2,beq2,lb,ub,opts);

% ---------------------- unpack ------------------------------------------
sol.exitflag = exitflag;
sol.output   = output;
sol.fval_USD = fval;

% Sizes
sol.NumPV        = round(xx(nPV));
sol.NumBatt      = round(xx(nBatt));
sol.NumCSP_Turb  = round(xx(nTurb));
sol.NumCSP_Store = round(xx(nStore));
sol.NumCSP_Field = round(xx(nField));

% Hourly series
get = @(k) reshape(xx(k:nh:nh*T*S), T, S);

sol.PV_to_Load_kWh   = get(1);
sol.PV_to_Batt_kWh   = get(2);
sol.PV_Dump_kWh      = get(3);

sol.Field_to_Turb_kWh= get(4);
sol.Field_to_TES_kWh = get(5);
sol.TES_to_Turb_kWh  = get(6);
sol.Turb_kWh         = get(7);

sol.CSP_to_Load_kWh  = get(8);
sol.CSP_to_Batt_kWh  = get(9);
sol.CSP_Dump_kWh     = get(10);

sol.Ebatt_kWh        = get(11);
sol.Batt_to_Load_kWh = get(12);
sol.ETES_kWh         = get(13);
sol.Unserved_Site_kWh= get(14);

sol.Export_or_Dump_kWh = sol.PV_Dump_kWh + sol.CSP_Dump_kWh;

% Per-scenario totals (expected with p_s)
sol.Unserved_total_kWh = sum(sol.Unserved_Site_kWh,1);
sol.Export_total_kWh   = sum(sol.Export_or_Dump_kWh,1);

% Economics (compact)
Cap_turb_kW = sol.NumCSP_Turb  * unit.CSP_turb_unit_kW_e;
Cap_TES_kWh = sol.NumCSP_Store * unit.CSP_store_unit_kWh_e;
A_Field_m2  = sol.NumCSP_Field * unit.CSP_field_unit_m2;

Capex_USD = ...
    sol.NumPV   * unit.PV_kWp   * unit.PV_Capex + ...
    sol.NumBatt * unit.Batt_kWh * unit.Batt_Capex + ...
    Cap_turb_kW * unit.CSP_Capex_turb_USD_per_kWe + ...
    Cap_TES_kWh * unit.CSP_Capex_store_USD_per_kWh_e + ...
    A_Field_m2  * unit.CSP_Capex_field_USD_per_m2;

Fixed_Opex_USD = ...
    sol.NumPV   * unit.PV_kWp   * unit.PV_Opex + ...
    sol.NumBatt * unit.Batt_kWh * unit.Batt_Opex + ...
    Cap_turb_kW * unit.CSP_Opex_USD_per_kWe_yr;

AnnualisedCAPEX_USD = econ.cost_of_CAPEX * Capex_USD;

% Expected export revenue
ExpRev_USD = econ.Price_export_energy_kWh_USD * sum(sum(sol.Export_or_Dump_kWh .* p_s'),1);

% Expected unmet penalty (same coefficient used in f_loss)
Unserved_USD = econ.unserved_cost_USD_per_kWh * sum(sum(sol.Unserved_Site_kWh .* p_s'),1);

sol.cost_breakdown = struct( ...
  'AnnualisedCAPEX_USD', AnnualisedCAPEX_USD, ...
  'Fixed_Opex_USD',      Fixed_Opex_USD, ...
  'ExportRevenue_USD',   ExpRev_USD, ...
  'Unserved_USD',        Unserved_USD );

% ---------------------- post-hoc H2 from surplus (does NOT affect sizing)
kWh2kgH2 = unit.elec_eta / 33.3;
sol.H2_surplus_kg = kWh2kgH2 .* sol.Export_or_Dump_kWh;  % T×S series
sol.H2_surplus_kg_total = sum(sol.H2_surplus_kg, 'all');

end
