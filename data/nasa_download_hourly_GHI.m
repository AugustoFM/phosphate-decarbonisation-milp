% nasa_download_hourly_GHI.m
% -------------------------------------------------------------------------
% Downloads NASA‑POWER ALLSKY_SFC_SW_DWN (hourly GHI) for one calendar year
% and stores it as *energy* per square metre for each hour (kWh m‑2).
%
% Output file : irradiance_hourly.mat
%               • GHI_kWhm2   – 8 760 × 1  numeric   (kWh m‑2 per hour)
%               • TimeUTC     – 8 760 × 1  datetime  (UTC, hourly)
% -------------------------------------------------------------------------
clear; clc;

% ------------ location & year -------------------------------------------
lat  = 31.200;                  % ° N
lon  = 39.850;                  % ° E
year = 2022;                    % calendar year

% ------------ NASA POWER query ------------------------------------------
startDate = sprintf('%d0101', year);
endDate   = sprintf('%d1231', year);

url = sprintf([ ...
   'https://power.larc.nasa.gov/api/temporal/hourly/point?' ...
   'parameters=ALLSKY_SFC_SW_DWN&community=RE&latitude=%.4f&longitude=%.4f&' ...
   'start=%s&end=%s&format=JSON'], ...
   lat, lon, startDate, endDate);

fprintf('\nRequesting NASA POWER hourly GHI …\n%s\n\n', url);

raw  = webread(url, weboptions('Timeout',60,'ContentType','text'));
json = jsondecode(raw);

% ------------ pick the correct branch (depends on POWER version) --------
if isfield(json,'features')
    ghiStruct = json.features(1).properties.parameter.ALLSKY_SFC_SW_DWN;
else
    ghiStruct = json.properties.parameter.ALLSKY_SFC_SW_DWN;
end

keys   = sort(fieldnames(ghiStruct));   % 'x2022010100', …
nH     = numel(keys);

GHI_Wm2 = zeros(nH,1);
TimeUTC = NaT(nH,1,'TimeZone','UTC');

for k = 1:nH
    tag          = keys{k};
    GHI_Wm2(k)   = ghiStruct.(tag);                % mean W m‑2 in that hour
    stamp        = tag(regexp(tag,'\d'));          % strip leading 'x'
    TimeUTC(k)   = datetime(stamp,'InputFormat','yyyyMMddHH','TimeZone','UTC');
end

% --- convert "mean power" to "hourly energy" -----------------------------
% mean W m‑2 × 1 h → kWh m‑2  (divide by 1000 then implicitly ×1 h)
GHI_kWhm2 = GHI_Wm2 / 1000;

save('irradiance_hourly.mat','GHI_kWhm2','TimeUTC');
fprintf('✓  irradiance_hourly.mat saved (%d hourly points).\n', nH);
