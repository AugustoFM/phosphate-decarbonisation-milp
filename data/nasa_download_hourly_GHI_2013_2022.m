% nasa_download_hourly_GHI_2013_2022.m
% -------------------------------------------------------------------------
% Downloads NASA-POWER ALLSKY_SFC_SW_DWN (hourly GHI) for 2013-2022, Al-Jalamid
% and stores irradiance_hourly_2013-2022.mat:
%     • GHI_kWhm2  –  n×1  double     (kWh m-2 per hour, UTC)
%     • TimeUTC    –  n×1  datetime   (UTC, hourly)
% -------------------------------------------------------------------------
clear;  clc;

% ------------ site coordinates ------------------------------------------
lat  = 31.200;         % ° N  (Al-Jalamid)
lon  = 39.850;         % ° E
yrFirst = 2013;        % first calendar year
yrLast  = 2022;        % last  calendar year  (inclusive)
timeStd = 'LST';       % NASA time-standard: 'LST' or 'UTC'

% ------------ build master timeline -------------------------------------
TimeUTC = (datetime(yrFirst,1,1,0,0,0,'TimeZone','UTC') : ...
           hours(1) : ...
           datetime(yrLast,12,31,23,0,0,'TimeZone','UTC')).';
nH = numel(TimeUTC);

TimeUTC( month(TimeUTC)==2 & day(TimeUTC)==29 ) = [];
nH = numel(TimeUTC);                    % recompute length

fprintf('Downloading hourly GHI for %d–%d …\n', yrFirst, yrLast);
opts = weboptions('Timeout',60);
kv   = containers.Map('KeyType','char','ValueType','double');

for y = yrFirst:yrLast
    url = sprintf(['https://power.larc.nasa.gov/api/temporal/hourly/point?', ...
        'parameters=ALLSKY_SFC_SW_DWN&community=RE&latitude=%.4f&longitude=%.4f&', ...
        'start=%4d0101&end=%4d1231&time-standard=%s&format=JSON'], ...
        lat, lon, y, y, upper(timeStd));

    raw = webread(url, opts);

    % Robust decode – webread returns struct in newer MATLAB, char in older
    if ischar(raw) || isstring(raw)
        J = jsondecode(raw);
    else
        J = raw;
    end

    P = J.properties.parameter.ALLSKY_SFC_SW_DWN;
    fns = fieldnames(P);

    for k = 1:numel(fns)
        key  = regexprep(fns{k}, '^\D+', '');    % strip leading non-digits
        if numel(key) == 10                      % expect yyyymmddHH
            kv(key) = double(P.(fns{k})) / 1000; % W m-2 → kWh m-2 per h
        end
    end
    fprintf('  ✓ %d complete (%d values)\n', y, numel(fns));
end

% ------------ align & interpolate gaps ----------------------------------
GHI_kWhm2 = nan(nH,1);
for i = 1:nH
    tag = datestr(TimeUTC(i),'yyyymmddHH');
    if isKey(kv, tag);  GHI_kWhm2(i) = kv(tag);  end
end

missing = sum(isnan(GHI_kWhm2));
if missing
    warning('%d / %d hours missing – filling by linear interp.', missing, nH);
    GHI_kWhm2 = fillmissing(GHI_kWhm2,'linear','MaxGap',3,'EndValues','nearest');
end

% ------------ save ------------------------------------------------------
save('irradiance_hourly_2013_2022.mat','GHI_kWhm2','TimeUTC');
fprintf('✓ Saved irradiance_hourly_2013-2022.mat  (%d hours)\n', nH);

% ------------ optional sanity plot (first 168 h) ------------------------
try
    figure;
    plot(TimeUTC(1:168), GHI_kWhm2(1:168));
    grid on;  ylabel('kWh m^{-2} h^{-1}');
    title(sprintf('GHI – first week of %d', yrFirst));
catch
end
