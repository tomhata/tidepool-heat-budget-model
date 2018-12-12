function [poolTemp, totalQvals, poolData, airTemp, oceanTemp, serialDate, totalTemp] = thermalModel1_2_4_4
%% Notes
% |thermalModel1_2_4_4| updated Dec 10 2018 by Tom Hata.
% This function is designed to generate a heat and mass budget model of a
% tidepool in the intertidal. Long-term weather data must be imported.
% Currently used data set is 'tenyear_final_serialdate.txt', which contains
% continuous environmental data for approx 10 years at 10 min intervals.
% See IMPORT DATA SETTINGS section to determine file name and rows to 
% import. Length of output data = (endRow - startRow + 1) = n.
%  
% This function outputs the following:
% |poolTemp|   : 1 x n array of temperature data (C) of the tidepool
% |totalQvals| : 5 x n matrix of heat flux data for each point in time (W):
%       1. shortwave radiation
%       2. net longwave radiation
%       3. forced convection along air/water or air/rock (if dry) interface
%       4. free convection along water/rock interface
%       5. evaporation
% |poolData|: 5 x n matrix of tide pool properties as follows:
%       1. radius of surface of pool in contact with air (m)
%       2. pool depth (m)
%       3. pool volume (m^3)
%       4. salinity (ppt)
%       5. mass of water (kg) (NOTE: This does not include dissolved salt
%           mass. If salt water mass is desired, it can be calculated as
%           mass of water x (1 + salinity/1000).
% |airTemp|     : n x 1 array of air temperature (C)
% |oceanTemp|   : n x 1 array of ocean temperature (C)
% |serialDate|  : serialized time and date information. Use datevec to read
% |totalTemp|   : OPTIONAL. m x n array of temperature data (C) of entire
%       model through time. Used for debugging. Use with large data set
%       runs the risk of running out of memory.
%
%% Version History
% v1.2.4.4 Minor editing of variables.
% v1.2.4.3 Changed output variable |totalRVolSalMass| to
% |totalRHVolSalMass| to extract current depth data.
% v1.2.4.2 Remove some conflicting use of global and subfunction variables.
% adjusted tidepool update calculation for case with vertical walls
% v1.2.4.1 Corrected error, surface convection calculated during empty pool now.
% Added emissivity for dry granite.
% v1.2.3 Merged access to total temperature of rock column based on number
% of output arguments. fixed salinity being forced to imposed max value during
% physical property calculation
% v1.2.2 fixing issues in previous version with longwave calculation.
% adding emissivitySky = 1 for evening foggy skies. Using original
% calculation by Diana.
% v1.2.1 it appears unnecessary to look at data within subloop, because
% temperature data look stable. Data collected in 1.2.1 only looks at data
% output of main time loop to allow larger time spans.
% Pool area was not being recalculated in v1
% Check with Mark, longwave emissivity of granite should be similar to
% water (for when pool is dried).
% Changed calculation of beta constant of water to (mass water * (1 +
% salinity/1000)) to reflect that some salt precipitates out.
% Calculations become unstable when remaining water mass becomes very
% small. setting mass to 0 if mass is <1% of starting mass appears to fix
% instability.
% Updated temperature distribution of TPool to rock once pool is dry.
% Updated coefficient matrix to apply longwave, shortwave, and forced
% convection to rock surface when pool is dry.

%% Initialize constants and variables
% IMPORT DATA SETTINGS
fileName = 'tenyear_final_serialdate.txt'; % Weather data file name
%|startRow| Starting row number in data set
%|endRow| Ending row number in data set
startRow = 127297; % Row 127297 is 1/1/2002 0:00
%endRow = 526032; % Last row of data set is 526032, 7/31/2009 23:50
endRow = 227296;

% MODEL SETTINGS - Pool is modeled as a truncated cone.
% To check heat flux values, radius can be set to 0.564 (area = 1m^2)
% and wall angle = pi/2.
rPoolBottom = 0.25; % Radius of bottom surface of pool (meters)
wallPoolAngle = pi/4; % Angle of tidepool walls relative to horizontal (rad)
depthPoolMax = 0.1; % Maximum pool depth (meters), also initial value
areaPoolBottom = pi * rPoolBottom ^ 2; % surface area of pool bottom (m^2)
salinityIni = 35; % Initial salinity of pool (ppt)
humidRel = 0.75; % Relative humidity. Replace with real data if available.

% ON/OFF AND LIMIT SETTINGS
viewFactorCalcOn = true; % Calculate view factor for longwave radiation
% |salinityLimit| (ppt). For use in subfunction |calcSeawaterProp|.
% If >= 0, calculate properties of tidepool using this value. If < 0, use
% current salinity.
salinityLimit = 120;
maxSalinity = 360; % max salinity of tide pool (ppt)
% |massWaterCutoff| minimum mass of WATER in pool (distinct from water + 
% dissolved salt) relative to initial value (fraction). Below this value, 
% mass is set to 0 and pool is considered dry.
massWaterCutoff = 0.02;

% SPATIAL AND TEMPORAL DISTRIBUTION
numRockNodes = 100; % Number of nodes in bedrock. Suggested 100
nodeDistance = 0.01; % Distance between nodes in rock (m). Suggested 0.01
dt = 10; % Time between steps (s). Suggested 10 or less
timeInterval = 600; % Interval between data points (s). Code can be changed
% to automatically calculate between points in data set.
subLoopSteps = ceil(timeInterval / dt); % Number of steps for subloop to step through

% PHYSICAL PROPERTIES
absorptivityRock = 0.60; % Shortwave absorptivity of granite
conductivityRock = 3; % W/(m*K)
cpRock = 790; % J/(kg*K)
densityRock = 2601; % granite, kg/m^3
emissivityWater = 0.96;
emissivityGranite = 0.86;
SB = 5.67e-8; % Stefan-Boltzmann constant (W / (m^2 * K^-4))
% Constant for coefficient matrix used in finite element model
betaConstantRock = 1 / (densityRock * cpRock * nodeDistance * areaPoolBottom);



%% POOL REPLENISHMENT - CHANGE AS NEEDED

%replenishPool = false;
replenishPool = true; % True to turn on tidepool replenishment subfunction
% |swashFactor| multiplier for sig wave height based on exposure. 0.5 for
% protected area, 0.65 for exposed area.
swashFactor = 0.65;
% |replenishCriteria| generic variable to trigger pool replenishment.
% Trial value is amplifiedTide > 2.7m.
replenishCriteria = 2.7;

%% Load data and prepare model
% Load weather data from longterm data set. Note that this data set does
% not contain relative humidity.
weatherData = loadWeatherData();

% |parseWeatherData|Parse weather data into separate arrays and matricies.
[day, irradiancePyr, airTemp, tidalHeight, waveHeight, windSpeed, ...
    oceanTemp, serialDate, timeInMin] = parseWeatherData();

clear weatherData

% CHANGE THIS AS NEEDED. SIMPLE CRITERIA TO DETERMINE WHEN TO REFILL POOL.
amplifiedTide = tidalHeight + (swashFactor * waveHeight);
fillPool = (amplifiedTide > replenishCriteria) .* replenishPool;

%% Initialize values
% Initialize temperature distribution of model, |tempDistrib|. N-by-1 array
% Set all temperatures to ocean temperature.
% tempDistrib(1) = tidepool temperature
% tempDistrib(N-1) should be set to ocean temperature each time step.
% tempDistrib(N) = 1 (for matrix multiplication)
tempDistrib = [ oceanTemp(1) * ones(numRockNodes + 1,1); 1];

% Test data dump
poolTemp = zeros( 1, length(day));
totalQvals = zeros(5, length(day));
poolData = zeros( 5, length(day));

% Only create totalTemp if nargout == 7.
if nargout == 7
    totalTemp = zeros( length(tempDistrib), length(day));
end



%% Initialize  main values for main loop
TPool = tempDistrib(1); % Update tide pool temperature
TRock = tempDistrib(2); % Update rock surface temperature
salinityCurr = salinityIni; % Update salinity

% Calculate tidepool water properties
[densityWater, isobExpWater, latentHeatWater, conductWater, kinemViscWater, cpWater, thermDiffWater] = calcSeawaterProp();

% Initialize tidepool parameters.
[rPoolIni, volPoolIni, massSalt, massWaterIni] = initializePool();
depthPoolCurr = depthPoolMax;
rPoolCurr = rPoolIni;
areaPoolCurr = pi * rPoolCurr^2;
massWaterCurr = massWaterIni;
volPoolCurr = volPoolIni;

% Prepare coefficient matrix for finite element model
coeffMatrix = createCoeffMatrix();


%% Calculate shortwave irradiation
% |solarTime| Reconcile solar and clock time.
solarTime = calcSolarTime();

% |calcSolarIrradiance| Calculate solar altitude and total shortwave
% irradiance.
[irradianceSolar, solarZenith] = calcSolarIrradiance();

% |calcSurfaceTransmit| Calculate transmittance of light through water
% surface based on sun angle.
transmitTotal = calcSurfaceTransmit();

% Calculate shortwave energy input onto upper rock surface. This done all
% at once, as it is independent of tidepool and rock temperature.
QswAll = calcShortWave();


tic
%% Begin main loop
for i = 1: length(day)
    
    z = waitbar(i/length(day));
    % Set current temperatures and wind speed
    TAir = airTemp(i);
    windSpeedCurr = windSpeed(i);
    tempDistrib(end-1) = oceanTemp(i);
    
    % Current shortwave irradiance
    Qsw = QswAll(i);
    
    % Calculate evaporative energy and mass flux if pool is not empty and
    % |fillPool| is false.
    % Use longwave emissivity of water.
    if massWaterCurr > 0 && ~fillPool(i)
        [Qevap, massFlux] = calcEvap();
        Qlw = calcLongWave(emissivityWater);
        
    else
        % If pool is dry or if pool is being refilled, evaporative mass flux = 0
        % If pool is dry, longwave emissivity of granite. If the pool is
        % not dry but being refilled, longwave heat flux is set to 0
        % anyway.
        Qevap = 0;
        massFlux = 0;
        Qlw = calcLongWave(emissivityGranite);
    end
    
    % Refill pool with ocean temp seawater if |fillPool| is true
    if fillPool(i)
                
        tempDistrib(1) = oceanTemp(i);
        Qlw = 0;
        depthPoolCurr = depthPoolMax;
        rPoolCurr = rPoolIni;
        areaPoolCurr = pi * rPoolCurr^2;
        massWaterCurr = massWaterIni;
        volPoolCurr = volPoolIni;
        salinityCurr = salinityIni;
        
    end
    
    %% Subloop for finite element model using Runge-Kutta
    for j = 1:subLoopSteps
        
        % Tested updating with subloop. Pool cools too quickly. It's either an
        % issue with Q calculations or coefficients used in matrix. Confirmed this
        % was the issue. Pushed to main branch.
        TPool = tempDistrib(1);
        TRock = tempDistrib(2);
        
        if massWaterCurr > 0 && ~fillPool(i)
            % If |fillpool| is true, set surface convection to 0, but calculate
            % bottom free convection.
            QcvSurf = calcConvectionSurf();
            QcvBtm = calcConvectionBottom();
                      
        elseif fillPool(i)
            % If |fillpool| is true, set surface convection to 0, but calculate
            % bottom free convection.
            TPool = oceanTemp(i);
            QcvSurf = 0;
            QcvBtm = calcConvectionBottom();

        else         
            % Pool is empty and not refilled. Forced surface convection of rock
            % layer but no free convection.
            TPool = TRock;
            tempDistrib(1) = TRock;
            QcvSurf = calcConvectionSurf(); % |QcvSurf| must be applied to rock layer
            QcvBtm = 0;
        end
        
        % Update coefficient matrix and calculate using Runge-Kutta method
        updateCoeffMatrix();        
        tempDistrib = rk4();
        
    end
    
    %% Update values for next cycle
    TPool = tempDistrib(1);
    TRock = tempDistrib(2);
    
    % Update water mass and salinity
    massWaterCurr = massWaterCurr + massFlux * timeInterval;
    if massWaterCurr < massWaterCutoff * massWaterIni
        massWaterCurr = 0;
    end
    
    salinityCurr =  massSalt / massWaterCurr * 1000;
    if salinityCurr > maxSalinity
        salinityCurr = maxSalinity;
        
    end
    
    % Update water properties
    [densityWater, isobExpWater, latentHeatWater, conductWater, kinemViscWater, cpWater, thermDiffWater] = calcSeawaterProp();
    
    % Update tidepool properties
    
    [volPoolCurr, areaPoolCurr, depthPoolCurr, rPoolCurr] = updatePool();
    
    
    totalQvals(:,i) = [Qsw; Qlw; QcvSurf; QcvBtm; Qevap];
    poolData(:,i) = [rPoolCurr; depthPoolCurr; volPoolCurr; salinityCurr; massWaterCurr];
    poolTemp(i) = TPool;
    
    if exist('totalTemp','var')
        totalTemp(:,i) = tempDistrib;
    end
end

%% End Loop, prepare data for export
close(z)

poolTemp = poolTemp - 273; % tide pool temperature is 1st row of totalTemp
airTemp = airTemp - 273;
oceanTemp = oceanTemp - 273;

if exist('totalTemp','var')
    totalTemp(end,:) = [];
    totalTemp = totalTemp - 273; % Convert K to C
end


toc

%% Subfunctions
    function weatherData = loadWeatherData()
        
        % |loadWeatherData| imports longterm historical weather data.
        % Standard input file used is 'tenyear_final_serialdate.txt', see README
        % file for details. This file is a tab-delimited text file with 13 columns:
        % 1  - Year
        % 2  - Day of year
        % 3  - Hour:min without colon (1:10AM = 110)
        % 4  - Solar irradiance, shortwave (W/m^2)
        % 5  - Air Temp (C)
        % 6  - Predicted tide ht (m) from Monterey station (Station ID = 9413450)
        % 7  - Measured tide ht (m) from Monterey station (Station ID = 9413450)
        % 8  - Significant Wave Height Hs (cm)
        % 9  - Wind speed (m/s)
        % 10 - Wind direction, degrees east of north
        % 11 - Sea surface temperature (C)
        % 12 - Empty
        % 13 - Matlab serial date number, readable with datestr command
        %
        % startRow and endRow are desired start and end row of data.
        
        % Load data
        weatherData = load(fileName);
        
        % Remove rows after desired end row
        if endRow < length(weatherData)
            weatherData( endRow + 1 : end , : ) = [];
        end
        
        % Remove rows preceding desired start row
        if startRow > 1
            weatherData( 1 : startRow - 1  , : ) = [];
        end
        
    end

    function [day, irradiancePyr, airTemp, tidalHeight, waveHeight, windSpeed,oceanTemp, serialDate, timeInMin] = parseWeatherData()
        
        %|parseWeatherData| Parse weather data into data categories to provide
        %greater flexibility in subsequent functions. Convert serialized date to
        %n-by-6 matrix of time and date data [yr, month, day, hr, min, sec]. For
        %temperature data, add 273 to convert to Kelvin.
        
        day           = weatherData(:, 2);
        irradiancePyr = weatherData(:, 4);
        airTemp       = weatherData(:, 5) + 273;
        tidalHeight   = weatherData(:, 7);
        waveHeight    = weatherData(:, 8)/100; % cm to m
        windSpeed     = weatherData(:, 9);
        oceanTemp     = weatherData(:,11) + 273;
        serialDate    = weatherData(:,13);
        dateTime      = datevec(weatherData(:,13));
        timeInMin     = dateTime(:,4) * 60 + dateTime(:,5);
    end

    function [rPoolIni, volPoolIni, massSalt, massWaterIni] = initializePool()
        
        % |initializePool| calculates initial state of a tidepool modeled as an
        % inverted truncated circular cone, with smooth, flat bottom of radius
        % |rPoolBottom|, total depth |depthPoolIni|, and sloped walls of angle
        % |wallPoolAngle| relative to horizontal (radians). Radius at surface of
        % the pool |rPoolSurfIni| > |rPoolBottom|.
        % Pool contains volume |volPoolIni|.
        
        rPoolIni = rPoolBottom + depthPoolMax * cot(wallPoolAngle);
        
        volPoolIni = 1/3 * pi * (rPoolBottom^2 + rPoolBottom * ...
            rPoolIni + rPoolIni^2) * depthPoolMax;
        
        massPoolIni = densityWater * volPoolIni;
        massSalt = massPoolIni * salinityCurr / (1000 + salinityCurr);
        massWaterIni = massPoolIni - massSalt;
        
    end

    function [rho, betaP, hfg, k, nu, cp, alpha] = calcSeawaterProp()
        % Calculate and update seawater properties
        % Inputs: |Temp| (K), |Salinity| (ppt)
        % Outputs:
        %  |rho|   - density (kg / m^3)
        %  |betaP| - isobaric expansivity (K^-1)
        %  |hfg|   - latent heat of evaporation (J / kg)
        %  |k|     - thermal conductivity (W / (m * K))
        %  |nu|    - kinematic viscosity (m^2 / s)
        %  |cp|    - specific heat (J / (kg * K))
        %  |alpha| - thermal diffusivity (m^2 / s)
        
        % Calculations are taken from Matlab scripts available at
        % http://web.mit.edu/seawater/. For these calculations, salinity ranges are
        % 0 - 150 ppt. An evaporating pool can exceed these values, as the maximum
        % solubility of salt in water is approximate 360 ppt. To accomodate for
        % calculation errors, calculations using this function can potentially use
        % a limited salinity range.
        
        %% Adjust units
        T = TPool - 273.15; %Temp in C
        
        % Set salinity to |salinityLimit| if a limit is set and is below current
        % salinity value.
        
        if salinityLimit >= 0 && salinityLimit < salinityCurr
            S = salinityLimit; %Salt in ppt
        else
            S = salinityCurr;
        end
        
        
        s = salinityCurr/1000; %Salt in fraction
        
        %% |rho| Density
        
        a = [
            9.9992293295E+02
            2.0341179217E-02
            -6.1624591598E-03
            2.2614664708E-05
            -4.6570659168E-08
            ];
        
        b = [
            8.0200240891E+02
            -2.0005183488E+00
            1.6771024982E-02
            -3.0600536746E-05
            -1.6132224742E-05
            ];
        
        rho_w = a(1) + a(2)*T + a(3)*T^2 + a(4)*T^3 + a(5)*T^4;
        D_rho = b(1)*s + b(2)*s*T + b(3)*s*T^2 + b(4)*s*T^3 + b(5)*s^2*T^2;
        rho   = rho_w + D_rho;
        
        %% |betaP| Isobaric expansivity
        drho_wdT = a(2) + 2*a(3)*T + 3*a(4)*T^2 + 4*a(5)*T^3;
        dD_rhodT = b(2)*s + 2*b(3)*s*T + 3*b(4)*s*T^2 + 2*b(5)*s^2*T;
        drho_sw_sharqdT = drho_wdT + dD_rhodT;
        
        betaP = - drho_sw_sharqdT / rho;
        
        %% |hfg| latent heat of evaporation
        c = [
            2.5008991412E+06
            -2.3691806479E+03
            2.6776439436E-01
            -8.1027544602E-03
            -2.0799346624E-05
            ];
        
        hfg_w = c(1) + c(2)*T + c(3)*T^2 + c(4)*T^3 + c(5)*T^4;
        hfg   = hfg_w * (1 - 0.001*S);
        
        %% |k| thermal conductivity
        T_90 = 1.00024*T;      %convert from T_90 to T_68
        S_P = S / 1.00472;    %convert from S to S_P
        
        k = 10^( log10(240 + 0.0002 * S_P) + ...
            0.434 * ( 2.3 - ( 343.5 + 0.037 * S_P)/  ...
            (T_90 + 273.15)) * (1 - ( T_90 + 273.15) ...
            /(647.3 + 0.03 * S_P))^(1/3) - 3);
        
        %% |nu| kinematic viscosity
        d = [
            1.5700386464E-01
            6.4992620050E+01
            -9.1296496657E+01
            4.2844324477E-05
            1.5409136040E+00
            1.9981117208E-02
            -9.5203865864E-05
            7.9739318223E+00
            -7.5614568881E-02
            4.7237011074E-04
            ];
        
        mu_w = d(4) + 1./(d(1)*(T+d(2)).^2+d(3));
        
        viscA  = d(5) + d(6) * T + d(7) * T.^2;
        viscB  = d(8) + d(9) * T + d(10)* T.^2;
        
        % |mu| dynamic viscosity
        mu = mu_w * (1 + viscA*s + viscB*s^2);
        nu  = mu / rho;
        
        %% |cp| specific heat
        
        T68 = 1.00024*(T+273.15);      %convert from T_90 to T_68
        
        A = 5.328 - 9.76 * 10 ^ (-2) * S + 4.04*10^(-4)*(S).^ 2;
        B = -6.913 * 10 ^ (-3) + 7.351 * 10 ^ (-4) * (S) - 3.15*10^(-6)*(S).^2;
        C = 9.6 * 10 ^ (-6) - 1.927 * 10 ^ (-6) * (S) + 8.23 * 10^(-9) *(S).^2;
        D = 2.5 * 10 ^ (-9) + 1.666 * 10 ^ (-9) * (S) - 7.125 * 10^(-12)*(S).^2;
        
        cp = 1000.*(A + B.*T68 + C.*(T68.^2) + D.*(T68.^3));
        
        %% |alpha| thermal diffusivity
        
        alpha = k./(rho.*cp);
    end

    function solarTime = calcSolarTime()
        
        % |calcSolarTime| Reconcile solar and clock time. solarTime is an n-by-1
        % array in minutes of the current day. Inputs are n-by-1 array of day of
        % year |day| and n-by-6 matrix of date and time (year, month, day, hour,
        % minute, second).
        
        % Longitude time correction |tLong| should be passed through as a global
        % constant in the main function.
        
        longitude = 121.88; % Longitude of HMS (degrees west)
        stdLongitude = 120; % Time zone longitude
        tLong = (stdLongitude - longitude) * 4; % longitude time correction (min)
        
        EOT = 0.1237 * cos(((2*pi) / 365.242) .* (day + 88.289)) + ...
            0.1654 * cos(((2*pi) / 182.621) .* (day - 127.029));
        
        solarTime = timeInMin + EOT + tLong;
        
    end

    function [irradianceSolar, solarZenith] = calcSolarIrradiance()
        
        % |calcSolarIrradiance| calculates the total solar irradiance
        % |irradianceSolar| (n-by-1 array) from zenith angle of the sun
        % |solarZenith|, latitude in radians |latitude|, hour angle |omega| and
        % declination angle |delta|. Hour angle |omega| is calculated using
        % solar time calculated by function |calcSolarTime|.
        
        % |latitude| should  be a passed through as a global constant in main
        % function.
        
        %% Calculate solar altitude from location, date, and time data.
        latitude = 36.62 * pi/180 ; % Latitude of HMS (radians)
        
        % |omega| Calculate hour angle (degrees east from north)
        omega = 2 * pi * ((solarTime - (12 * 60)) / (24 * 60));
        
        % |delta| Calculate the declination angle (degrees up from horizontal)
        delta = asin( 0.4093 * cos( 0.0172 .* (day - 173)));
        
        % |solarZenith| Calculate zenith angle (degrees from sun's culmination)
        % From Gates 1921
        solarZenith = acos( sin(latitude) .* sin(delta) + ...
            cos(latitude) .* cos(delta) .* cos(omega));
        
        % |solarAltitude. Back calculate the direct solar irradiance from the
        % cosinusoidal pyranometer data, taking into account that pyranometer
        % readings < 50 are questionable.
        
        % Solar altitude in radians.
        solarAltitude = (pi/2) - solarZenith;
        
        %% Adjust pyranometer irradiance reading to calculate total solar irradiance
        
        % Temporary vector to hold solar irradiance values
        irradianceTemp = irradiancePyr ./ sin(solarAltitude);
        
        irradianceSolar = zeros(size(irradianceTemp));
        
        % Keep only readings where solar altitude was > 0.087 and solar wattage was
        % > 50. All other values in solar remain zero.
        
        irradianceSolar( solarAltitude > 0.087 & irradiancePyr > 50) = ...
            irradianceTemp(solarAltitude > 0.087 & irradiancePyr > 50);
        
        % At extremely low angles, the calculations could give spuriously large
        % irradiances. Limit them to some reasonable value.
        
        irradianceSolar( irradianceSolar > 1200 ) = 1200;
        irradianceSolar( irradianceSolar < 0) = 0;
        
        % Assume daylight is between 5am and 7pm.
        irradianceSolar( timeInMin <= 300 | timeInMin > 1140) = 0;
        solarZenith( timeInMin <= 300 | timeInMin > 1140) = pi/2;
        %solarAltitude( timeInMin <= 300 | timeInMin > 1140) = 0;
    end

    function transmitTotal = calcSurfaceTransmit()
        
        % |calcSurfaceTransmit| calculates the fraction of light transmitted
        % through water's surface to the rock substrate at the bottom of the pool.
        % Angle of incidence is defined by n-by-1 array |solarZenith|. Reflectance
        % of perpendicular and parallel components of polarized light are
        % calculated, and total reflectance is calculated as the average the two
        % values. Transmittance of light |transmitTotal| is an n-by-1 array
        % calculated as 1 - total reflectance, and is a fraction between 0 and 1.
        
        refracIndex = 1.33;
        
        % Calculate reflectance of perpendicular component of light on water
        % surface.
        reflecPerpendicular = ( (sin(solarZenith) - tan(solarZenith)/refracIndex) ./ ...
            (sin(solarZenith) + tan(solarZenith)/refracIndex)).^2;
        
        % Calculate reflectance of parallel component of light on water surface.
        reflecParallel = ( (tan(solarZenith) - sin(solarZenith)/refracIndex) ./ ...
            (tan(solarZenith) + sin(solarZenith)/refracIndex)).^2;
        
        % Total reflectance is avererage of perpendicular and parallel reflectance.
        reflecTotal = 0.5 * (reflecPerpendicular + reflecParallel);
        
        % Transmittance is 1 - total reflectance.
        transmitTotal = 1 - reflecTotal;
        
    end

    function QswAll = calcShortWave()
        
        % |calcShortWave| calculates energy input into system through shortwave
        % irradiation for all points in time|QswAll|. Inputs are n-by-1 arrays
        % solar irradiance |irradianceSolar|, solar zenith angle |solarZenith|,
        % and transmittance of solar energy |transmitTotal|. Output is n-by-1
        % array QswAll, in units W. Qsw should be applied to surface layer of
        % bedrock, not water layer.
        
        % Calculate projected area of bottom of the pool (where energy is absorbed
        % by rock substrate) relative to sun.
        % May need to add additional correction depending on wall geometry.
        areaProj =  areaPoolBottom * cos(solarZenith);
        
        % |QswAll| calculate shortwave energy input for all time periods onto
        % rock substrate at the bottom of the pool, based on projected area
        QswAll = areaProj .* irradianceSolar .* transmitTotal .* absorptivityRock;
        
    end

    function Qlw = calcLongWave(emissivitySurf)
        
        % |calcLongWave| calculates energy input into system through longwave
        % radiation for current point in time |Qsw|, in units W. Qsw should be
        % applied to water layer of model. Inputs are:
        % current air temperature |TAir|(Kelvin)
        % current tidepool temperature |TPool| (Kelvin)
        % true/false signal to calculate view factor |viewFactorCalcOn| (0 or 1)
        % initial tidepool depth |depthPoolIni| (m)
        % current tidepool depth |depthPoolCurr| (m)
        % initial tidepool radius |radiusPoolIni| (m)
        % current tidepool radius |radiusPoolCurr| (m)
        % Also required are Stefan-Boltzmann Constant (5.67e-8 W / (m^2 * K^-4))
        % and emissivity of water 0.96.
        % In main program, all variables and constants can be passed through as
        % global variables.
        % CHECK WITH MARK. LONGWAVE EMISSIVITY OF GRANITE SHOULD BE CLOSE
        % TO WATER.
        
        %% Calculated variables
        
        % |emissivityAir| Air emissivity also depends on relative humidity, but
        % this effect is not significant. The full equation is:
        % emissivityAir = 7.37e-7 * Tair^2.46 * humidRel^0.143.
        % Decreasing relative humidity from 1 to 0.75 only changes air emissivity
        % by 5%. Therefore it is currently omitted from calculation. Relative
        % humidity data was not available during the creation of v1 of this code.
        
        emissivityAir = 7.37e-7 * TAir^2.46;
        
        if timeInMin(i) <= 300 || timeInMin(i) > 11400
            emissivityAir = 1;
        end
        
        %% Calculate view factor
        % Calculate view factor |viewFactor| of pool if |viewFactorCalcOn| is set
        % to 1 or true. Else, assume full view of sky (viewFactor = 1).
        
        if viewFactorCalcOn
            
            viewFactor = calcViewFactor;
        else
            
            viewFactor = 1;
        end
        
        %% Calculate longwave radiation
        %Qlw = emissivityWater * SB * areaPoolCurr * viewFactor * ...
        %    (emissivityAir * TAir^4 - TPool^4);
        %Qlw = emissivityWater * SB * areaPoolCurr * viewFactor * ...
        %    (TAir^4 * (emissivityAir - 1) + 4 * TAir^3 * (TAir - TPool));
        
        % Alternative calculation based on linear approximation. Used in Diana's
        % code. Present for confirmation that above equation works.
        Qlw = emissivitySurf * SB * areaPoolCurr * viewFactor * TAir^4 * (emissivityAir -1) ...
            + 4 * emissivitySurf * SB * areaPoolCurr * viewFactor * TAir^3 * (TAir - TPool);
        
        %% Subfunction |calcViewFactor|
        
        function viewFactor = calcViewFactor()
            
            % |seeSky| calculates the view factor |viewFactor| from the surface of the
            % pool, which is the fraction of sky that is "viewable" the by the surface.
            % View factor is used to calculate long wave emissivity.
            %
            % The calculation is taken from the center of the surface of the pool, and
            % requires inputs initial pool depth in meters |depthPoolIni|,
            % currentPoolDepth |depthPoolCurr|, and radius at the top of the pool
            % |radiusPoolIni|.
            %
            % This calculation may not be necessary if view factor is assumed to be 1.
            % This function is meant to be a subfunction within function |calcLongWave|
            
            depthFromTop = depthPoolMax - depthPoolCurr;
            skyAngle = atan( rPoolIni / depthFromTop );
            viewFactor = skyAngle / (pi/2);
            
        end
    end

    function coeffMatrix = createCoeffMatrix()
        % |createCoeffMatrix| creates matrix of coefficients to solve
        % finite element model. Adds only conduction constants. Matrix must
        % be updated with coefficients from other heat flux sources during
        % each step through loop.
        
        conductionConstantRock = conductivityRock / ...
            (densityRock * cpRock * nodeDistance^2);
        
        coeffMatrix = zeros( length( tempDistrib + 1));
        
        % Conduction between surface rock node and adjacent node.
        coeffMatrix(2,2) = - conductionConstantRock;
        coeffMatrix(2,3) = conductionConstantRock;
        
        % Populate remaining matrix with conduction constants for bedrock.
        for u = 3 : length(tempDistrib) - 2
            
            coeffMatrix(u, u-1) = conductionConstantRock;
            coeffMatrix(u, u  ) = -2 * conductionConstantRock;
            coeffMatrix(u, u+1) = conductionConstantRock;
        end
    end

    function [Qevap, massFlux] = calcEvap()
        %% Variables for evaporation calculation.
        % |TFilm| Temperature of air film directly above water
        TFilm = 0.5 * (TAir + TPool);
        
        % |kinemViscAir| Temperature-adjusted kinematic viscosity of air film
        kinemViscAir = -1.25e-5 + 9.28e-8 * TFilm;
        
        % |diffaw| Binary diffusion coefficient between air and water.Regression
        % of data from Bolz and Tuve 1976. Data found second-hand through web
        % search. Need to find a source. (m^2/s)
        diffaw = -2.775e-6 + 4.479e-8 * TFilm + 1.656e-10 * TFilm^2;
        
        % |reynolds| number of air blowing over tidepool. Characteristic length is
        % diameter of pool.
        reynolds = 2 * rPoolCurr * windSpeedCurr / kinemViscAir;
        
        % |schmidt| number
        schmidt = kinemViscAir / diffaw;
        
        % |Sherwood| number, convection over flat plate.
        sherwood = 0.664 * reynolds^(1/2) * schmidt^(1/3);
        
        % |hDiff| mass transfer coefficient (m/s)
        hDiff = sherwood * diffaw / (2 * rPoolCurr);
        
        % |vaporDensityAir| water vapor density in air in g / m^3. Saturation vapor
        %density at temperature TAir multiplied by relative humidity |humidRel|.
        vaporDensityAir = humidRel * calcSaturationD(TAir);
        
        % |vaporDensitySurf| vapor density at tide pool surface. Saturation vapor
        % density depends on surface's relative water content (see Luke Hunt thesis)
        vaporDensitySurf = calcHumidRelSurf() * calcSaturationD(TPool);
        
        %% Mass and energy flux calculation
        % |massFlux| is the rate of evaporation of water (g/s).
        % |massFlux| = (mass transfer coeff) * (diff in vapor density) * surf area
        % If surface vapor density exceeds air vapor density, mass and energy
        % flux will be negative.
        
        % !!! Written without conditional statements. Results may be thrown off if
        % air vapor density exceeds water vapor density as pool dries.
        
        % |massFlux| evaporative mass flux (g/s)
        massFlux = hDiff * (vaporDensityAir - vaporDensitySurf) * areaPoolCurr / 1000;
        
        % |Qevap| Latent heat of evaporation. Latent heat is in J/kg, so is
        % divided by 1e3 to convert to J/g
        Qevap = latentHeatWater * massFlux;
        
        
        %% |calcEvap| Subfunctions
        function satDensity = calcSaturationD(T)
            % |calcSaturationD| 4th order polynomial to calculate saturation vapor
            % density of water in air for temperature |T| in Kelvin. Output is
            % saturation vapor pressure |satDensity| in units g / m^3.
            
            coeffs = [3.32e-6 9.17e-5 0.0121 0.315 4.92];
            satDensity = polyval( coeffs, (T-273));
        end
        
        function humidRelSurf = calcHumidRelSurf()
            % |CalcHumidRelSurf| calculates relative humidity above a surface
            % |humidRelSurf| as a function of the surface's water content, as
            % described in Luke Hunt's thesis. |relWaterContent| Relative water
            % content is calculated from |salinity| (ppt).
            
            relWaterContent =  1 - salinityCurr / 1000;
            humidRelSurf = 0.21 * log(relWaterContent) + 1.04;
            
            % Maintain relative humidity within physical boundaries.
            if humidRelSurf > 1
                humidRelSurf = 1;
            elseif humidRelSurf < 0
                humidRelSurf = 0;
            end
        end
    end

    function  QcvSurf = calcConvectionSurf()
        % |calcConvectionSurf| calculates sensible heat convection across the
        % surface of the tidepool due to forced convection of air flow across
        % surface. |QcvSurf| is the rate of energy input or removed from the
        % tidepool (W).
        
        % |TFilm| Temperature of air film directly above water
        TFilm = 0.5 * (TAir + TPool);
        
        % |kinemViscAir| kinematic viscosity of air film (m^2 / s)
        kinemViscAir = -1.25e-5 + 9.28e-8 * TFilm;
        
        %|conductAir| Thermal conductivity of air [W / (m*k)]
        conductAir = 0.005013 + 7.2e-5 * TFilm;
        
        % |prandtl| Prandtl number of air at approx 300K. Prandtl number depends
        % on temperature but varies little in our operating temperature range.
        prandtl = 0.71;
        
        % |reynolds| number of air blowing over tidepool. Characteristic length is
        % diameter of pool.
        reynolds = 2 * rPoolCurr * windSpeedCurr / kinemViscAir;
        
        % |Nusselt| number of laminar flow over flat plate.
        nusselt = 0.664 * reynolds^(1/2) * prandtl^(1/3);
        
        % |hConv| convective heat transfer coefficient
        hConv = nusselt * conductAir / (2 * rPoolCurr);
        
        % |QcvSurf| calculate energy transfer due to forced convection
        QcvSurf = hConv * areaPoolCurr * (TAir - TPool);
    end

    function QcvBtm = calcConvectionBottom()
        % |calcConvectionBottom| calculate heat transfer between tidepool and
        % rock surface at the bottom of the pool via free (natural) convection.
        % |QcvBtm| is energy transfer rate (W) to the tidepool from the rock
        % surface. Conversly, -|Qcvbtm| is the energy transferred to the rock
        % surface from the tidepool.
        
        gravity = 9.8; % m/s^2
        
        % |Lc| characteristic length of system, which can be subjective. For a
        % horizontal flat plate, Lc = area / perimeter = r/2 (for circular disk).
        % Another value can be tidepool depth.
        Lc = rPoolBottom / 2;
        
        % |rayleigh| Rayleigh Number. Positive if TPool > TRock.
        rayleigh = gravity * isobExpWater * abs(TPool - TRock) * Lc^3 / ...
            (thermDiffWater * kinemViscWater);
        
        % |nusselt| number calculation depends on whether
        % TPool is > or < TRock. Incropera and DeWitt p.551.
        if TPool == TRock
            
            nusselt = 0;
            
        elseif  TPool > TRock
            
            % Top surface of cooled plate.
            nusselt = 0.27 * rayleigh^(1/4);
            
            % Top surface of warmed plate. Nusselt number calculation depends on
            % Rayleigh number.
        elseif rayleigh < 1e7
            
            nusselt = 0.54 * rayleigh^(1/4);
        else
            
            nusselt = 0.15 * rayleigh^(1/3);
        end
        
        % |hConvFree| heat transfer coefficient calculate using Nusselt number.
        hConvFree = conductWater * nusselt / Lc;
        
        % |QcvBtm| heat transfer from tidepool to rock. Conversely, heat transfer
        % from rock to pool is -QcvBtm
        QcvBtm = - hConvFree * areaPoolBottom * (TPool - TRock);
    end

    function updateCoeffMatrix()
        
        % For use in matrix to solve rate of temperature change. Divide each heat
        % flux source (excluding conduction) by this value to calculate total
        % temperature change for Tidepool and layer of rock beneath. This needs to
        % be recalculated each iteration, as water volume will change.
        
        
        if massWaterCurr > 0
            
            betaConstantWater = 1/ (cpWater * (massWaterCurr * (1 + salinityCurr/1000)));
            
            % Add heat fluxes to tidepool node. Longwave radiation, evaporative
            % cooling, forced convection between air and tidepool, free convection
            % between tidepool and rock surface.
            coeffMatrix(1,end) = betaConstantWater * (Qlw + Qevap + QcvSurf + QcvBtm);
            
            % Add heat fluxes to surface rock node. Shortwave radiation and free
            % convection between tidepool and rock surface (-QcvBtm).
            coeffMatrix(2,end) = betaConstantRock * (Qsw - QcvBtm);
            
        else
            
            % betaConstantWater = 0;
            % If current water mass is 0, pool is dry. longwave radiation
            % and forced convection should therefore be removed from rock
            % surface.
            %coeffMatrix(1,end) = betaConstantWater * (Qevap + QcvSurf + QcvBtm);
            coeffMatrix(1,end) =0;
            
            % Add heat fluxes to surface rock node. Shortwave radiation,
            % longwave radiation, and forced convection.
            coeffMatrix(2,end) = betaConstantRock * (Qlw + Qsw + QcvSurf);
            
        end
    end

    function temperatureData = rk4()
        % Resolve using 4th order Runge-Kutta
        % y' = F(x,y)
        % y0 = f(x0) -> yn = f(xn), initial values known
        % time step h
        % k1 = F( xn       , yn        )
        % k2 = F( xn + h/2 , yn + k1/2 )
        % k3 = F( xn + h/2 , yn + k2/2 )
        % k4 = F( xn + h   , yn + k3   )
        % y(n+1) = 1/6 * h * (k1 + 2*k2 + 2*k3 + k4)
        
        % Calculate k1 = slope at initial temp
        k1 = coeffMatrix * tempDistrib;
        
        % Calculate k2 at half time step using k1
        k2 = coeffMatrix * (tempDistrib + 0.5 * dt * k1);
        
        % Calculate k3 at half time step using k2
        k3 = coeffMatrix * (tempDistrib + 0.5 * dt * k2);
        
        % Calculate k4 at full time step using k3
        k4 = coeffMatrix * (tempDistrib + dt * k3);
        
        % Resolve temperature distibution using k values
        temperatureData = tempDistrib + 1/6 * (k1 + 2*k2 + 2*k3 + k4) * dt;
    end

    function [vol, apool, h, r2] = updatePool()
        % Calculate current depth (h) and surface radius from base radius, wall
        % angle, and volume.
        %h = 1/(tan(pi/2 - wallPoolAngle)) * (((3 * volPoolCurr * ...
        %    tan(pi/2 - wallPoolAngle) / pi + rPoolBottom^3))^(1/3) - rPoolBottom);
        %r2 = rPoolBottom + h / tan(wallPoolAngle);
        
        if massWaterCurr > 0
            
            vol = (massWaterCurr * (1 + salinityCurr/1000)) / densityWater;
            
            if wallPoolAngle < pi/2
                % Calculate current depth and surface radius from base radius, wall
                % angle, and volume.
                h = 1/(tan(pi/2 - wallPoolAngle)) *  ...
                    (((3 * vol * tan(pi/2 - wallPoolAngle) / ...
                    pi + rPoolBottom^3))^(1/3) - rPoolBottom);
                r2 = rPoolBottom + h / tan(wallPoolAngle);
            
            elseif wallPoolAngle == pi/2
                
                h = vol / (pi * rPoolBottom^2);
                r2 = rPoolBottom;                
            end
                        
        else
            
            % If the pool is dry, set assoc. parameters to empty.
            vol = 0;
            h = 0;
            r2 = rPoolBottom;
        end
        
        apool = pi * r2^2;
        
    end

end
