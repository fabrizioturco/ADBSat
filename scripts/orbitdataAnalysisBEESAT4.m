function tmData = orbitdataAnalysisBEESAT4(tmData, activityLevel)
% "low", "moderate", "high", or "actual"

    %% SETTINGS
    modName = 'beesat';
    % modName = 'full_simplified';
    % modName = 'flp_compact';
    % Path to model file
    ADBSat_path = ADBSat_dynpath;
    modIn = fullfile(ADBSat_path,'inou','obj_files',[modName,'.obj']);
    modOut = fullfile(ADBSat_path,'inou','models');
    resOut = fullfile(ADBSat_path,'inou','results');
    
    % tmData = readtable(inFilePath);
    % fprintf('TM data loaded.\n')
    
    fprintf('Selected level of solar and geomagnetic activity: ' + activityLevel + '\n')
    
    % Satellite data
    mass = 1;
    
    % Model parameters
    AnO =  1;
    shadow = 1;
    inparam.gsi_model = ['sentman'];
    inparam.Tw = 300; % Wall Temperature [K]
    solar = 0;
    inparam.sol_cR = 0.15; % Specular Reflectivity
    inparam.sol_cD = 0.25; % Diffuse Reflectivity
    
    verb = 0;
    del = 0;
    
    % Import model
    [modOut] = ADBSatImport(modIn, modOut, verb);
    
    lineLength = 0;
    fprintf("Currently analysing: ")
    tic
    for i=1:height(tmData)
        %% INPUT
        %Input conditions
        alt = tmData.alt(i);
        lat = tmData.lat(i);
        lon = tmData.lon(i);
        timestamp = char(tmData.timestamp(i));
        if i~=0 % do not delete line if first step
            fprintf(repmat('\b',1,lineLength)); % delete previous print output
        end
        lineLength  = fprintf(timestamp); % print current timestep
        y           = year(tmData.timestamp(i));
        dayofyear   = day(tmData.timestamp(i),'dayofyear');
        UTseconds   = second(tmData.timestamp(i),'secondofday');
    
        % solar and geomagnetic activity indexes
        if activityLevel == "low"
            f107Average = 65;
            f107Daily = 65;
            magneticIndex = ones(1,7) * 0;
        elseif activityLevel == "moderate"
            f107Average = 140;
            f107Daily = 140;
            magneticIndex = ones(1,7) * 15;
        elseif activityLevel == "high"
            f107Average = 250;
            f107Daily = 250;
            magneticIndex = ones(1,7) * 45;
        else
            f107Average = tmData.f107A(i);
            f107Daily = tmData.f107(i);
            magneticIndex = str2double(split(erase(cell2mat(tmData.aph{i}),["[","]"])))';
        end
     
        env = [alt, lat, lon, y, dayofyear, UTseconds, f107Average, f107Daily, magneticIndex, AnO]; % Environment variables
        
        if AnO
            Oflag = 'Oxygen';
        else
            Oflag = 'NoOxygen';
        end
        
        % SESAM: Semiempirical Model for Satellite Energy-Accommodation Coefficients
        [T, rho] = atmosnrlmsise00(alt, lat, lon, y, dayofyear, UTseconds, f107Average, f107Daily, magneticIndex, Oflag);
    
        inparam.alpha = 7.5E-17*rho(2)*T(2) / (1+7.5E-17*rho(2)*T(2)); % SESAM for accommodation coefficient (altitude dependent)
        if inparam.alpha <= 0.85
                inparam.alpha = 0.85;
        elseif inparam.alpha >= 1
            inparam.alpha = 1;
        end
        
        %% ANALYSIS
       
        % Calculate
        [ADBout] = ADBSatFcn(modOut, resOut, inparam, 0, 0, shadow, solar, env, del, verb);
        result_min = load(ADBout);
        [ADBout2] = ADBSatFcn(modOut, resOut, inparam, 45, 35.2644, shadow, solar, env, del, verb);
        result_max = load(ADBout2);
        beta_min = -result_min.Cf_w(1)*result_min.Aref/mass;
        beta_max = -result_max.Cf_w(1)*result_max.Aref/mass;
        if tmData.mode{i} == "EPM"
            beta = beta_min;
        else
            beta = (beta_min+beta_max)/2;
        end
        
        tmData.beta_min(i) = beta_min;
        tmData.beta_max(i) = beta_max;        
        tmData.alpha(i) = inparam.alpha;
        tmData.rho(i) = rho(6);
        tmData.T(i) = T(1);
        tmData.beta(i) = beta;
    
    end 
    time = toc;
    fprintf('\n');
    fprintf('Mean Results:\n')
    fprintf('beta \t= %.4g\n', mean(tmData.beta));
    fprintf('rho \t= %.4g kg/m^3\n', mean(tmData.rho));
    hour = floor(time/3600);
    min = floor((time - hour*3600)/60);
    s = time - hour*3600 - min*60;
    fprintf('\nComputation time: %.4gh %.4gmin %.4gs\n', hour, min, s);
    fprintf('\n');
    fprintf('\n');
    
    % tmdata_path_res = char(inFilePath);
    % outFilePath = "adbsat_analysis_res.csv";
    % writetable(tmData,outFilePath);
end