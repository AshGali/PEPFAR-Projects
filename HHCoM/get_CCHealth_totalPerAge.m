% DESCRIPTION: Gets ccHealthState and totalPerAge calculations from
% vaxCEA_multSims_processResults_PEPFAR.m, for use with
% MATLAB Engine in Python.
%
% INPUTS:
%   fileInd      - index for vaxSimResult filename
%   sceString    - scenario string (e.g. '1')
%   baseFileName - scenario results folder name. If omitted,
%                  built like processResults (line 79).
%
% OUTPUTS:
%   ccHealthState   - [nTimepoints x (age+1) x endpoints]
%   totalPerAge     - [nTimepoints x gender x (age+1)]
%   monthlyTimespan - [1 x nTimepoints] years


function [ccHealthStateReshape, totalPerAgeReshape, ccHealthStateCols, totalPerAgeCols, monthlyTimespan] = get_CCHealth_totalPerAge(fileInd, sceString, baseFileName)

if nargin < 3 || isempty(baseFileName)
    baseFileName = ['VaccineKenyaPrEPCea_Mar14_PEPFARstop_S', sceString]; % ***SET ME***
end
fileIndStr = num2str(fileInd);
n = fileIndStr;

%% Load only the params we actually need out of loadUp2_S1
%  into a cell array and pull out just the ones this function uses.
% NOTE: confirm these index positions still match loadUp2_S1's
% output order (see vaxCEA_multSims_processResults_PEPFAR.m
% unpacking list) 
%   CHECK: does nargout('loadUp2_S1') == numel(c) ?
nOut = abs(nargout('loadUp2_S1'));
c = cell(1, nOut);
[c{:}] = loadUp2_S1(1, 0, [], [], [], 1);

% TEMPORARY
j = 1;
k = 1;

timeStep = c{2};  startYear = c{3};
disease  = c{7};  viral = c{8};  hpvVaxStates = c{9};  hpvNonVaxStates = c{10};
endpoints = c{11}; intervens = c{12}; gender = c{13}; age = c{14}; risk = c{15};
toInd = c{19};

lastYear = 2026; % manually set in futureSim, see processResults.m line 59

%% Timespan
monthlyTimespan = startYear : timeStep : lastYear;
monthlyTimespan = monthlyTimespan(1 : end-1);
nTimepoints = length(monthlyTimespan);

%% Load results, concatenate historical + future (see lines 173-209)
% ***SET ME***: update file patterns
curr = load([pwd , '/HHCoM_Results/toNow_sk1822_stochMod_baseline_2dose_nowanePEPFARHistoricalBaseTest' , fileIndStr]);
result = load([pwd , '/HHCoM_Results/' , baseFileName , '/vaxSimResult' , fileIndStr , '.mat']);
popVec = [curr.popVec(1:end, :); result.popVec(2:end, :)]; % only popVec needed here

% load results from vaccine run into cell array
vaxResult{n} = result;
disp(size(vaxResult{n}));
disp(size(curr));

% concatenate vectors/matrices of population up to current year to population
% matrices for years past current year
% curr is historical model results 
% vaxResult is future model results
% this section of code combines the historical results with future
% results
% notice for vaxResult you start at row 2. likely because of
% 2023 being double counted in both. 
vaxResult{n}.popVec = [curr.popVec(1 : end  , :); vaxResult{n}.popVec(2 : end , :)]; % consolidating historical population numbers with future

% Indices of calib runs to plot
% Temporarily commenting out to only run one scenario first to test out
% code
%***SET ME***: this will likely change once we spit out all of the sim results
fileInds = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', ...
                '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', ...
                '22', '23', '24', '25'};    % 22Apr20Ph2V11 ***************SET ME****************
%fileInds = {'18'}; % FORTESTING
nRuns = 1;
endYear = 2025;

%lastYear = 2026; % manually set in futureSim
monthlyTimespan = [startYear : timeStep : lastYear]; % list all the timespans in a vector
monthlyTimespan = monthlyTimespan(1 : end-1); % remove the very last date
monthlyTimespanFut = [endYear : timeStep : lastYear]; % screening time span starts at 2021
monthlyTimespanFut = monthlyTimespanFut(1 : end-1); 
nTimepoints = length(monthlyTimespan);
nTimepointsFut = length(monthlyTimespanFut); 
fivYrAgeGrpsOn = 1;
diseaseVec_vax = {[1:2], [3:7], 8}; % HIV negative grouped together, and then all the HIV positive states 

% scenarios = {'1.1', '1.2', '2.1', '2.2', '3.1'}; ***SET ME***: specify the scenarios to loop through
scenarios = {'0'}; 

% parallelizing the for loop
loopSegments = {0 , round(length(scenarios)/2) , length(scenarios)}; % running 10 scenarios ***SET ME***: the number of scenarios will be different
loopSegmentsLength = length(loopSegments); 

%for k = 1 : loopSegmentsLength-1 
%   parfor j = loopSegments{k}+1 : loopSegments{k+1} % for testing (parfor)
    
 for j = [1] % FORTESTING

    sceNum = j - 1; 
    sceString = scenarios{j}; % turn sceNum into string sceString
    sce = sceNum + 1; % add one since indices start at 1 (so scenarios will be 1-10 in this case) 

    % Initialize result matrices 
    ccHealthState = zeros(nTimepoints, age+1, endpoints, nRuns); 
    totalPerAge = zeros(nTimepoints, gender, age+1, nRuns); 


    % Feeding in the zeroed result matrix, spitting out the same matrix but with all the counts added in for that scenario
    %[ccHealthState, totalPerAge] = ...
    %   vaxCEA_multSims_processResults_PEPFAR(1 , sceString , {'0'}, fileInds, vax, ccDeaths, hivDeaths,  ccHealthState, hpvHealthState, newCC,  nonDisabHealthState, totalPerAge, screenTreat, screenSympCCTreat, hivHealthState, newHiv, prepCov, newCirc);  
     
    %% CC HEALTH STATES (PREVALENCE) - see lines 284-294
    %ccHealthState = zeros(nTimepoints, age+1, endpoints);
    for a = 1 : age
        for x = 1 : endpoints
            vaxInds1 = toInd(allcomb(1:disease, 1:viral, 1:hpvVaxStates, 6, x, 1:intervens, 2, a, 1:risk));
            vaxInds2 = toInd(allcomb(1:disease, 1:viral, 6, [1:5 7], x, 1:intervens, 2, a, 1:risk));
            vaxInds = [vaxInds1; vaxInds2]; 
            ccHealthState(1:end, a, x, j) = sum(vaxResult{n}.popVec(:, vaxInds), 2);
        end
    end

    %% TOTAL NUMBER OF PEOPLE PER AGE GROUP - see lines 341-355
    %totalPerAge = zeros(nTimepoints, gender, age+1);
    for a = 1 : age
        for g = 1 : gender
            inds = toInd(allcomb(1:disease, 1:viral, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:endpoints, 1:intervens, g, a, 1:risk));
            totalPerAge(:, g, a) = sum(popVec(:, inds), 2);
        end
    end
    inds1 = toInd(allcomb(1:disease, 1:viral, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:endpoints, 1:intervens, 1, 1:age, 1:risk));
    inds2 = toInd(allcomb(1:disease, 1:viral, 1:hpvVaxStates, 1:hpvNonVaxStates, 1:endpoints, 1:intervens, 2, 1:age, 1:risk));
    totalPerAge(:, 1, 17) = sum(popVec(:, inds1), 2);
    totalPerAge(:, 2, 17) = sum(popVec(:, inds2), 2);



% turn all the result matrices into 2D 
    for param = 1 : nRuns
        for a = 1 : (age + 1)
            for g = 1 : gender
            % turning total per age matrix into 2D 
                if (param == 1 && a == 1 && g == 1)
                    totalPerAgeReshape = [transpose(monthlyTimespan), g.*ones(nTimepoints,1), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ...
                                            totalPerAge(:, g, a, param)];
                                    else 
                    totalPerAgeReshape = [totalPerAgeReshape; 
                                            transpose(monthlyTimespan), g.*ones(nTimepoints,1), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ...
                                           totalPerAge(:, g, a, param)]; 
                end 
            end 

 
               for x = 1 : endpoints 
                    %dInd = 1;
                    if (param == 1 && a == 1 && x == 1)
                        %ccHealthStateReshape = [transpose(monthlyTimespan), dInd.*ones(nTimepoints,1) , a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                        %                    sce.*ones(nTimepoints,1), ccHealthState(:, a, x, param)];
                        % REMOVED DIND
                        ccHealthStateReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                            sce.*ones(nTimepoints,1), ccHealthState(:, a, x, param)];
                    
                    else 
                        %ccHealthStateReshape = [ccHealthStateReshape; 
                        %                    transpose(monthlyTimespan), dInd.*ones(nTimepoints,1) , a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                        %                    sce.*ones(nTimepoints,1), ccHealthState(:, a, x, param)];
               
                        ccHealthStateReshape = [ccHealthStateReshape; 
                                            transpose(monthlyTimespan), a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                            sce.*ones(nTimepoints,1), ccHealthState(:, a, x, param)];
                                            
                    end 
               end 

        end
        disp(['Complete Scenario ', num2str(sce), ', Parameter ', num2str(param)])
    end 
  

% turn into arrays

%ccHealthStateReshape1 = array2table(ccHealthStateReshape, 'VariableNames', {'year',  'hivState' , 'age', 'endpoint', 'paramNum', ...
%    'sceNum', 'count'});
ccHealthStateReshape1 = array2table(ccHealthStateReshape, 'VariableNames', {'year', 'age', 'endpoint', 'paramNum', ...
    'sceNum', 'count'});
totalPerAgeReshape1 = array2table(totalPerAgeReshape, 'VariableNames', {'year', 'gender', 'age', 'paramNum', 'sceNum', 'count'}); 

% instead of returning the table objects, return matrix + column names
totalPerAgeCols = {'year', 'gender', 'age', 'paramNum', 'sceNum', 'count'};
ccHealthStateCols = {'year', 'age', 'endpoint', 'paramNum', 'sceNum', 'count'}; 

end % closing the parfor loops 
%end

end