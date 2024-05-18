
initCobraToolbox(false);
load('C:\Users\Giulia Albeniz\cobratoolbox\newModels_giuliaproject.mat', 'myModel_rxnadded_16');
%% To add metabolites to metabolic model
clc
% Load the metabolic model
load('C:\Users\Giulia Albeniz\Downloads\giuliaproject.mat', 'myModel_rxnadded_16');
% Define the path to the CSV file
csvFilePath = 'C:\Users\Giulia Albeniz\Downloads\projectmetabolites10.csv';

% Read the metabolite information from the CSV file
metabolitesData = readtable(csvFilePath);
% Loop through each row in the CSV file
for i = 1:size(metabolitesData, 1)
   % Extract metabolite information from the current row
   metaboliteID = char(metabolitesData{i, 8});
   metName = char(metabolitesData{i, 7});
   metFormula = char(metabolitesData{i, 9});
   metCharge = metabolitesData{i, 10};  % Extract charge directly as numeric;  % Convert charge to numeric
  
   % Add the metabolite to the model if it doesn't already exist
   if ~any(strcmp(myModel_rxnadded_16.mets, metaboliteID))
       myModel_rxnadded_16 = addMetabolite(myModel_rxnadded_16, metaboliteID, 'metName', metName, 'metFormula', metFormula, 'Charge', metCharge);
      
       % Display the information of the added metabolite
       disp(['Added metabolite: ', metaboliteID]);
       disp(['  Name: ', metName]);
       disp(['  Formula: ', metFormula]);
       disp(['  Charge: ', num2str(metCharge)]);
   else
       disp(['Metabolite ', metaboliteID, ' already exists in the model. Skipping...']);
   end
end
% save updated model
save('C:\Users\Giulia Albeniz\Downloads\giuliaproject.mat', 'myModel_rxnadded_16');

%% to add reactions

clc
% Define the paths to the CSV file and the original metabolic model
csvFilePath = 'C:\Users\Giulia Albeniz\Downloads\matlabdetails9.csv';
originalModelFilePath = 'C:\Users\Giulia Albeniz\Downloads\giuliaproject.mat';
% Load the original metabolic model
load(originalModelFilePath, 'myModel_rxnadded_16');

% Before adding reactions
if isfield(myModel_rxnadded_16, 'rxns')
   numReactionsBefore = numel(myModel_rxnadded_16.rxns);
   disp(['Number of reactions before adding: ', num2str(numReactionsBefore)]);
else
   disp('Warning: The loaded model does not contain reaction information.');
   numReactionsBefore = 0;
end
% Read the CSV file
reactionsData = readtable(csvFilePath);
% Iterate through each row in the CSV file
for i = 1:size(reactionsData, 1)
   % Extract reaction information from the current row
   reactionID = char(reactionsData{i, 1});
   reactionName = char(reactionsData{i, 5});
   reactionFormula = char(reactionsData{i, 3});
  
   % Check if the reaction already exists in the model
   if isfield(myModel_rxnadded_16, 'rxns') && any(strcmp(myModel_rxnadded_16.rxns, reactionID))
       disp(['Reaction ', reactionID, ' already exists in the model. Skipping...']);
   else
       % Add the new reaction to the existing model
       myModel_rxnadded_16 = addReaction(myModel_rxnadded_16, reactionID, 'reactionName', reactionName, 'reactionFormula', reactionFormula);
   end
end
% After adding reactions
if isfield(myModel_rxnadded_16, 'rxns')
   numReactionsAfter = numel(myModel_rxnadded_16.rxns);
   disp(['Number of reactions after adding: ', num2str(numReactionsAfter)]);
else
   disp('Warning: The loaded model does not contain reaction information.');
end

% Save the updated model
save('C:\Users\Giulia Albeniz\Downloads\giuliaproject.mat', 'myModel_rxnadded_16');

%% add gene rules
% Load the model
model_path = 'C:\Users\Giulia Albeniz\Downloads\giuliaproject.mat';
model_data = load(model_path);
myModel_rxnadded_16 = model_data.myModel_rxnadded_16;

% Read the data from the CSV file
data = readtable('C:\Users\Giulia Albeniz\Downloads\matlabdetails10.csv');

% Extract reaction IDs and ensemble IDs from the table
reaction_ids = data{:, 7};  % Assuming reaction IDs are in the 7th column
ensemble_ids = data{:, 6};  % Assuming ensemble IDs are in the 6th column

% Iterate over each reaction ID and its corresponding ensemble ID
for i = 1:length(reaction_ids)
    % Get the reaction ID and ensemble ID
    reaction_id = reaction_ids{i};
    ensemble_id = ensemble_ids{i};
    
    % Check if the reaction ID exists in the model
    if any(strcmp(myModel_rxnadded_16.rxns, reaction_id))
        % Add the ensemble ID to the grRule field
        myModel_rxnadded_16.grRules{strcmp(myModel_rxnadded_16.rxns, reaction_id)} = ensemble_id;
    else
        disp(['Reaction with ID ', reaction_id, ' not found in the model.']);
    end
end

% Save the updated model
save('C:\Users\Giulia Albeniz\Downloads\giuliaproject.mat', 'myModel_rxnadded_16');


%%
initCobraToolbox(false);

%%
clc
load('C:\Users\Giulia Albeniz\Downloads\giuliaproject.mat', 'myModel_rxnadded_16');
[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, Elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(myModel_rxnadded_16);


%% Add demand reactions
metaboliteList = {'MET00020[n]', 'MET00020[r]', 'MET00021[c]', 'MET00022[c]', 'MET00023[c]', 'MET00026[n]', 'MET00026[r]', 'MET00028[c]', 'MET00030[c]', 'MET00030[e]', 'MET00033[g]', 'MET00034[r]', 'MET00035[c]', 'MET00036[c]', 'MET00037[m]', 'MET00037[l]', 'MET00037[e]', 'MET00038[c]', 'MET00041[c]', 'MET00042[c]', 'MET00042[m]', 'MET00043[c]', 'MET00043[m]', 'MET00043[n]', 'MET00044[c]', 'MET00044[m]', 'MET00044[n]', 'MET00045[c]', 'MET00047[c]', 'MET00047[g]', 'MET00047[n]', 'MET00048[c]', 'MET00049[r]', 'MET00049[n]', 'MET00050[r]', 'MET00053[n]', 'MET00053[r]', 'MET00055[n]', 'MET00055[r]', 'MET00056[c]', 'MET00058[m]', 'MET00065[g]', 'MET00066[c]', 'MET00069[e]', 'MET00069[c]', 'MET00070[r]', 'MAM02786[n]', 'MAM02789[x]', 'MAM02789[r]', 'MAM02555[e]', 'MET00049[g]', 'MET00049[c]', 'MAM02786[g]', 'MET00042[r]', 'MAM01983[r]', 'MAM01802[n]'};

myModel_rxnadded_16_obj_exc = addDemandReaction(myModel_rxnadded_16_obj, metaboliteList);
%%
myModel_rxnadded_16_obj_exc_max = optimizeCbModel(myModel_rxnadded_16_obj_exc);

%%
clc
% Load the model
load('C:\Users\Giulia Albeniz\Downloads\giuliaproject.mat', 'myModel_rxnadded_16');

% Define the list of reactions
reactionsList = {
    'RXN00017', 'RXN00018', 'RXN00019', 'RXN00020', 'RXN00021', 'RXN00022', 'RXN00023', 'RXN00024', 'RXN00025', 'RXN00026', ...
    'RXN00027', 'RXN00028', 'RXN00029', 'RXN00030', 'RXN00031', 'RXN00032', 'RXN00033', 'RXN00034', 'RXN00035', 'RXN00036', ...
    'RXN00037', 'RXN00038', 'RXN00041', 'RXN00042', 'RXN00043', 'RXN00044', 'RXN00045', 'RXN00046', ...
    'RXN00047', 'RXN00048', 'RXN00049', 'RXN00050', 'RXN00051', 'RXN00052', 'RXN00053', 'RXN00054', 'RXN00055', 'RXN00056', ...
    'RXN00057', 'RXN00058', 'RXN00059', 'RXN00060', 'RXN00061', 'RXN00062', 'RXN00063', 'RXN00064', 'RXN00065', 'RXN00066', ...
    'RXN00067', 'RXN00068', 'RXN00069', 'RXN00070', 'RXN00071', 'RXN00072', 'RXN00073', 'RXN00074', 'RXN00075', 'RXN00076', ...
    'RXN00077', 'RXN00078', 'RXN00079', 'RXN00080', 'RXN00081', 'RXN00082', 'RXN00083', 'RXN00084', 'RXN00085', 'RXN00086', ...
    'RXN00087', 'RXN00088', 'RXN00089', 'RXN00090', 'RXN00091', 'RXN00092', 'RXN00093', 'RXN00094', 'RXN00095', 'RXN00096', ...
    'RXN00097', 'RXN00098', 'RXN00099', 'RXN00100', 'RXN00101', 'RXN00102', 'RXN00103', 'RXN00104', 'RXN00105', 'RXN00106', ...
    'RXN00107', 'RXN00108', 'RXN00109', 'RXN00110', 'RXN00111', 'RXN00112', 'RXN00113', 'RXN00114', 'RXN00115', 'RXN00116', ...
    'RXN00117', 'RXN00118', 'RXN00119', 'RXN00120', 'RXN00121'
};

% Flatten reactionsList to ensure it's a 1-dimensional cell array
reactionsListFlat = reshape(reactionsList.', 1, []);

% Calculate flux variability
[minFlux, maxFlux] = fluxVariability(myModel_rxnadded_16_obj_exc, 0, 'max', reactionsListFlat, 0, 1, 'FBA');

% Check the number of elements in each array
numReactions = numel(reactionsListFlat);
numMinFlux = numel(minFlux);
numMaxFlux = numel(maxFlux);

% Ensure all arrays have the same number of elements
if numReactions == numMinFlux && numReactions == numMaxFlux
    % Convert minFlux and maxFlux to cell arrays
    minFluxCell = num2cell(minFlux);
    maxFluxCell = num2cell(maxFlux);

    % Concatenate reactionsListFlat, minFluxCell, and maxFluxCell
    fluxData = [reactionsListFlat', minFluxCell, maxFluxCell];

    % Define the filename including the full path
    filename = 'C:\Users\Giulia Albeniz\Downloads\flux_data.csv';

    % Write the cell array to a CSV file
    writecell(fluxData, filename);
    disp('CSV file saved successfully.');
else
    disp('Error: Number of elements in arrays are not consistent.');
end





