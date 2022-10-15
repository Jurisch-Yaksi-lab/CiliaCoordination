function analysis_loop(CBF)

% List all mat files in the folder
stk_files = dir([CBF.sourceP,'*.mat']);


%% Set up file logging
diary(fullfile(CBF.folderP, sprintf('log_%s.txt',datetime("now", 'Format','dMMMy_HH-mm-ss'))))

% Initialize 
fprintf('%s\n',datetime("now"));
fprintf('%s\n',version);

%% Run through all recordings 

progress()
for iRec = 1:length(stk_files)
    progress(iRec,length(stk_files),1)
      
    % load the recording
    [~, CBF.name, ~] = fileparts(stk_files(iRec).name); % Give a names
    var = who(matfile(fullfile(CBF.sourceP,CBF.name)));
    data = double(load(fullfile(CBF.sourceP,CBF.name),var{1}).(var{1}));
    fprintf('%s\n\n',fullfile(CBF.sourceP,CBF.name));
    
    % Set some metadata
    CBF.x = size(data,1); 
    CBF.y = size(data,2); 
    
    % Retrieve the frame rate of acquisition
    try
        [~,value]=import_json([CBF.sourceP, erase(CBF.name,CBF.metadata_ID), '.json']);
        CBF.Fs = value(8); % Frequency of acquisition
    catch
        fprintf('there is probably no associated .json metadata file\n')
    end

   % Define the target path
    CBF.targetP = fullfile(CBF.folderP, CBF.name, filesep);
    [~, ~] = mkdir(CBF.targetP);
    
    
    Master_analysis_for_loop
  
    fprintf('%s was saved to %s\n\n',CBF.name, fullfile(CBF.targetP));
end

diary off