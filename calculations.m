%% LOAD IN .MAT FILES FROM PROCESSING
clc;
clear all;
folderPath = "Saved_Process";

matFiles = dir(fullfile(folderPath, '*.mat'));

for file = 1:length(matFiles)

    % Get the full file name
    fileName = matFiles(file).name;
    fullFilePath = fullfile(folderPath, fileName);


    load(fullFilePath);

    % define eventframes with labels 
    RHS(:, 1) = rmintab; %col 1 is the frames from RHS
    RHS(:, 2) = repmat(1, length(rmintab), 1); %col 2 is the label to categorize this as RHS
    LTO(:, 1) = ltoeloc;
    LTO(:, 2) = repmat(2, length(ltoeloc), 1);
    LHS(: ,1) = lmintab;
    LHS(:, 2) = repmat(3, length(lmintab), 1);
    RTO(:, 1) = rtoeloc;
    RTO(:, 2) = repmat(4, length(rtoeloc), 1);
    allevents = vertcat(RHS, LTO, LHS, RTO);
    allevents = sortrows(allevents, 1);
    
    %%
    
    P_output = struct();
    step_counter = 1;
    
    for hs_frame = 1:length(allevents)
        % Try catch to end the loop when we get to the end
        try 
            allevents(hs_frame + 2, 2);
        catch
            disp("End of Events")
            break
        end
        
        % Skip over toe offs
        if allevents(hs_frame,2) == 2 || allevents(hs_frame,2) == 4
            continue
        end
    
    
        % Skip over turns
        if abs(allevents(hs_frame,1) - allevents(hs_frame+2,1)) > 200
            continue
        end
        
        % Assign current and previous frame
        current_frame = allevents(hs_frame, 1);
        next_frame = allevents(hs_frame + 2, 1);
        
    
        %%% Calculations
        
       % Determine leg
        leg_side = allevents(hs_frame, 2);
        
        % We do this to set the correct data for indexing based on which foot
        % is the starting foot
        if leg_side == 1
            side = "Right";
            foot_final = rfoot_final;
            other_foot_final = lfoot_final;
        elseif leg_side == 3
            side = "Left";
            foot_final = lfoot_final;
            other_foot_final = rfoot_final;
        end
    
        
        % Determing parameters based on index of our current heel strike and
        % finding the difference using the next opposite heel strike
        step_length = abs(foot_final(current_frame, 1) - other_foot_final(next_frame, 1));
        step_width = abs(foot_final(current_frame, 3) - other_foot_final(next_frame, 3));
        gait_speed = (step_length) / ((.01)*(next_frame - current_frame));
    
        steps(step_counter, 1) = step_length;
        steps(step_counter, 2) = step_width;
        steps(step_counter, 3) = gait_speed;
        
        
        subject = strcat("Subject_", headerInfo.Subject);  % Get the subject from header info
        trial = strcat("Trial_", headerInfo.Trial);     % Get the trial name from header info
    
        % get subject number as single digit
        % Convert the cell value to a string
        str = string(headerInfo.Subject);
        % Replace all non-numeric characters with an empty string
        str = regexprep(str, '[^0-9\.]', '');
        % Convert the string back to a number and store it in the cell array
        headerInfo.Subject = str2double(str);
        % Convert the number back to a string and store it in the cell array
        headerInfo.Subject = num2str(headerInfo.Subject);
    
        Output(file).Subject = headerInfo.Subject;
        Output(file).Trial = headerInfo.Trial;
        Output(file).Task = headerInfo.Task;
        Output(file).Elevation = headerInfo.Elevation;
        Output(file).Steplength = steps(:,1);
        Output(file).Stepwidth = steps(:,2);
        Output(file).gaitspeed =  steps(:,3);
        FullOutput = Output;
        step_counter = step_counter + 1;
    
    end

    clear RHS LHS LTO RTO rmintab rmaxtab lmintab lmaxtab ltoeloc ltoemin rtoeloc rtoemin allevents

    FullOutput(file).meanStepLength = mean(FullOutput(file).Steplength);
    FullOutput(file).medianStepLength = median(FullOutput(file).Steplength);
    FullOutput(file).sdStepLength = std(FullOutput(file).Steplength);
    FullOutput(file).maxStepLength = max(FullOutput(file).Steplength);
    FullOutput(file).minStepLength = min(FullOutput(file).Steplength);
    FullOutput(file).meanStepWidth = mean(FullOutput(file).Stepwidth);
    FullOutput(file).medianStepWidth = median(FullOutput(file).Stepwidth);
    FullOutput(file).sdStepWidth = std(FullOutput(file).Stepwidth);
    FullOutput(file).maxStepWidth = max(FullOutput(file).Stepwidth);
    FullOutput(file).minStepWidth = min(FullOutput(file).Stepwidth);
    FullOutput(file).meanGaitSpeed = mean(FullOutput(file).gaitspeed);
    FullOutput(file).medianGaitSpeed = median(FullOutput(file).gaitspeed);
    FullOutput(file).sdGaitSpeed = std(FullOutput(file).gaitspeed);
    FullOutput(file).maxGaitSpeed = max(FullOutput(file).gaitspeed);
    FullOutput(file).minGaitSpeed = min(FullOutput(file).gaitspeed);



end


savedir = strcat(pwd, '/PDTVR');
fullfile(savedir, 'FullOutput.mat')
save(fullfile(savedir, 'FullOutput.mat'), 'FullOutput');
% save(fullfile(savedir, 'fullframesoutput.mat'), 'fullframesoutput')

%%
% save for future opening
 savedir = strcat(pwd, '/PDTVR');
  save(fullfile(savedir, 'FullOutput.mat'), 'FullOutput');
  load(fullfile(savedir, 'FullOutput.mat'))
  % load(fullfile(savedir, 'fullframesoutput.mat'))

FullOutput = load(strcat(pwd, '/PDTVR/FullOutput.mat'));
FullOutput = FullOutput.FullOutput;
FullOutput = struct2table(FullOutput, 'AsArray', true);
% fullframesoutput = struct2table(fullframesoutput);

%% Exporting the Cadence into it's own matrix/csv

% 
% CadenceOnly = FullOutput.cadence;
% 
% writematrix(CadenceOnly, 'CadenceOnly.csv');

%%
%FullOutput = struct2table(FullOutput);
sz = [1, 7];
varNames = {'Subject', 'Task', 'Elevation', 'StepNum', 'Steplength', 'Stepwidth', 'gaitspeed'};
varTypes = {'string', 'string', 'string', 'double', 'double', 'double', 'double'};
TableOutput = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
TableOutput(1, :) = {"Subject", "Task", "Elevation", "StepNum", "Steplength", "Stepwidth", "gaitspeed"};

for k = 1:size(FullOutput,1)

    % Subject ID 
    Subject = FullOutput(k, 1);
    % Task 
    Task = FullOutput(k, 3);
    % Elevation 
    Elevation = FullOutput(k, 4);

    NumbSteps = size(FullOutput.Steplength{k})
    StepNum = table([1 : NumbSteps]');
    StepNum.Properties.VariableNames = "StepNum";
    % Steplength
    Steplength = array2table(FullOutput.Steplength{k});
    Steplength.Properties.VariableNames = "Steplength";
    % Stepwidth
    Stepwidth = array2table(FullOutput.Stepwidth{k});
    Stepwidth.Properties.VariableNames = "Stepwidth";
    % Gaitspeed
    gaitspeed = array2table(FullOutput.gaitspeed{k})
    gaitspeed.Properties.VariableNames = "gaitspeed";
    % Cadence
    % Cadence = FullOutput(k, 8);
    

    Subject = repmat(Subject,NumbSteps);
    Task = repmat(Task,NumbSteps); 
    Elevation = repmat(Elevation, NumbSteps);
    % Cadence = repmat(Cadence, NumbSteps)
    
%     datamatrix(k, 1) = Subject;
%     datamatrix(k, 2) = Task;
%     datamatrix(k, 3) = Elevation;
%     datamatrix(k, 4) = StepNum
%     datamatrix(k, 5) = vertcat(Steplength_start, Steplength_end);
    
    % Combine them into a new table
    % datamatrix = horzcat(Subject, Task, Elevation, StepNum, Steplength, Stepwidth, gaitspeed, Cadence);
    datamatrix = horzcat(Subject, Task, Elevation, StepNum, Steplength, Stepwidth, gaitspeed);
    TableOutput = vertcat(TableOutput, datamatrix);


% % Display the new table
% disp(newTable);

end 
writetable(TableOutput, 'HMD_Output78.csv');
