clc
clear

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

    


    mos_matrix = [];
    start_leg = allevents(1, 2);

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

        %% MOS

        % MOS Inputs
        COM_AP = mlumbar(:, 1);
        COM_ML = mlumbar(:, 3);
        COM_UP = mlumbar(:, 2);
        
        rank_ap = rfoot_final(: , 1);
        rank_ml = rfoot_final(: , 3);
        rank_up = rfoot_final(: , 2);
        
        lank_ap = lfoot_final(: , 1);
        lank_ml = lfoot_final(: , 3);
        lank_up = lfoot_final(: , 2);
        
        rtoe_ap = rfoot_final(: , 1);
        
        ltoe_ap = lfoot_final(: , 1);

        hs = allevents(hs_frame, 1);
        opp_to = allevents(hs_frame + 1, 1);
        opp_hs = allevents(hs_frame + 2, 1);
        to = allevents(hs_frame + 3, 1);

        hs_foot = allevents(hs_frame ,2);

        % create CoM vector and  Ankle vector
        % at heelstrike
        CoM_Vec_at_hs = [COM_AP(hs);COM_ML(hs);COM_UP(hs)];
        RAnk_Vec_at_hs = [rank_ap(hs);rank_ml(hs);rank_up(hs)];
        LAnk_Vec_at_hs = [lank_ap(hs);lank_ml(hs);lank_up(hs)];
        % at opposite heelstrike
        CoM_Vec_at_ohs = [COM_AP(opp_hs);COM_ML(opp_hs);COM_UP(opp_hs)];
        RAnk_Vec_at_ohs = [rank_ap(opp_hs);rank_ml(opp_hs);rank_up(opp_hs)];
        LAnk_Vec_at_ohs = [lank_ap(opp_hs);lank_ml(opp_hs);lank_up(opp_hs)];
        % at toe off
        CoM_Vec_at_to = [COM_AP(to);COM_ML(to);COM_UP(to)];
        RAnk_Vec_at_to = [rank_ap(to);rank_ml(to);rank_up(to)];
        LAnk_Vec_at_to = [lank_ap(to);lank_ml(to);lank_up(to)];
        % at opposite toe off
        CoM_Vec_at_oto = [COM_AP(opp_to);COM_ML(opp_to);COM_UP(opp_to)];
        RAnk_Vec_at_oto = [rank_ap(opp_to);rank_ml(opp_to);rank_up(opp_to)];
        LAnk_Vec_at_oto = [lank_ap(opp_to);lank_ml(opp_to);lank_up(opp_to)];

        % calculate center of mass velocity
        fr = 100;
        dt = 1/fr;
        
        for i = 1:length(COM_AP)-1
            CoM_vel_AP(i,1) = (COM_AP(i+1)-COM_AP(i))/dt;
            CoM_vel_ML(i,1) = (COM_ML(i+1)-COM_ML(i))/dt;
        end
        
    
        if leg_side ~= start_leg
            continue
        end

        if leg_side == 3% left heel strike/gc
        LBoS_AP_hs = ltoe_ap(hs);
        LBoS_ML_hs = lank_ml(hs);
        RBoS_AP_hs = rtoe_ap(opp_hs);
        RBoS_ML_hs = rank_ml(opp_hs);
        LBoS_AP_to = ltoe_ap(to);
        LBoS_ML_to = lank_ml(to);
        RBoS_AP_to = rtoe_ap(opp_to);
        RBoS_ML_to = rank_ml(opp_to);
    
        %acl_MoS(CoM_Vec,ank_vec, CoM, CoM_vel, BoS)
    
        L_MoS_AP_hs = calc_MoS(CoM_Vec_at_hs,LAnk_Vec_at_hs,COM_AP(hs),CoM_vel_AP(hs),LBoS_AP_hs);
        R_MoS_AP_hs = calc_MoS(CoM_Vec_at_ohs,RAnk_Vec_at_ohs,COM_AP(opp_hs),CoM_vel_AP(opp_hs),RBoS_AP_hs);
        L_MoS_ML_hs = calc_MoS(CoM_Vec_at_hs,LAnk_Vec_at_hs,COM_ML(hs),CoM_vel_ML(hs),LBoS_ML_hs);
        R_MoS_ML_hs = calc_MoS(CoM_Vec_at_ohs,RAnk_Vec_at_ohs,COM_ML(opp_hs),CoM_vel_ML(opp_hs),RBoS_ML_hs);
        L_MoS_AP_to = calc_MoS(CoM_Vec_at_to,LAnk_Vec_at_to,COM_AP(to),CoM_vel_AP(to),LBoS_AP_to);
        R_MoS_AP_to = calc_MoS(CoM_Vec_at_oto,RAnk_Vec_at_oto,COM_AP(opp_to),CoM_vel_AP(opp_to),RBoS_AP_to);
        L_MoS_ML_to = calc_MoS(CoM_Vec_at_to,LAnk_Vec_at_to,COM_ML(to),CoM_vel_ML(to),LBoS_ML_to);
        R_MoS_ML_to = calc_MoS(CoM_Vec_at_oto,RAnk_Vec_at_oto,COM_ML(opp_to),CoM_vel_ML(opp_to),RBoS_ML_to);
    
        
        else % right leg heelstrike
           LBoS_AP_hs = ltoe_ap(opp_hs);
            LBoS_ML_hs = lank_ml(opp_hs);
            RBoS_AP_hs = rtoe_ap(hs);
            RBoS_ML_hs = rank_ml(hs);
            LBoS_AP_to = ltoe_ap(opp_to);
            LBoS_ML_to = lank_ml(opp_to);
            RBoS_AP_to = rtoe_ap(to);
            RBoS_ML_to = rank_ml(to);
        
            L_MoS_AP_hs = calc_MoS(CoM_Vec_at_ohs,LAnk_Vec_at_ohs,COM_AP(opp_hs),CoM_vel_AP(opp_hs),LBoS_AP_hs);
            R_MoS_AP_hs = calc_MoS(CoM_Vec_at_hs,RAnk_Vec_at_hs,COM_AP(hs),CoM_vel_AP(hs),RBoS_AP_hs);
            L_MoS_ML_hs = calc_MoS(CoM_Vec_at_ohs,LAnk_Vec_at_ohs,COM_ML(opp_hs),CoM_vel_ML(opp_hs),LBoS_ML_hs);
            R_MoS_ML_hs = calc_MoS(CoM_Vec_at_hs,RAnk_Vec_at_hs,COM_ML(hs),CoM_vel_ML(hs),RBoS_ML_hs);
            L_MoS_AP_to = calc_MoS(CoM_Vec_at_oto,LAnk_Vec_at_oto,COM_AP(opp_to),CoM_vel_AP(opp_to),LBoS_AP_to);
            R_MoS_AP_to = calc_MoS(CoM_Vec_at_to,RAnk_Vec_at_to,COM_AP(to),CoM_vel_AP(to),RBoS_AP_to);
            L_MoS_ML_to = calc_MoS(CoM_Vec_at_oto,LAnk_Vec_at_oto,COM_ML(opp_to),CoM_vel_ML(opp_to),LBoS_ML_to);
            R_MoS_ML_to = calc_MoS(CoM_Vec_at_to,RAnk_Vec_at_to,COM_ML(to),CoM_vel_ML(to),RBoS_ML_to);
        
            
        end
    
        
        
        
        mos_gait = [L_MoS_AP_hs R_MoS_AP_hs L_MoS_ML_hs R_MoS_ML_hs L_MoS_AP_to R_MoS_AP_to L_MoS_ML_to R_MoS_ML_to];

        mos_matrix = [mos_gait; mos_matrix];
        
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


end