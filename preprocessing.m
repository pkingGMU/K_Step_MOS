%%% Preprocessing focused on adjusting global to local coordinate system
%%% and filtering

clear all
close all
clc

% Local Imports
addpath(genpath('Functions'));

% Data directory and file list
directory = strcat(pwd, '/Data');
fileList = dir(directory);
fileList([fileList.isdir],:)= [];
fileList = fileList(contains({fileList.name}, '.csv'));

for file = 1:height(fileList)
    % Refresh variables
    clear RHS LHS LTO RTO rmintab rmaxtab lmintab lmaxtab ltoeloc ltoemin rtoeloc rtoemin allevents

    filename = fileList(file).name;
    %% Set up the Import Options and import the data
        opts = delimitedTextImportOptions("NumVariables", 16);
        opts.ExtraColumnsRule='ignore';

        % Specify range and delimiter
        opts.DataLines = [1, 2];
        opts.Delimiter = ",";

        %Import two data tables
        headerInfo = readtable(strcat(directory,"/",filename), opts);
        clear opts;

        headerInfo.Properties.VariableNames = ["Subject","Trial", "Task", ...
            "Elevation", "Corner_BottomLeft_PositionX", "Corner_BottomLeft_PositioY", ...
            "Corner_BottomLeft_PositionZ",...
            "Corner_BottomRight_PositionX", "Corner_BottomRight_PositioY", ...
            "Corner_BottomRight_PositionZ", "Corner_TopRight_PositionX","Corner_TopRight_PositioY", ...
            "Corner_TopRight_PositionZ", "Corner_TopLeft_PositionX","Corner_TopLeft_PositioY",...
            "Corner_TopLeft_PositionZ"];
        headerInfo = headerInfo(2,:);
        data = readtable(strcat(directory,"/",filename), 'HeaderLines', 2);
        %% define plank coordinates
        
        BottomLeft = str2double([headerInfo.Corner_BottomLeft_PositionX, headerInfo.Corner_BottomLeft_PositioY,...
            headerInfo.Corner_BottomLeft_PositionZ]);

        BottomRight= str2double([headerInfo.Corner_BottomRight_PositionX, headerInfo.Corner_BottomRight_PositioY,...
            headerInfo.Corner_BottomRight_PositionZ]);

        TopLeft = str2double([headerInfo.Corner_TopLeft_PositionX, headerInfo.Corner_TopLeft_PositioY,...
            headerInfo.Corner_TopLeft_PositionZ]);

        TopRight= str2double([headerInfo.Corner_TopRight_PositionX, headerInfo.Corner_TopRight_PositioY,...
            headerInfo.Corner_TopRight_PositionZ]);
        
        %%%Getting trajectories for sensors in linear matrix
        LFootlin = [data.leftFootPositionX data.leftFootPositionY data.leftFootPositionZ];
        RFootlin = [data.rightFootPositionX data.rightFootPositionY data.rightFootPositionZ];
        HMDlin = [data.HeadSetPositionX data.HeadSetPositionY data.HeadSetPositionZ];
        FLumbarlin = [data.FrontLumbarPositionX data.FrontLumbarPositionY data.FrontLumbarPositionZ];
        BLumbarlin = [data.BackLumbarPositionX data.BackLumbarPositionY data.BackLumbarPositionZ];

        %create average vectors from existing coordinates
        origin = mean([BottomLeft; BottomRight],1);
        endpoint = mean([TopLeft; TopRight],1);

        %set direction of vectors
        x = (endpoint - origin);
        z = (BottomRight - BottomLeft);
        y = cross(z,x);

        %normalize z using x because its most reliable, using x,y to orient vector in
        %proper direction
        z = cross (x,y);

        %logic statement to check the vectors are orthogonal
        dot(z,y) == 0 & dot(z,x) == 0;

        %determine magnitude of vectors
        i = x ./ norm(x);
        j = y ./ norm(y);
        k = z ./ norm(z);

        G2A = [i; j; k;];
        %A2G = G2A';

        xDir = [origin; i+origin];
        yDir = [origin; j+origin];
        zDir = [origin; k+origin];

        %%Redefine plank positions
        %subtract origin value from all marker data
        %transposed first for matrix math and transposed again for graphing
        BLeft = [G2A*[BottomLeft-origin]']';
        BRight = [G2A*[BottomRight-origin]']';
        TLeft = [G2A*[TopLeft-origin]']';
        TRight = [G2A*[TopRight-origin]']';

        %plot coordinate system back in local frame, not global
        i2 = [G2A*i']';
        j2 = [G2A*j']';
        k2 = [G2A*k']';

        %bringing back to local 000 origin not the newly defined origin for plot
        xDir = [[0,0,0]; i2];
        yDir = [[0,0,0]; j2];
        zDir = [[0,0,0]; k2];

        %% Creating corrected trajectories
        %linear
        lfoot = (G2A * (LFootlin - origin)')';
        rfoot = (G2A * (RFootlin - origin)')';
        HMD = (G2A * (HMDlin - origin)')';
        flumbar = (G2A * (FLumbarlin - origin)')';
        blumbar = (G2A * (BLumbarlin - origin)')';

        %%flipping the axes for so that vertical is positive
        HMD = (HMD .* [1, -1, 1]);
        lfoot = (lfoot .* [1, -1, 1]);
        rfoot = (rfoot .* [1, -1, 1]);
        flumbar = (flumbar .* [1, -1, 1]);
        blumbar = (blumbar .* [1, -1, 1]);

        %% interpolate signal for equal time sample
        %%HMD
        HMD_time = timeseries(HMD, data.ElapsedTime);
        new_interp_time = (0:0.01:data.ElapsedTime(end,:));
        HMD_int = resample(HMD_time, new_interp_time);

        %get rid of NAN in data set (first few rows) replace with 0
        HMD_int.Data(isnan(HMD_int.Data)) = 0;

        %changes time series to a matrix for later use
        HMD_corrected = HMD_int.Data;

        %%rfoot
        [rfoot_corrected, rfoot_final, rfoot_smooth8, rfoot_int] = interpolation(rfoot, data.ElapsedTime, new_interp_time);
        %%lfoot
        [lfoot_corrected, lfoot_final, lfoot_smooth8, lfoot_int] = interpolation(lfoot, data.ElapsedTime, new_interp_time);
        %%flumbar
        [flumbar_corrected, flumbar_final, flumbar_smooth8, flumbar_int] = interpolation(flumbar, data.ElapsedTime, new_interp_time);
        %%blumbar
        [blumbar_corrected, blumbar_final, blumbar_smooth8, blumbar_int] = interpolation(blumbar, data.ElapsedTime, new_interp_time);

        %% Gap Filling
   
        % left foot
        [lfoot_final] = gap_fill(lfoot_final, lfoot_smooth8, lfoot_int, origin);
        % right foot
        [rfoot_final] = gap_fill(rfoot_final, rfoot_smooth8, rfoot_int, origin);
        % front lumbar
        [flumbar_final] = gap_fill(flumbar_final, flumbar_smooth8, flumbar_int, origin);
        % front lumbar
        [blumbar_final] = gap_fill(blumbar_final, blumbar_smooth8, blumbar_int, origin);

        % Find the rows in rfoot_final that have at least one non-zero element
        non_empty_rows_rfoot_final = any(rfoot_final,2);

        % Find the rows in lfoot_final that have at least one non-zero element
        non_empty_rows_lfoot_final = any(lfoot_final,2);

        % Find the rows in flumbar_final that have at least one non-zero element
        non_empty_rows_flumbar_final = any(flumbar_final,2);

        % Find the rows in blumbar_final that have at least one non-zero element
        non_empty_rows_blumbar_final = any(blumbar_final,2);

        % Find the rows in rfoot_final that have at least one non-zero element
        non_empty_rows_HMD_corrected = any(HMD_corrected,2);

        % Find the common rows that have at least one non-zero element in all three matrices
        common_non_empty_rows = non_empty_rows_rfoot_final & non_empty_rows_lfoot_final & non_empty_rows_HMD_corrected;

        % Extract the common rows from all three matrices
        rfoot_final = rfoot_final(common_non_empty_rows,:);
        lfoot_final = lfoot_final(common_non_empty_rows,:);
        flumbar_final = flumbar_final(common_non_empty_rows,:);
        blumbar_final = blumbar_final(common_non_empty_rows,:);
        HMD_corrected = HMD_corrected(common_non_empty_rows,:);

        %filtering parameters
        fc = 6;
        fs = 100;
        [b,a] = butter(2,fc/(fs/2));

        %filter corrected trajectory
        HMD_filt= filtfilt(b,a,HMD_corrected);

        %%rfoot
        rfoot_filt= filtfilt(b,a,rfoot_final);
        %%lfoot
        lfoot_filt= filtfilt(b,a,lfoot_final);
        %%rfoot
        flumbar_filt= filtfilt(b,a,flumbar_final);
        %%lfoot
        blumbar_filt= filtfilt(b,a,blumbar_final);


        %%use velocity of HMD loop calculator here
        velHead = zeros(length(HMD_filt),3);
        velRFoot = zeros(length(rfoot_filt),3);
        velLFoot = zeros(length(lfoot_filt),3);
        velflumbar = zeros(length(flumbar_filt),3);
        velblumbar = zeros(length(blumbar_filt),3);

        for n=2:length(velHead)-1
            velHead(n,:) = (HMD_filt(n+1,:) - HMD_filt(n-1,:)) ./(new_interp_time(n+1) - new_interp_time(n-1));
        end

        for n=2:length(velRFoot)-1
            velRFoot(n,:) = (rfoot_filt(n+1,:) - rfoot_filt(n-1,:)) ./(new_interp_time(n+1) - new_interp_time(n-1));
        end

        for n=2:length(velLFoot)-1
            velLFoot(n,:) = (lfoot_filt(n+1,:) - lfoot_filt(n-1,:)) ./(new_interp_time(n+1) - new_interp_time(n-1));
        end

        for n=2:length(velflumbar)-1
            velflumbar(n,:) = (flumbar_filt(n+1,:) - flumbar_filt(n-1,:)) ./(new_interp_time(n+1) - new_interp_time(n-1));
        end

        for n=2:length(velblumbar)-1
            velblumbar(n,:) = (blumbar_filt(n+1,:) - blumbar_filt(n-1,:)) ./(new_interp_time(n+1) - new_interp_time(n-1));
        end

        %%smooth the filtered data to calculate still period better
        lfoot_velocitysmooth = movmean(vecnorm(velLFoot(2:end,1:2)'), 20);
        rfoot_velocitysmooth = movmean(vecnorm(velRFoot(2:end,1:2)'), 20);
        
        flumbar_velocitysmooth = movmean(vecnorm(velflumbar(2:end,1:2)'), 20);
        blumbar_velocitysmooth = movmean(vecnorm(velblumbar(2:end,1:2)'), 20);
       
        %% creating a distance vector between the HMD and the foot tracker
        rGE_vector = (rfoot_filt(:,1)-HMD_filt(:,1)).*sign(velHead(:,1));
        lGE_vector = (lfoot_filt(:,1)-HMD_filt(:,1)).*sign(velHead(:,1));

        %% Midpoint lumbar
        mlumbar = (flumbar_filt + blumbar_filt) / 2;

        figure;
        hold on
        plot(mlumbar(:,1), 'b'); 
        plot(mlumbar(:,2), 'r'); 
        plot(mlumbar(:,3), 'g'); 
        hold off

        %find peaks that represent heel contact-max distance between HMD and
        %tracker

        [rmaxtab, rmintab] = findpeaks(rGE_vector,'MinPeakProminence',0.1);
        [lmaxtab,lmintab] = findpeaks(lGE_vector,'MinPeakProminence',0.1);

        [rtoemin, rtoeloc] = findpeaks(-rGE_vector, 'MinPeakProminence', 0.1);
        [ltoemin, ltoeloc] = findpeaks(-lGE_vector, 'MinPeakProminence', 0.1);


        %% Right Plot Selection 
        
        % Step 3: Plot the original signal and detected peaks
        figure;
        plot(rGE_vector, 'b');  % Plot the right signal in blue
        hold on;
        plot(rmintab, rmaxtab, 'mo', 'MarkerFaceColor', 'm');  % Mark right peaks with red circles
        
        plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
        plot(lmintab, lmaxtab, 'co', 'MarkerFaceColor', 'c');  % Mark left peaks with magenta circles
        
        plot(rGE_vector, 'b');  % Plot the right signal in blue
        hold on;
        plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark right peaks with red circles
        
        plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
        plot(ltoeloc, -ltoemin, 'go', 'MarkerFaceColor', 'g');  % Mark left peaks with magenta circles
        
        title('Detected Peaks (Click on peaks to remove from PINK DOTS)');
        xlabel('Sample Index');
        ylabel('Signal Amplitude');
        
        % Step 4: Allow user to manually select peaks (click on graph)
        disp('Click on peaks you want to remove. Press Enter when done.');
        [x_manualr, ~] = ginput;  % User clicks on graph to select peaks
        
        
        % Step 5: Match manually selected points to the closest peaks
        tolerance = 60;  % Increased tolerance to allow for better matching
        selected_peaksr = false(size(rmintab));  % Logical array to track selected peaks
        
        for i = 1:length(x_manualr)
            % Find the index of the closest peak (in terms of x-values/indices) to the selected point
            [~, idxr] = min(abs(rmintab - round(x_manualr(i))));
            
            % Check if the selected point is within the tolerance range of a detected peak index
            if abs(rmintab(idxr) - round(x_manualr(i))) <= tolerance
                selected_peaksr(idxr) = true;  % Mark peak as selected
            end
        end
        
        % Step 6: Filter the peaks to delete the selected ones
        rmintab(selected_peaksr) = [];
        rmaxtab(selected_peaksr) = [];
        
        
        %% Left Plot Selection 
        
        % Step 3: Plot the original signal and detected peaks
        figure;
        plot(rGE_vector, 'b');  % Plot the right signal in blue
        hold on;
        plot(rmintab, rmaxtab, 'mo', 'MarkerFaceColor', 'm');  % Mark right peaks with red circles
        
        plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
        plot(lmintab, lmaxtab, 'co', 'MarkerFaceColor', 'c');  % Mark left peaks with magenta circles
        
        plot(rGE_vector, 'b');  % Plot the right signal in blue
        hold on;
        plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark right peaks with red circles
        
        plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
        plot(ltoeloc, -ltoemin, 'go', 'MarkerFaceColor', 'g');  % Mark left peaks with magenta circles
        
        title('Detected Peaks (Click on peaks to remove from CYAN DOTS)');
        xlabel('Sample Index');
        ylabel('Signal Amplitude');
        
        
        % Step 4: Allow user to manually select peaks (click on graph)
        disp('Click on peaks you want to remove. Press Enter when done.');
        [x_manuall, ~] = ginput;  % User clicks on graph to select peaks
        
        
        % Step 5: Match manually selected points to the closest peaks
        tolerance = 60;  % Increased tolerance to allow for better matching
        selected_peaksl = false(size(lmintab));  % Logical array to track selected peaks
        
        for j = 1:length(x_manuall)
            % Find the index of the closest peak (in terms of x-values/indices) to the selected point
            [~, idxl] = min(abs(lmintab - round(x_manuall(j))));
            
            % Check if the selected point is within the tolerance range of a detected peak index
            if abs(lmintab(idxl) - round(x_manuall(j))) <= tolerance
                selected_peaksl(idxl) = true;  % Mark peak as selected
            end
        end
        
        % Step 6: Filter the peaks to delete the selected ones
        lmintab(selected_peaksl) = [];
        lmaxtab(selected_peaksl) = [];
        
        
        % Step 7: Update the plot with the selected peaks
        figure;
        plot(lGE_vector, 'b');  % Plot the signal
        hold on;
        plot(lmintab, lmaxtab, 'ro', 'MarkerFaceColor', 'r');  % Mark selected peaks with green circles
        title('Manually Selected Peaks');
        hold on
        plot(rGE_vector, 'b');  % Plot the signal
        hold on;
        plot(rmintab, rmaxtab, 'bo', 'MarkerFaceColor', 'b');  % Mark selected peaks with green circles
        title('Manually Selected Peaks');
        plot(rGE_vector, 'b');  % Plot the right signal in blue
        hold on;
        plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark right peaks with red circles
        
        plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
        plot(ltoeloc, -ltoemin, 'go', 'MarkerFaceColor', 'g');  % Mark left peaks with magenta circles
        
        
        disp('Selected peaks are now displayed in red.');
        disp('Selected peaks are now displayed in blue.');
        
        
        
        %% Right Toe Plot Selection 
        
        % Step 3: Plot the original signal and detected peaks
        figure;
        plot(rGE_vector, 'b');  % Plot the right signal in blue
        hold on;
        plot(rmintab, rmaxtab, 'mo', 'MarkerFaceColor', 'm');  % Mark right peaks with red circles
        
        plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
        plot(lmintab, lmaxtab, 'co', 'MarkerFaceColor', 'c');  % Mark left peaks with magenta circles
        
        plot(rGE_vector, 'b');  % Plot the right signal in blue
        hold on;
        plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark right peaks with red circles
        
        plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
        plot(ltoeloc, -ltoemin, 'go', 'MarkerFaceColor', 'g');  % Mark left peaks with magenta circles
        
        title('Detected Peaks (Click on peaks to remove from BLUE DOTS)');
        xlabel('Sample Index');
        ylabel('Signal Amplitude');
        
        % Step 4: Allow user to manually select peaks (click on graph)
        disp('Click on peaks you want to remove. Press Enter when done.');
        [x_mantoer, ~] = ginput;  % User clicks on graph to select peaks
        
        
        % Step 5: Match manually selected points to the closest peaks
        tolerance = 60;  % Increased tolerance to allow for better matching
        selected_peakstoer = false(size(rtoeloc));  % Logical array to track selected peaks
        
        for i = 1:length(x_mantoer)
            % Find the index of the closest peak (in terms of x-values/indices) to the selected point
            [~, idxr] = min(abs(rtoeloc - round(x_mantoer(i))));
            
            % Check if the selected point is within the tolerance range of a detected peak index
            if abs(rtoeloc(idxr) - round(x_mantoer(i))) <= tolerance
                selected_peakstoer(idxr) = true;  % Mark peak as selected
            end
        end
        
        % Step 6: Filter the peaks to delete the selected ones
        rtoeloc(selected_peakstoer) = [];
        rtoemin(selected_peakstoer) = [];
        
        
        %% Left Plot Selection 
        
        figure;
        plot(rGE_vector, 'b');  % Plot the right signal in blue
        hold on;
        plot(rmintab, rmaxtab, 'mo', 'MarkerFaceColor', 'm');  % Mark right peaks with red circles
        
        plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
        plot(lmintab, lmaxtab, 'co', 'MarkerFaceColor', 'c');  % Mark left peaks with magenta circles
        
        plot(rGE_vector, 'b');  % Plot the right signal in blue
        hold on;
        plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark right peaks with red circles
        
        plot(lGE_vector, 'Color', [0, 0.5, 0]);  % Plot the left signal in darker green
        plot(ltoeloc, -ltoemin, 'go', 'MarkerFaceColor', 'g');  % Mark left peaks with magenta circles
        
        
        title('Detected Peaks (Click on peaks to remove from GREEN DOTS)');
        xlabel('Sample Index');
        ylabel('Signal Amplitude');
        
        % Step 4: Allow user to manually select peaks (click on graph)
        disp('Click on peaks you want to remove. Press Enter when done.');
        [x_mantoel, ~] = ginput;  % User clicks on graph to select peaks
        
        
        % Step 5: Match manually selected points to the closest peaks
        tolerance = 60;  % Increased tolerance to allow for better matching
        selected_peakstoel = false(size(ltoeloc));  % Logical array to track selected peaks
        
        for j = 1:length(x_mantoel)
            % Find the index of the closest peak (in terms of x-values/indices) to the selected point
            [~, idxl] = min(abs(ltoeloc - round(x_mantoel(j))));
            
            % Check if the selected point is within the tolerance range of a detected peak index
            if abs(ltoeloc(idxl) - round(x_mantoel(j))) <= tolerance
                selected_peakstoel(idxl) = true;  % Mark peak as selected
            end
        end
        
        % Step 6: Filter the peaks to delete the selected ones
        ltoeloc(selected_peakstoel) = [];
        ltoemin(selected_peakstoel) = [];
        
        
        % Step 7: Update the plot with the selected peaks
        figure;
        plot(lGE_vector, 'r');  % Plot the signal
        hold on;
        plot(ltoeloc, -ltoemin, 'ro', 'MarkerFaceColor', 'r');  % Mark selected peaks with green circles
        title('Manually Selected Peaks');
        hold on
        plot(rGE_vector, 'b');  % Plot the signal
        hold on;
        plot(rtoeloc, -rtoemin, 'bo', 'MarkerFaceColor', 'b');  % Mark selected peaks with green circles
        title('Manually Selected Peaks');
        plot(lGE_vector, 'r');  % Plot the signal
        hold on;
        plot(lmintab, lmaxtab, 'ro', 'MarkerFaceColor', 'r');  % Mark selected peaks with green circles
        title('Manually Selected Peaks');
        hold on
        plot(rGE_vector, 'b');  % Plot the signal
        hold on;
        plot(rmintab, rmaxtab, 'bo', 'MarkerFaceColor', 'b');  % Mark selected peaks with green circles
        title('Manually Selected Peaks');
        disp('Selected peaks are now displayed in red.');
        disp('Selected peaks are now displayed in blue.');
        
        clearvars -except mlumbar rmintab rmaxtab ltoeloc rtoeloc ltoemin rtoemin lmintab lmaxtab lfoot_final rfoot_final headerInfo ...
            fileList directory
        cd("Saved_Process")
        save(strcat(string(headerInfo.Subject), string(headerInfo.Trial)));
        cd ..



end




