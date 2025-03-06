function [corrected, final, smooth, int] = interpolation(data, elapsed_time, new_interp_time)
%INTERPOLATION S
    
    time = timeseries(data, elapsed_time);
    int = resample(time, new_interp_time);

    %get rid of NAN in data set (first few rows) replace with 0
    int.Data(isnan(int.Data)) = 0;

    %changes time series to a matrix for later use
    corrected = int.Data;
    final = corrected;

    %for gap filling
    smooth = movmean(corrected, 8);
end

