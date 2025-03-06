function [data_final] = gap_fill(data_final, data_smooth8, data_int, origin)
%GAP_FILL 
        fs = 100;
        
        if any(data_smooth8(:,2)>15)
            errorThresholdL= -origin(:,2)+15;
        else
            errorThresholdL= -origin(:,2);
        end

        k = 1;
        vrTime = data_int.Time;
        [badL] = find(data_smooth8(:,2)<errorThresholdL);
        if isempty(badL)
            gap_fill = false;
        else
            gap = find(diff(badL)>k*1);

            if isempty(gap)
                gap = length(badL);
            end

            badL_wGap = [];
            for n=1:length(gap)
                if n==1
                    badL_wGap(n,:) = [max([badL(1)-k; 1]) badL(gap(n))+k];
                else
                    badL_wGap(n,:) = [badL(gap(n-1)+1)-k badL(gap(n))+k];
                end
            end

            badL_complete = [];
            for n=1:length(badL_wGap(:,1))
                badL_complete = [badL_complete;[badL_wGap(n,1):badL_wGap(n,2)]'];
            end

            data_interpolated = [vrTime, data_smooth8];
            data_interpolated(badL_complete,:) = [];


            for n=1:3
                s(:,n) = interp1(data_interpolated(:,1), data_interpolated(:,n+1), vrTime(badL_complete), 'l');
            end
            data_final(badL_complete,:) = s;
            gap_fill = true;
        end
end

