function [MoS] = calc_MoS(CoM_vec,ank_vec, CoM, CoM_vel, BoS)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% distance between CoM and ankle of heelstriking foot:
    L = norm(CoM_vec-ank_vec);

    % calculate w0 = sqrt(gravity/l);
    w0 = sqrt(9.81/L);

    % compue xCoM: CoM position + velocity/w0
    xCoM = CoM + (CoM_vel/w0);

    % BoSAP was defined by the toe marker or anterior boundary of the leading foot in seven studies [11,12,13, 18, 27, 39, 43],
    %BoSML was defined using the lateral malleolar marker [14, 23,
    %31, 40] ref numbers from: https://bmcmusculoskeletdisord.biomedcentral.com/articles/10.1186/s12891-021-04466-4

% Base of Support - xCoM (see
% https://www.sciencedirect.com/science/article/pii/S0966636224006222) 
    MoS = BoS-xCoM; 
  
end