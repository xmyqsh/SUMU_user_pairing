    MCS(1,:) = [2 0.07617 0.1523]; 
    MCS(2,:) = [2 0.11719 0.2344]; 
    MCS(3,:) = [2 0.18848 0.3770]; 
    MCS(4,:) = [2 0.30078 0.6016]; 
    MCS(5,:) = [2 0.43848 0.8770]; 
    MCS(6,:) = [2 0.58789 1.1758];       
    MCS(7,:) = [4 0.36914 1.4766]; 
    MCS(8,:) = [4 0.47852 1.9141]; 
    MCS(9,:) = [4 0.60156 2.4063]; 
    MCS(10,:) = [6 0.45508 2.7305]; 
    MCS(11,:) = [6 0.55371 3.3223]; 
    MCS(12,:) = [6 0.65039 3.9023]; 
    MCS(13,:) = [6 0.75391 4.5234]; 
    MCS(14,:) = [6 0.85254 5.1152];   
    MCS(15,:) = [6 0.92578 5.5547]; 
    
 SINR_CQI(1,:) = [-8.0536   -7.8676   -7.6769   -7.6200   -7.5176   -7.4846]; 
 BLER_CQI(1,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(2,:) = [-6.2136   -6.0276   -5.8369   -5.7800   -5.6776   -5.6446]; 
 BLER_CQI(2,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(3,:) = [-4.0836   -3.8976   -3.7069   -3.6500   -3.5476   -3.5146]; 
 BLER_CQI(3,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(4,:) = [-2.0436   -1.8576   -1.6669   -1.6100   -1.5076   -1.4746]; 
 BLER_CQI(4,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(5,:) = [-0.0436    0.1424    0.3331    0.3900    0.4924    0.5254]; 
 BLER_CQI(5,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(6,:) = [1.8564    2.0424    2.2331    2.2900    2.3924    2.4254]; 
 BLER_CQI(6,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(7,:) = [3.8584    4.0307    4.2364    4.2876    4.3972    4.4286]; 
 BLER_CQI(7,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(8,:) = [5.6884    5.8607    6.0664    6.1176    6.2272    6.2586]; 
 BLER_CQI(8,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(9,:) = [7.6584    7.8307    8.0364    8.0876    8.1972    8.2286]; 
 BLER_CQI(9,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(10,:) = [9.5424    9.7133    9.9361   10.0154   10.1191   10.1425]; 
 BLER_CQI(10,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(11,:) = [11.4424   11.6133   11.8361   11.9154   12.0191   12.0425]; 
 BLER_CQI(11,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(12,:) = [13.3924   13.5633   13.7861   13.8654   13.9691   13.9925]; 
 BLER_CQI(12,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(13,:) = [15.2424   15.4133   15.6361   15.7154   15.8191   15.8425]; 
 BLER_CQI(13,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(14,:) = [17.1724   17.3433   17.5661   17.6454   17.7491   17.7725]; 
 BLER_CQI(14,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
 SINR_CQI(15,:) = [19.0624   19.2333   19.4561   19.5354   19.6391   19.6625]; 
 BLER_CQI(15,:) = [0.9       0.5       0.1       0.05      0.01      0.001]; 
    
 TargetBLER = 0.01;
    MCSLevelNum =  size(BLER_CQI,1); 
    TargetSINR = zeros(MCSLevelNum,1); 
    for iMCS = MCSLevelNum : -1 : 1 %??10%???????? 
        TargetSINR(iMCS)   = interp1(BLER_CQI(iMCS,:),SINR_CQI(iMCS,:),TargetBLER,'linear','extrap'); 
    end 
    tLinkLevelCurve.SINR_dB       = single(SINR_CQI); 
    tLinkLevelCurve.BLER          = single(BLER_CQI); 
    tLinkLevelCurve.TargetSINR_dB = single(TargetSINR); 
    tLinkLevelCurve.TargetBLER    = single(TargetBLER); 
    
    tLinkLevelCurve.MCS           = single(MCS); 
    
    tLinkLevelCurve.R{1}              = MCS(1:6,2);  %1-6 CQI????????? 
    tLinkLevelCurve.ModuSINR{1}       = SINR_CQI(1:6,:); 
    tLinkLevelCurve.ModuBLER{1}       = BLER_CQI(1:6,:); 
    tLinkLevelCurve.R{2}              = MCS(7:9,2);  %1-6 CQI????????? 
    tLinkLevelCurve.ModuSINR{2}       = SINR_CQI(7:9,:); 
    tLinkLevelCurve.ModuBLER{2}       = BLER_CQI(7:9,:); 
    tLinkLevelCurve.R{3}              = MCS(10:15,2);  %1-6 CQI????????? 
    tLinkLevelCurve.ModuSINR{3}       = SINR_CQI(10:15,:); 
    tLinkLevelCurve.ModuBLER{3}       = BLER_CQI(10:15,:); 
    
    tLinkLevelCurve.ModuIndex         = [1 7 10 16]; 
