function C= SinrC(sinr)
%  容量曲线
ALL = [-6.75	-5.5	-5	-4.4	-3	-2	-1	0	0.25	1.25....
    2.5	3.25	4.5	4.9	5.25	6.5	8.1	9.2	10.5	11.25	11.75....
    12.25	12.75	13	14.25	15.8	16.75	17.3	18.6	18.8;...
    0.15	0.21	0.25	0.28	0.4	0.5	0.65	0.8	0.85	1 ....
    1.2	1.35	1.5	1.6	1.7	2	2.4	2.65	3	3.2	3.35	3.45...
    3.55	3.6	4	4.5	4.8	5	5.35	5.4];
Sinr = ALL(1,:); % 信噪比
Sbe =  ALL(2,:); % 谱效率

% plot(SNR,SE,'--ro');
% grid on;

LayerNum = max(size(sinr));

C = 0;
for n = 1:LayerNum
SINR = 10*log10(sinr(n,1));
if SINR < -6.75
    C = C+ 0;
elseif SINR > 18.8
    C = C+ 5.5;
else
C = C+interp1(Sinr,Sbe,SINR);
end
end