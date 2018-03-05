% Matlab implementation of the "Geometric Mean Decomposition"
% version of Hager, December 3, 2003
% slightly modified by Yi, April 19, 2004
% Copyright 2003, University of Florida, Gainesville, Florida
% 
%A = U*S*V' is the singular value decomposition of A
%           U, V unitary, S diagonal matrix with nonnegative
%           diagonal entries in decreasing order
%  = Q*R*P' is the geometric mean decomposition of A
%           P, Q unitary, R real upper triangular with r_ii =
%           geometric mean of the positive singular values of A,
%           1 <= i <= p, p = number of positive singular values
% All singular values smaller than tol treated as zero

function [Q1, R1, P1] = gmd_zcy_streamreduce (h)
%"User order and subchannel selection for power minimization in mimo broadcast channels using BD-GMD"

TH_contribution_rate = 0.1;

[m n]=size(h);

[U,S,V]=svd(h);
T1=[eye(m,m);zeros(n-m,m)];
S=S*T1;
V=V*T1;

RI = sum(diag(S)./(sum(S)*ones(m,1)) >= TH_contribution_rate);   %   rank adaptation
streamnum = RI;

[u s v]=gmd_zcy(S(1:streamnum,1:streamnum));
Q1=U(:,1:streamnum)*u;
R1=s;
P1=V(:,1:streamnum)*v;

end