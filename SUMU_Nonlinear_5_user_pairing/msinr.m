function sinr=msinr(hs,hi,snr)
% 计算MMSE接收机的SINR
% hs为信号的等效信道,hi为干扰的等效信道,snr为信噪比
% 不失一般性,假设信号的功率为1,噪声的功率为1/snr(即噪信比)

% if size(hs,1)~=size(hi,1)
%     error('matrix dimensions must agree.');
% end
% if ~isscalar(snr)
%     error('snr must be scalar');
% end
% if ~isreal(snr)
%     error('snr must be real');
% end
R=hi*hi'+(1/snr)*eye(size(hi,1));
t=diag(inv(hs'*inv(R)*hs+eye(size(hs,2))));
sinr=1./real(t)-1;