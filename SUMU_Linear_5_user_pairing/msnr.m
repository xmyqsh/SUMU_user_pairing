function snr=msnr(hs,ebno)
% 计算MMSE接收机的SNR
% hs为信号的等效信道,snr为信噪比
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
R = (1/ebno)*eye(max(size(hs)));
t=diag(inv(hs'*inv(R)*hs+eye(size(hs,2))));
snr=1./real(t)-1;