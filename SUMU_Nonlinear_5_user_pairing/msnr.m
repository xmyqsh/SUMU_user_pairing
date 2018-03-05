function snr=msnr(hs,ebno)
% ����MMSE���ջ���SNR
% hsΪ�źŵĵ�Ч�ŵ�,snrΪ�����
% ��ʧһ����,�����źŵĹ���Ϊ1,�����Ĺ���Ϊ1/snr(�����ű�)

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