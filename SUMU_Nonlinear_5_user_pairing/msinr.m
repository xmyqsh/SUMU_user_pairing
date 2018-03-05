function sinr=msinr(hs,hi,snr)
% ����MMSE���ջ���SINR
% hsΪ�źŵĵ�Ч�ŵ�,hiΪ���ŵĵ�Ч�ŵ�,snrΪ�����
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
R=hi*hi'+(1/snr)*eye(size(hi,1));
t=diag(inv(hs'*inv(R)*hs+eye(size(hs,2))));
sinr=1./real(t)-1;