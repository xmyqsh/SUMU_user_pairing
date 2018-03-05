function [Hf chan_par] = calc_channel(syspar,scmpar,linkpar,antpar,Frame_time,chan_par_backup)

global Frame_time

ni=nargin;
linkpar.MsVelocity = 0.8333;
speed_of_light=2.99792458e8;
wavelength=speed_of_light/scmpar.CenterFrequency;%2.99792458e8/2e9=0.1499
scmpar.SampleDensity=((wavelength ./ linkpar.MsVelocity.')./2*1000);%�����ܶ�=(��0.1499./0.8333.'��./2*1000)

if ni == 5
    [chan, delays, chan_par] = scm(scmpar, linkpar, antpar);
else
    [chan, delays, chan_par] = scm(scmpar, linkpar, antpar, chan_par_backup);
end

% delta_t=(wavelength./linkpar.MsVelocity.')./2/scmpar.SampleDensity
% delaysSampleNumbers=delays/delta_t;
% pause;

Fk  = ( (1:syspar.ntones) - ceil(syspar.ntones/2) ) * 180e3; 
projMatIdx=repmat(reshape(repmat(1:syspar.ntones,syspar.Nr,1),[],1),1,syspar.Nt);
chanMat=repmat(chan,syspar.ntones,1);%chanMat5��4��6ҳ
for pathIdx = 1:size(chan,3)
    projMat=exp(-1i*2*pi*delays(pathIdx)*Fk(projMatIdx));%ÿһ������λʱ��
    for tIdx=1:size(chan,4)
        chanMat(:,:,pathIdx,tIdx)=chanMat(:,:,pathIdx,tIdx).*projMat;
    end
end
Hf=sum(chanMat,3);
Hf= squeeze(Hf);
for TTI = 1:size(Hf,3)%���ڵ�HF�ĵ���ά�����������
    Hftmp(size(Hf,1)*(TTI-1)+(1:size(Hf,1)),:) = Hf(:,:,TTI);%%������������봫��ʱ϶��ʲô��ϵ
end
Hf = Hftmp;

% Hf=reshape(permute(squeeze(Hf),[1 3 2]),[],size(Hf,2));