function [syspar,scmpar,antpar]=argdeliv(ctrlpar)
if ~isscalar(ctrlpar.maxlayer)
    error('');
end
if ctrlpar.maxlayer>min(ctrlpar.Nr,ctrlpar.Nt)
    error('��֧�ֵĴ������Ŀ');
end
if ctrlpar.codesize<4 || ctrlpar.codesize>9
    error('�����������ȴ���');
end
if (any(strcmpi(ctrlpar.precoding,{'PEV','SLNR','BD','RZF'}))==0)
    error('��֧�ֵ�Ԥ����ʽ');
end
syspar=ctrlpar;%����ϵͳ��������
scmpar=scmparset;%SCM�ŵ���������
antpar=antparset;%���߲�������

scmpar.NumBsElements=syspar.Nt;
scmpar.NumMsElements=syspar.Nr;
scmpar.NumTimeSamples=syspar.ntimesamples;
scmpar.Scenario=syspar.scenario;
scmpar.ScmOptions=syspar.ScmOptions;

if strcmp(syspar.antconfig,'diversity-ant')
    antpar.BsElementPosition=4;%���߼�ľ���4������
end

if strcmp(scmpar.ScmOptions,'polarized')
    d=antpar.BsElementPosition;
    antpar.BsElementPosition=reshape(repmat(0:d:d*(syspar.Nt/2-1),2,1),1,[]);
end

az=antpar.BsGainAnglesAz;
ag=-1*min(12*(az/70).^2,20) + 14 + 0;
ag=10.^(ag/20); % In antparset.m, it says  The complex field patterns are
%given in linear scale. The antenna gain is
%20*log10(abs(BsGainPattern)).��λ��������Ĺ�ϵ��������

if strcmp(scmpar.ScmOptions,'polarized')
    bs_array=zeros(syspar.Nt,2,1,length(az));
    for t=1:2:(syspar.Nt)
        bs_array(t,:,:,:)=dipole(ag,45);
        bs_array(t+1,:,:,:)=dipole(ag,-45);
    end
else
    bs_array=zeros(syspar.Nt,1,1,length(az));
    for t=1:(syspar.Nt)
        bs_array(t,:,:,:)=ag;
    end
end
 antpar.BsGainPattern=bs_array;