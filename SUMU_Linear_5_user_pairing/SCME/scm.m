function [H, delays, full_output]=scm(scmpar,linkpar,antpar,initvalues)
%SCM 3GPP Spatial Channel Model (extended)
%   H=SCM(SCMPAR,LINKPAR,ANTPAR) is a 5D-array of channel coefficients. For
%   explanation of the input parameter structs, see SCMPARSET, LINKPARSET,
%   and ANTPARSET. SIZE(H)=[U S N T K], where U is the number of MS (RX) 
%   elements, S is the number of BS (TX) elements, N is the number of paths, 
%   T is the number of time samples, and K is the number of links. If K=1,
%   the final dimension will be dropped, i.e. H is a 4D-array.
%
%   [H DELAYS]=SCM(...) outputs also a [KxN] matrix of path delays (in seconds). 
%
%   [H DELAYS BULKPAR]=SCM(...) outputs also the struct BULKPAR, whose fields
%   are as follows:
%
%   When scmpar.ScmOptions is 'none' or 'urban_canyon':
%
%   delays          - path delays in seconds [KxN]
%   path_powers     - relative path powers [KxN]
%   aods            - angles of departure in degrees over (-180,180) [KxNxM]
%   aoas            - angles of arrival in degrees over (-180,180) [KxNxM]
%   subpath_phases  - final phases for subpaths in degrees over (0,360) [KxNxM]
%   path_losses     - path losses in linear scale [Kx1]
%   shadow_fading   - shadow fading losses in linear scale [Kx1]
%   delta_t         - time sampling intervals for all links [Kx1]
%
%   In addition, when scmpar.ScmOptions is 'los' (in addition to the above):
%
%   K_factors       - K factors for all links [Kx1]
%   Phi_LOS         - final phases for LOS paths in degrees over (-180,180) [Kx1]
%
%   In addition, when scmpar.ScmOptions='polarized':
%
%   subpath_phases  - final phases for subpaths in degrees over (0,360)
%                     [Kx4xNxM], where the second dimension are the [VV VH HV HH]
%                     components (iid).
%   xpd             - cross-polarization ratios in linear scale [Kx2xN],
%                     where the (:,1,:)th dimension is the V-to-H power coupling, 
%                     and (:,2,:)th dimension is the H-to-V power coupling.
%
%   In addition, when scmpar.IntraClusterDsUsed='yes', the following fields
%   are included [L=LENGTH(NumSubPathsPerMidpath)]:
%   NumSubPathsPerMidpath   - sum over this vector is M  [1xL]
%   MidPathOrder            - order in which midpaths are summed [1xM]
%   MidPathPowers           - powers of midpaths relative to path power [1xL]
%   MidPathDelays           - delays of midpaths relative to path delay [1xL]
%   path_powers_all         - powers of all midpaths for all links [Kx(L*N)]
%
%   When scmpar.DriftShadowFading='yes', 
%   sf_init                 - final values of shadow fading [Kx1] (can be 
%                             used as initial values in recursive use of SCM)
%   
%   When scmpar.DriftDelaysAngles='yes', the delays vector is a KxNxT
%   array. In this case the delays are given wrt first time instant, 
%   i.e., they can be negative.
%
%   Also delays has size [Kx(L*N)]. The number of midpaths per path and the
%   delays and powers of midpaths are constant for all links but may vary in
%   different scenarios. 
%
%   [H ...]=SCM(...,INIT_VALUES) uses initial values given in the struct
%   INIT_VALUES, instead of random parameter generation. INIT_VALUES has
%   the same format as BULKPAR, except that SUBPATH_PHASES are now the
%   initial phases. Also, time sampling intervals (delta_t) are not used
%   (they are recalculated for every call of SCM).
%
%   The 'far scatterer clusters' option [1, Sec. 5.5.2] is not currently
%   supported. The SCM options are mutually exclusive, i.e. one cannot, for
%   instance, choose 'polarized' and 'los' simultaneously.
%
%   Example:
%       % to generate matrices for 10 links with default parameters
%       H=scm(scmparset,linkparset(10),antparset);
%       % to generate matrices for 'urban_macro' scenario
%       scmpar=scmparset;scmpar.Scenario='urban_macro';
%       H=scm(scmpar,linkparset(10),antparset);
%
%   Ref. [1]: 3GPP TR 25.996 v6.1.0 (2003-09)
%        [2]: D. Baum et al, IEEE VTC'05 (spring)
%
%   See also SCMPARSET, LINKPARSET, ANTPARSET

%   Authors: Jari Salo (HUT), Giovanni Del Galdo (TUI), Pekka Ky鰏ti (EBIT), 
%   Daniela Laselva (EBIT), Marko Milojevic (TUI), Christian Schneider (TUI)
%   $Revision: 0.33$  $Date: May 20, 2005$


% Note: all units are in degrees, meters, Hertz (1/s) and meters/second (m/s)

global Frame_time


ni=nargin;
if (ni<3 || ni>4)
    error('SCM requires three or four input arguments !')
end



if isfield(linkpar,'BsNumber')==true && isfield(linkpar,'MsNumber')==false
    error('The parameter linkpar.BsNumber has been renamed as linkpar.MsNumber!')
end


% SCM parameters, common to all links
Scenario=scmpar.Scenario;
SampleDensity=scmpar.SampleDensity;
NumTimeSamples=scmpar.NumTimeSamples;
N=scmpar.NumPaths;
M=scmpar.NumSubPathsPerPath;
CenterFrequency=scmpar.CenterFrequency;
ScmOptions=scmpar.ScmOptions;
DelaySamplingInterval=scmpar.DelaySamplingInterval;
PathLossModel=scmpar.PathLossModel;
RandomSeed=scmpar.RandomSeed;
UniformTimeSampling=scmpar.UniformTimeSampling;
PathLossModelUsed=scmpar.PathLossModelUsed;
ShadowingModelUsed=scmpar.ShadowingModelUsed;
AnsiC_core=scmpar.AnsiC_core;
LookUpTable=scmpar.LookUpTable;
FixedPdpUsed = scmpar.FixedPdpUsed;
IntraClusterDsUsed = scmpar.IntraClusterDsUsed;
FixedAnglesUsed = scmpar.FixedAnglesUsed;
DriftShadowFading=scmpar.DriftShadowFading;

% antenna parameters 
BsGainPattern=antpar.BsGainPattern;
BsGainAnglesAz=antpar.BsGainAnglesAz;
BsElementPosition=antpar.BsElementPosition;
MsGainPattern=antpar.MsGainPattern;
MsGainAnglesAz=antpar.MsGainAnglesAz;
MsElementPosition=antpar.MsElementPosition;
InterpFunction=antpar.InterpFunction;
InterpMethod=antpar.InterpMethod;   

% link parameters
MsBsDistance=linkpar.MsBsDistance;
ThetaBs=linkpar.ThetaBs;
ThetaMs=linkpar.ThetaMs;
OmegaMs=linkpar.OmegaMs; 
MsVelocity=linkpar.MsVelocity; 
MsDirection=linkpar.MsDirection;
MsHeight=linkpar.MsHeight; 
BsHeight=linkpar.BsHeight; 
MsNumber=linkpar.MsNumber;
DriftDelaysAngles=scmpar.DriftDelaysAngles;



% check that the scenario is a valid string检查场景是不是有效的字符串
if(any(strcmpi(Scenario,{'suburban_macro','urban_macro','urban_micro'}))==0) 
    error('scmpar.Scenario must be ''suburban_macro'', ''urban_macro'', or ''urban_micro''')
end

% check that the ScmOptions is a valid string检查SCM选项是不是有效的字符串
if(any(strcmpi(ScmOptions,{'none','polarized','los','urban_canyon'}))==0) 
    error('scmpar.ScmOptions must be ''none'', ''polarized'', ''los'', or ''urban_canyon'' ')
end



% check that SCM options comply with the selected scenario
if (strcmpi(ScmOptions,'urban_canyon')==1 && strcmpi(Scenario,'suburban_macro')==1 )
    scmpar.Scenario='urban_macro';
    warning('MATLAB:UrbanCanyonWrongScenario','Urban canyon option cannot be selected with "suburban_macro" -> scenario changed to "urban_macro"')
end


% check that delay/AoA drifting is not used with ANSI-C core
if (strcmpi(AnsiC_core,'yes')==1 && strcmpi(DriftDelaysAngles,'yes')==1 )
    scmpar.AnsiC_core='no';
    warning('MATLAB:NotSupportedByAnsiC_core','This version does not support drifting of delays/AoAs with ANSI-C core -> changed to scmpar.AnsiC_core=''no''.')
end



% extract the number of links
NumLinks=length(MsBsDistance);

% Check that the struct linkpar has the same number of parameters in
% each of its fields. This is also the number of links/users.
if (    NumLinks ~= length(ThetaBs)     ||... 
        NumLinks ~= length(ThetaMs)     ||...
        NumLinks ~= length(OmegaMs)     ||...
        NumLinks ~= length(MsVelocity)  ||...
        NumLinks ~= length(MsDirection) ||...
        NumLinks ~= length(MsHeight)    ||...
        NumLinks ~= length(BsHeight)    ||...
        NumLinks ~= length(MsNumber))
    error('All fields in input struct LINKPAR must be of same size!')
end



% Set random seeds if given 
if (isempty(RandomSeed)==0)
    rand('state',RandomSeed);
    randn('state',RandomSeed);
end



% determine the size of the MIMO system
% S - number of BS array antenna elements
if (numel(BsGainPattern)==1)
    S=scmpar.NumBsElements;
else
    S=size(BsGainPattern,1);
end

% U - number of MS array antenna elements
if (numel(MsGainPattern)==1)
    U=scmpar.NumMsElements;
else
    U=size(MsGainPattern,1);
end

% check that element displacement（位移） vector is of right size
if (length(BsElementPosition)~=S && length(BsElementPosition)~=1)
    error('antpar.BsElementPosition has wrong size!')
end

if (length(MsElementPosition)~=U && length(MsElementPosition)~=1)
    error('antpar.MsElementPosition has wrong size!')
end


% check that LUT size is a power-of-two没有用到
% this check is now also in scm_mex_core.c
if (strcmpi(AnsiC_core,'yes')==1)
    if (LookUpTable>0)
        if (2^nextpow2(LookUpTable)-LookUpTable~=0)
            scmpar.LookUpTable=2^nextpow2(LookUpTable);
            warning('MATLAB:LUTSizeChanged',['scmpar.LookUpTable is not a power-of-2: size changed to ' num2str(scmpar.LookUpTable) '.'])
        end
    end
end


% if fixed PDP is used, NumPaths must be six没有用到
if (strcmpi(FixedPdpUsed,'yes')==1 && N~=6)
    scmpar.NumPaths=6; N=6;
    warning('MATLAB:NumPathsChangedPdp',['Using fixed PDP, scmpar.NumPaths changed to ' num2str(scmpar.NumPaths) '.'])
elseif (strcmpi(FixedAnglesUsed,'yes')==1 && N~=6)  % if fixed AoD/AoAs are used, NumPaths must be six
    scmpar.NumPaths=6; N=6;
    warning('MATLAB:NumPathsChangedAoa',['Using fixed AoD/AoAs, scmpar.NumPaths changed to ' num2str(scmpar.NumPaths) '.'])
end



% GENERATION OF RANDOM "BULK" PARAMETERS FOR ALL LINKS
switch (ni)
        
    case (3)    % do the basic thing
        
        % check that M=20
        if (M ~= 20)
            scmpar.NumSubPathsPerPath=20; M=20;
            warning('MATLAB:NumSubPathsChanged','NumSubPathsPerPath is not 20! Using NumSubPathsPerPath=20 instead.')
        end
        
    
        % generate bulk parameters for all links
        bulkpar=generate_bulk_par(scmpar,linkpar,antpar);
        

        
        % for interpolation
        aods=bulkpar.aods;
        aoas=bulkpar.aoas;
                
        
    case (4)    % do not generate random link parameters, use initial values
        
        % take bulk parameters from input struct
        bulkpar=initvalues;
        
        % for interpolation
        aods=bulkpar.aods;
        aoas=bulkpar.aoas;
        
        
end     



% ANTENNA FIELD PATTERN INTERPOLATION
% Interpolation is computationally intensive（运算密集）, so avoid it if possible.
% Since SCM does not support elevation, dismiss the elevation dimension (for now)
% NOTE: aods/aoas should be given in degrees.
BsGainIsScalar=0;
MsGainIsScalar=0;
if numel(BsGainPattern)>1
    if (strcmp(ScmOptions,'polarized')==1)
        BsGainPatternInterpolated = zeros([2 S size(aods)]); % [polarizations(2) elements links N(6) M(20)]
        BsGainPatternInterpolated(1,:,:,:,:)=feval(InterpFunction,squeeze(BsGainPattern(:,1,1,:)),BsGainAnglesAz,aods, InterpMethod); % V
        BsGainPatternInterpolated(2,:,:,:,:)=feval(InterpFunction,squeeze(BsGainPattern(:,2,1,:)),BsGainAnglesAz,aods, InterpMethod); % H
        BsGainPatternInterpolated=permute(BsGainPatternInterpolated,[3 2 1 4 5]); % [link rx_element polarization path subpath]
    else
        BsGainPatternInterpolated=feval(InterpFunction,squeeze(BsGainPattern(:,1,1,:)),BsGainAnglesAz,aods, InterpMethod); % V only
        BsGainPatternInterpolated=permute(BsGainPatternInterpolated,[2 1 3 4]);
    end
else    % if BsGainPattern is scalar
    if (strcmp(ScmOptions,'polarized')==1)
        BsGainPatternInterpolated=repmat(BsGainPattern, [NumLinks S 2 N M]);    % [link rx_element polarization path subpath]
        BsGainIsScalar=1;        
    else
        BsGainPatternInterpolated=repmat(BsGainPattern, [NumLinks S N M]);
        BsGainIsScalar=1;
    end
end

if numel(MsGainPattern)>1
    if (strcmp(ScmOptions,'polarized')==1)
        MsGainPatternInterpolated=zeros([2 U size(aoas)]);% [polarizations(2) elements links N(6) M(20)]
        MsGainPatternInterpolated(1,:,:,:,:)=feval(InterpFunction,squeeze(MsGainPattern(:,1,1,:)),MsGainAnglesAz,aoas, InterpMethod); % V
        MsGainPatternInterpolated(2,:,:,:,:)=feval(InterpFunction,squeeze(MsGainPattern(:,2,1,:)),MsGainAnglesAz,aoas, InterpMethod); % H
        MsGainPatternInterpolated=permute(MsGainPatternInterpolated,[3 2 1 4 5]); % [link Ms_element polarization path subpath]
    else
        MsGainPatternInterpolated=feval(InterpFunction,squeeze(MsGainPattern(:,1,1,:)),MsGainAnglesAz,aoas, InterpMethod); % V only
        MsGainPatternInterpolated=permute(MsGainPatternInterpolated,[2 1 3 4]);
    end
else    % if MsGainPattern is scalar
    if (strcmp(ScmOptions,'polarized')==1)
        MsGainPatternInterpolated=repmat(MsGainPattern, [NumLinks U 2 N M]);    % [link rx_element polarization path subpath]
        MsGainIsScalar=1;
    else
        MsGainPatternInterpolated=repmat(MsGainPattern, [NumLinks U N M]);
        MsGainIsScalar=1;
    end
end

% Note: The gain patterns at this point have size(MsGainPatternInterpolated) = [link rx_element path subpath]
%  OR
% size(MsGainPatternInterpolated) = [link rx_element polarization path subpath]
% (the same for BsGainPatternInterpolated)

% Do antenna field pattern interpolation for the LOS path
if (strcmp(ScmOptions,'los')==1)
    if numel(BsGainPattern)>1
        BsGain_Theta_BS= feval(InterpFunction,squeeze(BsGainPattern(:,1,1,:)),BsGainAnglesAz,ThetaBs(:), InterpMethod); % V only
        BsGain_Theta_BS= BsGain_Theta_BS.'; % size()= [NumLinks S]
    else
        BsGain_Theta_BS=repmat(BsGainPattern,[NumLinks S]);
    end
    
    if numel(MsGainPattern)>1
        MsGain_Theta_MS= feval(InterpFunction,squeeze(MsGainPattern(:,1,1,:)),MsGainAnglesAz,ThetaMs(:), InterpMethod); % V only 
        MsGain_Theta_MS= MsGain_Theta_MS.'; % size()= [NumLinks U]
    else
        MsGain_Theta_MS= repmat(MsGainPattern,[NumLinks U]);
    end
else 
    % Set dummy values in case LOS option is not used
    BsGain_Theta_BS=NaN;
    MsGain_Theta_MS=NaN;
    
end





% CHANNEL MATRIX GENERATION
[H delta_t FinalPhases FinalPhases_LOS drift_delays] = scm_core( scmpar,...
                                                    linkpar,...
                                                    antpar,...
                                                    bulkpar,...
                                                    BsGainPatternInterpolated,...
                                                    BsGain_Theta_BS,...             % gain of LOS path
                                                    MsGainPatternInterpolated,...
                                                    MsGain_Theta_MS,...             % gain of LOS path
                                                    Frame_time,...                            % offset time (not used typically)
                                                    BsGainIsScalar,...      
                                                    MsGainIsScalar);

% final phases
bulkpar.subpath_phases=FinalPhases;

if strcmpi(scmpar.DriftDelaysAngles,'yes')==true
    bulkpar.delays=drift_delays;
end

% time sampling grid
bulkpar.delta_t=delta_t;

[U S N T K]=size(H);

% If path loss and shadowing are to be multiplied into the output
if ( (strcmp(PathLossModelUsed,'yes')==1) || strcmp(ShadowingModelUsed,'yes')==1 )
    
    if (K==1) % only one link
        if (strcmpi(PathLossModelUsed,'yes')==1)
            H=sqrt(bulkpar.path_losses).*H;   % path loss in linear scale
        end
        
        if strcmpi(ShadowingModelUsed,'yes')==1
            if strcmpi(DriftShadowFading,'yes')==0
                H=H*sqrt(bulkpar.shadow_fading);           % shadow fading in linear scale
            else    % time-variant shadow fading
                Hmat=reshape(H,U*S*N,T);  % a matrix with NumLinks cols
                if nargin>3     % initial values given
                    [sf_proc, sf_final]=sfad(scmpar,linkpar, initvalues.sf_init);
                    sf_mat=diag(sparse(sqrt(sf_proc)));    % shadow fading is in linear scale 
                    Hmat=Hmat*sf_mat;           % multiply shadow fading into each link
                    
                else
                    [sf_proc, sf_final]=sfad(scmpar,linkpar);
                    sf_mat=diag(sparse(sqrt(sf_proc)));    % shadow fading is in linear scale 
                    Hmat=Hmat*sf_mat;           % multiply shadow fading into each link
                    
                end
                
                H=reshape(Hmat,U,S,N,T);      % put back to original size
                bulkpar.sf_init=sf_final;    
            end
            
            
        end
        
    else    % if more than one link
        
        if (strcmpi(PathLossModelUsed,'yes')==1)
            Hmat=reshape(H,U*S*N*T,K);  % a matrix with NumLinks cols
            pl_mat=diag(sparse(sqrt(bulkpar.path_losses))); 
            Hmat=Hmat*pl_mat;           % multiply path loss into each link
            H=reshape(Hmat,U,S,N,T,K);      % put back to original size

        end
        
        if (strcmpi(ShadowingModelUsed,'yes')==1)
            if strcmpi(DriftShadowFading,'yes')==0
                Hmat=reshape(H,U*S*N*T,K);  % a matrix with NumLinks cols
                sf_mat=diag(sparse(sqrt(bulkpar.shadow_fading)));    % shadow fading is in linear scale 
                Hmat=Hmat*sf_mat;           % multiply shadow fading into each link
                H=reshape(Hmat,U,S,N,T,K);      % put back to original size

            else
                if nargin>3     % initial values given
                    
                    [sf_proc, sf_final]=sfad(scmpar,linkpar, initvalues.sf_init);
                    for q=1:K
                        Hmat=reshape(H(:,:,:,:,q),U*S*N,T);  % a matrix with NumLinks cols
                        sf_mat=diag(sparse(sqrt(sf_proc(q,:))));    % shadow fading is in linear scale 
                        Hmat=Hmat*sf_mat;           % multiply shadow fading into each link
                        H(:,:,:,:,q)=reshape(Hmat,U,S,N,T);      % put back to original size
                    end
                else
                    [sf_proc, sf_final]=sfad(scmpar,linkpar);
                    for q=1:K
                        Hmat=reshape(H(:,:,:,:,q),U*S*N,T);  % a matrix with NumLinks cols
                        sf_mat=diag(sparse(sqrt(sf_proc(q,:))));    % shadow fading is in linear scale 
                        Hmat=Hmat*sf_mat;           % multiply shadow fading into each link
                        H(:,:,:,:,q)=reshape(Hmat,U,S,N,T);      % put back to original size
                    end
                    
                end
                bulkpar.sf_init=sf_final;
                
            end
        end
        
        
    end
    
end
    
    
% GENERATE OUTPUT
no=nargout;
if (no>1)
    delays=bulkpar.delays;
end

if (no>2)
    switch lower(ScmOptions)
        case {'none','urban_canyon','polarized'}
            full_output=bulkpar;
            
        case ('los')
            bulkpar.Phi_LOS=FinalPhases_LOS;
            full_output=bulkpar;
            
    end
end
    






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function that maps inputs from (-inf,inf) to (-180,180)
function y=prin_value(x)
y=mod(x,360);
y=y-360*floor(y/180);



function [sf_proc, Zfinal]=sfad(scmpar,linkpar, varargin)
% SFAD Generate lognormal shadow fading process
% P=SFAD(SCMPAR,LINKPAR) returns lognormal shadow fading processes with
% standard deviation and correlation distance defined in [2]. P is a K x N
% matrix, where K is the number of links and N=scmpar.NumTimeSamples.
%
% [P FINAL]=SFAD(SCMPAR,LINKPAR,INIT) uses the initial state INIT for generating
% shadow fading. INIT is a Kx1 vector. FINAL is a Kx1 vector of final
% conditions.
%
%   Ref. [2]:   D. Baum et al, "An Interim Channel Model for Beyond-3G
%               Systems", VTC'2005 spring

%   Authors: Pekka Ky鰏ti, Jari Salo
%   $Revision= 0.11$, $Date: May 19, 2005$


K=length(linkpar.MsBsDistance);
MsVelocity = linkpar.MsVelocity;  % This is needed in UniformTimeSampling case
UniformTimeSampling = scmpar.UniformTimeSampling;
SampleDensity = scmpar.SampleDensity;
CenterFrequency = scmpar.CenterFrequency;
N = scmpar.NumTimeSamples;
Scenario=scmpar.Scenario;
ScmOptions=scmpar.ScmOptions;

SpeedOfLight = 2.99792458e8;

% Calculate channel sampling distance in meters
if strcmpi(UniformTimeSampling,'no')
    % Distance per channel sample: ds = v*1/fs = c/(2*SD*fc)
    ds = SpeedOfLight/(2*SampleDensity*CenterFrequency);
    ds=repmat(ds,K,1);
else    % Channel time sampling according to highest MSvelocity
    SD = SampleDensity*max(MsVelocity)./MsVelocity';
    % Distance per channel sample: ds = v*1/fs = c/(2*SD*fc)
    ds = SpeedOfLight./(2*SD*CenterFrequency);
end



switch lower(Scenario)
    % See [2].
    case 'urban_micro'
        CorrDist=5;     % in meters
        
        switch ScmOptions
            case {'none', 'polarized', 'urban_canyon'}
                SFstd=10;       % in dB
            case 'los'
                SFstd=4;        % in dB
        end
        
    case 'urban_macro'
        CorrDist=50;     % in meters
        switch ScmOptions
            case {'none', 'polarized', 'urban_canyon'}
                SFstd=10;       % in dB
            case 'los'
                SFstd=4;        % in dB
        end
        
    case 'suburban_macro'
        CorrDist=250;     % in meters
        switch ScmOptions
            case {'none', 'polarized', 'urban_canyon'}
                SFstd=10;       % in dB
            case 'los'
                SFstd=4;        % in dB
        end
end


% Generate iid Gaussian SF sequences (KxN)
uncorr = randn(K,N);

if length(varargin)==0      % if initial state not given
    Zinit = randn(K,1);      % IIR filter initial state
else                        % if initial state given
    Zinit=10*log10(varargin{1});  
    if ~isfinite(Zinit)
        warning('Non-finite shadow fading initialization!')
    end
end

corr = zeros(K,N);

% Get proper IIR coefficients for K links
lambda = exp(-ds./CorrDist*log(2));
forwcoeff = SFstd.*sqrt(1-lambda.^2);
backwcoeff = [ones(K,1) -lambda];

for k=1:K
    [corr(k,:),Zfinal(k)] = filter(forwcoeff(k),backwcoeff(k,:),uncorr(k,:),Zinit(k));
end
    

Zfinal=10.^(0.1*Zfinal(:));     % to lognormal 
sf_proc=10.^(0.1*corr);      % to lognormal




% %%%%%%%%%%%%%%%%
% %%%%%%%%%%%
% %%%%%%%%
% %%%%%        That's all folks !!!
% %%
% %