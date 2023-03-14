function [ seq,save_acs_values_count,  varargout ] = FLASH_VERSE_ACS_FS(seq, ...
    MB,sys, fovr,fovp, Nx,Ny, blip, ...
    li_blip_pattern, ni_acqu_pattern,B1max, nLacs,Nstartup, varargin)
%% [ seq,save_acs_values_count,  varargout ] = FLASH_SMS3_VERSE_ACSfinal_FS_4github(seq, ...
%     MB,sys, fovr,fovp, Nx,Ny, blip, ...
%     li_blip_pattern, ni_acqu_pattern,B1max, nLacs,Nstartup, varargin)
% 
%% ##########################################################################
%                         FLASH SMS readout
%                               with
%                   Auto-calibration blip patterns
%
%
% % this FLASH readout inlcudes:
%          - VERSE-MB
%                     (for minimizing the achievable TR and the SAR impact
%                      of SMS),
%
%          - blipped-bSSFP
%                     (to induce CAIPI-shifts without disturbing the bSSFP
%                     frequency response),
%
%          - Auto-calibration blip patterns
%                     (to ensure consistency between the data and coil
%                      sensitivity information).
%
% Andreia S Gaspar @ ISR, IST  March 2021
% andreia.gaspar@tecnico.ulisboa.pt
%%
% Slice gap in units of slice-thickness [dimensionless].
if nargin ==13
    RF_parameter = varargin{1};
    
    % RF parameters
  
    sliceThickness   = RF_parameter.thick*1e-3;
    rf_dur = RF_parameter.rf_dur; % us
    rf_apo = RF_parameter.rf_apo;
    rf_bwt = RF_parameter.rf_bwt;
    bs = RF_parameter.bs;
    
elseif nargin==11
    RF_parameter = varargin{1};
    
    % RF parameters
 
    sliceThickness   = RF_parameter.thick*1e-3;
    rf_dur = RF_parameter.rf_dur; % us
    rf_apo = RF_parameter.rf_apo;
    rf_bwt = RF_parameter.rf_bwt;
    bs = RF_parameter.bs;
    
elseif nargin==14
    RF_parameter = varargin{1};
    
    % RF parameters
  
    sliceThickness   = RF_parameter.thick*1e-3;
    rf_dur = RF_parameter.rf_dur; % us
    rf_apo = RF_parameter.rf_apo;
    rf_bwt = RF_parameter.rf_bwt; % was 1.5
    bs = RF_parameter.bs;
    
    rfinput=1;
    info_MB_grad_rf =  varargin{2};
    
else
    alpha=10; % deg
    inv_angle = 180; %deg
    sliceThickness =6*1e-3; %mm
    rf_dur=490; % us
    rf_apo=0.5;
    rf_bwt=1.5; % was 1.5
    bs=4;
    
end

alpha = 10;                        % flip angle
TE=3e-3;
TR=5e-3;                        % only a single value for now


% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
adc_duration =2560/2*1e-6; %us 2560; 3.2e-3;
pre_duration = 0.5e-3; %1e-3;


% Create 'alpha' degree slice selection pulse and gradient


[rf, gz, gzReph] = make_MBVERSE_Pulse(alpha*pi/180,'Duration',rf_dur*1e-6,...
    'SliceThickness',sliceThickness ,'apodization',rf_apo,'timeBwProduct',rf_bwt,...
    'gzrduration', pre_duration,...
    'system',sys, 'MBfactor', MB, 'B1max', B1max, 'SliceGap', bs);

rf_fisrt = rf;
gz_first=gz;


% Define other gradients and ADC events
deltakr=1/fovr;
deltakp=1/fovp;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltakr,'FlatTime',adc_duration,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',pre_duration,'system',sys);


phaseAreas = ((0:nLacs-1)-nLacs/2)*deltakp;

% gradient spoiling
gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltakr,'system',sys);
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys);

% Calculate timing
delayTE=ceil((TE - mr.calcDuration(gxPre) - mr.calcDuration(gz)/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR=ceil((TR - mr.calcDuration(gxPre) - mr.calcDuration(gz) ...
    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;

rf_phase=0;
rf_inc=0;

rf.delay=rf.delay;
gz.delay = rf.delay;

gamma_mT = 2*pi*4.257*1e4; %<--- Gyromagnetric ratio [rad/mT/s]
% bs = 4;                   %<--- Slice gap in units of slice-thickness [dimensionless].
fshift = MB;
gz_blip_area = pi/(2*pi*bs*sliceThickness*fshift);
gzReph_1st = gzReph;
count_lines = 0;

for l=1:nLacs
    for iMB = 1:MB
        count_lines = count_lines+1;
        %seq.addBlock(rf_fs,gz_fs); % fat-sat
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
        
        
        % added to have correct pattern for MB>2 and with different patterns
        ki_factor=2*mod(l+iMB, MB)-MB+1;
        save_acs_values_count(count_lines,1) = l ;
        save_acs_values_count(count_lines,2) = sign(ki_factor);
        gzr_blip = mr.makeTrapezoid('z' , sys, 'Area', gzReph_1st.area + gz_blip_area*ki_factor , ...
            'Duration', mr.calcDuration(gzReph_1st));
        
        
        seq.addBlock(mr.makeDelay(delayTE));
        seq.addBlock(rf,gz);
        gyPre = mr.makeTrapezoid('y','Area',phaseAreas(l),'Duration',pre_duration,'system',sys);
        seq.addBlock(gxPre,gyPre,gzr_blip);
        
        if delayTE>0
            seq.addBlock(mr.makeDelay(delayTE));
        end
        
        seq.addBlock(gx,adc); % usual readout

        
        gyPre.amplitude=-gyPre.amplitude;
        seq.addBlock(mr.makeDelay(delayTR),gxSpoil,gyPre,gzSpoil)
    end
end




%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

