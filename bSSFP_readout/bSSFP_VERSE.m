function [ seq, varargout ] = bSSFP_VERSE(seq,  MB,sys, fovr,fovp, Nx,Ny, blip, ...
    li_blip_pattern, ni_acqu_pattern,B1max, z0, varargin)
% bSSFP_VERSE_4github(seq,  MB,sys, fovr,fovp, Nx,Ny, blip, ...
%     li_blip_pattern, ni_acqu_pattern,B1max, z0, varargin)
%
%##########################################################################
%                            bSSFP SMS readout
%
%
% % this bSSFP readout inlcudes:
%          - VERSE-MB
%                     (for minimizing the achievable TR and the SAR impact
%                      of SMS),
%
%          - blipped-bSSFP
%                     (to induce CAIPI-shifts without disturbing the bSSFP
%                     frequency response),
%
% Andreia S Gaspar @ ISR, IST  March 2021
% andreia.gaspar@tecnico.ulisboa.pt
%%

Nstartup = 11;


if isempty(seq)
    seq=mr.Sequence();              % Create a new sequence object
end

rfinput=0;
% Slice gap in units of slice-thickness [dimensionless].
if nargin ==13
    RF_parameter = varargin{1};
    
    % RF parameters
    alpha  = RF_parameter.alpha;
    thick  = RF_parameter.thick;
    rf_dur = RF_parameter.rf_dur; % us
    rf_apo = RF_parameter.rf_apo;
    rf_bwt = RF_parameter.rf_bwt;
    bs = RF_parameter.bs;
    
elseif nargin==11
    RF_parameter = varargin{1};
    
    % RF parameters
    alpha  = RF_parameter.alpha;
    thick  = RF_parameter.thick;
    rf_dur = RF_parameter.rf_dur; % us
    rf_apo = RF_parameter.rf_apo;
    rf_bwt = RF_parameter.rf_bwt;
    bs = RF_parameter.bs;
    
elseif nargin==14
    RF_parameter = varargin{1};
    
    % RF parameters
    alpha  = RF_parameter.alpha;
    thick  = RF_parameter.thick;
    rf_dur = RF_parameter.rf_dur; % us
    rf_apo = RF_parameter.rf_apo;
    rf_bwt = RF_parameter.rf_bwt; % was 1.5
    bs = RF_parameter.bs;
    
    rfinput=1;
    info_MB_grad_rf =  varargin{2};
    
else
    alpha=35; % deg
    inv_angle = 180; %deg
    thick=6; %mm
    rf_dur=490; % us
    rf_apo=0.5;
    rf_bwt=1.5; % was 1.5
    bs=4;
    
end


% ADC duration (controls TR/TE)
adc_dur = 2560/2; %us 2560;



% RF parameters
% alpha=35; % deg
% inv_angle = 180; %deg
% thick=6; %mm
% rf_dur=490; % us
% rf_apo=0.5;
% rf_bwt=1.5; % was 1.5


% Define other gradients and ADC events
deltak=1/fovr;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6,'system',sys);
gz_readout = mr.makeTrapezoid('z','FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6,'system',sys); % to measure the inversion
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys);


gx_parts=mr.splitGradientAt(gx,ceil((adc.dwell*adc.numSamples+adc.delay+adc.deadTime)/sys.gradRasterTime)*sys.gradRasterTime);
gx_parts(1).delay=mr.calcDuration(gxPre);
gx_1=mr.addGradients({gxPre,gx_parts(1)},'system',sys);
adc.delay=adc.delay+mr.calcDuration(gxPre);
gx_parts(2).delay=0;
gxPre.delay=mr.calcDuration(gx_parts(2));
gx_2=mr.addGradients({gx_parts(2),gxPre},'system',sys);

% Create 'alpha' degree slice selection pulse and gradient
if  rfinput==0
    [rf_VERSE, gz, gzReph] = make_MBVERSE_Pulse(alpha*pi/180,'Duration',rf_dur*1e-6,...
        'SliceThickness',thick*1e-3,'apodization',rf_apo,'timeBwProduct',rf_bwt,...
        'gzrduration', mr.calcDuration(gx_2),...
        'system',sys, 'MBfactor', MB, 'B1max', B1max, 'SliceGap', bs);
    
    %     nlocaltion = z0*1e-3; % add z0 offset
    %     rf.freqOffset = gz.amplitude*nlocaltion;
    rf_dur=(rf_VERSE.t(end)-sys.rfRingdownTime)*1e6; % us Calculate the same duration as VERSE version
    
    [rf, ~, ~] = make_MB_Pulse(alpha*pi/180,'Duration',rf_dur*1e-6,...
        'SliceThickness',thick*1e-3,'apodization',rf_apo,'timeBwProduct',rf_bwt,...
        'gzrduration', mr.calcDuration(gx_2),...
        'system',sys, 'MBfactor', MB, 'SliceGap', bs);
    
    rf_fisrt = rf;
    gz_first=gz;
else
    rf = info_MB_grad_rf.rf;
    rf_fisrt=rf;
    gz = info_MB_grad_rf.gz;
    gzReph=info_MB_grad_rf.gz_reph;
end

gamma_mT = 2*pi*4.257*1e4; %<--- Gyromagnetric ratio [rad/mT/s]
% bs = 4;                   %<--- Slice gap in units of slice-thickness [dimensionless].
fshift = MB;
gz_blip_area = pi/(2*pi*bs*thick*1e-3*fshift);
gzReph_1st = gzReph;

% blip = 0;
if blip==1
    gzr_blip_plus = mr.makeTrapezoid('z' , sys, 'Area', gzReph.area + gz_blip_area , ...
        'Duration', mr.calcDuration(gzReph));
    gzr_blip_minus = mr.makeTrapezoid('z' , sys, 'Area', gzReph.area - gz_blip_area , ...
        'Duration', mr.calcDuration(gzReph));
    %      li_blip_pattern = mod(nACQ,2); % define odd or ev
else
    gzr_blip_plus = gzReph;
    gzr_blip_minus =gzReph;
    %     li_blip_pattern = 0; % define odd or even lines pattern
end

deltakp = 1/fovp;

phaseAreas = ((0:Ny)-Ny/2)*deltakp; %############################This has to be confimed I used Ny instead of Ny-1 to be able to loop i=2:Ny+1 but I think it works fine just add one element to the vector


gz.delay=mr.calcDuration(gzReph);
gz_1=mr.addGradients({gzReph,gz},'system',sys);



[rf]=mr.align('right',rf,gz);

gz_2=gzReph;




% Added to correct ky trajectory 30062021
if mr.calcDuration(gz_2) > mr.calcDuration(gxPre)
    gx_1.delay = gx_1.delay + mr.calcDuration(gz_2) - mr.calcDuration(gxPre);
    adc.delay=adc.delay+gx_1.delay;
end


% Calculate timing
gxPre.delay=0; % otherwise duration below is misreported

% Added to correct ky trajectory 30062021
pe_dur=mr.calcDuration(gx_2)+ gx_1.delay;

gz_1.delay=max(mr.calcDuration(gx_2)-rf.delay,0);
gz_1_Delay = max(mr.calcDuration(gx_2)-rf.delay,0);
rf.delay=rf.delay+gz_1.delay+3e-5;%0.49e-5;


% finish timing calculation
add_to_TR = mr.calcDuration(rf)-mr.calcDuration(gz);
if add_to_TR>0
    TR=mr.calcDuration(gz_1)+mr.calcDuration(gx_1)+add_to_TR;
else
    TR=mr.calcDuration(gz_1)+mr.calcDuration(gx_1);
end
TE=TR/2;

% undersampling pattern
Ny_aq =  length(ni_acqu_pattern);

fprintf('Acquisition window is: %g ms\n', TR*1e3*(Ny_aq+Nstartup));


for nramp = 1:Nstartup
    if  mod(nramp,2)
        rf.phaseOffset =0;
        adc.phaseOffset=0;
    else
        rf.phaseOffset = -pi;
        adc.phaseOffset= -pi;
    end
    
    if blip==1
        if nramp ==1
            gz_1 = gz_1;
        else
            gz_1 = gz_1_previous;
        end
        
        if MB >1
            % added to have correct pattern for MB>2 and with different patterns
            ki_factor=2*mod(nramp+li_blip_pattern, MB)-MB+1;
            
            gzr_blip = mr.makeTrapezoid('z' , sys, 'Area', gzReph_1st.area + gz_blip_area*ki_factor , ...
                'Duration', mr.calcDuration(gzReph_1st));
            
            gzr_blip_rewind = mr.makeTrapezoid('z' , sys, 'Area', gzReph_1st.area - gz_blip_area*ki_factor , ...
                'Duration', mr.calcDuration(gzReph_1st));
            
            [ ~, ~, gz_2] = creat_GzbSSFP_VERSE(gz,rf_fisrt, gzr_blip ,gz_1_Delay, sys);
            
            [ gz_1_previous, ~, ~] = creat_GzbSSFP_VERSE(gz,rf_fisrt, gzr_blip_rewind ,gz_1_Delay, sys);
            
        else
            blip=0; % at the moment blip is not perform in SB
            [gz_1, rf, gz_2] = creat_GzbSSFP(gz_first,rf_fisrt, gzReph_1st , sys);
            if  mod(nramp,2)
                rf.phaseOffset =0;
                adc.phaseOffset=0;
            else
                rf.phaseOffset = -pi;
                adc.phaseOffset= -pi;
            end
        end
    else
        if MB >1
            [gz_1, ~, gz_2] = creat_GzbSSFP_VERSE(gz,rf_fisrt, gzReph_1st ,gz_1_Delay, sys);
        else
            [gz_1, rf, gz_2] = creat_GzbSSFP(gz_first,rf_fisrt, gzReph_1st , sys);
            if  mod(nramp,2)
                rf.phaseOffset =0;
                adc.phaseOffset=0;
            else
                rf.phaseOffset = -pi;
                adc.phaseOffset= -pi;
            end
        end
    end
    
    rf05=rf;
    rf05.signal=nramp*(1/Nstartup)*rf.signal;
    
    gyPre_2 = mr.makeTrapezoid('y','Area',phaseAreas(1),'Duration',pe_dur,'system',sys); %  PE step 1 just to have correct gradient
    gyPre_1 = mr.makeTrapezoid('y','Area',-phaseAreas(mod(2+Ny-2,Ny)+1),'Duration',pe_dur,'system',sys); % previous PE step
    if nramp==1
        seq.addBlock(rf05,gz_1); % without gyPre_1, gx_2
        
        seq.addBlock(gx_1,gyPre_2, gz_2); % readout without measurement
    else
        seq.addBlock(rf05,gz_1, gyPre_1, gx_2);
        
        seq.addBlock(gx_1,gyPre_2, gz_2); % readout without measurement
    end
end


% Loop over phase encodes and define sequence blocks
for ni=1:length(ni_acqu_pattern) %Ny
    i = ni_acqu_pattern(ni);
    if  mod(ni+nramp,2)
        rf.phaseOffset =0;
        adc.phaseOffset=0;
    else
        rf.phaseOffset = -pi;
        adc.phaseOffset= -pi;
    end
    
    if blip==1
        if nramp ==1
            gz_1 = gz_1;
        else
            gz_1 = gz_1_previous;
        end
        
        if MB >1
            % added to have correct pattern for MB>2 and with different patterns
            ki_factor=2*mod(ni+nramp+li_blip_pattern, MB)-MB+1;
            
            gzr_blip = mr.makeTrapezoid('z' , sys, 'Area', gzReph_1st.area + gz_blip_area*ki_factor, ...
                'Duration', mr.calcDuration(gzReph_1st));
            
            gzr_blip_rewind = mr.makeTrapezoid('z' , sys, 'Area', gzReph_1st.area - gz_blip_area*ki_factor , ...
                'Duration', mr.calcDuration(gzReph_1st));
            
            [ ~, ~, gz_2] = creat_GzbSSFP_VERSE(gz,rf_fisrt, gzr_blip ,gz_1_Delay, sys);
            
            [ gz_1_previous, ~, ~] = creat_GzbSSFP_VERSE(gz,rf_fisrt, gzr_blip_rewind ,gz_1_Delay, sys);
            
            
            
        else
            blip=0;
            [gz_1, rf, gz_2] = creat_GzbSSFP(gz_first,rf_fisrt, gzReph_1st , sys);
            if  mod(ni+nramp,2)
                rf.phaseOffset =0;
                adc.phaseOffset=0;
            else
                rf.phaseOffset = -pi;
                adc.phaseOffset= -pi;
            end
        end
    else
        if MB >1
            [ gz_1, ~, gz_2] = creat_GzbSSFP_VERSE(gz,rf_fisrt,  gzReph_1st ,gz_1_Delay, sys);
        else
            [gz_1, rf, gz_2] = creat_GzbSSFP(gz_first,rf_fisrt, gzReph_1st , sys);
            if  mod(ni+nramp,2)
                rf.phaseOffset =0;
                adc.phaseOffset=0;
            else
                rf.phaseOffset = -pi;
                adc.phaseOffset= -pi;
            end
        end
    end
    
    gyPre_2 = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',pe_dur,'system',sys); % current PE step
    if i>1
        n_step = 1 + i-ni_acqu_pattern(ni-1);
        gyPre_1 = mr.makeTrapezoid('y','Area',-phaseAreas(mod(i+Ny-n_step,Ny)+1),'Duration',pe_dur,'system',sys); % previous PE step
        seq.addBlock(rf,gz_1, gyPre_1, gx_2);
    else
        gyPre_1 = mr.makeTrapezoid('y','Area',-phaseAreas(mod(i+Ny-1,Ny)+1),'Duration',pe_dur,'system',sys); % previous PE step
        seq.addBlock(rf,gz_1, gyPre_1, gx_2);
    end
    
    
    seq.addBlock(gx_1,gyPre_2, gz_2,adc); % readout
    
end
% finish the x-grad shape
gyPre_final = mr.makeTrapezoid('y','Area',-phaseAreas(i),'Duration',pe_dur,'system',sys); % current PE step

seq.addBlock(gx_2,gyPre_final); % changes added gyPre_final 16/06/2020


if nargout==2
    varargout{1} = TR;
elseif nargout==3
    info_MB_grad_rf.rf = rf_fisrt;
    info_MB_grad_rf.gz = gz_first;
    info_MB_grad_rf.gz_reph = gzReph_1st;
    varargout{2} =  info_MB_grad_rf;
end
