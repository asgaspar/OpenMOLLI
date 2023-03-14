function [ seq,save_acs_values_count,  varargout ] = bSSFP_VERSE_ACSinterleave_FS(seq,  MB,sys, fovr,fovp, Nx,Ny, blip, ...
    li_blip_pattern, ni_acqu_pattern,B1max, nLacs,Nstartup, varargin)
% [ seq,save_acs_values_count,  varargout ] = bSSFP_SMS3_VERSE_ACSinterleave_nL_FS_noRUp_4github(seq,  MB,sys, fovr,fovp, Nx,Ny, blip, ...
%     li_blip_pattern, ni_acqu_pattern,B1max, nLacs,Nstartup, varargin)
%%                       bSSFP readout
%                           LINEAR

%##########################################################################
%                            bSSFP SMS readout
%                                   with
%                           Auto-calibration blip patterns
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
%          - Auto-calibration blip patterns
%                     (to ensure consistency between the data and coil
%                      sensitivity information).
%
% Andreia S Gaspar @ ISR, IST  March 2021
% andreia.gaspar@tecnico.ulisboa.pt
%%
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
    rf_bwt = RF_parameter.rf_bwt; % was 1.5
    bs = RF_parameter.bs;
    
elseif nargin==11
    RF_parameter = varargin{1};
    
    % RF parameters
    alpha  = RF_parameter.alpha;
    thick  = RF_parameter.thick;
    rf_dur = RF_parameter.rf_dur; % us
    rf_apo = RF_parameter.rf_apo;
    rf_bwt = RF_parameter.rf_bwt; % was 1.5
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



% Define other gradients and ADC events
deltak=1/fovr;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys);



gx_parts=mr.splitGradientAt(gx,ceil((adc.dwell*adc.numSamples+adc.delay+adc.deadTime)/sys.gradRasterTime)*sys.gradRasterTime);
gx_parts(1).delay=mr.calcDuration(gxPre);
gx_1=mr.addGradients({gxPre,gx_parts(1)},'system',sys);
adc.delay=adc.delay+mr.calcDuration(gxPre);
gx_parts(2).delay=0;
gxPre.delay=mr.calcDuration(gx_parts(2));
gx_2=mr.addGradients({gx_parts(2),gxPre},'system',sys);


if  rfinput==0
    [rf, gz, gzReph] = make_MBVERSE_Pulse(alpha*pi/180,'Duration',rf_dur*1e-6,...
        'SliceThickness',thick*1e-3,'apodization',rf_apo,'timeBwProduct',rf_bwt,...
        'gzrduration', mr.calcDuration(gx_2),...
        'system',sys, 'MBfactor', MB, 'B1max', B1max, 'SliceGap', bs);
    
    %     nlocaltion = z0*1e-3; % add z0 offset
    %     rf.freqOffset = gz.amplitude*nlocaltion;
    
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


gzr_blip_plus = mr.makeTrapezoid('z' , sys, 'Area', gzReph.area + gz_blip_area , ...
    'Duration', mr.calcDuration(gzReph));
gzr_blip_minus = mr.makeTrapezoid('z' , sys, 'Area', gzReph.area - gz_blip_area , ...
    'Duration', mr.calcDuration(gzReph));
%      li_blip_pattern = mod(nACQ,2); % define odd or ev


deltakp = 1/fovp;

phaseAreas = ((0:Ny)-Ny/2)*deltakp; %############################This has to be confimed I used Ny instead of Ny-1 to be able to loop i=2:Ny+1 but I think it works fine just add one element to the vector

gz_parts=mr.splitGradientAt(gz,rf.delay+rf.t(end));

gz.delay=mr.calcDuration(gzReph);
gz_1=mr.addGradients({gzReph,gz},'system',sys);

[rf]=mr.align('right',rf,gz);

gz_2=gzReph; %mr.addGradients({gz_parts(2),gzReph},'system',sys);




% Added to correct ky trajectory 30062021
if mr.calcDuration(gz_2) > mr.calcDuration(gxPre)
    gx_1.delay = gx_1.delay + mr.calcDuration(gz_2) - mr.calcDuration(gxPre);
    adc.delay=adc.delay+gx_1.delay;
end

gzSpoil_INV=mr.makeTrapezoid('z','Area',-3.5e3,'Duration',9e-3,'system',sys);

% Calculate timing
gxPre.delay=0; % otherwise duration below is misreported

% Added to correct ky trajectory 30062021
pe_dur=mr.calcDuration(gx_2)+ gx_1.delay;

% adjust delays to align objects
gz_1.delay=max(mr.calcDuration(gx_2)-rf.delay,0);
gz_1_Delay = max(mr.calcDuration(gx_2)-rf.delay,0);
rf.delay=rf.delay+gz_1.delay+3e-5;%0.49e-5;
% rf.delay=rf.delay+1e-5;

% finish timing calculation
add_to_TR = mr.calcDuration(rf)-mr.calcDuration(gz);
if add_to_TR>0
    TR=mr.calcDuration(gz_1)+mr.calcDuration(gx_1)+add_to_TR;
else
    TR=mr.calcDuration(gz_1)+mr.calcDuration(gx_1);
end
TE=TR/2;

% undersampling pattern
% ni_acqu_pattern = round([1:2:(Ny/2-13) (Ny/2-12):(Ny/2+12) (Ny/2+13):2:Ny]);
% ni_acqu_pattern = 1:Ny;
Ny_aq =  length(ni_acqu_pattern);

fprintf('Acquisition window is: %g ms\n', TR*1e3*(Ny_aq+Nstartup));

% alpha / x preparation: ramp up of 10 RF same TR with gradients PE, RO
% create 0.5*alpha prep pulse
% create 0.5*alpha prep pulse
if Nstartup~=0
    for nramp = 1:Nstartup
        if  mod(nramp,2)
            rf.phaseOffset =-pi;
            adc.phaseOffset=-pi;
        else
            rf.phaseOffset = 0;
            adc.phaseOffset= 0;
        end
        
        
        if nramp ==1
            gz_1 = gz_1;
        else
            gz_1 = gz_1_previous;
        end
        
        
        % added to have correct pattern for MB>2 and with different patterns
        ki_factor=2*mod(nramp+li_blip_pattern, MB)-MB+1;
        
        gzr_blip = mr.makeTrapezoid('z' , sys, 'Area', gzReph_1st.area + gz_blip_area*ki_factor , ...
            'Duration', mr.calcDuration(gzReph_1st));
        
        gzr_blip_rewind = mr.makeTrapezoid('z' , sys, 'Area', gzReph_1st.area - gz_blip_area*ki_factor , ...
            'Duration', mr.calcDuration(gzReph_1st));
        
        [ ~, ~, gz_2] = creat_GzbSSFP_VERSE(gz,rf_fisrt, gzr_blip ,gz_1_Delay, sys);
        
        [ gz_1_previous, ~, ~] = creat_GzbSSFP_VERSE(gz,rf_fisrt, gzr_blip_rewind ,gz_1_Delay, sys);
        
        
        
        
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
else
    nramp =0;
end
iMB = 0;
count_rf = 0;
count_lines = 0;
% Loop over phase encodes and define sequence blocks
for ni=1:length(ni_acqu_pattern) %Ny
    i = ni_acqu_pattern(ni);
    if  mod(ni+nramp+count_rf,2)
        rf.phaseOffset =-pi;
        adc.phaseOffset=-pi;
    else
        rf.phaseOffset = 0;
        adc.phaseOffset= 0;
    end
    
    
    gz_1 = gz_1_previous;
    
    
    
    % added to have correct pattern for MB>2 and with different patterns
    ki_factor=2*mod(ni+nramp+li_blip_pattern, MB)-MB+1;
    
    gzr_blip = mr.makeTrapezoid('z' , sys, 'Area', gzReph_1st.area + gz_blip_area*ki_factor, ...
        'Duration', mr.calcDuration(gzReph_1st));
    
    gzr_blip_rewind = mr.makeTrapezoid('z' , sys, 'Area', gzReph_1st.area - gz_blip_area*ki_factor , ...
        'Duration', mr.calcDuration(gzReph_1st));
    
    [ ~, ~, gz_2] = creat_GzbSSFP_VERSE(gz,rf_fisrt, gzr_blip ,gz_1_Delay, sys);
    
    [ gz_1_previous, ~, ~] = creat_GzbSSFP_VERSE(gz,rf_fisrt, gzr_blip_rewind ,gz_1_Delay, sys);
    
    
    
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
    
    if i==(Ny/2-nLacs/2+1)
        saveinit_ni=ni;
    end
    
    count_lines = count_lines+1;
    save_acs_values(ni, 1) = i;
    save_acs_values_count(count_lines,1) = i ;
    save_acs_values_count(count_lines,2) = sign(ki_factor);
    if  i>(Ny/2-nLacs/2) && i<(Ny/2+nLacs/2+1)
        savefi_ni=ni;
        %                  keyboard
        ni_vector_ACS = saveinit_ni:savefi_ni;
        ACSLines_vector=(Ny/2-11):(Ny/2+12);
        
        for iMB = 1:MB-1
            count_rf = count_rf+1;
            
            if  mod(ni+count_rf+nramp,2)
                rf.phaseOffset =-pi;
                adc.phaseOffset=-pi;
            else
                rf.phaseOffset = 0;
                adc.phaseOffset= 0;
            end
            
            save_acs_values(ni, iMB+1) = i;
            
            gz_1 = gz_1_previous;
            
            % added to have correct pattern for MB>2 and with different patterns
            ki_factor=2*mod(ni+nramp+li_blip_pattern+iMB, MB)-MB+1;
            
            count_lines = count_lines+1;
            save_acs_values_count(count_lines,1) = i ;
            save_acs_values_count(count_lines,2) = sign(ki_factor);
            
            gzr_blip = mr.makeTrapezoid('z' , sys, 'Area', gzReph_1st.area + gz_blip_area*ki_factor, ...
                'Duration', mr.calcDuration(gzReph_1st));
            
            gzr_blip_rewind = mr.makeTrapezoid('z' , sys, 'Area', gzReph_1st.area - gz_blip_area*ki_factor , ...
                'Duration', mr.calcDuration(gzReph_1st));
            
            [ ~, ~, gz_2] = creat_GzbSSFP_VERSE(gz,rf_fisrt, gzr_blip ,gz_1_Delay, sys);
            
            [ gz_1_previous, ~, ~] = creat_GzbSSFP_VERSE(gz,rf_fisrt, gzr_blip_rewind ,gz_1_Delay, sys);
            
            
            gyPre_2 = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',pe_dur,'system',sys); % current PE step
            
            i_next = ni_acqu_pattern(ni+1);
            n_step = 1 + i-ni_acqu_pattern(ni-1);
            gyPre_1 = mr.makeTrapezoid('y','Area',-phaseAreas(mod(i_next+Ny- n_step ,Ny)+1),'Duration',pe_dur,'system',sys); % previous PE step
            seq.addBlock(rf,gz_1, gyPre_1, gx_2);
            
            
            
            seq.addBlock(gx_1,gyPre_2, gz_2,adc);
        end
        
    end
end
% finish the x-grad shape
gyPre_final = mr.makeTrapezoid('y','Area',-phaseAreas(i),'Duration',pe_dur,'system',sys); % current PE step

seq.addBlock(gx_2,gyPre_final); % changes added gyPre_final 16/06/2020


%%


if nargout==2
    varargout{1} = TR;
elseif nargout==3
    info_MB_grad_rf.rf = rf_fisrt;
    info_MB_grad_rf.gz = gz_first;
    info_MB_grad_rf.gz_reph = gzReph_1st;
    varargout{2} =  info_MB_grad_rf;
end



