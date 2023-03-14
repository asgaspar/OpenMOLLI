%%                          Open-MOLLI-SMS
%                               af2

% Open-MOLLI-SMS, integrating simultaneous multi-slice (SMS) with
% open-source ProMyoT1, to obtain all slices in a fast single-shot
% auto-calibrated sequence.


% Three methods were combined to achieve fast T1-mapping with Open-MOLLI-SMS:
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


%This demo will generate Open-MOLLI-SMS for MB 2 and 3 with in-plane
%acceleration (af=2)

% Andreia S Gaspar @ ISR, IST  Oct 2021
% andreia.gaspar@tecnico.ulisboa.pt
%
% Preliminary work at: 
% AS Gaspar, NA Silva, RG Nunes. “Open-Source Myocardial T1 mapping 
% accelerated with SMS: combining an auto-calibrated blip-bSSFP readout
% with VERSE-MB pulses” Proc. of Annual Meeting ISMRM 2022, London, 2022.

%% Add paths for Pulseq and VERSE oprimization 
path4pulseq ='../pulseq';  % define pulseq path 
path4VERSEoptimization ='../Multiband-RF-master';  % define pulseq path 

addpath(genpath(path4pulseq))
addpath(genpath(path4VERSEoptimization ))

%% Define system Parameters
write_flag = 1;

% set system limits
% had to slow down ramps and increase adc_duration to avoid stimulation
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);


%%
fovr=384e-3; Nx = 256;
fovp=327e-3; Ny = 144; % fovp = Ny*2.7e-3   % Define FOV and resolution  % should be 2.7 mm Use even number
Nstartup = 11;

opts.Interpreter = 'tex';
OpenMOLLIParameters = inputdlg({'FOV_{readout} [mm]','FOV_{phase} [mm]','N_{readout}','N_{phase}' , 'Trigger Delay [ms]'},...
    'OpenMOLLIparameters', [1 20; 1 20; 1 20 ; 1 20; 1 30], {'384','327','256','144', '550'}, opts);

% seq.plot('timeRange', [6.1865 6.2]);
try
    fovr = str2num(OpenMOLLIParameters{1})*1e-3; %[m]
    fovp = str2num(OpenMOLLIParameters{2})*1e-3; %[m]
    Nx =  str2num(OpenMOLLIParameters{3});
    sms    Ny =  str2num(OpenMOLLIParameters{4});
    Tdelay_trig = str2num(OpenMOLLIParameters{5})*1e-3;
catch
    f = msgbox(sprintf('     Default Parameters used:    \n FOVr:   %3.0f mm (%3.0f) \n FOVp:  %3.0f mm (%3.0f)', fovr*1e3, Nx , fovp*1e3, Ny));
    fovr=384e-3; Nx = 256;
    fovp=327e-3; Ny = 144;
    Tdelay_trig = 550*1e-3; % [s]
    
end


%%
% ***********************************
% ADC duration (controls TR/TE)
adc_dur = 2560/2; %us 2560;

% ***********************************
% RF parameters
RF_parameter.alpha =35;
RF_parameter.rf_dur =490; %us
RF_parameter.rf_apo =0.5;
RF_parameter.rf_bwt =1.5;
RF_parameter.thick = 8; %mm


% ***********************************
% Autocalibration Parameters
autocali_TI8 = 1;
nLacs = 48; % number of lines for autocalibration signal
nLacs_FS = 48;
ACS_factor = 0;

% ***********************************
% Other parameters
blip =1;         %<--- Blip flag (0 - no blip, 1 - gradient blip)
z0 = 0;         %<--- Centre location position which is not in use


% ***********************************
% Create inversion pulse
% Inversion RF pulse
validTypeInv = 'HS1_Kellman2013_10';
rf_inv = mr.makeHyperSecPulse(inv_angle*pi/180,'Duration',10.24e-3,...
   'system',sys,  'use', 'inversion', 'TypeInv',validTypeInv);
InvDur =  rf_inv.t(end); %s


gzSpoil_INV=mr.makeTrapezoid('z','Area',-3.5e3,'Duration',9e-3,'system',sys);


% ***********************************
% Crusher between images
gzCrusher_im=mr.makeTrapezoid('z','FlatArea',-0.4e3,'FlatTime', 5e-4,'RiseTime', 5e-4,'Duration',4e-3,'system',sys);%'FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6

FlatArea_crusher = [0  0.4e3  0.32e3   0.27e3   +0.32e3  0  0.53e3  -0.32e3 ];
FlatTime_crusher = [0   5e-4  6.3e-4   6.3e-4     6.3e-4  0   9.3e-4  6.3e-4  ];
crushers_in_between = 1;

%%
for   Tdelay_trig = 0.35%:0.05:0.55 % [s]
    for MB = 3
        delay_ACS=1;
        Nstartup_acs = 11;
        bSSFP_acs_flag=1;
        B1max = 5e-3;   %<--- Peak B1 amplitude [mT]
        for zgap = 24 %16:8:60
            
            bs = zgap/RF_parameter.thick;
            
            gzSpoil_INV=mr.makeTrapezoid('z','Area',-3.5e3,'Duration',9e-3,'system',sys);
            
            
            li_blip_pattern_bssfp =0 ;
            [~, TR] = bSSFP_VERSE([], MB,sys, fovr,fovp, Nx,Ny, blip, li_blip_pattern_bssfp, 1:2:Ny,B1max);
            
            gzCrusher_im=mr.makeTrapezoid('z','FlatArea',-0.4e3,'FlatTime', 5e-4,'RiseTime', 5e-4,'Duration',4e-3,'system',sys);%'FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6
            FlatArea_crusher = [0  0.4e3  0.32e3   0.27e3   +0.32e3  0  0.53e3  -0.32e3 ];
            FlatTime_crusher = [0   5e-4  6.3e-4   6.3e-4     6.3e-4  0   9.3e-4  6.3e-4  ];
            crushers_in_between = 1;
            
            % Tdelay_trig = 550e-3; % seconds
            
            % define the trigger to play out
             trig_BeweenInversion=mr.makeTrigger('physio1','duration', Tdelay_trig, 'delay', 100e-6); %  duration after trig to acquire in diastole 500 ms
            
             fprintf('cardiac acquisition window is: %g ms\n', TR*(Ny/2+nLacs*(MB-1))*1e3);
            
            
            
            
            % undersampling pattern
            n_fs= 2.5;
            % ni_acqu_pattern_ACS = round([1:2:(Ny/2-12) (Ny/2-11):(Ny/2+12) (Ny/2+13):2:Ny]); % 1:Ny; %1:2:Ny; % round([1:2:(Ny/2-12) (Ny/2-11):(Ny/2+12) (Ny/2+13):2:Ny]); %[1:2:Ny 128]; %round([1:2:(Ny/2-12) (Ny/2-11):(Ny/2+11) (Ny/2+12):2:Ny]); %1:Ny; %round([1:2:(Ny/2-12) (Ny/2-11):(Ny/2+11) (Ny/2+12):2:Ny]);:2:Ny;
            ni_acqu_pattern_ACS = round([1:2:(Ny/2-nLacs_FS/2) (Ny/2-nLacs_FS/2+1):(Ny/2+nLacs_FS/2) (Ny/2+nLacs_FS/2+1):2:Ny]);
            ni_acqu_pattern_ACS_R2 = round(1:2:Ny);
            Nyacquired_matrix{1} = 1:Ny;
            
            
            Nyacquired_matrix{2} = 1:2:Ny;                % R=2
            Nyacquired_matrix{3} = ni_acqu_pattern_ACS;
        
            
            for ipattern = 2%2%size(Nyacquired_matrix,2)
                seq=mr.Sequence();              % Create a new sequence object
                ni_acqu_pattern = Nyacquired_matrix{ipattern};
                Ny_aq =  length(ni_acqu_pattern);
                disp(['#### Pattern ' num2str(ipattern) ' af = ', num2str(Ny/Ny_aq)])
                
                NImage_perInv = [5 3];
                TI_Vector = [0.100 0.100+0.08 ];% [s] TI1 minimum TI of 100 msec, TI increment of 80 msec,  Messroghli 2007
                TI_Vector_real = zeros(size(TI_Vector));
                
                % blip pattern
                %       0            default pattern
                %       1            blip pattern 2
                %       2            blip pattern 3 (only for MB > 2)
                blip_pattern_vector = ACS_factor*ones(1,sum( NImage_perInv));  % zero
                
                
                %##########################################################################
                %
                %                            Inversions
                %
                %##########################################################################
            
                
                for nInv=1:(length(NImage_perInv))
                    if  nInv>1
                        seq.addBlock(trig_BeweenInversion);
                        seq.addBlock(trig_BeweenInversion);
                        seq.addBlock(trig_BeweenInversion); % 2 cardiac cycles between inversions. One  triger of then is already at the end of the nACQ loop
                    end
                    
                    
                    
                    try
                        assert(all(TI_Vector(nInv)>= (TR*(Nstartup+Ny_aq/2)+mr.calcDuration(gzSpoil_INV))));
                        delayINV = TI_Vector(nInv)-TR*(Nstartup+Ny_aq /2)-mr.calcDuration(gzSpoil_INV);
                        TI_Vector_real(nInv) = TI_Vector(nInv);
                        seq.addBlock(mr.makeDelay(delayINV)); %TI TI_Vector
                        
                        if nInv>1 && TI_adjustFlag
                            TI_Vector_real(nInv) = TI_Vector_real(nInv-1) + 0.08;
                            delayINV = TI_Vector_real(nInv)-TR*(Nstartup+Ny_aq /2)-mr.calcDuration(gzSpoil_INV);
                            seq.addBlock(mr.makeDelay(delayINV)); %TI TI_Vector
                            disp(['TI changed from ' num2str(TI_Vector(nInv)*1e3) ' ms to ' num2str(TI_Vector_real(nInv-1)*1e3+80) ' ms'])
                        end
                        
                    catch
                        delayINV = 0;
                        if nInv==1
                            TI_Vector_real(nInv) = TR*(Ny_aq/2+Nstartup) +mr.calcDuration(gzSpoil_INV); % minimum TI allowed
                            TI_adjustFlag = 1;
                            
                            
                        else
                            TI_Vector_real(nInv) = TR*(Ny_aq/2+Nstartup) +mr.calcDuration(gzSpoil_INV) + 0.08;
                            delayINV = 0.08; % add 80 ms in the second inversion even if TI is limited so is different from first one
                            %             seq.addBlock(mr.makeDelay(0.08));
                        end
                        disp(['TI changed from ' num2str(TI_Vector(nInv)*1e3) ' ms to ' num2str(TI_Vector_real(nInv)*1e3) ' ms'])
                        
                    end
                    
                    trig_inv=mr.makeTrigger('physio1','duration', Tdelay_trig - delayINV- InvDur -mr.calcDuration(gzSpoil_INV) - (Nstartup+1)*TR   , 'delay', 100e-6); %  duration after trig to acquire in diastole 500 ms
                    seq.addBlock(trig_inv); % wait for the cardiac trigger
                    
                    %##########################################################################################
                    seq.addBlock(rf_inv);
                    %##########################################################################################

                    seq.addBlock(gzSpoil_INV); %spoiler after inversion
                    
                    if delayINV>0
                        seq.addBlock(mr.makeDelay(delayINV)); %TI TI_Vector
                    end
                    
                    % Loop over phase encodes and define sequence blocks
                    for nACQ=1:NImage_perInv(nInv)
                        if nACQ >1
                            if crushers_in_between
                                gzCrusher_im=mr.makeTrapezoid('z','FlatArea',FlatArea_crusher(nACQ),'FlatTime', FlatTime_crusher(nACQ),'system',sys);%'FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6
                                if  autocali_TI8==1 && nACQ==NImage_perInv(nInv)&& nInv ==length(NImage_perInv)
                                    trig_crush=mr.makeTrigger('physio1','duration', Tdelay_trig -(Nstartup+1)*TR - mr.calcDuration( gzCrusher_im) - TR*nLacs*(MB-1), 'delay', 100e-6); %  duration after trig to acquire in diastole 550 ms
                                    
                                else
                                    trig_crush=mr.makeTrigger('physio1','duration', Tdelay_trig -(Nstartup+1)*TR - mr.calcDuration( gzCrusher_im), 'delay', 100e-6); %  duration after trig to acquire in diastole 550 ms
                                end
                                seq.addBlock(trig_crush); % wait for the cardiac trigger
                                seq.addBlock(gzCrusher_im); %spoiler before image acquisition
                            else
                                gzCrusher_im=mr.makeTrapezoid('z','FlatArea',FlatArea_crusher(nACQ),'FlatTime', FlatTime_crusher(nACQ),'system',sys);%'FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6
                                if  autocali_TI8==1 && nACQ==NImage_perInv(nInv)&& nInv ==length(NImage_perInv)
                                    trig_crush=mr.makeTrigger('physio1','duration', Tdelay_trig -(Nstartup+1)*TR - mr.calcDuration( gzCrusher_im) - TR*nLacs*(MB-1), 'delay', 100e-6); %  duration after trig to acquire in diastole 550 ms
                                    
                                else
                                    trigtrig_crush=mr.makeTrigger('physio1','duration', Tdelay_trig -(Nstartup+1)*TR - mr.calcDuration( gzCrusher_im), 'delay', 100e-6); %  duration after trig to acquire in diastole 550 ms
                                end
                                seq.addBlock(trig_crush); % wait for the cardiac trigger
                                seq.addBlock(mr.makeDelay(mr.calcDuration( gzCrusher_im)))
                            end
                        end
                        
                        li_blip_pattern = blip_pattern_vector(nACQ*nInv);
                        if autocali_TI8==1 && nACQ==NImage_perInv(nInv)&& nInv ==length(NImage_perInv)
                            if ipattern==1
                                ni_acqu_pattern_ACSi =  ni_acqu_pattern; % FS
                                R_inplane=1;
                            else
                                ni_acqu_pattern_ACSi = ni_acqu_pattern_ACS_R2; % ni_acqu_pattern_ACS;
                                R_inplane=2;
                            end
                            
                            [ seq ,save_acs_values_count] = bSSFP_VERSE_ACSinterleave(seq, MB,sys, fovr,fovp, Nx,Ny, blip, li_blip_pattern, ni_acqu_pattern_ACSi,B1max,nLacs*R_inplane, RF_parameter, info_MB_grad_rf);
                          
                            if delay_ACS
                                delayACS_s = Tdelay_trig -(Nstartup+1)*TR - TR*(nLacs/2)*(MB-1); %s
                                 trig_crushACS=mr.makeTrigger('physio1','duration', Tdelay_trig -(Nstartup+1)*TR - TR*(nLacs/2)*(MB-1), 'delay', 100e-6); %  duration after trig to acquire in diastole 550 ms
                                seq.addBlock( trig_crushACS); % wait for the cardiac trigger
                            end
                            
                            if bSSFP_acs_flag
                                
                                [ seq ,save_acs_values_count_2] = bSSFP_VERSE_ACSinterleave_FS(seq, MB,sys, fovr,fovp,...
                                    Nx,nLacs_FS+2, blip, li_blip_pattern, 1:(nLacs_FS+2),B1max,nLacs_FS, Nstartup_acs,RF_parameter, info_MB_grad_rf);
                            else
                                
                                [seq ,save_acs_values_count_2] = ...
                                    FLASH_VERSE_ACS_FS(seq, MB,sys, fovr,fovp,...
                                    Nx,nLacs_FS+2, blip, li_blip_pattern, 1:(nLacs_FS+2),B1max,nLacs_FS,...
                                    Nstartup,RF_parameter, info_MB_grad_rf);
                              
                            end
                            
                        else
                            if  nACQ==1 && nInv ==1
                                [ seq, ~, info_MB_grad_rf ] = bSSFP_VERSE(seq, MB,sys, fovr,fovp, Nx,Ny, blip, li_blip_pattern, ni_acqu_pattern,B1max,z0,RF_parameter);
                            else
                                [ seq ] = bSSFP_VERSE(seq, MB,sys, fovr,fovp, Nx,Ny, blip, li_blip_pattern, ni_acqu_pattern,B1max,z0, RF_parameter,  info_MB_grad_rf);
                                
                            end
                        end
                    end
                end
                
                if write_flag
                    afactor =  Ny/Ny_aq;
                    
                    nameseq =  sprintf('OpenMOLLI_SMS_MB%1.0f_af%1.1f_NyAcq%.0f_fovr%.0f_fovp%.0f_TD%3.0f_B1%.0f_TR%01.01f_zg%2.0f_st%1.0f_nLacs%2.0fFS%2.0f_TI8i_forkacs_Delay%1.0f_Ns%02.0f_bssfp%01.0fFLASH',  ...
                        MB,afactor, Ny_aq , fovr*1e3, fovp*1e3,Tdelay_trig*1e3, B1max*1e3, TR*1e3,bs*thick, thick, nLacs,nLacs_FS, delay_ACS*delayACS_s*1e3,Nstartup_acs,bSSFP_acs_flag);
                    
                    
                    path_base = ['./' date '_nl/ACS_af2/'];
                    if delay_ACS
                        path_save = sprintf('OpenMOLLI_SMS/MB%1.0f_B1%01.0f/st%1.0f/zg%2.0f/af%1.1f/Delay_ACS', ...
                            MB,B1max*1e3,thick, bs*thick, afactor);
                    else
                        path_save = sprintf('OpenMOLLI_SMS_af2/MB%1.0f_B1%01.0f/st%1.0f/zg%2.0f/af%1.1f/NoDelay_ACS', ...
                            MB,B1max*1e3,thick, bs*thick, afactor);
                    end
                    
                    
                    try
                        mkdir([path_base path_save])
                        cd([path_base path_save])
                    catch
                        cd([path_base path_save])
                    end
                    
                    %% prepare sequence export
                    seq.setDefinition('Name',  nameseq );
                    seq.setDefinition('FOV', [fovr fovp thick*1e-3*(bs*(MB-1)+1)]);
                    seq.setDefinition('Inversion times', num2str(TI_Vector_real(:)'));
                    seq.setDefinition('Blip Pattern', num2str(blip_pattern_vector(:)'));
                    seq.setDefinition('Date', datestr(datetime));
                    seq.write([ nameseq  '.seq'])  % Write to pulseq
                    
                    % create readout pulseq sequence
                    seq_bssfp=mr.Sequence();
                    [seq_bssfp] = bSSFP_VERSE(seq_bssfp, MB,sys, fovr,fovp, Nx,Ny, blip, li_blip_pattern, ni_acqu_pattern,B1max, z0, RF_parameter);
                    seq_bssfp.setDefinition('FOV', [fovr fovp thick*1e-3*(bs*(MB-1)+1)]);
                    seq_bssfp.write([ 'bSSFP_' nameseq(17:end-8)  '.seq'])  % Write to pulseq

                    
                   
                end
                
                
            end
            
            HR = 60;
            RR = 60e3/HR;
            delay = 20/24*RR;
            ACquisition_interval = [delay-Ny_aq/2*TR*1e3 delay+Ny_aq/2*TR*1e3];
            disp(fprintf('if HR= %2.0f bpm (RR %.0f), ACquisition within RR [%3.0f %3.0f] ms', HR,RR,ACquisition_interval(1),ACquisition_interval(2)))
            if ACquisition_interval(2)>RR
                warning(['warning: Tdelay should be <' num2str(ACquisition_interval(1)) 'ms'])
            end
        end
    end
    close all
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


