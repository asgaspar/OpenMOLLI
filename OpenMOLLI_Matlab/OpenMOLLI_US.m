%%#########################################################################
%
%                         Open-MOLLI
%
%                             US
%                       (af1 af2 af1.7)
%
%%#########################################################################
% Open-MOLLI: Open-source Inversion recovery myocardial T1
%           mapping sequence for fast prototyping
%           Andreia S Gaspar, Nuno A. da Silva, Rita G Nunes
%
%
% Andreia S Gaspar @ ISR, IST  Dec 2020
% andreia.gaspar@tecnico.ulisboa.pt
%
% Parameters used as input:
%   FOVr FOVp - field-of-view
%   Nx Ny - Matrix size
%   Tdelay_trig - Trigger delay
%
%
% The acquisition pattern is defined with variable ni_acqu_pattern
% in this implementation a sequence *.seq is build for each one of three
% pattern acceleration factor (af) of 1, 2 and 1.7 (fully sampled centre)
% are created, e.g.:
%       OpenMOLLI_invivo_af1.0_NyAcq144_fovr384_fovp327_TD550.seq
%       OpenMOLLI_invivo_af1.7_NyAcq84_fovr384_fovp327_TD550.seq
%       OpenMOLLI_invivo_af2.0_NyAcq72_fovr384_fovp327_TD550.seq
%
%%

write_flag = 1;


opts.Interpreter = 'tex';
OpenMOLLIParameters = inputdlg({'FOV_{readout} [mm]','FOV_{phase} [mm]','N_{readout}','N_{phase}' , 'Trigger Delay [ms]'},...
    'Open-MOLLI parameters', [1 20; 1 20; 1 20 ; 1 20; 1 30], {'384','327','256','144', '550'}, opts);


try
    fovr = str2num(OpenMOLLIParameters{1})*1e-3; %[m]
    fovp = str2num(OpenMOLLIParameters{2})*1e-3; %[m]
    Nx =  str2num(OpenMOLLIParameters{3});
    Ny =  str2num(OpenMOLLIParameters{4});
    Tdelay_trig = str2num(OpenMOLLIParameters{5})*1e-3;
catch
    f = msgbox(sprintf('     Default Parameters used:    \n FOVr:   %3.0f mm (%3.0f) \n FOVp:  %3.0f mm (%3.0f)', fovr*1e3, Nx , fovp*1e3, Ny));
    fovr=384e-3; Nx = 256;
    fovp=327e-3; Ny = 144;
    Tdelay_trig = 550*1e-3; % [s]
    
end



% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 125, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);


% ADC duration (controls TR/TE)
adc_dur = 2560/2; %us 2560;

% RF parameters
Nstartup = 11;
alpha=35; % deg
inv_angle = 180; %deg
thick=6; %mm
rf_dur=490; % us
rf_apo=0.5;
rf_bwt=1.5;


% Create 'alpha' degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rf_dur*1e-6,...
    'SliceThickness',thick*1e-3,'apodization',rf_apo,'timeBwProduct',rf_bwt,'system',sys);


% Inversion RF pulse
validTypeInv =  'HS1_Kellman2013_10';
rf_inv = mr.makeHyperSecPulse(inv_angle*pi/180,'Duration',4.75e-3,...
    'system',sys,  'use', 'inversion', 'TypeInv',validTypeInv);


InvDur =  rf_inv.t(end); %s



% Define other gradients and ADC events
deltak=1/fovr;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6,'system',sys);
gz_readout = mr.makeTrapezoid('z','FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6,'system',sys); % to measure the inversion
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys);

deltak_p = 1/fovp;
phaseAreas = ((0:Ny)-Ny/2)*deltak_p;

% gradient gz reshape
gz_parts=mr.splitGradientAt(gz,rf.delay+rf.t(end));
gz_parts(1).delay=mr.calcDuration(gzReph);
gz_1=mr.addGradients({gzReph,gz_parts(1)},'system',sys);
[rf, temp]=mr.align('right',rf,gz_1);
gz_parts(2).delay=0;
gzReph.delay=mr.calcDuration(gz_parts(2));
gz_2=mr.addGradients({gz_parts(2),gzReph},'system',sys);

% gradient gx reshape
gx_parts=mr.splitGradientAt(gx,ceil((adc.dwell*adc.numSamples+adc.delay+adc.deadTime)/sys.gradRasterTime)*sys.gradRasterTime);
gx_parts(1).delay=mr.calcDuration(gxPre);
gx_1=mr.addGradients({gxPre,gx_parts(1)},'system',sys);
adc.delay=adc.delay+mr.calcDuration(gxPre);
gx_parts(2).delay=0;
gxPre.delay=mr.calcDuration(gx_parts(2));
gx_2=mr.addGradients({gx_parts(2),gxPre},'system',sys);


gzSpoil_INV=mr.makeTrapezoid('z','Area',-3.5e3,'Duration',9e-3,'system',sys);

% Calculate timing
gxPre.delay=0; % otherwise duration below is misreported
pe_dur=mr.calcDuration(gx_2);


gz_1.delay=max(mr.calcDuration(gx_2)-rf.delay,0);
rf.delay =  rf.delay+gz_1.delay; % -0.5e-5;

% finish timing calculation
TR=mr.calcDuration(gz_1)+mr.calcDuration(gx_1);
TE=TR/2;

gzCrusher_im=mr.makeTrapezoid('z','FlatArea',-0.4e3,'FlatTime', 5e-4,'RiseTime', 5e-4,'Duration',4e-3,'system',sys);%'FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6
FlatArea_crusher = [0  0.4e3  0.32e3   0.27e3   +0.32e3  0  0.53e3  -0.32e3 ];
FlatTime_crusher = [0   5e-4  6.3e-4   6.3e-4     6.3e-4  0   9.3e-4  6.3e-4  ];
crushers_in_between = 1;


% define the trigger to play out
trig=mr.makeTrigger('physio1','duration', Tdelay_trig -(Nstartup+1)*TR,...
    'delay', 100e-6); %  duration after trig to acquire in diastole 500 ms
trig_BeweenInversion=mr.makeTrigger('physio1','duration', Tdelay_trig, ...
    'delay', 100e-6); %  duration after trig to acquire in diastole 500 ms




% undersampling pattern
n_fs= 2.5;
ni_acqu_pattern_ACS = round([1:2:(Ny/2-12) (Ny/2-11):(Ny/2+12) (Ny/2+13):2:Ny]);
Nyacquired_matrix{1} = 1:Ny;
Nyacquired_matrix{2} = 1:2:Ny;
Nyacquired_matrix{3} = ni_acqu_pattern_ACS;

for ipattern = 1:size(Nyacquired_matrix,2)
    seq=mr.Sequence();              % Create a new sequence object
    ni_acqu_pattern = Nyacquired_matrix{ipattern};
    Ny_aq =  length(ni_acqu_pattern);
    disp(['#### Pattern ' num2str(ipattern) ' af = ', num2str(Ny/Ny_aq)])
    fprintf('cardiac acquisition window is: %g ms\n', TR*Ny_aq*1e3);
    
    NImage_perInv = [5 3];
    TI_Vector = [0.100 0.100+0.08 ];% [s] TI1 minimum TI of 100 msec, TI increment of 80 msec,  Messroghli 2007
    TI_Vector_real = zeros(size(TI_Vector));
    
    
    %##########################################################################
    %
    %                             bSSFP no INV 1 image
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
                    trig_crush=mr.makeTrigger('physio1','duration', Tdelay_trig -(Nstartup+1)*TR - mr.calcDuration( gzCrusher_im), 'delay', 100e-6); %  duration after trig to acquire in diastole 550 ms
                    seq.addBlock(trig_crush); % wait for the cardiac trigger
                    seq.addBlock(gzCrusher_im); %spoiler before image acquisition
                else
                    gzCrusher_im=mr.makeTrapezoid('z','FlatArea',FlatArea_crusher(nACQ),'FlatTime', FlatTime_crusher(nACQ),'system',sys);%'FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6
                    trigtrig_crush=mr.makeTrigger('physio1','duration', Tdelay_trig -(Nstartup+1)*TR - mr.calcDuration( gzCrusher_im), 'delay', 100e-6); %  duration after trig to acquire in diastole 550 ms
                    seq.addBlock(trig_crush); % wait for the cardiac trigger
                    seq.addBlock(mr.makeDelay(mr.calcDuration( gzCrusher_im)))
                end
            end
            
            %-----------------------------------------------------------------------
            % Loop over the start-up RF pulses
            % alpha / x preparation: ramp up of 10 RF same TR with gradients PE, RO
            for nramp = 1:Nstartup
                if  mod(nramp,2)
                    rf.phaseOffset =0;
                    adc.phaseOffset=0;
                else
                    rf.phaseOffset = -pi;
                    adc.phaseOffset= -pi;
                end
                
                rf05=rf;
                rf05.signal=nramp*(1/Nstartup)*rf.signal;
                
                %  PE step 1 just to have correct gradient
                gyPre_2 = mr.makeTrapezoid('y','Area',phaseAreas(1),...
                    'Duration',pe_dur,'system',sys); 
                gyPre_1 = mr.makeTrapezoid('y','Area',...
                    -phaseAreas(mod(2+Ny-2,Ny)+1),'Duration',...
                    pe_dur,'system',sys); % previous PE step
                if nramp==1
                    rf05_1 = rf;
                    rf05_1.signal=nramp*(1/Nstartup)*rf.signal;
                    
                    seq.addBlock(rf05_1,gz_1)% without gyPre_1, gx_2
                    
                    seq.addBlock(gx_1,gyPre_2, gz_2); % readout without measurement
                else
                    seq.addBlock(rf05,gz_1, gyPre_1, gx_2);
                    
                    seq.addBlock(gx_1,gyPre_2, gz_2); % readout without measurement
                end
            end
            
            
            
            
            %-----------------------------------------------------------------------
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
                rf_phase = -pi;
                gyPre_2 = mr.makeTrapezoid('y','Area',phaseAreas(i),...
                    'Duration',pe_dur,'system',sys); % current PE step
                if i>1
                    
                    n_step = 1 + i-ni_acqu_pattern(ni-1);
                    gyPre_1 = mr.makeTrapezoid('y','Area',...
                        -phaseAreas(mod(i+Ny-n_step,Ny)+1),'Duration',...
                        pe_dur,'system',sys); % previous PE step
                    seq.addBlock(rf,gz_1, gyPre_1, gx_2);
                    
                else
                    gyPre_1 = mr.makeTrapezoid('y','Area',...
                        -phaseAreas(mod(i+Ny-1,Ny)+1),'Duration',...
                        pe_dur,'system',sys); % previous PE step
                    seq.addBlock(rf,gz_1, gyPre_1, gx_2);
                end
                
                
                
                seq.addBlock(gx_1,gyPre_2, gz_2,adc); % readout
                
                
            end
            % finish the x-grad shape
            gyPre_final = mr.makeTrapezoid('y','Area',-phaseAreas(i),...
                'Duration',pe_dur,'system',sys); % current PE step
            
            seq.addBlock(gx_2,gyPre_final);
            
            
            
            
        end
    end
    
    if write_flag
        afactor =  Ny/Ny_aq;
        nameseq =  sprintf('OpenMOLLI_invivo_af%1.1f_NyAcq%.0f_fovr%.0f_fovp%.0f_TD%3.0f', afactor, Ny_aq , fovr*1e3, fovp*1e3,Tdelay_trig*1e3);
        seq.setDefinition('Name',  nameseq );
        seq.setDefinition('FOV', [fovr fovp thick*1e-3]);
        seq.setDefinition('Inversion times', num2str(TI_Vector_real(:)'));
        seq.write([ nameseq  '.seq'])  % Write to pulseq
        if ispc
            try
                seq.write(['E:\Documents\PULSEQ\' nameseq '.seq']) % Write to pulseq
            catch
                warning('Not saved for external drive!')
                msgbox('Not saved for external drive!')
            end
        end
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

