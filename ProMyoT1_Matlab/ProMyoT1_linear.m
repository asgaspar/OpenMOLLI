%%#########################################################################
%
%                           ProMyoT1
%
%                           LINEAR
%%#########################################################################
%
% ProMyoT1: Open-source Inversion recovery myocardial T1
%           mapping sequence for fast prototyping
%           Andreia S Gaspar, Nuno A. da Silva, Rita G Nunes
%
% Andreia S Gaspar @ ISR, IST  Dec 2020
% andreia.gaspar@tecnico.ulisboa.pt
%
%
%% Sequence definition
seq=mr.Sequence();              % Create a new sequence object
fov=200e-3; Nx = 128; Ny =128;     % Define FOV and resolution
Nstartup = 11;

% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 125, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% ADC duration (controls TR/TE)
adc_dur = 2560/2; %us

% RF parameters
alpha=35; % deg
inv_angle = 180; %deg
thick=6; %mm
rf_dur=490; % us
rf_apo=0.5;
rf_bwt=1.5; % was 1.5

% define the trigger to play out
trig_out=mr.makeDigitalOutputPulse('ext1','duration', 100e-6,'delay',100e-6); % possible channels: 'osc0','osc1','ext1'
trig=mr.makeTrigger('physio1','duration', 500000e-6, 'delay', 100e-6); %  duration after trig to acquire in diastole 500 ms


% Create 'alpha' degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rf_dur*1e-6,...
    'SliceThickness',thick*1e-3,'apodization',rf_apo,'timeBwProduct',rf_bwt,'system',sys);


hearbeats = 1; % per image

% Inversion RF pulse
validTypeInv =  'HS1_Kellman2013_10';
rf_inv = mr.makeHyperSecPulse(inv_angle*pi/180,'Duration',4.75e-3,...
    'system',sys,  'use', 'inversion', 'TypeInv',validTypeInv);



lines_per_segment = round(Ny/hearbeats);
Ns=ceil(Ny/lines_per_segment);
Ny=Ns*lines_per_segment; % it can be that because of the rounding above we measure few more k-space lines...


% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6,'system',sys);
gz_readout = mr.makeTrapezoid('z','FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6,'system',sys); % to measure the inversion
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys);

phaseAreas = ((0:Ny)-Ny/2)*deltak;

% gradient gz reshape
gz_parts=mr.splitGradientAt(gz,rf.delay+rf.t(end));
gz_parts(1).delay=mr.calcDuration(gzReph);
gz_1=mr.addGradients({gzReph,gz_parts(1)},'system',sys);
[rf]=mr.align('right',rf,gz_1);
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
gxPre.delay=0;
pe_dur=mr.calcDuration(gx_2);


gz_1.delay=max(mr.calcDuration(gx_2)-rf.delay,0);
rf.delay=rf.delay+gz_1.delay;

% finish timing calculation
TR=mr.calcDuration(gz_1)+mr.calcDuration(gx_1);
TE=TR/2;

gzCrusher_im=mr.makeTrapezoid('z','FlatArea',-0.4e3,'FlatTime', 5e-4,'RiseTime', 5e-4,'Duration',4e-3,'system',sys);%'FlatArea',Nx*deltak,'FlatTime',adc_dur*1e-6


fprintf('the sequence will acquire %d lines per segment resulting in a temporal resolution of %g ms per phase\n', lines_per_segment, TR*lines_per_segment*1e3);
fprintf('cardiac acquisition window is: %g ms\n', TR*lines_per_segment*1e3);

NImage_perInv = [5 3];
TI_Vector = [0.100 0.100+0.08 ];% [s] TI1 minimum TI of 100 msec, TI increment of 80 msec,  Messroghli 2007
TI_Vector_real = zeros(size(TI_Vector));
for nInv=1:(length(NImage_perInv))
    if  nInv>1
        seq.addBlock(trig);
        seq.addBlock(trig); % 2 cardiac cycles between inversions. One  triger of then is already at the end of the nACQ loop
    end
    
    trig_inv=mr.makeTrigger('physio1','duration', 0.550-TI_Vector(nInv) , 'delay', 100e-6); %  duration after trig to acquire in diastole 500 ms
    seq.addBlock(trig_inv); % wait for the cardiac trigger
    
    
    %##########################################################################################
    seq.addBlock(rf_inv);
    %##########################################################################################
    delay_inv=ceil(TI_Vector(nInv) - TR*(Nstartup+Ny/2));
    
    seq.addBlock(gzSpoil_INV); %spoiler after inversion
    
    try
        assert(all(TI_Vector(nInv)>= (TR*(Nstartup+Ny/2)+mr.calcDuration(gzSpoil_INV))));
        delayINV = TI_Vector(nInv)-TR*(Nstartup+Ny/2)-mr.calcDuration(gzSpoil_INV);
        TI_Vector_real(nInv) = TI_Vector(nInv);
        seq.addBlock(mr.makeDelay(delayINV)); %TI TI_Vector
        
    catch
        delayINV = 0;
        if nInv==1
            TI_Vector_real(nInv) = TR*(Ny/2+Nstartup) +mr.calcDuration(gzSpoil_INV); % minimum TI allowed
        else
            TI_Vector_real(nInv) = TR*(Ny/2+Nstartup) +mr.calcDuration(gzSpoil_INV) + 0.08;
            seq.addBlock(mr.makeDelay(0.08));
        end
        disp(['TI changed from ' num2str(TI_Vector(nInv)*1e3) ' ms to ' num2str(TI_Vector_real(nInv)*1e3) ' ms'])
    end
    
    
    
    for nACQ=1:NImage_perInv(nInv)
        seq.addBlock(gzCrusher_im); %spoiler before image acquisition
        for s=1:Ns
            
            % alpha / x preparation: ramp up of 11 RF same TR with gradients PE, RO
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
            for i=1:Ny
                
                if  mod(i+nramp,2)
                    rf.phaseOffset =0;
                    adc.phaseOffset=0;
                else
                    rf.phaseOffset = -pi;
                    adc.phaseOffset= -pi;
                end
                rf_phase = -pi;
                gyPre_2 = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',pe_dur,'system',sys); % current PE step
                if i>1
                    gyPre_1 = mr.makeTrapezoid('y','Area',-phaseAreas(mod(i+Ny-2,Ny)+1),'Duration',pe_dur,'system',sys); % previous PE step
                    seq.addBlock(rf,gz_1, gyPre_1, gx_2);
                else
                    gyPre_1 = mr.makeTrapezoid('y','Area',-phaseAreas(mod(i+Ny-1,Ny)+1),'Duration',pe_dur,'system',sys); % previous PE step
                    seq.addBlock(rf,gz_1, gyPre_1, gx_2);
                end
                
                
                
                seq.addBlock(gx_1,gyPre_2, gz_2,adc); % readout
                
                
            end
            % finish the x-grad shape
            gyPre_final = mr.makeTrapezoid('y','Area',-phaseAreas(i+1),'Duration',pe_dur,'system',sys); % current PE step
            
            seq.addBlock(gx_2,gyPre_final); % changes added gyPre_final 16/06/2020
            
            
        end
        seq.addBlock(trig); % wait for the cardiac trigger
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

%% prepare sequence export
seq.setDefinition('FOV', [fov fov thick*1e-3]);
seq.setDefinition('Name', 'ProMyoT1_Linear_11RF_FOV200x200_res1.56x1.56x6');
seq.setDefinition('Inversion times', num2str(TI_Vector_real(:)'));

seq.write('ProMyoT1_Linear.seq')  % Write to pulseq file


%% plot sequence and k-space diagrams

seq.plot();


%%
% new single-function call for trajectory calculation
[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();

% plot k-spaces
time_axis=(1:(size(ktraj,2)))*sys.gradRasterTime;
figure; plot(time_axis, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points


