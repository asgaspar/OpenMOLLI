function [rf, gz, gzr] = make_MBVERSE_Pulse(flip,varargin)
%makeSincPulse Create a slice selective since pulse.


validPulseUses = {'excitation','refocusing','inversion'};
validDirection = {'x','y','z'};

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'make_MB_Pulse';
    
    % RF params
    addRequired(parser, 'flipAngle', @isnumeric);
    addOptional(parser, 'system', mr.opts(), @isstruct);
    addParamValue(parser, 'duration', 0, @isnumeric);
    addParamValue(parser, 'freqOffset', 0, @isnumeric);
    addParamValue(parser, 'phaseOffset', 0, @isnumeric);
    addParamValue(parser, 'timeBwProduct', 4, @isnumeric);
    addParamValue(parser, 'apodization', 0, @isnumeric);
    addParamValue(parser, 'centerpos', 0.5, @isnumeric);
    addParamValue(parser, 'MBfactor', 2, @isnumeric);
    addParamValue(parser, 'B1max', 10e-3, @isnumeric);
    addParamValue(parser, 'SliceGap', 4, @isnumeric); %<--- Slice gap in units of slice-thickness [dimensionless].
    % Slice params
    addParamValue(parser, 'maxGrad', 0, @isnumeric);
    addParamValue(parser, 'maxSlew', 0, @isnumeric);
    addParamValue(parser, 'sliceThickness', 0, @isnumeric);
    addParamValue(parser, 'delay', 0, @isnumeric);
    addParamValue(parser, 'gzrduration', 0, @isnumeric); %% added AG 11/03/2021
     % whether it is a refocusing pulse (for k-space calculation)
    addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));
end
parse(parser, flip, varargin{:});
opt = parser.Results;


BW = opt.timeBwProduct/opt.duration;
alpha = opt.apodization;
N = round(opt.duration/1e-6);
t = (1:N)*opt.system.rfRasterTime;
tt = t - opt.duration*opt.centerpos;
window = (1.0-alpha+alpha*cos(2*pi*tt/opt.duration)); % Hamming window for apodization 
signal = window.*sinc(BW*tt);
flip = sum(signal)*opt.system.rfRasterTime*2*pi; % normalize
signal = signal*opt.flipAngle/flip;

rf.type = 'rf';
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;
rf.deadTime = opt.system.rfDeadTime;
rf.ringdownTime = opt.system.rfRingdownTime;
rf.delay = opt.delay;


% Design Phase-optimized MB pulse
mb = opt.MBfactor;  
bs = opt.SliceGap;                   % Slice gap in units of slice-thickness
AM_only = 0; 

% Scale input single-band RF to mT:
verse_singleband = 1;
dt_os = 1;

maxgslew = floor(opt.system.maxSlew/(42.57*1e6)*0.999)*1e3;        % Maximum gradient amplitude [mT/m/s]
maxg = (opt.system.maxGrad/(opt.system.gamma*1e-3));               % Maximum gradient amplitude [mT/m]
b1max =  opt.B1max;         % Peak B1 [mT]


gamma_mT = 2*pi*4.257*1e4; %<--- Gyromagnetric ratio [rad/mT/s]
Gsel = 2*pi*BW/(gamma_mT*opt.sliceThickness);
Gz = Gsel*ones(N,1);
dt = opt.system.rfRasterTime;
girf = [42 42]*1e-6; 

rfsb_mT = signal/(42.57*1e3);

%#########################################################################
%                       VERSE optimization 
[rfvMB,gvMB,gvMB_actual]= dz_MBverse(rfsb_mT,Gz,dt,maxg,...
    maxgslew,b1max,verse_singleband,mb,bs*opt.sliceThickness,...
    dt_os,AM_only,girf);
%#########################################################################

Nvmb = length(rfvMB);
dtvmb = dt/dt_os;
tvmb = 0:dtvmb:(Nvmb-1)*dtvmb;

N = round(tvmb(end)/1e-6);
t = (1:N)*opt.system.rfRasterTime;
tt = t - opt.duration*opt.centerpos;

rf.t = [tvmb(1:dt_os:end)+opt.system.rfRasterTime ];
rf.signal = rfvMB(1:dt_os:end)'*(42.57*1e3); % rfmb'; % [Hz];

%Correc predited duration (correct error in IDEA)
numadd =  round((ceil(rf.t(end)*1e5)-rf.t(end)*1e5)*10); 
tFill = (1:round(numadd))*1e-6;  % Round to microsecond
rf.t = [rf.t rf.t(end)+tFill];
rf.signal = [rf.signal, zeros(size(tFill))];

if ~isempty(opt.use)
    rf.use=opt.use;
end
if rf.deadTime > rf.delay
    rf.delay = rf.deadTime;
end

if nargout > 1
    assert(opt.sliceThickness > 0,'SliceThickness must be provided');
    if opt.maxGrad > 0
        opt.system.maxGrad = opt.maxGrad;
    end
    if opt.maxSlew > 0
        opt.system.maxSlew = opt.maxSlew;
    end
    
    % gz waveform conversion from mT/m to Hz/m
    waveform =(42.57*1e3)*gvMB(:,3)' ; % [Hz/m]
    
    %Downsample acording with raster time (in the function make gradient 
    %     the gradient time is performed according with the sys.gradRasterTime )
    downsampling_factor = opt.system.gradRasterTime/opt.system.rfRasterTime;
    gwaveformdown = round([downsample(waveform,round(downsampling_factor*dt_os))]*1e2)/1e2;  

    area = sum(gwaveformdown)*opt.system.gradRasterTime; %amplitude*opt.duration;
    area_round= round(area*1e6)/1e6; 
    g_toround =  area_round - area;
    gwaveformdown(end) = gwaveformdown(end)+g_toround/opt.system.gradRasterTime;
    gz = mr.makeArbitraryGrad('z', opt.system, 'waveform',gwaveformdown);

    
    if opt.gzrduration ~= 0
        gzr= mr.makeTrapezoid('z' , opt.system, 'Area', -area*(1-opt.centerpos)-0.5*(area-area), ...
            'Duration', opt.gzrduration);
    else
        gzr= mr.makeTrapezoid('z' , opt.system, 'area', -area*(1-opt.centerpos)-0.5*(sum(gz.waveform)-area));
    end

    if rf.delay < (gz.delay)
        rf.delay = gz.delay; % these are on the grad raster already which is coarser 
    end
end


if rf.ringdownTime > 0
    tFill = (1:round(rf.ringdownTime/1e-6))*1e-6;  % Round to microsecond
    rf.t = [rf.t rf.t(end)+tFill];
    rf.signal = [rf.signal, zeros(size(tFill))];
end
    

function y = sinc(x)
    % sinc Calculate the sinc function:
    %   sinc(x) = sin(pi*x)/(pi*x)
    
    i = find(x == 0);                                                              
    x(i) = 1;
    y = sin(pi*x)./(pi*x);                                                     
    y(i) = 1;   
end

end