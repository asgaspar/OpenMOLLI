function rf = makeHyperSecPulse(flip,varargin)
%makeHyperSecPulse Create a slice selective HyperSec pulse.
%   rf=makeHyperSecPulse(flip, 'Duration', dur) Create sinc pulse
%   with given flip angle (rad) and duration (s). But these will be also be
%   defined by the TypeInv which in this case is defined as the HS1 pulse
%   described by Kellman in Kellman P, Herzka DA, Hansen MS. Adiabatic
%   inversion pulses for myocardial T1 mapping. Magn Reson Med.
%   2014;71(4):1428-1434. doi:10.1002/mrm.24793
%   And described by
%
%   rf=makeHyperSecPulse(..., 'FreqOffset', f,'PhaseOffset',p)
%   Create sinc pulse with frequency offset (Hz) and phase offset (rad).
%
%Andreia S Gaspar @ IST/ISR 2021 Adapted from makeSincPulse.m
%%
validPulseUses = {'excitation','refocusing','inversion'};
validTypeInv = {'HS1_Kellman2013_10'};

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeSincPulse';
    
    % RF params
    addRequired(parser, 'flipAngle', @isnumeric);
    addOptional(parser, 'system', mr.opts(), @isstruct);
    addParamValue(parser, 'duration', 0, @isnumeric);
    addParamValue(parser, 'freqOffset', 0, @isnumeric);
    addParamValue(parser, 'phaseOffset', 0, @isnumeric);
    addParamValue(parser, 'centerpos', 0.5, @isnumeric);
    addParamValue(parser, 'delay', 0, @isnumeric);
    % whether it is a refocusing pulse (for k-space calculation)
    addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));
    addOptional(parser, 'TypeInv', '', @(x) any(validatestring(x,validTypeInv)));
end
parse(parser, flip, varargin{:});
opt = parser.Results;


if strncmp(opt.TypeInv, validTypeInv{1}, 18) %Kellman 2013
    Beta = 3.42; %rad/s
    miu = 5.2;
    opt.duration = 0.01024; % [s]
    A0_adiab = sqrt(miu)*Beta*2/opt.duration/pi; % [Hz] adiabatic level
    A0 = 2*A0_adiab;
    N = round(opt.duration/1e-6);
    tt = linspace(-1, 1, N );
    signal0 = (sech(Beta*tt)).^(1-1i*miu); % HyperSec function
    signal_normal = signal0;
    signal = A0*signal_normal; %*1524
end

N = round(opt.duration/1e-6);
t = (1:N)*opt.system.rfRasterTime;


rf.type = 'rf';
rf.signal = signal;
rf.t = t;
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;
rf.deadTime = opt.system.rfDeadTime;
rf.ringdownTime = opt.system.rfRingdownTime;
rf.delay = opt.delay;
if ~isempty(opt.use)
    rf.use=opt.use;
end
if rf.deadTime > rf.delay
    rf.delay = rf.deadTime;
end


if rf.ringdownTime > 0
    tFill = (1:round(rf.ringdownTime/1e-6))*1e-6;  % Round to microsecond
    rf.t = [rf.t rf.t(end)+tFill];
    rf.signal = [rf.signal, zeros(size(tFill))];
end

end