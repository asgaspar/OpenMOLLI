function [x_final]= adjust_T1_3parameters_PSIR( xdata, ydata,x0)
% The function adjust_T1_3parameters_PSIR( xdata, ydata,x0) estimates the
% T1 value for a pixel according with the model: 
% Three parameter fit
%             S_PSIR(t) = A - B exp(-t/T1*)
%
% Receives the following inputs: 
%
%   xdata - TI vector [ms]
%   ydata - MR signal evolution, works best if using PSIR reconstruction
%   x0    - Initialization of A, B and T1*
%
%Outputs: 
%   x_final -  Estimated T1* 
%
% Note: 
%       Bound parameters: A   = 0 - 10 
%                         B   = 0 - 20
%                         T1* = 0 - 6000 [ms]
%
%
% Andreia S Gaspar @ ISR, IST  Dec 2020
%                  andreia.gaspar@tecnico.ulisboa.pt
%%
% Model
fun=@(param,time) (param(1)-param(2)*exp(-time./param(3)));

%Definition of optimization parameters
optionsLSQ = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt',...
    'MaxFunctionEvaluations',1500);

% Initialization
x0 = double(x0);

try
    % Different initializations
    [x_tmp(1,:),RESNORM(1) Residual]= lsqcurvefit(fun,x0,xdata,double(ydata./abs(max(ydata(:))))', [0 0 0],[10 20 6000], optionsLSQ);
    [x_tmp(2,:),RESNORM(2) Residual]= lsqcurvefit(fun,[1 2 2500],xdata,double(ydata./abs(max(ydata(:))))', [0 0 0],[10 20 6000], optionsLSQ);
    [x_tmp(3,:),RESNORM(3) Residual]= lsqcurvefit(fun,[1 2 1000],xdata,double(ydata./abs(max(ydata(:))))', [0 0 0],[10 20 6000], optionsLSQ);
    [x_tmp(4,:),RESNORM(4) Residual]= lsqcurvefit(fun,[1 2 100],xdata,abs(ydata./abs(max(ydata(:))))', [0 0 0],[10 20 6000], optionsLSQ);
    
    ide = find(min(RESNORM)==RESNORM);
    x = x_tmp(ide(1),:);
catch
    x = [0; 0; 0] ;
end

if x(2)<x(1)  %% should be B>A
    [x]= lsqcurvefit(fun,[x(1) 2 1000],xdata,double(ydata./max(ydata(:)))', [0 x(1) 0],[10 20 6000], optionsLSQ);
end


x_final = x;