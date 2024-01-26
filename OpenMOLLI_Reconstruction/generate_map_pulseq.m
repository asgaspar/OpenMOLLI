function  [Param_map, varargout] = generate_map_pulseq(image_series,  Time_vector)
%

% recon_scanner 
[nl, nc, ntime] = size(squeeze(image_series));

%**************************Generate signal evolutions%*********************
recon_scanner_permuted = permute(squeeze(image_series), [3 1 2]);
recon_scanner_p_vector = reshape(recon_scanner_permuted,[ntime nl*nc] );

xdata = Time_vector;

%************************** start par for *********************************
mypool=parpool(14);
parfor npixels = 1:nl*nc
    
    ydata_noise = recon_scanner_p_vector(:, npixels);
    
    %######################## MOLLI T1 ####################################
    [T1corr(npixels,:)] = adjust_T1_3parameters_PSIR(xdata, ydata_noise,[1 2 300]);
    %######################## MOLLI T1 ####################################
end
delete(gcp('nocreate'))


A_recon =  reshape(squeeze(T1corr(:,1)),[nl nc] );
B_recon =  reshape(squeeze(T1corr(:,2)),[nl nc] );
T1star_recon =   reshape(squeeze(T1corr(:,3)),[nl nc] );

%######################## LL correction ###################################

varargout{1} = T1star_recon ;
Parameter = T1star_recon.* (B_recon./A_recon -1 );

% reshape to map
Param_map = Parameter;


end

