%##########################################################################
%
%                   Open-MOLLI Reconstruction 
%
%                         with GRAPPA
%##########################################################################
% Recon_data_OpenMOLLI(data_file_path,nContrast, sampling,  Ny_fs)

% Andreia S Gaspar @ ISR, IST  Dec 2023
% andreia.gaspar@tecnico.ulisboa.pt

% The following script receives a file *.dat acquired in Siemens MR machine 
% with Open-MOLLI.seq sequence. 

% The default acquisition parameters with Open-MOLLI were matrix of 256x144, 
% 8 contrasts and acceleration factor of 2 with 24 fully sampled lines in
% the centre of k-space  acquisition. 
% Changes in the acquisition will require changes in the reconstruction
% parameters. 

%% ************************************************************************
%                            Select data
if ispc
    path = [cd '\'];
else
    path = [cd '/'];
end
pattern='*.dat';

D=dir([path pattern]);
[names,I]=sort([D(:).datenum]);
data_file_path=[D(I(1)).name];

%% ************************************************************************
%                           Read *.dat
twix_obj = mapVBVD(data_file_path);

twix_obj{2}.image.flagRemoveOS =0;
try
    data = twix_obj{2}.image.unsorted();
catch
    data = twix_obj{2}.image();
end

%% ************************************************************************
%              Acquisition and Reconstruction parameters 
nContrast = 8;
Ny_fs= 144; % Number of lines when fully sampled 

Ny_recon=218; % Reconstructed Ny
Nx_recon=256; % Reconstructed Nx

Nx = size(data,1); % Acquired Nx
Ny = size(data,3)./nContrast; % Acquired Ny

nCoils = size(data, 2); % Numer of coils

%#####################################################################
% sampling:
%           1 - fully sampled k-space which is undersampled in this step
%           2 - data with fully sampled k-space centre in all images
%           3 - data with only 1 image with fully sampled k-space centre
%           4 - fully sampled k-space




% Same pattern applied when generating the OpenMOLLI.seq file
ni_acqu_pattern = round([1:2:(Ny_fs/2-12) (Ny_fs/2-11):(Ny_fs/2+12) (Ny_fs/2+13):2:Ny_fs]);
fully_sampled_centre_lines = ni_acqu_pattern((diff(ni_acqu_pattern)==1)); 

data_p = permute(data, [1, 3, 2]);
data_final = reshape(data_p, [Nx Ny nContrast nCoils]);
data_final2 = zeros(Nx,Ny_fs, nContrast, nCoils);

%%#####################################################################
% sampling:
af = 2; %acceleration factor
if af==2
    data_final2(:,ni_acqu_pattern,:,:) = data_final(:,:,:,:);
    data_final2_zeros = zeros(Nx,Ny_fs, nContrast, nCoils);
    data_final2_zeros(:,1:2:Ny_fs,:,:) = data_final2(:,1:2:Ny_fs,:,:);
else
    data_final2 = data_final(:,:,:,:);
end

%% ************************************************************************
%               Reconstruction loop over each contrast
acs_respective = 1;

% Initializa variables
image_recon_coils_sets = zeros(nCoils ,Nx_recon,Ny_recon,  nContrast);
sos_grappa = zeros(Nx_recon, Ny_recon, nContrast);
sos_rotated = zeros( Nx_recon,Ny_recon,  nContrast);

for iInv=1:nContrast
    if acs_respective
        Nima_acs = iInv;
    else
        Nima_acs = 1;
    end
    
    
    acs_dim = squeeze(data_final2(:,fully_sampled_centre_lines ,Nima_acs,:));
    acs_dimP = permute(acs_dim, [3 2 1]);
    
    %######################################################################
    %                       GRAPPA reconstruction 
   
    
    if af==2
        %##################################################################
        % Kernel size and calibration data
        kSize = [5,5];
        
        kCalib = acs_dim;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %  Calculate and Apply GRAPPA weights
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        sig_dim_US_temp = squeeze(data_final2(:,:,iInv,:));
        sig_dim_US_temp(:,2:2:Ny_fs,:) = 0;
        DATA =sig_dim_US_temp;
        
        %%GRAPPA
        [res] = GRAPPA(DATA,kCalib, kSize, 0.01,1);
        
        recon_or = res;
    else
        recon_or = ifft2s(squeeze(data_final2(:,:,iInv,:)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply Tukey filter in k space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kspace_down = recon_or;
    kspace.tukey = 0.5;
    filterWindow = window3_tukeytaper([size(kspace_down,1) size(kspace_down,2) 1], @tukeywin, kspace.tukey);
    kspace_down_filt      = kspace_down .* filterWindow;
    
   
    % zeropadded for PSIR recon
    recon_zeropaded   = padarray( kspace_down_filt , [0 37 0], 0, 'pre'); % pre or post?
    recon_zp = padarray(  recon_zeropaded  , [0 37 0], 0, 'pos'); % pre or post?
    
    imagem_recon1 = zeros(size(recon_zp));
    
    for ii = 1:nCoils
        imagem_recon1(:,:, ii) = fftshift(fft2(fftshift(( recon_zp(:,:, ii)))));
        image_recon_coils_sets(ii,:,:,iInv) = imagem_recon1(:,:, ii) ;
    end
    
    sos_grappa(:,:,iInv)=abs(sum(imagem_recon1.*conj(imagem_recon1),ndims(imagem_recon1))).^(1/2);
    sos_rotated(:,:,iInv) = flip(flip(sos_grappa(:,:,iInv), 2), 1);
    
    
end

phase = angle(sum( image_recon_coils_sets, 1)/nCoils);
imspace_PSIR = sign(squeeze(phase)).*squeeze(sos_grappa); 

sos_rotated_one = cat(1, zeros(1,Ny_fs+74, nContrast),sos_rotated);
sos_rotated_two = cat(2, zeros(Nx,1, nContrast), sos_rotated_one(1:Nx, :,:));
sos_final = sos_rotated_two(1:Nx, 1:Ny_fs+74, :);


%% ************************************************************************
%                              Display                
figure('units','normalized','outerposition',[0 0 1 1])
try
    noCols = min(5,ceil(size(sos_final,3)/2));
    outIm = makemosaic(sos_final,noCols);
    
    imagesc([0 twix_obj{2}.hdr.Config.PhaseFoV*noCols],...
        [0 twix_obj{2}.hdr.Config.ReadFoV*( nContrast/noCols)],...
        double((double(outIm)./max(double(outIm(:))))), [0 1])
catch
    imab(sos_final),
end

xlabel('Phase FOV [mm]')
ylabel('Readout FOV [mm]')
colormap gray
if acs_respective
    name_figure = ['recon_GRAPPA_respective_2x'   data_file_path];
else
    name_figure = ['recon_GRAPPA_acs' num2str(Nima_acs) '_2x'   data_file_path];
end
pl = title(name_figure);
set(pl, 'interpreter', 'none')

%% ************************************************************************
%                              TI values   
%
% Time_vector is obtained from the *.seq file
seq = mr.Sequence();              % Create a new sequence object
seq.read(seq_file_path,'detectRFuse');
[ktraj_adc, ktraj, t_excitation, t_refocusing] = seq.calculateKspace();
TI_Vector_real = seq.getDefinition('Inversion'); 
% If phantom the following TI will obtained according with fixed RR 
% If in vivo the following TI will obtained according with squence timing

%% ************************************************************************
%                              T1 estimation   
Param_map_PSI = generate_map_pulseq(imspace_PSIR, Time_vector);


