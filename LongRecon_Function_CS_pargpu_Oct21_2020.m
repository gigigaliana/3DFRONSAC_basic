%This code starts to use functions for recons
%Also, orientation is fixed and no longer permutable as of July 2020


%%%% 10/21/2020: I'm suspicious because fnlcg seems to change the output even when xfm is zero and tv is teeny, yet larger TV makes little difference.
%%%%  Also, it introduces a huge scale factor around 33, and the data term starts huuuuuuuge.  weird weird weird.
%%%% One possibility is that the weird res=XFM*res before input to fnlcg is not implemented correctly, so I rewrite fnl and the calling script to 
%%%% expect the input s.t. res = ReconObject (straight matrix inversion.
close all, clear all


for data_idx=[2] %which dataset
    for Rz = [2,1]
        if Rz==2
            Ry=4;
        elseif Rz==1
            Ry=1;
        end;
        
        clearvars -except data_idx Ry Rz
        warning off all;
        addpath('ScanData/');
        addpath(genpath('Scripts/'));

        
        OutputSuffix='Subject3_1123';
        B1_map_Filename=cat(2,'FieldsAndProfiles/GG_April_Poly_B1Maps_shiftedFit_',OutputSuffix,'.mat');
        NLG_map_Filename=cat(2,'FieldsAndProfiles/GG_July_FRONField_MatsAndTrajs_UnwrappedPhase_NetImgFTlikeData.mat');
        NLG         = 'on ';
        
        
        %% Prepare the    initial parameters:
        
        LinearMode='Cartesian ';%'EchoPlanar';%%'Radial    ';%%% % %%%%'Radial   %%'Spiral    ';
        CaipiOn = 0;
        R_RO        = 1;
        ncoils=8;
        ReconSingleCoil=0; %0=multicoil
        
        %%%%%%% Nyq_xyz sets k_max for fully sampled data
        Nyq_x       = 128;  %this is used for setting trajectory
        Nyq_z       = 128;  %it's also the size of fully sampled PE directions (but NOT 1st dim Ns)
        Nyq_y       = 128;
        
        Nx          = 184;  % Recon matrix size, should match 2xmap size, or else maps will be trimmed to that size assuming equal/double mm/pix
        Ny          = 128;
        Nz          = 128;  % Recon matrix size in z, should also match map size
        FronmatN    = 96;   % should allow map size =! N, but note this is usually pre-trimmed to a cube = FOV
        FronmatNz   = 64;   % should allow map size =! Nz
        
        %%%%%% Ns, NumPE_y/z set spacing/skipping in k vecs
        NumPE_y     = Nyq_y/Ry; %error if NumPE>Nyq_y, this is actual NumPE after undersampling
        NumPE_z     = Nyq_z/Rz;
        Ns          = 1024;      % Samples during the Readout
        
        EncodingFOVx        = 0.25;      % Field of View [m] used to set trajectory
        EncodingFOVy        = 0.25;
        EncodingSlabThk     = 0.25;      % Field of View in z used to set trajectory
        
        ReconFOVx        = 0.25*(Nx/128);% Field of View [m] used to set Linear MATs for recon
        ReconFOVy        = 0.25;
        ReconSlabThk     = 0.25;         % Field of View in z used to set Linear MATs for recon
        
        SaveVars=0;
        curr_path = fileparts(mfilename('fullpath'));
        cd(curr_path);
        
        
        
        %% READ IN SIGNAL
        
        if data_idx==1
            if NLG=='on '
                DataFile='ScanData/191123_Subject3/meas_MID281_FR_3D_GRE_8ch_64cyc_5_1amp_FID85885.dat';
                dat=ReadSiemensMeasVB17(DataFile,'off');
                DataName= cat(2,'Fron_GRE_',OutputSuffix);
            else
                DataFile='ScanData/191123_Subject3/meas_MID283_Cart_3D_GRE_8ch_FID85887.dat';
                dat=ReadSiemensMeasVB17(DataFile,'off');
                DataName= cat(2,'Cart_GRE',OutputSuffix);
            end;
        elseif data_idx==2
            if NLG=='on '
                DataFile='ScanData/191123_Subject3/meas_MID285_FR_3D_MPRAGE_8ch_64cyc_5_1amp_FID85889.dat';
                dat=ReadSiemensMeasVB17(DataFile,'off');
                DataName= cat(2,'Fron_MPRAGE',OutputSuffix);
            else
                DataFile='ScanData/191123_Subject3/meas_MID287_Cart_3D_MPRAGE_8ch_FID85891.dat';
                dat=ReadSiemensMeasVB17(DataFile,'off');
                DataName= cat(2,'Cart_MPRAGE',OutputSuffix);
            end;
        end;
        
        disp(cat(2,'Trimming ',num2str(size(dat,1)-ncoils*Nyq_y*Nyq_z),' lines of data.  Check pars if this is a lot'));
        dat=dat(1:ncoils*Nyq_y*Nyq_z,:);
        dat=reshape(dat,ncoils, Nyq_y, Nyq_z, Ns);
        dat=permute(dat, [4 2 3 1]); %[read phase slab/partition coil]
        
        %datus = undersampled data
        datus=dat(:,1:Ry:end,1:Rz:end,:);
        datus=datus(1:R_RO:end,:,:,:);
        
        %% LOAD B1 MAPS AND TRIM
        
        if ncoils > 1
            %This B1_2DxCoil map should be have same mm/pix
            load(B1_map_Filename,'B1_2DxCoil','OutputSuffix');
        elseif ncoils ==1
            B1_2DxCoil=ones(Nx,Ny,Nz,1);
        end;
        
        [B1Nx B1Ny B1Nz B1Nc]=size(B1_2DxCoil);
        
        %%%%%%%%%%%%Resize or trim coils to N and Nz%%%%%%%%%%%%%%%%%%%
        
        if sum([B1Nx B1Ny B1Nz] ~= [Nx Ny Nz])>0
            
            %%%% Routine to resize coil profiles taken at
            %%%% different resolution, same FOV
            %                     for ci=1:ncoils
            %                         B1_2DxCoil_res(:,:,:,ci)=imresize3(B1_2DxCoil(:,:,:,ci), [Nx Ny Nz]); %8 channel coil
            %                     end
            %                     B1_2DxCoil=B1_2DxCoil_res;
            %                     clear B1_2DxCoil_res;
            %                     disp('looks like coil profiles needed resizing; equal FOV assumed');
            
            
            %%%% Routine to trim coil profiles to the recon FOV
            %%%% taken at same mm/pixel resolution
            B1trim=(B1Nx-Nx)/2;
            B1_2DxCoil=B1_2DxCoil(B1trim+1:B1Nx-B1trim,:,:,:);
            clear B1_2DxCoil_res;
            disp('Looks like B1 maps need trimming... equal voxel size being assumed');
            
        end;
        
        B1_2DxCoil=flipdim(B1_2DxCoil,1);
        %%%%%%%%%%%% END MULTICHANNEL CASE %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% DEFINE/LOAD ENCODING FIELDS
        xxx=(((-1*Nx/2):1/1:(Nx/2)-.5/1)) * (ReconFOVx/Nx); % [m]
        yyy=(((-1*Ny/2):1/1:(Ny/2)-.5/1)) * (ReconFOVy/Ny); % [m]
        zzz=(((-1*Nz/2):1/1:(Nz/2)-.5/1)) * (ReconSlabThk/Nz);
        [XXX,YYY,ZZZ]=ndgrid(xxx,yyy,zzz);
        
        
        
        %%%% Defining pure/unmeasured fields (kx/ky)
        SEMmat(:,:,:,1) = XXX.*(XXX.^2 - 3*(YYY.^2));     % C3 %ideal nonlinears
        SEMmat(:,:,:,2) = YYY.*(3*(XXX.^2) - YYY.^2);     % S3 %ideal nonlinears
        SEMmat(:,:,:,3) = -0.5*(XXX.^2 + YYY.^2)+ZZZ.^2;  % Z2 %ideal nonlinears
        SEMmat(:,:,:,4) = XXX; %linears
        SEMmat(:,:,:,5) = YYY; %linears
        SEMmat(:,:,:,6) = ZZZ; %linears
        
        %%%%%%%%%% adding SEMS from FRON field
        %%%%%%%%%% % true nonlinears from 7 to NumComponents
        
        if NLG == 'on '
            % These FRONSAC maps may cover a larger FOV than FOVx,
            % but they will be trimmed to an equal FOV assuming
            % mm/pix of map = 2*(mm/pix of recon)
            
            load(NLG_map_Filename,'FronsacMAT','FronsacTraj');
            NumComponents=size(FronsacMAT,2);
            FronsacMAT=reshape(FronsacMAT,[FronmatN FronmatN FronmatNz NumComponents]);
            
            %%%This flip instituted July 2020
            FronsacMAT=flipdim(FronsacMAT,1);
            
            %%% This assumes mm/pix in NLG map = 2*mm/pix in data
            NLGTrim=(FronmatN-(Nx/2))/2;
            FronsacMAT=FronsacMAT(NLGTrim+1:FronmatN-NLGTrim,17:80,:,:);
            disp(cat(2,'FronsacMAT trimmed to ',num2str(size(FronsacMAT))));
            
            for n=1:NumComponents
                SEMmat(:,:,:,n+6)=imresize3(FronsacMAT(:,:,:,n),[Nx Ny Nz]);
            end
            
            disp(cat(2,'Trimmed FronsacMAT resampled to matrix size ',num2str(size(SEMmat))));
            
        else
            NumComponents=0;
        end
        
        
        %                 %%%% Mask for cylindrical FOV, long in the HF direction,
        %                 %%%% *********ASSUMES READ IS HF DIRECTION**********
        %                 mask=zeros([Nx Ny Nz]);
        %                 mask( (YYY.^2 + ZZZ.^2) < FOVy.^2) = 1;
        %
        %                 B1_2DxCoil = B1_2DxCoil.*repmat(mask,[1 1 1 ncoils]);
        
        %%%%%%%%%%%%%%%%%%%%%% from ND code:
        
        %% DEFINE ENCODING TRAJECTORIES
        %IMPORTANT: THIS IS NOT THE TRAJECTORY FOR Nx,Ny,Nz, it
        %is a Nyquist trajectory for Nyq_x,Nyq_y,Nyq_z res over
        %FOVx,FOVx,SlabThk, which is then
        %over/undersampled to Ns,NumPE_y,NumPE_z
        [traj_mat_t] = 2*pi*ReturnTrajectory_Fronsac_3D_caipi_C3S3Z2( 'Cartesian ',Nyq_x,Nyq_z,Ns,NumPE_y,NumPE_z,EncodingFOVx,EncodingSlabThk,CaipiOn); %[rad]*[cyc/m]
        
        %For PSF recon, we only need the readout encoding
        traj_mat_1y1z     =  traj_mat_t(1:Ns,:); %includes component 4 which defines kx
        traj_mat_1y1z(:,5)=0; %ky=0 for PSF recon
        traj_mat_1y1z(:,6)=0; %kz=0 for PSF recon
        
        %NLG field trajectories get appended
        if NLG == 'on '
            traj_mat_1y1z(:,7:(NumComponents+6))=FronsacTraj;
        end;
        
        %% FT DATA FOR PSF RECON
        %To match B1 maps, necks should point down and
        %noses to the right
        for c=1:ncoils
            fftShiftedData=fftshift(fftshift(fftshift(datus(:,:,:,c),1),2),3); %undersampled data
            FT_data(:,:,:,c)=fftshift(fftshift(fftshift(fftn(fftShiftedData,[Ns Nyq_y/Ry Nyq_z/Rz]),1),2),3); %folded image
            ShiftedData_full=ifftshift(ifftshift(ifftshift(dat(:,:,:,c),1),2),3); %original data
            Full_FT_data(:,:,:,c)=fftshift(fftshift(fftshift(fftn(ShiftedData_full,[Ns Nyq_y Nyq_z]),1),2),3); %full FOV img
        end;
        FT_data_flip=FT_data;
        FT_data_flip=flipdim(FT_data_flip,1);
        FT_data_flip=flipdim(FT_data_flip,2);
        
        
        
        
        %% SET DATA NAME
        %This line adds acceleration info to data choice
        DataName=cat(2,'R',num2str(Ry),'x',num2str(Rz),'_',DataName);
        
        
        
        
        %% CODE FOR SINGLE COIL RECON MODE
        if ReconSingleCoil==1
            FT_data_flip=FT_data_flip(:,:,:,1);
            B1_2DxCoil=ones(Nx,Ny,Nz);
            ncoils=1;
            DataName=cat(2,DataName,'_1ch');
        end;
        
        
        %% GENERATE PARS FOR ENCODING MATRIX
        
        A_vars.SEMmat          = SEMmat; ; %dimensions Nx,Ny,Nz,NumComponents
        A_vars.B1_2DxCoil      = B1_2DxCoil; %dimensions Nx,Ny,Nz,NumCoils
        A_vars.traj_mat_1y1z   = traj_mat_1y1z; %dimensions a_res,NumComponents
        A_vars.Ry              = Ry; %assumes Nyq_y sampling is NumPE_y*Ry
        A_vars.Rz              = Rz; %assumes Nyq_z sampling is NumPE_z*Rz
        A_vars.CentralSliceOnlyFlag = 0; %1 means only the central slice is used/reconstructed
        
        tic
        [ReconObject] = Ainv_timesb_psf_gpu(A_vars,FT_data_flip);
        toc
 
                
            end %TVWeight loop
        end %xfmWeight loop
        
    end; %Rz loop
end; %data_idx loop


