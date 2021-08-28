function [ReconObject] = Ainv_timesb_psf_gpu(A_vars,FT_data)
%This generates PSF based recons for FT(kspacedata) in a FRONSAC
%acquisition
%A is the encoding matrix that converts an image X into coil-weighted FRONSAC
%data, so Ainv*X takes coil-weighted FRONSAC data and turns it into an
%image

addpath(genpath('Scripts'));
disp('CALCULATING PRODUCT WITH Atranspose');

SEMmat          =A_vars.SEMmat; %dimensions Nx,Ny,Nz,NumComponents
B1_2DxCoil      =A_vars.B1_2DxCoil; %dimensions Nx,Ny,Nz,NumCoils
traj_mat_1y1z   =A_vars.traj_mat_1y1z; %dimensions a_res,NumComponents
Ry              =A_vars.Ry; %assumes Nyq_y sampling is NumPE_y*Ry
Rz              =A_vars.Rz; %assumes Nyq_z sampling is NumPE_z*Rz

[Nx Ny Nz NumComponents]= size(SEMmat);
[Ns NumPE_y NumPE_z ncoils]=size(FT_data);

ReconObject     =zeros(Nx,Ny,Nz);
ReconObject_gpu =gpuArray(ReconObject);

if A_vars.CentralSliceOnlyFlag==0
    z_start=1; z_end=Nz/Rz;
elseif A_vars.CentralSliceOnlyFlag==1
    z_start=Nz/Rz/2; z_end=Nz/Rz/2;
    disp('Calculating a single z-slice (or slice set)');
end;

for  k=z_start:z_end
   % tic
    
    if Rz >1
        kz_set=sort(mod([k+Nz/2/Rz-1:Nz/Rz:(k+Nz)],Nz)+1);
    else
        kz_set=k;
    end;
    
    for j=1:Ny/Ry
        
        if Ry >1
            ky_set=sort(mod([j+Ny/2/Ry-1:Ny/Ry:(j+Ny)],Ny)+1);
        else
            ky_set=j;
        end;
        
        %Pulling out relevant encoding
        %fields for this/these lines
        SEM_mat_1y1z=reshape(squeeze(SEMmat(:,ky_set,kz_set,:)),[Nx Ry Rz NumComponents]);
        Cmat_1y1z=gpuArray(reshape(B1_2DxCoil(:,ky_set,kz_set,:),[Nx Ry Rz ncoils])); %8 ch
        %
        %                                           %%%%%%Gather up encodings for these
        %%%%%%folded lines of ky and kz
        for ry_idx1=1:Ry
            for rz_idx1=1:Rz
                DirectEncoding(:,ry_idx1,rz_idx1,:)=gpuArray(squeeze(SEM_mat_1y1z(:,ry_idx1,rz_idx1,:))*traj_mat_1y1z');
            end;
        end;
        
        [A_gpu] = MakeAmatrix_1y1zset(DirectEncoding,Cmat_1y1z,ncoils,Ns,Ry,Rz);
        
        FT_data_reshaped_gpu=gpuArray(reshape(FT_data(:,j,k,:),[Ns*ncoils 1]));
        ReconObject_gpu(:,ky_set,kz_set)=reshape( A_gpu'*FT_data_reshaped_gpu, [Nx Ry Rz]);
    end;
    
%     if mod(k,32)==1
%         disp(cat(2,'Recon: PE_z step ', num2str(k),' out of ',num2str(Nz/Rz),' Duration: ',num2str(toc)));
%     end;
end;
ReconObject=gather(ReconObject_gpu);
end

