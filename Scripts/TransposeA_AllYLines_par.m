function [ReconObject_gpu] = TransposeA_AllYLines_par(j,k,kz_set,SEMmat,traj_mat_1y1z,B1_2DxCoil,FT_data,LittleVarContainer)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
t=num2cell(LittleVarContainer);
 [Nx,Ny,Nz,Ry,Rz,NumComponents,ncoils,Ns]=deal(t{:});
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
  ReconObject_gpu=ReconObject_gpu(:,:,kz_set);
end

