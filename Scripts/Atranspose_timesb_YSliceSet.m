function  [ReconObject_SomeSlices] = Atranspose_timesb_YSliceSet(k,U)


Ry=U.Ry;
Rz=U.Rz;

[Nx Ny Nz NumComponents]= size(U.SEMmat);
[Ns NumPE_y NumPE_z ncoils]=size(U.FT_data);

ReconObject_gpu =gpuArray(zeros(Nx,Ny,Nz));

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
        SEM_mat_1y1z=reshape(squeeze(U.SEMmat(:,ky_set,kz_set,:)),[Nx Ry Rz NumComponents]);
        Cmat_1y1z=gpuArray(reshape(U.B1_2DxCoil(:,ky_set,kz_set,:),[Nx Ry Rz ncoils])); %8 ch
        %
        %                                           %%%%%%Gather up encodings for these
        %%%%%%folded lines of ky and kz
        for ry_idx1=1:Ry
            for rz_idx1=1:Rz
                DirectEncoding(:,ry_idx1,rz_idx1,:)=gpuArray(squeeze(SEM_mat_1y1z(:,ry_idx1,rz_idx1,:))*U.traj_mat_1y1z');
            end;
        end;
        
        [A_gpu] = MakeAmatrix_1y1zset(DirectEncoding,Cmat_1y1z,ncoils,Ns,Ry,Rz);
        
        FT_data_reshaped_gpu=gpuArray(reshape(U.FT_data(:,j,k,:),[Ns*ncoils 1]));
        ReconObject_gpu(:,ky_set,kz_set)=reshape( A_gpu'*FT_data_reshaped_gpu, [Nx Ry Rz]);
    end;
    ReconObject_SomeSlices(:,:,1:Rz)=gather(ReconObject_gpu(:,:,kz_set));