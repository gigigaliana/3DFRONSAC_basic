function [A_gpu] = MakeAmatrix_1y1zset(DirectEncoding,Cmat_1y1z,ncoils,Ns,Ry,Rz)
%This function used by A_timesX_psf_gpu, Atranspose_timesb_psf_gpu, and Ainv_timesb_psf_gpu
%it expects to receive gpu variables for DirectEncoding and Cmat, and returns a gpu variable A_gpu

Nx=size(Cmat_1y1z,1);
        A_gpu=gpuArray(zeros(Ns*ncoils,Nx*Ry*Rz));        
        iP=1:Nx; %For every x-location in the recon
            for ry_idx=1:Ry
                for rz_idx=1:Rz
                    Phasors_allX_gpu=exp(1i*squeeze(DirectEncoding(iP,ry_idx,rz_idx,:))); %this is Nx x Ns
                    SinglePointData(:,ry_idx,rz_idx,iP)=reshape(bsxfun(@times,Phasors_allX_gpu,Cmat_1y1z(iP,ry_idx,rz_idx,:)),[Nx Ns*ncoils]).';
                    for cc=1:ncoils
                        A_gpu((cc-1)*Ns+1:cc*Ns,iP+(ry_idx-1)*Nx +(rz_idx-1)*(Nx*Ry))=squeeze(flipud(fftshift(fft(fftshift(SinglePointData((cc-1)*Ns+1:cc*Ns,ry_idx,rz_idx,:),1),[],1),1)));
                    end
                end
            end 

end

