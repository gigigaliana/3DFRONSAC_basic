function [Ax] = MakeSignal(a_res,b_res,c_res,ChunkSize_k,ChunkSize_p,DirectEncoding,Cmat,NP,object)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
nC=length(Cmat(1,:));

if size(nonzeros(object))>1
for F = 1:(a_res*b_res*c_res/ChunkSize_k) %scan through k-space
    
   % rem_time=tloop*(a_res*b_res*c_res/ChunkSize_k - F)/(60*60);
   if ChunkSize_k>a_res
       disp(['Chunk ' num2str(F) ' out of ' num2str(a_res*b_res*c_res/ChunkSize_k) ' remaining ' num2str(0) ' hrs'])
   end;
    iE_k= (1:ChunkSize_k) + (F-1)*ChunkSize_k;
    iP = 1:NP; %doing all the pixels at once
    iE_p = iP;
GradField=zeros(numel(iE_k),numel(iE_p));
%     for iF=1:size(traj_mat,2);
%         GradField= GradField+ traj_mat(iE_k,iF)*SEM_mat(iE_p,iF).' ;
% % % % % % %         figure; subplot(2,1,1); imagesc(reshape(GradField,[64 64]))
% % % % % % %                 subplot(2,1,2); imagesc(reshape(traj_mat(iE_k,iF)*SEM_mat(iE_p,iF).',[64 64])); title(num2str(traj_mat(iEk,iF)));
%     end
     GradField=DirectEncoding(iE_p,iE_k).';
     thisLine = exp( 1i*( GradField  ) ) ;
     
    for iC = 1:nC
        iK = (iC-1)*(a_res*b_res*c_res) + iE_k;
        
        Ax(iK) = ( thisLine .* squeeze(repmat(Cmat(iP,iC),[1 1 ChunkSize_k])).' ) * object(iP);
    end
    
end;

else

    TargetP=find(object,1);
for F = 1:(a_res*b_res*c_res/ChunkSize_k) %scan through k-space
    
   % rem_time=tloop*(a_res*b_res*c_res/ChunkSize_k - F)/(60*60);
   if ChunkSize_k>a_res
       disp(['Chunk ' num2str(F) ' out of ' num2str(a_res*b_res*c_res/ChunkSize_k) ' remaining ' num2str(0) ' hrs'])
   end;
    iE_k= (1:ChunkSize_k) + (F-1)*ChunkSize_k;
    iP = TargetP; %doing all the pixels at once
    iE_p = iP;
GradField=zeros(numel(iE_k),numel(iE_p));
%     for iF=1:size(traj_mat,2);
%         GradField= GradField+ traj_mat(iE_k,iF)*SEM_mat(iE_p,iF).' ;
% % % % % % %         figure; subplot(2,1,1); imagesc(reshape(GradField,[64 64]))
% % % % % % %                 subplot(2,1,2); imagesc(reshape(traj_mat(iE_k,iF)*SEM_mat(iE_p,iF).',[64 64])); title(num2str(traj_mat(iEk,iF)));
%     end
     GradField=DirectEncoding(iE_p,iE_k).';
     thisLine = exp( 1i*( GradField  ) ) ;
     
    for iC = 1:nC
        iK = (iC-1)*(a_res*b_res*c_res) + iE_k;
        
        Ax(iK) = ( thisLine .* squeeze(repmat(Cmat(iP,iC),[1 ChunkSize_k])).' );
    end
    
end;    
    
    
    
    
end
end

