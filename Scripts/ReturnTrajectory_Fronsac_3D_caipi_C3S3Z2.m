function  [K_list] = ReturnTrajectory_Fronsac_3D_caipi_SimRecon( FullySampledTrajectory, Nyq_x,Nyq_z,Ns,NumPE_y,NumPE_z,FOV,SlabThk, CaipiOn);
                                                                               %[set full sampled kmax] %[over/undersampling   %real units  %caipi flag             
             
% Version 1
Ry=Nyq_x/NumPE_y; 
Rz=Nyq_z/NumPE_z;
 if Ry<1
     Ry=1;
 end;
 

%     Units here are Cycles/meter
 kr_OS          = (-.5:          (1/(Ns)):        (.5-(1/Ns)))          *   (Nyq_x/FOV);      %Nyquist kr range but finer steps (oversampling)
 ky             = (-.5:          (1/(NumPE_y)):   (.5-(1/NumPE_y)))     *   (Nyq_x/FOV);      %Units are cycles/m
 kz             = (-.5:          (1/(NumPE_z)):   (.5-(1/NumPE_z)))     *   (Nyq_z/SlabThk);  %In ky and kz, these are already undersampled

     kz2_constant =0;

 
 if CaipiOn == 1
     delta_kz=kz(2)-kz(1);
     for p=1:NumPE_y
         shift=mod(p-1,2);
     CaipiShift(p)=shift*delta_kz/2;
     end;
 else
     CaipiShift(1:NumPE_y)=0;
 end;

%% Create initial (fully sampled) trajectory


if FullySampledTrajectory == 'Cartesian '
for zz=1:NumPE_z
    
    for jj=1:NumPE_y
        kx_vec_full(:,jj,zz)=kr_OS(:);%kr.*cos((jj-1)*pi/(sqrt(2)*NumRead));
        ky_vec_full(:,jj,zz)=ky(jj).*ones(size(kr_OS));%kr.*sin((jj-1)*pi/(sqrt(2)*NumRead));
        kz_vec_full(:,jj,zz)=(kz(zz)+ CaipiShift(jj)).*ones(size(kr_OS));%kr.*sin((jj-1)*pi/(sqrt(2)*NumRead));
   
%         % Patloc    
%         kz2_vec_full(:,jj,zz)   = 0*kx_vec_full(:,jj,zz);
%         kc3_vec_full(:,jj,zz)   = kz2_constant*kx_vec_full(:,jj,zz);
%         ks3_vec_full(:,jj,zz) = kz2_constant*ky_vec_full(:,jj,zz);
%         
%         % Wave
%         
%         ky_vec_full(:,jj,zz)=ky(jj).*ones(size(kr_OS));%kr.*sin((jj-1)*pi/(sqrt(2)*NumRead));
%         kz_vec_full(:,jj,zz)=(kz(zz)+ CaipiShift(jj)).*ones(size(kr_OS));%kr.*sin((jj-1)*pi/(sqrt(2)*NumRead));

        
        kz2_vec_full(:,jj,zz)   = 0*kr_OS(:);
        kc3_vec_full(:,jj,zz)   = 0*kr_OS(:);
        ks3_vec_full(:,jj,zz)   = 0*kr_OS(:);  


    end; %end loop over jj=NumPE
end; %end loop over Nz
end;%end elseif branch for Spiral/Cartesian/Radial
   

kx_vec_full=kx_vec_full(:);
ky_vec_full=ky_vec_full(:);
kz_vec_full=kz_vec_full(:);
if NumPE_z==1
    kz_vec_full=zeros(size(kz_vec_full(:)));
end
kz2_vec_full=kz2_vec_full(:);
kc3_vec_full=kc3_vec_full(:);
ks3_vec_full=ks3_vec_full(:);
size(kx_vec_full);
size(kz2_vec_full);


K_list=cat(2,-kc3_vec_full,-ks3_vec_full,-kz2_vec_full,kx_vec_full,ky_vec_full,kz_vec_full);  % rads/m




end