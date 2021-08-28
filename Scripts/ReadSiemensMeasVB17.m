function [rawdata,loopcounters,sMDH] = ReadSiemensMeasVB17(measfile, dispopt)
%function [rawdata,loopcounters,sMDH] = ReadSiemensMeasVB15(measfile, dispopt)
% Reads the 2D measurement data got from Siemens scanners following 
% (mstart mate->3->10->2(Rawdata2Disk ON)); 
% I assume the file extesion .out.This mfunction extracts imaging information from the
% mini MDH in the .out file. No .asc file is needed.
%
% The .out file has 
%     1) an uint32 indicating hdr_len (only from VB15);
%     2) an common header of hdr_len including the uint32 followed by all ADC data (only from VB15);
%     3) each ADC acquisition is oversampled by 2 and has a 128 byte header;
%     4) each pixel consists of the read and imaginary parts, which are float (4 bytes);
%     5) 384/1024? bytes are put at the end. 
% 
% dispopt, - 'on' or 'off'. The default is 'on'. The FFTed images for slices of the first 3D 
% is displayed only for your reference. if your data were segmented, PF, ReadOutOffcentre, or 
% arranged in the special way, you might not see the correct images. You have to turn this 
% option off, and work on the rawdata output from this function for the right raw line order.
%
% Output:
% rawdata, - raw data in matrix; oversmapled (x2) in complex (A+iB), indexed by loopcounters and
% ulChannelId;
% loopcounters,- all the sMDH counters and ChannelId corresponding to all the ADC data lines; 
% All dimension are in the same order as saved in sMDH.sLC;
% sMDH, - the MDH of the last ADC acquision except that extended sMDH.sLC is set to maximum 
% loop values.
%
% EXAMPLE:
%           [rawdata loopcounters sMDH] = ReadSiemensMeasVB15('Z:\n4\pkg\MeasurementData\FLASH\meas.out');
%
% Maolin Qiu YALE 5-28-2010 (maolin.qiu(at)yale(dot)edu)

if nargin < 1
    help ReadSiemensMeasVB17;
    [filename pathname] = uigetfile( ...
       {'*.out';'*.dat';'*.*'}, ...
        'Pick a Siemens MESUREMENT file');
    if ~filename & ~pathname
        disp(['You selected no file.']);
        return;
    else
        measfile = fullfile(pathname, filename);
    end
end

if nargin < 2
    dispopt = 'on';
end

% measfile

[pa na ex] = fileparts(measfile);
measfile = fullfile(pa, [na ex]);
if exist(measfile, 'file')
    disp(['Measurement data file: ' measfile]);
else
    disp(['Measurement data file does not exist: ' measfile]);
    return;
end

%% init output

data = []; dimensions = []; loopcounters = []; sMDH = [];

%% Constants used in sMDH

MDH_NUMBEROFEVALINFOMASK   = 2;
MDH_NUMBEROFICEPROGRAMPARA = 4;
MDH_FREEHDRPARA            = 4;

IDX_DIM                    = 20;

%%--------------------------------------------------------------------------%%
%% Definition of loop counter structure                                     %%
%% Note: any changes of this structure affect the corresponding swapping    %%
%%       method of the measurement data header proxy class (MdhProxy)       %%
%%--------------------------------------------------------------------------%%

sLoopCounter = struct( ...
  'ushLine',0,...                  %% unsigned short  line index                   %%
  'ushAcquisition',0,...           %% unsigned short  acquisition index            %%
  'ushSlice',0,...                 %% unsigned short  slice index                  %%
  'ushPartition',0,...             %% unsigned short  partition index              %%
  'ushEcho',0,...                  %% unsigned short  echo index                   %%	
  'ushPhase',0,...                 %% unsigned short  phase index                  %%
  'ushRepetition',0,...            %% unsigned short  measurement repeat index     %%
  'ushSet',0,...                   %% unsigned short  set index                    %%
  'ushSeg',0,...                   %% unsigned short  segment index  (for TSE)     %%
  'ushIda',0,...                   %% unsigned short  IceDimension a index         %%
  'ushIdb',0,...                   %% unsigned short  IceDimension b index         %%
  'ushIdc',0,...                   %% unsigned short  IceDimension c index         %%
  'ushIdd',0,...                   %% unsigned short  IceDimension d index         %%
  'ushIde',0 ...                   %% unsigned short  IceDimension e index         %%
);                                 %% sizeof : 28 byte             %%

%%--------------------------------------------------------------------------%%
%%  Definition of slice vectors                                             %%
%%--------------------------------------------------------------------------%%

sVector = struct( ...
  'flSag',0.0,...       %% float
  'flCor',0.0,...       %% float
  'flTra',0.0 ...       %% float
);

sSliceData = struct( ...
  'sSlicePosVec',sVector,...                   %% slice position vector               %%
  'aflQuaternion',zeros(1,4) ...               %% float rotation matrix as quaternion %%
);                                              %% sizeof : 28 byte                    %%

%%--------------------------------------------------------------------------%%
%%  Definition of cut-off data                                              %%
%%--------------------------------------------------------------------------%%

sCutOffData = struct( ...
  'ushPre',0,...               %% unsigned short  write ushPre zeros at line start %%
  'ushPost',0 ...              %% unsigned short  write ushPost zeros at line end  %%
);

%%--------------------------------------------------------------------------%%
%%  Definition of measurement data header                                   %%
%%--------------------------------------------------------------------------%%

sMDH = struct( ...
  'ulDMALength',0,...                                       %% unsigned long  DMA length [bytes] must be                        4 bytes %% first parameter                        
  'lMeasUID',0,...                                          %% long           measurement user ID                               4     
  'ulScanCounter',0,...                                     %% unsigned long  scan counter [1...]                               4
  'ulTimeStamp',0,...                                       %% unsigned long  time stamp [2.5 ms ticks since 00:00]             4
  'ulPMUTimeStamp',0,...                                    %% unsigned long  PMU time stamp [2.5 ms ticks since last trigger]  4
  'aulEvalInfoMask',zeros(1,MDH_NUMBEROFEVALINFOMASK),...   %% unsigned long  evaluation info mask field                        8
  'ushSamplesInScan',0,...                                  %% unsigned short # of samples acquired in scan                     2
  'ushUsedChannels',0,...                                   %% unsigned short # of channels used in scan                        2   =32
  'sLC',sLoopCounter,...                                    %% loop counters                                                    28  =60
  'sCutOff',sCutOffData,...                                 %% cut-off values                                                   4           
  'ushKSpaceCentreColumn',0,...                             %% unsigned short centre of echo                                    2
  'ushDummy',0,...                                          %% unsigned short for swapping                                      2
  'fReadOutOffcentre',0.0,...                               %% float          ReadOut offcenter value                           4
  'ulTimeSinceLastRF',0,...                                 %% unsigned long  Sequence time stamp since last RF pulse           4
  'ushKSpaceCentreLineNo',0,...                             %% unsigned short number of K-space centre line                     2
  'ushKSpaceCentrePartitionNo',0,...                        %% unsigned short number of K-space centre partition                2
  'aushIceProgramPara',zeros(1,MDH_NUMBEROFICEPROGRAMPARA),... %% unsigned short free parameter for IceProgram                  8  =88
  'aushFreePara',zeros(1,MDH_FREEHDRPARA),...               %% unsigned short free parameter                          4 * 2 =   8   
  'sSD',sSliceData,...                                      %% Slice Data                                                       28 =124
  'ulChannelId',0 ...                                       %% unsigned long	 channel Id must be the last parameter          4
);                                                          %% total length: 32 * 32 Bit (128 Byte)                             128
%% MDH_H %%

%% read header information

fid = fopen(measfile,'r');
hdr_len = fread(fid, 1, 'uint32'); %%..maolin modified according to Florian.Knoll at IDEA Board 4/30/2009
disp(['Leangth of the Header : ' num2str(hdr_len)]);
fseek(fid, hdr_len, 'bof'); %%..maolin modified according to Florian.Knoll at IDEA Board 4/30/2009
[tmp Nbytes] = fread(fid,inf,'uchar'); tmp = []; FileSize = ftell(fid); disp(['Original File Size : ' num2str(FileSize)]);
disp(['The size of the measurement data without the header is (bytes): ' num2str(Nbytes) '. This will take long if it is HUGE ... ...']);
%% frewind(fid);
%% fseek(fid, 32, -1); %%..maolin modified according to Florian.Knoll at IDEA Board 4/30/2009
fseek(fid, hdr_len, 'bof'); %%..maolin modified according to Florian.Knoll at IDEA Board 4/30/2009
LastADCData = 0;
lcounter = 1;

while ~feof(fid) & ~LastADCData
    
   sMDH.ulDMALength                             = fread(fid, 1, 'uint32');      % 4
   sMDH.lMeasUID                                = fread(fid, 1,  'int32');      % 8
   sMDH.ulScanCounter                           = fread(fid, 1, 'uint32');      % 12
   sMDH.ulTimeStamp                             = fread(fid, 1, 'uint32');      % 16
   sMDH.ulPMUTimeStamp                          = fread(fid, 1, 'uint32');      % 20
   for i = 1:MDH_NUMBEROFEVALINFOMASK       % 2
        sMDH.aulEvalInfoMask(i)                 = fread(fid, 1, 'uint32');      % 20 + 2 * 4 = 28
   end
   sMDH.ushSamplesInScan                        = fread(fid, 1, 'uint16');      % 30
   sMDH.ushUsedChannels                         = fread(fid, 1, 'uint16');      % 32
   sMDH.sLC.ushLine                             = fread(fid, 1, 'uint16');
   sMDH.sLC.ushAcquisition                      = fread(fid, 1, 'uint16');
   sMDH.sLC.ushSlice                            = fread(fid, 1, 'uint16');
   sMDH.sLC.ushPartition                        = fread(fid, 1, 'uint16');
   sMDH.sLC.ushEcho                             = fread(fid, 1, 'uint16');
   sMDH.sLC.ushPhase                            = fread(fid, 1, 'uint16');
   sMDH.sLC.ushRepetition                       = fread(fid, 1, 'uint16');
   sMDH.sLC.ushSet                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushSeg                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushIda                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushIdb                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushIdc                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushIdd                              = fread(fid, 1, 'uint16');
   sMDH.sLC.ushIde                              = fread(fid, 1, 'uint16');      % 32 + 14 * 2 = 60
   sMDH.sCutOff.ushPre                          = fread(fid, 1, 'uint16');
   sMDH.sCutOff.ushPost                         = fread(fid, 1, 'uint16');      % 60 + 2 * 2 = 64
   sMDH.ushKSpaceCentreColumn                   = fread(fid, 1, 'uint16');
   sMDH.ushDummy                                = fread(fid, 1, 'uint16');      % 64 + 2 * 2 = 68
   sMDH.fReadOutOffcentre                       = fread(fid, 1, 'float');       % 68 + 4 = 72
   sMDH.ulTimeSinceLastRF                       = fread(fid, 1, 'uint32');
   sMDH.ushKSpaceCentreLineNo                   = fread(fid, 1, 'uint16');
   sMDH.ushKSpaceCentrePartitionNo              = fread(fid, 1, 'uint16');      % 72 + 4 + 2 + 2 = 80
   for i = 1:MDH_NUMBEROFICEPROGRAMPARA    % 4
        sMDH.aushIceProgramPara(i)              = fread(fid, 1, 'uint16');      % 80 + 4 * 2 = 88
   end
   for i = 1:MDH_FREEHDRPARA  % 4
        sMDH.aushFreePara                       = fread(fid, 1, 'uint16');      % 88 + 4 * 2 = 96
   end
   sMDH.sSD.sVector.flSag                       = fread(fid, 1, 'float');
   sMDH.sSD.sVector.flCor                       = fread(fid, 1, 'float');
   sMDH.sSD.sVector.flTra                       = fread(fid, 1, 'float');       % 96 + 3 * 4 = 108
   for i = 1:4
        sMDH.aflQuaternion(i)                   = fread(fid, 1, 'float');       % 108 + 4 * 4 = 124
   end
   sMDH.ulChannelId                             = fread(fid, 1, 'uint32');      % 124 + 4 = 128 OK!
   
   if lcounter == 1
       aADC = -2*ones(1,IDX_DIM+sMDH.ushSamplesInScan*2); % the first 14 entries are the indices of this ADC in terms of sLoopCounter
       aADC(1)  = sMDH.sLC.ushLine;
       aADC(2)  = sMDH.sLC.ushAcquisition;
       aADC(3)  = sMDH.sLC.ushSlice;
       aADC(4)  = sMDH.sLC.ushPartition;
       aADC(5)  = sMDH.sLC.ushEcho;
       aADC(6)  = sMDH.sLC.ushPhase;
       aADC(7)  = sMDH.sLC.ushRepetition;
       aADC(8)  = sMDH.sLC.ushSet;
       aADC(9)  = sMDH.sLC.ushSeg;
       aADC(10) = sMDH.sLC.ushIda;
       aADC(11) = sMDH.sLC.ushIdb;
       aADC(12) = sMDH.sLC.ushIdc;
       aADC(13) = sMDH.sLC.ushIdd;
       aADC(14) = sMDH.sLC.ushIde;
       aADC(15) = sMDH.ulChannelId;
       % ... ...
       aADC(IDX_DIM+1:end) = fread(fid, sMDH.ushSamplesInScan*2, 'float'); % 4 bytes in each float
       LADC = 128 + sMDH.ushSamplesInScan*2*4;   % in bytes now
       NADC = round((Nbytes-384-3*52*4)/LADC); % discard last 3 acquisitions of 52 floats added by VB
       data = zeros(NADC, size(aADC,2)); % Optimize to speed up!
       data(lcounter, :) = aADC;
   else
       aADC = -2*ones(1,IDX_DIM+sMDH.ushSamplesInScan*2); % the first 14 entries are the indices of this ADC in terms of sLoopCounter
       aADC(1)  = sMDH.sLC.ushLine;
       aADC(2)  = sMDH.sLC.ushAcquisition;
       aADC(3)  = sMDH.sLC.ushSlice;
       aADC(4)  = sMDH.sLC.ushPartition;
       aADC(5)  = sMDH.sLC.ushEcho;
       aADC(6)  = sMDH.sLC.ushPhase;
       aADC(7)  = sMDH.sLC.ushRepetition;
       aADC(8)  = sMDH.sLC.ushSet;
       aADC(9)  = sMDH.sLC.ushSeg;
       aADC(10) = sMDH.sLC.ushIda;
       aADC(11) = sMDH.sLC.ushIdb;
       aADC(12) = sMDH.sLC.ushIdc;
       aADC(13) = sMDH.sLC.ushIdd;
       aADC(14) = sMDH.sLC.ushIde;
       aADC(15) = sMDH.ulChannelId;
       % ... ...
       aADC(IDX_DIM+1:end) = fread(fid, sMDH.ushSamplesInScan*2, 'float');
       if size(aADC,2) == size(data,2)
           data(lcounter,:) = aADC;
       else
           disp(['WARNING: Discarded scan lcounter = ' num2str(lcounter) ' size(aADC,2) = ' num2str(size(aADC,2)) ' that is different from the previously detected ' num2str(size(data,2)-IDX_DIM)]);
       end
   end

   curr_pos = ftell(fid);
   if curr_pos >= hdr_len + Nbytes - 384, LastADCData = 1; end % exclude the last 384 bytes of tail
   lcounter = lcounter + 1;
   
end % of while

fclose(fid);

loopcounters = round(data(:,1:IDX_DIM)+1); % all start from 1, unused -1
dimensions = max(loopcounters);

% Let sMDH have the maximum values of all loop counters

sMDH.sLC.ushLine = dimensions(1)  ; 
sMDH.sLC.ushAcquisition = dimensions(2)  ; 
sMDH.sLC.ushSlice = dimensions(3)  ; 
sMDH.sLC.ushPartition = dimensions(4)  ; 
sMDH.sLC.ushEcho = dimensions(5)  ; 
sMDH.sLC.ushPhase = dimensions(6)  ; 
sMDH.sLC.ushRepetition = dimensions(7)  ; 
sMDH.sLC.ushSet = dimensions(8)  ; 
sMDH.sLC.ushSeg = dimensions(9)  ; 
sMDH.sLC.ushIda = dimensions(10) ; 
sMDH.sLC.ushIdb = dimensions(11) ; 
sMDH.sLC.ushIdc = dimensions(12) ; 
sMDH.sLC.ushIdd = dimensions(13) ; 
sMDH.sLC.ushIde = dimensions(14) ; 
sMDH.sLC.ulChannelId = dimensions(15) ;
% ... ...

rawdata = complex(data(:,IDX_DIM+1:2:end),data(:,IDX_DIM+2:2:end));

if strcmp(dispopt, 'on')
    for slice = 1:sMDH.sLC.ushSlice
        % select raw lines -
        lines = find(loopcounters(:,2)  == 1 ...        % sMDH.sLC.ushAcquisition
                   & loopcounters(:,3)  == slice ...    % sMDH.sLC.ushSlice
                   & loopcounters(:,4)  == 1 ...        % sMDH.sLC.ushPartition
                   & loopcounters(:,5)  == 1 ...        % sMDH.sLC.ushEcho
                   & loopcounters(:,6)  == 1 ...        % sMDH.sLC.ushPhase
                   & loopcounters(:,7)  == 1 ...        % sMDH.sLC.ushRepetition
                   & loopcounters(:,8)  == 1 ...        % sMDH.sLC.ushSet %  & loopcounters(:,9)  == 1 ...        % sMDH.sLC.ushSeg
                   & loopcounters(:,10) == 1 ...        % sMDH.sLC.ushIdb
                   & loopcounters(:,11) == 1 ...        % sMDH.sLC.ushIda
                   & loopcounters(:,12) == 1 ...        % sMDH.sLC.ushIda
                   & loopcounters(:,13) == 1 ...        % sMDH.sLC.ushIda
                   & loopcounters(:,14) == 1 ...        % sMDH.sLC.ushIda
                   & loopcounters(:,15) == sMDH.ulChannelId + 1 ...        % sMDH.ulChannelId
                   );
        araw = rawdata(lines,:);
        lc = loopcounters(lines,:);
        
        % get rid of extra lines, e.g., the phase correction lines
        istart = find(lc(:,1) == 1);
        istop  = find(lc(:,1) == sMDH.sLC.ushLine);
        araw = araw(istart(1):istop(1),:);
        % sMDH.fReadOutOffcentre
        % ... ... to add
        % upper PF 
        rows = sMDH.sLC.ushLine; cols = sMDH.ushSamplesInScan/2; % overampling fator of 2
        if rows < cols
            A = zeros(cols, sMDH.ushSamplesInScan); % each row is one complex ADC data
           for i = 1:cols-rows
               A(i,:) = araw(rows-i+1,:);
           end
          A(cols-rows+1:cols,:) = araw;
          araw = A;
        end
        tmp = fftshift(fft2(fftshift(araw)));
        aimg = tmp(:,1+sMDH.ushSamplesInScan/4:3*sMDH.ushSamplesInScan/4); % remove oversampling

        mag =   abs(aimg); % magnitude image
        pha = angle(aimg); % phase image
        atitle = ['Slice:' num2str(slice)];
        figure('Name', ['PHASE ' atitle]); imagesc(pha); axis image; colormap(gray);
        figure('Name', ['MAGNITUDE ' atitle]); imagesc(mag); axis image; colormap(gray);
    end
end


