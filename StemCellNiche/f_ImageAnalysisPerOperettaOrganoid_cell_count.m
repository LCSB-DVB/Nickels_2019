function  [ObjectsThisOrganoid] = FINAL_f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, ch1, ch2, ch3, ch4, ChannelNames, PreviewPath);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % vol(ch1, 0, 500) % Alexa 488 >>> Ki67
    % vol(ch2, 0, 4000) % Alexa 647 >>> SOX2
    % vol(ch3, 0, 5000) % HOECHST 33342 >>> Hoechst imtool(max(ch2, [], 3))
    % vol(ch4, 0, 5000) % TRITC >>> LMX1a

    %% Initialize variables
    NucleiMask = [];
    Ki67Mask = [];
    SOX2Mask = [];
    LMX1aMask = [];
    
    %% Segment nuclei
    %vol(ch3, 0, 1000)
    ch3BlurSmall = imfilter(double(ch3), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch3BlurSmall)
    ch3BlurBig = imfilter(double(ch3), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch3BlurBig) %%Ki67nd of flatfield corraction, to account for different bk in the pic
    ch3DoG = ch3BlurSmall - ch3BlurBig; %vol(ch3DoG, 0, 200, 'hot')
    NucleiMask = ch3DoG > 75; %vol(NucleiMask)
    NucleiMask = bwareaopen(NucleiMask, 20);%vol(NucleiMask)
    ch3LP = imfilter(ch3, fspecial('gaussian', 11, 1), 'symmetric');%vol(ch3LP, 0, 4000, 'hot')
    NucMaskHigh =  (ch3LP > 2000) .* NucleiMask; %vol(NucMaskHigh, 0, 1) % tried 3000 (includes lot of death nuclei)
    NucMaskAlive = NucleiMask & ~NucMaskHigh; % vol(NucMaskAlive)
    
    %% Ki67 (ch1)
    
     %vol(ch2, 0, 2000)
    ch1MedFilt = [];
    SizeZ = size(ch1, 3);
    parfor p = 1:size(ch1, 3)
        ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p));
    end
    %vol(ch1MedFilt, 0, 1500, 'hot')    
    Ki67Mask = ch1MedFilt > 200; 
    Ki67Mask = bwareaopen(Ki67Mask, 200);
    
    Ki67NucleiMask = Ki67Mask & NucleiMask;
                    
    %% SOX2 (ch2)
    
    %vol(ch2, 0, 500)
    ch2MedFilt = [];
    SizeZ = size(ch2, 3);
    parfor p = 1:size(ch2, 3)
        ch2MedFilt(:,:,p) = medfilt2(ch2(:,:,p));
    end
    %vol(ch1MedFilt, 0, 1500, 'hot')    
    SOX2Mask = ch2MedFilt > 200; 
    SOX2Mask = bwareaopen(SOX2Mask, 200);
    
    SOX2NucleiMask = SOX2Mask & NucleiMask;
    %vol(SOX2Mask)   

    %% LMX1a (ch4)         
    
    ch4MedFilt = [];
    SizeZ = size(ch4, 3);
    for p = 1:SizeZ
        ch4MedFilt(:,:,p) = medfilt2(ch4(:,:,p));
    end
    %vol(ch4MedFilt, 0, 5000)
    ch4Blur = imfilter(ch4MedFilt, fspecial('gaussian', 10, 1), 'symmetric');
    %vol(ch4Blur, 0, 5000, 'hot')
    LMX1aMask = ch4Blur > 700;
    
    LMX1aNucleiMask = LMX1aMask & NucleiMask;
    %vol(LMX1aMask)
    
    %% Double Mask Ki6767 and SOX22
    DoubleMasks = SOX2Mask & Ki67Mask; %vol(DoubleMasks)
       
    %% Previews 
    
    % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    %it(BarMask)

    PreviewKi67 = imoverlay2(imadjust(max(ch1,[],3),[0 0.07]), bwperim(max(Ki67Mask,[],3)), [1 0 0]);
    PreviewKi67 = imoverlay2(PreviewKi67, BarMask, [1 1 1]);
    %imtool(PreviewKi67)
    
    PreviewHoechst = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewHoechst = imoverlay2(PreviewHoechst, BarMask, [1 1 1]);
    % imtool(PreviewHoechst)
        
    PreviewSOX2 = imoverlay2(imadjust(max(ch2, [], 3), [0 0.02]), bwperim(max(SOX2Mask,[],3)), [0 0 1]);
    PreviewSOX2 = imoverlay2(PreviewSOX2, BarMask, [1 1 1]);
    %imtool(PreviewSOX2)
    
    PreviewLMX1a = imoverlay2(imadjust(max(ch4, [], 3), [0 0.07]), bwperim(max(LMX1aMask,[],3)), [0 0 1]);
    PreviewLMX1a = imoverlay2(PreviewLMX1a, BarMask, [1 1 1]);
    %imtool(PreviewLMX1a)
    
    PreviewDoubleMask = imoverlay2(imadjust(max(ch1, [], 3), [0 0.07]), bwperim(max(DoubleMasks,[],3)), [0 0 1]);
    PreviewDoubleMask = imoverlay2(PreviewDoubleMask, BarMask, [1 1 1]);
    %imtool(PreviewLMX1a)
     
    PreviewNucMaskAlive = imoverlay2(imadjust(max(ch3, [], 3), [0 0.07]), bwperim(max(NucMaskAlive,[],3)), [0 0 1]);
    PreviewNucMaskAlive = imoverlay2(PreviewNucMaskAlive, BarMask, [1 1 1]);
    %imtool(PreviewLMX1a)
    
    PreviewNucMaskHigh = imoverlay2(imadjust(max(ch3, [], 3), [0 0.07]), bwperim(max(NucMaskHigh,[],3)), [0 0 1]);
    PreviewNucMaskHigh = imoverlay2(PreviewNucMaskHigh, BarMask, [1 1 1]);
    %imtool(PreviewLMX1a)
    
    PreviewSox2Nuclei = imoverlay2(imadjust(max(ch2, [], 3), [0 0.07]), bwperim(max(SOX2NucleiMask,[],3)), [0 0 1]);
    PreviewSox2Nuclei = imoverlay2(PreviewSox2Nuclei, BarMask, [1 1 1]);
%     %imtool(PreviewLMX1a)
%     
    PreviewKi67Nuclei = imoverlay2(imadjust(max(ch1, [], 3), [0 0.07]), bwperim(max(Ki67NucleiMask,[],3)), [0 0 1]);
    PreviewKi67Nuclei = imoverlay2(PreviewKi67Nuclei, BarMask, [1 1 1]);
%     %imtool(PreviewLMX1a)
%     
    PreviewLMX1aNuclei = imoverlay2(imadjust(max(ch4, [], 3), [0 0.07]), bwperim(max(LMX1aNucleiMask,[],3)), [0 0 1]);
    PreviewLMX1aNuclei = imoverlay2(PreviewLMX1aNuclei, BarMask, [1 1 1]);
%     %imtool(PreviewLMX1a)

         
    %imwrite(PreviewLMX1a, [PreviewPath, filesep, CellLine, '_', 'LMX1a', '_', num2str(s), '.png'])
    IdentityString = [Label.AreaName{:}, '_Idx_', num2str(Label.Idx)];
    imwrite(PreviewKi67, [PreviewPath, filesep, IdentityString, '_', 'Ki67', '.png'])
    imwrite(PreviewHoechst, [PreviewPath, filesep, IdentityString, '_', 'Hoechst', '.png'])
    imwrite(PreviewSOX2, [PreviewPath, filesep, IdentityString, '_', 'SOX2', '.png'])
    imwrite(PreviewLMX1a, [PreviewPath, filesep, IdentityString, '_', 'LMX1a', '.png'])
    imwrite(PreviewDoubleMask, [PreviewPath, filesep, IdentityString, '_', 'PreviewDoubleMaskKi67SOX2', '.png'])
    imwrite(PreviewNucMaskAlive, [PreviewPath, filesep, IdentityString, '_', 'PreviewNucMaskAlive', '.png'])
    imwrite(PreviewNucMaskHigh, [PreviewPath, filesep, IdentityString, '_', 'PreviewNucMaskHigh', '.png'])
    
    imwrite(PreviewSox2Nuclei, [PreviewPath, filesep, IdentityString, '_', 'PreviewSox2Nuclei', '.png'])
    imwrite(PreviewKi67Nuclei, [PreviewPath, filesep, IdentityString, '_', 'PreviewKi67Nuclei', '.png'])
    imwrite(PreviewLMX1aNuclei, [PreviewPath, filesep, IdentityString, '_', 'PreviewLMX1aNuclei', '.png'])
   
    
    %% Feature extraction
    
    ObjectsThisOrganoid = table();
    ObjectsThisOrganoid.LabelIdx = {Label.Idx};
    ObjectsThisOrganoid.AreaName = {Label.AreaName};
    ObjectsThisOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsThisOrganoid.LMX1aMaskSum = sum(LMX1aMask(:));
    ObjectsThisOrganoid.Ki67MaskSum = sum(Ki67Mask(:));
    ObjectsThisOrganoid.SOX2MaskSum = sum(SOX2Mask(:));
    ObjectsThisOrganoid.NucMaskHigh = sum(NucMaskHigh(:));
    ObjectsThisOrganoid.NucMaskAlive = sum(NucMaskAlive(:));
    
    ObjectsThisOrganoid.SOX2NucleiMask = sum(SOX2NucleiMask(:));
    ObjectsThisOrganoid.Ki67NucleiMask = sum(Ki67NucleiMask(:));
    ObjectsThisOrganoid.LMX1aNucleiMask = sum(LMX1aNucleiMask(:));
    
    
end

