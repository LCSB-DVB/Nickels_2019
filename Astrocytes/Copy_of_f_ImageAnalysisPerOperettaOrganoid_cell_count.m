function  [ObjectsThisOrganoid] = FINAL_f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, ch1, ch2, ch3, ch4, ChannelNames, PreviewPath);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % vol(ch1, 0, 500) % Alexa 488 >>> S100b
    % vol(ch2, 0, 4000) % Alexa 647 >>> GFAP
    % vol(ch3, 0, 5000) % HOECHST 33342 >>> Hoechst imtool(max(ch2, [], 3))
    % vol(ch4, 0, 5000) % TRITC >>> TUJ1

    %% Initialize variables
    NucleiMask = [];
    GFAPMask = [];
    S100bMask = [];
    Tuj1Mask = [];
    
    %% Segment nuclei
    %vol(ch3, 0, 10000)
    ch3BlurSmall = imfilter(double(ch3), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch3BlurSmall)
    ch3BlurBig = imfilter(double(ch3), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch3BlurBig) %%kind of flatfield corraction, to account for different bk in the pic
    ch3DoG = ch3BlurSmall - ch3BlurBig; %vol(ch3DoG, 0, 200, 'hot')
    NucleiMask = ch3DoG > 75; %vol(NucleiMask)
    NucleiMask = bwareaopen(NucleiMask, 20);%vol(NucleiMask)
    ch3LP = imfilter(ch3, fspecial('gaussian', 11, 1), 'symmetric');%vol(ch3LP, 0, 4000, 'hot')
    NucMaskHigh =  (ch3LP > 1500) .* NucleiMask; %vol(NucMaskHigh, 0, 1) % tried 3000 (includes lot of death nuclei)
    NucMaskAlive = NucleiMask & ~NucMaskHigh; % vol(NucMaskAlive)
    
    %% GFAP (ch2)

    %vol(ch2, 0, 2000)
    ch2MedFilt = []; 
    SizeZ = size(ch2, 3);
    parfor p = 1:SizeZ
        ch2MedFilt(:,:,p) = medfilt2(ch2(:,:,p));
    end
    %vol(ch2MedFilt, 0, 4000, 'hot')
    GFAPMask = ch2MedFilt > 1500; %previously 2000
    GFAPDoG = imfilter(ch2, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 31, 10), 'symmetric');
    %vol(GFAPDoG, 0, 300, 'hot')
    GFAPDoGMask = GFAPDoG > 300;
    %vol(GFAPDoGMask, 0, 1)
    GFAPMask = GFAPMask & GFAPDoGMask;
    %vol(GFAPMask, 0, 1)
    GFAPMask = bwareaopen(GFAPMask, 300);

%                                     %%  skeleton3D GFAP
%                                     disp('Start skel')
%                                     tic
%                                     skelGFAP = Skeleton3D(GFAPMask);
%                                     toc
%                                     disp('Skel done')
%                                 %     vol(skelTH, 0, 1)
%                                     [AdjacencyMatrixGFAP, nodeGFAP, linkGFAP] = Skel2Graph3D(skelGFAP,0);                       
%                                     %imtool(AdjacencyMatrixTH, [])
%                                     NodeGFAP = zeros(size(GFAPMask), 'uint8');
%                                     NodeIdxs = vertcat(nodeGFAP(:).idx);
%                                     NodeGFAP(NodeIdxs) = 1;
%                                 %     vol(NodeTH)    
%                                     if size(NodeIdxs, 1) == 0
%                                         return
%                                     end
%                                     NodeGFAPPreview = uint8(skelGFAP) + NodeGFAP + uint8(GFAPMask); 
%                                     NodeGFAPPreview2D = max(NodeGFAPPreview, [], 3);
%                                     %it(NodeTHPreview2D)
%                                 %   vol(NodeTHPreview, 0, 3, 'jet')    
%                                     NodeDegreeVectorGFAP = sum(AdjacencyMatrixGFAP, 1);

%     ZeroNodeExplanationNeeded = 0;
%     if ZeroNodeExplanationNeeded
%         ZeroNodes = find(NodeDegreeVectorGFAP == 0);
%         ZeroNodesLinIdx = vertcat(nodeGFAP(ZeroNodes).idx);
%         ZeroNodeMask = zeros(size(GFAPMaskClipped), 'uint8');
%         ZeroNodeMask(ZeroNodesLinIdx) = 1; %vol(ZeroNodeMask)
%         NodePreviewZeroCase = uint8(skelGFAP) + NodeMaskGFAP + 10*uint8(ZeroNodeMask) + uint8(GFAPMask);
%     end  

%                         %% GFAP Fragmentation
% 
%                         % Define structuring element for surface detection
%                         Conn6 = strel('sphere', 1); % 6 connectivity
%                         % Detect surface
%                         SurfaceGFAP = GFAPMask & ~(imerode(GFAPMask, Conn6));
%                         %vol(SurfaceGFAP)

    %% S100b (ch1)
    
    %vol(ch1, 0, 2000)
    ch1MedFilt = [];
    SizeZ = size(ch1, 3);
    parfor p = 1:size(ch1, 3)
        ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p));
    end
    %vol(ch1MedFilt, 0, 1500, 'hot')    
    S100bMask = ch1MedFilt > 400; 
    S100bMask = bwareaopen(S100bMask, 200);
    %vol(S100bMask)
    
    %% Tuj1 (ch4)         
    ch4MedFilt = [];
    SizeZ = size(ch4, 3);
    for p = 1:SizeZ
        ch4MedFilt(:,:,p) = medfilt2(ch4(:,:,p));
    end
    %vol(ch4MedFilt, 0, 5000)
    ch4Blur = imfilter(ch4MedFilt, fspecial('gaussian', 10, 1), 'symmetric');
    %vol(ch4Blur, 0, 5000, 'hot')
    Tuj1Mask = ch4Blur > 1500;
    %Tuj1Mask = Tuj1Mask & ~ NucleiMask;
    %vol(Tuj1Mask)
    
%                              %% Tuj1 Fragmentation
% 
%                             % Define structuring element for surface detection
%                             Conn7 = strel('sphere', 1); % 6 connectivity
%                             % Detect surface
%                             SurfaceTUJ1 = Tuj1Mask & ~(imerode(Tuj1Mask, Conn7));
%                             %vol(SurfaceTUJ1)

%                          %%  skeleton3D TUJ1
% 
%                         disp('Start skel')
%                         tic
%                         skelTUJ1 = Skeleton3D(Tuj1Mask);
%                         toc
%                         disp('Skel done')
%                     %     vol(skelTH, 0, 1)
%                         [AdjacencyMatrixTUJ1, nodeTUJ1, linkTUJ1] = Skel2Graph3D(Tuj1Mask,0);                       
%                         %imtool(AdjacencyMatrixTH, [])
%                         nodeTUJ1 = zeros(size(Tuj1Mask), 'uint8');
%                         NodeIdxs2 = vertcat(nodeTUJ1(:).idx);
%                         nodeTUJ1(NodeIdxs2) = 1;
%                     %     vol(NodeTH)    
%                         if size(NodeIdxs2, 1) == 0
%                             return
%                         end
%                         NodeTUJ1Preview = uint8(skelTUJ1) + nodeTUJ1 + uint8(Tuj1Mask); 
%                         NodeTUJ1Preview2D = max(NodeTUJ1Preview, [], 3);
%                         %it(NodeTHPreview2D)
%                     %   vol(NodeTHPreview, 0, 3, 'jet')    
%                         NodeDegreeVectorTUJ1 = sum(AdjacencyMatrixTUJ1, 1);

%     clear ZeroNodesLinIdx  ZeroNodeExplanationNeeded  ZeroNodeMask  ZeroNodes  NodePreviewZeroCase
%     ZeroNodeExplanationNeeded = 0;
%     if ZeroNodeExplanationNeeded
%         ZeroNodes = find(NodeDegreeVectorTUJ1 == 0);
%         ZeroNodesLinIdx = vertcat(nodeTUJ1(ZeroNodes).idx);
%         ZeroNodeMask = zeros(size(TUJ1MaskClipped), 'uint8');
%         ZeroNodeMask(ZeroNodesLinIdx) = 1; %vol(ZeroNodeMask)
%         NodePreviewZeroCase = uint8(skelTUJ1) + NodeMaskTUJ1 + 10*uint8(ZeroNodeMask) + uint8(Tuj1Mask);
%     end  
    
    
    %% Double Mask S100b GFAP
    DoubleMasks = S100bMask & GFAPMask; %vol(DoubleMasks)

    %% Perinuclear Volume (to detect amount of cells positive for GFAP)
    
    %vol(NucleiMask)
    NucleiMaskSingleCells = f_RemoveBigObjects (NucleiMask, 10000); 
    NucDil = imdilate(imdilate(NucleiMaskSingleCells, strel('disk', 4)), strel('sphere',1));
    NucPerim = logical(NucDil) & ~logical(NucleiMaskSingleCells);
    %vol(NucPerim)
    GFAPMaskinNucPerim = GFAPMask & NucPerim;% vol(THMaskinNucPerim)
    
    
     %% Perinuclear Volume (to detect amount of cells positive for TUJ1)
    
    %vol(NucleiMask)
    NucleiMaskSingleCells = f_RemoveBigObjects (NucleiMask, 10000); 
    NucDil = imdilate(imdilate(NucleiMaskSingleCells, strel('disk', 4)), strel('sphere',1));
    NucPerim = logical(NucDil) & ~logical(NucleiMaskSingleCells);
    %vol(NucPerim)
    TUJ1MaskinNucPerim = Tuj1Mask & NucPerim;% vol(THMaskinNucPerim)
    
    
    %% Percent GFAP pos
    %split perinuc
    D = bwdist(NucleiMaskSingleCells);
    %vol(D, 0, 20, 'hot')
    %it(D(:,:,1))
    disp('start watershed')
    tic
    W = watershed(D);
    toc
    disp('watershed done')
    %vol(W)
    NucPerimStencil = uint16(W) .* uint16(imreconstruct(logical(imdilate(NucPerim, strel('disk', 1))), logical(NucleiMaskSingleCells))); % This line was causing the error 20171207 % Function imreconstruct expected MARKER and MASK to have the same class.
    %vol(NucPerimStencil)
    %vol(NucPerim)
    %vol(NucleiMaskSingleCells)
    %toto = imreconstruct(logical(NucPerim), logical(NucleiMaskSingleCells));
       
    PeriNucMask = logical(NucPerimStencil);
    PeriNucMask = bwareaopen(PeriNucMask, 500);
%     vol(PeriNucMask)
    PerinucLM = bwlabeln(PeriNucMask);
    
    % GFAP 
    PeriNucObjectsGFAP = regionprops('table', PerinucLM, double(GFAPMask), 'PixelValues');
    GFAPproportions = rowfun(@(x) sum(x{:})/length(x{:}), PeriNucObjectsGFAP, 'InputVariables', 'PixelValues');
    GFAPPos = array2table(table2array(GFAPproportions) > 0.01);
    GFAPPos.Properties.VariableNames(end) = {'GFAPpos'};
    PeriNucObjectsGFAP = [PeriNucObjectsGFAP, GFAPproportions, GFAPPos];
    PeriNucObjectsGFAP.Properties.VariableNames(end-1) = {'GFAPproportion'};
    PeriNucObjectsCompactGFAP = PeriNucObjectsGFAP(:, {'GFAPproportion','GFAPpos'});
    GFAPPercent = (sum(PeriNucObjectsCompactGFAP.GFAPpos)/height(PeriNucObjectsCompactGFAP))*100;
    % S100B 
    PeriNucObjectsS100b = regionprops('table', PerinucLM, double(S100bMask), 'PixelValues');
    S100bproportions = rowfun(@(x) sum(x{:})/length(x{:}), PeriNucObjectsS100b, 'InputVariables', 'PixelValues');
    S100bPos = array2table(table2array(S100bproportions) > 0.01);
    S100bPos.Properties.VariableNames(end) = {'S100bpos'};
    PeriNucObjectsS100b = [PeriNucObjectsS100b, S100bproportions, S100bPos];
    PeriNucObjectsS100b.Properties.VariableNames(end-1) = {'S100bproportion'};
    PeriNucObjectsCompactS100b = PeriNucObjectsS100b(:, {'S100bproportion','S100bpos'});
    S100bPercent = (sum(PeriNucObjectsCompactS100b.S100bpos)/height(PeriNucObjectsCompactS100b))*100;
    
%     PeriNucObjectsDoubleMasks = regionprops('table', PerinucLM, double(DoubleMasks), 'PixelValues');
%     DoubleMasksproportions = rowfun(@(x) sum(x{:})/length(x{:}), PeriNucObjectsDoubleMasks, 'InputVariables', 'PixelValues');
%     DoubleMasksPos = array2table(table2array(DoubleMasksproportions) > 0.01);
%     DoubleMasksPos.Properties.VariableNames(end) = {'DoubleMasks'};
%     PeriNucObjectsDoubleMasks = [PeriNucObjectsDoubleMasks, DoubleMasksproportions, DoubleMasksPos];
%     PeriNucObjectsDoubleMasks.Properties.VariableNames(end-1) = {'DoubleMasksproportion'};
%     PeriNucObjectsCompactDoubleMasks = PeriNucObjectsDoubleMasks(:, {'DoubleMasksproportion','DoubleMaskspos'});
%     DoubleMasksPercent = (sum(PeriNucObjectsCompactDoubleMasks.DoubleMaskspos)/height(PeriNucObjectsCompactDoubleMasks))*100;

    % TUJ1
    PeriNucObjectsTUJ1 = regionprops('table', PerinucLM, double(Tuj1Mask), 'PixelValues');
    TUJ1proportions = rowfun(@(x) sum(x{:})/length(x{:}), PeriNucObjectsTUJ1, 'InputVariables', 'PixelValues');
    TUJ1Pos = array2table(table2array(TUJ1proportions) > 0.01);
    TUJ1Pos.Properties.VariableNames(end) = {'TUJ1Pos'};
    PeriNucObjectsTUJ1 = [PeriNucObjectsTUJ1, TUJ1proportions, TUJ1Pos];
    PeriNucObjectsTUJ1.Properties.VariableNames(end-1) = {'TUJ1proportions'};
    PeriNucObjectsCompactTUJ1 = PeriNucObjectsTUJ1(:, {'TUJ1proportions','TUJ1Pos'});
    TUJ1Percent = (sum(PeriNucObjectsCompactTUJ1.TUJ1Pos)/height(PeriNucObjectsCompactTUJ1))*100;
    
     
    %% Previews 
    
    % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    %it(BarMask)

    PreviewGFAP = imoverlay2(imadjust(max(ch2,[],3),[0 0.07]), bwperim(max(GFAPMask,[],3)), [1 0 0]);
    PreviewGFAP = imoverlay2(PreviewGFAP, BarMask, [1 1 1]);
    %imtool(PreviewGFAP)
    
    PreviewHoechst = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewHoechst = imoverlay2(PreviewHoechst, BarMask, [1 1 1]);
    % imtool(PreviewHoechst)
        
    PreviewS100b = imoverlay2(imadjust(max(ch1, [], 3), [0 0.02]), bwperim(max(S100bMask,[],3)), [0 0 1]);
    PreviewS100b = imoverlay2(PreviewS100b, BarMask, [1 1 1]);
    %imtool(PreviewS100b)
    
    PreviewTuj1 = imoverlay2(imadjust(max(ch4, [], 3), [0 0.07]), bwperim(max(Tuj1Mask,[],3)), [0 0 1]);
    PreviewTuj1 = imoverlay2(PreviewTuj1, BarMask, [1 1 1]);
    %imtool(PreviewTuj1)
    
    PreviewDoubleMask = imoverlay2(imadjust(max(ch1, [], 3), [0 0.07]), bwperim(max(DoubleMasks,[],3)), [0 0 1]);
    PreviewDoubleMask = imoverlay2(PreviewDoubleMask, BarMask, [1 1 1]);
    %imtool(PreviewTuj1)
    
     
    PreviewNucMaskAlive = imoverlay2(imadjust(max(ch3, [], 3), [0 0.07]), bwperim(max(NucMaskAlive,[],3)), [0 0 1]);
    PreviewNucMaskAlive = imoverlay2(PreviewNucMaskAlive, BarMask, [1 1 1]);
    %imtool(PreviewTuj1)
    
    PreviewNucMaskHigh = imoverlay2(imadjust(max(ch3, [], 3), [0 0.07]), bwperim(max(NucMaskHigh,[],3)), [0 0 1]);
    PreviewNucMaskHigh = imoverlay2(PreviewNucMaskHigh, BarMask, [1 1 1]);
    %imtool(PreviewTuj1)
    
    
    
    %imwrite(PreviewTuj1, [PreviewPath, filesep, CellLine, '_', 'Tuj1', '_', num2str(s), '.png'])
    IdentityString = [Label.AreaName{:}, '_Idx_', num2str(Label.Idx)];
    imwrite(PreviewGFAP, [PreviewPath, filesep, IdentityString, '_', 'GFAP', '.png'])
    imwrite(PreviewHoechst, [PreviewPath, filesep, IdentityString, '_', 'Hoechst', '.png'])
    imwrite(PreviewS100b, [PreviewPath, filesep, IdentityString, '_', 'S100b', '.png'])
    imwrite(PreviewTuj1, [PreviewPath, filesep, IdentityString, '_', 'Tuj1', '.png'])
    imwrite(PreviewDoubleMask, [PreviewPath, filesep, IdentityString, '_', 'PreviewDoubleMaskGFAPS100b', '.png'])
    imwrite(PreviewNucMaskAlive, [PreviewPath, filesep, IdentityString, '_', 'PreviewNucMaskAlive', '.png'])
    imwrite(PreviewNucMaskHigh, [PreviewPath, filesep, IdentityString, '_', 'PreviewNucMaskHigh', '.png'])
   
    
    %% Feature extraction
    
    ObjectsThisOrganoid = table();
    ObjectsThisOrganoid.LabelIdx = {Label.Idx};
    ObjectsThisOrganoid.AreaName = {Label.AreaName};
    ObjectsThisOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsThisOrganoid.Tuj1MaskSum = sum(Tuj1Mask(:));
    ObjectsThisOrganoid.GFAPMaskSum = sum(GFAPMask(:));
    ObjectsThisOrganoid.S100bMaskSum = sum(S100bMask(:));
    ObjectsThisOrganoid.Tuj1ByNuc = sum(Tuj1Mask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.S100bMaskNuc = sum(S100bMask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.GFAPMaskByNuc = sum(GFAPMask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.NucMaskHigh = sum(NucMaskHigh(:));
    
                %     ObjectsThisOrganoid.GFAPFragmentation = sum(SurfaceGFAP(:)) / sum(GFAPMask(:));
                %     ObjectsThisOrganoid.SkelGFAP = sum(skelGFAP(:));
                %     ObjectsThisOrganoid.NodesGFAP = size(nodeGFAP, 2);
                %     ObjectsThisOrganoid.LinksGFAP = size(linkGFAP, 2);
                %     
                %     ObjectsThisOrganoid.TUJ1Fragmentation = sum(SurfaceTUJ1(:)) / sum(Tuj1Mask(:));
                %     ObjectsThisOrganoid.SkelTUJ1 = sum(skelTUJ1(:));
                %     ObjectsThisOrganoid.NodesTUJ1 = size(nodeTUJ1, 2);
                %     ObjectsThisOrganoid.LinksTUJ1 = size(linkTUJ1, 2);
                %     
    ObjectsThisOrganoid.GFAPPercent = GFAPPercent;
%%ObjectsThisOrganoid.DoubleMasksPercent = DoubleMasksPercent;
    ObjectsThisOrganoid.S100bPercent = S100bPercent; 
    ObjectsThisOrganoid.Tuj1Percent = TUJ1Percent;
    

end

