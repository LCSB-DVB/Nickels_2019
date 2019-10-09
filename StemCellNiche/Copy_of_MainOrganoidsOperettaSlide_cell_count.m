%% User inputs
%% addapted from Silvia 

clear
clc

SetupMode = 0; % 1 for creating numeric organoid labels OR 0 for linking the final analysis to human labels

%%%%%%% text file defining the samples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slide 1 NO
SlideLayout = 'JM_20181130_K7A6n1n2n3_T129A6n1n2n3_lmx1a_568_sox2_647_ki67_488.txt';
% slide 2 NO
%SlideLayout = 'JM_20181201_CORb9n1n2n3_CORC3n1n2n3_CORA6n1n2n3_lmx1a_568_sox2_6.txt';
% slide 3 
%SlideLayout = 'JM_20181130_K7C3n1n2n3_T129C3n1n2n3_lmx1a_568_sox2_647_ki67_488.txt';
% slide 4 
%SlideLayout = 'JM_20181201_K7b9n1n2n3_T129b9n1n2n3_lmx1a_568_sox2_647_ki67_488.txt';


% save path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% slide 1 
SavePath = 'S:\HCS_Platform\Data\JenniferModamio\optimisationProtocol\JM_20181130_K7A6n1n2n3_T129A6n1n2n3_lmx1a_568_sox2_647_ki67_488';
% slide 2 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\optimisationProtocol\JM_20181201_CORb9n1n2n3_CORC3n1n2n3_CORA6n1n2n3_lmx1a_568_sox2_6';
% slide 3 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\optimisationProtocol\JM_20181130_K7C3n1n2n3_T129C3n1n2n3_lmx1a_568_sox2_647_ki67_488';
% slide 4 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\optimisationProtocol\JM_20181201_K7b9n1n2n3_T129b9n1n2n3_lmx1a_568_sox2_647_ki67_488';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mkdir(SavePath)
PreviewSavePath = [SavePath, filesep, 'Previews'];
mkdir(PreviewSavePath)

%% Parallel pool control
delete(gcp('nocreate'))
myCluster = parcluster;
Workers = myCluster.NumWorkers;
% parpool(28) %for HURRICANE
parpool(Workers) % for MEGATRON

%% Run mode control
if SetupMode
    RunMode = 0;
else
    RunMode = 1;
end

%% Common part

channelID = 1; % channel to use for overview

% Input information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% slide 1
InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_20190115_K7A6n1n2n3_T129A6n1n2n3_lmx1a_568_sox2_647_ki67_488\a01ea58c-1b45-45f2-b4ad-fef88459a0a7\metadata.csv');
% slide 2
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_20190116_CORb9n1n2n3_CORC3n1n2n3_CORA6n1n2n3_lmx1a_568_sox2_6\ec6e13d0-29ed-4f41-9663-a21a9ece664f\metadata.csv');
% slide 3
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_20190116_K7C3n1n2n3_T129C3n1n2n3_lmx1a_568_sox2_647_ki67_488\b9355986-5626-4c30-ba86-b254ad698925\metadata.csv');
% slide 4
% OLD still 
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_20181201_K7b9n1n2n3_T129b9n1n2n3_lmx1a_568_sox2_647_ki67_488\e836648c-dd6e-45d1-a5c6-4a5f774d85da\metadata.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ChannelNames = unique(InfoTable.Channel);
    Channels = length(unique(InfoTable.Channel));
    Planes = unique(InfoTable.Plane)';
    Timepoints = unique(InfoTable.Timepoint)' + 1;
    [GroupsTable, GroupsIm5DCellArray] = FindGroups(InfoTable); % it(GroupsTable)
    
%% Setup mode

if SetupMode == 1
    
    Preview = CreateLabelHelpPreview(GroupsTable, PreviewSavePath);
    imwrite(Preview, [PreviewSavePath, filesep, 'layout.png'])
    Message = ['The plate Layout has been saved at ', [PreviewSavePath, filesep, 'layout.png'], '. Please save a text file without header, using tab separation, where the first column is the index number as shown in the preview and the second is the area name. Save the text file as SlideLayout_Date.txt in your working directory and set the variable SetupMode to 0 >>> Run'];
    h = msgbox(Message);
    
else    
%% Analysis mode
    
    % Load annotations
    Layout = readtable(SlideLayout)
    Layout.Properties.VariableNames = {'Idx', 'AreaName'};
    
    % Load images and organize in an XYC array
    Groups = unique(GroupsTable(GroupsTable > 0))';
    GroupPreviews = {};
    ObjectsAll = {};


    for g = Groups % Number of organoids
        XYMosaicCells = {};
        GroupZone = GroupsTable == g;
        [GroupIdxRowVec, GroupIdxColVec] = find(GroupZone); % linear indexes
        Elements = sum(GroupZone(:));
        InfoTablesThisGroup = {};
        for e = 1:Elements % Fields of a given organoid
            for c = 1:Channels
            InfoTableThisField = GroupsIm5DCellArray{GroupIdxRowVec(e), GroupIdxColVec(e)};
            InfoTablesThisGroup{e} = InfoTableThisField;
            InfoTableThisChannel = InfoTableThisField(strcmp(InfoTableThisField.Channel, ChannelNames{c}), :);
                clear Im4D
                for t = Timepoints
                    for p = Planes
                        InfoTableThisChannelThisPlane = InfoTableThisChannel(InfoTableThisChannel.Plane == p, :);
                        ImPathThisPlane = InfoTableThisChannelThisPlane.Path{:};   
                        Im4D(:,:,t,p) = imread(ImPathThisPlane); % it(Im4D(:,:,t,p))
                    end
                end
               XYMosaicCells{c}{GroupIdxRowVec(e), GroupIdxColVec(e)} = Im4D; % Operetta counterpart of XYmosaicCells for Opera
            end
        end

        InfoTableThisGroup = vertcat(InfoTablesThisGroup{:});

        %% Remove empty cells
        XYMosaicCells = cellfun(@(x) GroupClipper(x),  XYMosaicCells, 'UniformOutput', false);

        %% Stitch
        XYmosaicContourCell = cellfun(@(x) stdfilt(x, ones(3)), XYMosaicCells{1}, 'UniformOutput', false);
        XPositions = unique(InfoTableThisGroup.PositionX); % m
        YPositions = unique(InfoTableThisGroup.PositionY); % m
        ResolutionXY = 675 / 1360; % um per pixel
        MaxPixelDrift = 30;
        PreviewChannel = 1;
        ShowProgress = 0;
        [CroppedMosaic, StitchedIm] = f_stitching_operetta(XYMosaicCells, XYmosaicContourCell, XPositions, YPositions, ResolutionXY, MaxPixelDrift, PreviewChannel, ShowProgress);
        GroupPreviews{g} = max(CroppedMosaic{channelID},[],3); %it(GroupPreviews{g})
        
        %% Image analysis
        Label = Layout(g,:);
         try
            ObjectsThisOrganoid = Copy_of_f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, CroppedMosaic{1}, CroppedMosaic{2}, CroppedMosaic{3}, CroppedMosaic{4}, ChannelNames, PreviewSavePath);
            %ObjectsThisOrganoid = Copy_of_f_ImageAnalysisPerOperettaOrganoid_20171129(Label, CroppedMosaic{1}(1:500, 1:500, 10:30), CroppedMosaic{2}(1:500, 1:500, 10:30), CroppedMosaic{3}(1:500, 1:500, 10:30), ChannelNames, PreviewSavePath);
        catch
             Errors{g} = 'Image analysis failed';
             continue % next group g
         end
            ObjectsAll{g} = ObjectsThisOrganoid;

    end
    
    Objects = vertcat(ObjectsAll{:});
    save([SavePath, filesep, 'Objects.mat'], 'Objects');
    writetable(Objects, [SavePath, filesep, 'Objects.csv'])
    writetable(Objects, [SavePath, filesep, 'Objects.xlsx'])

    %% Preview of the whole slide

    SizeSingleIm = size(XYMosaicCells{1}{1,1});
    SizeSingleIm = SizeSingleIm(1:2);
    RescanGridSize = size(GroupsTable);
    GreatPreview = zeros(SizeSingleIm(1)*RescanGridSize(1), SizeSingleIm(2)*RescanGridSize(2), 'uint16');
    ImHeight = SizeSingleIm(1);
    ImWidth = SizeSingleIm(2);
    StartRCell = {};
    StartCCell = {};

    for g = Groups
        g
        StitchedGroupSize = size(GroupPreviews{g});
        ZoneNow = GroupsTable == g;
        [R,C] = find(ZoneNow)
        StartR = min(R);
        StartC = min(C);
        StartRPixel = ((StartR-1) * ImHeight) + 1;
        %EndRPixel = StartRPixel + (3 * ImHeight) - 1;
        EndRPixel = StartRPixel + StitchedGroupSize(1) - 1;
        StartCPixel = ((StartC-1) * ImWidth) + 1;
        %EndCPixel = StartCPixel + (3 * ImWidth) - 1;
        EndCPixel = StartCPixel + StitchedGroupSize(2) - 1;
        GreatPreview(StartRPixel:EndRPixel, StartCPixel:EndCPixel) = GroupPreviews{g};
        StartRCell{g} = StartRPixel;
        StartCCell{g} = StartCPixel;
    end

    Zoomfactor = 50;
    %GreatPreviewResized = imresize(imadjust(GreatPreview), 1/Zoomfactor);
    GreatPreviewResized = imresize(imadjust(GreatPreview, [0 0.02], [0 1]), 1/Zoomfactor);

    for g = Groups
        GreatPreviewResized = insertText(GreatPreviewResized, [round(StartCCell{g}/Zoomfactor), round(StartRCell{g}/Zoomfactor)], num2str(g), 'FontSize', 12, 'BoxColor', 'red', 'TextColor', 'white');
    end
    
    %imtool(GreatPreview)
    %imtool(GreatPreviewResized)
    imwrite(GreatPreviewResized, [SavePath, filesep, 'GreatPreview.png'])
    
    % it(GreatPreviewResized)
    % save([SavePath, filesep, 'WorkspaceIncludingObjects.mat'])
    
end




