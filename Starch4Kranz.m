% Final vein_density function, will involve splitting image in addition to
% enhancing the contrast of said image

% Will return the vein density, if the image was split or not and the blur
% factor applied

function [vein_density] = Starch4Kranz(filename, trim_factor, pixel_length_um, y_pixels, x_pixels, split_mode, shave, initial_blur_factor, branch_threshold, endpoint_threshold, trim_threshold, mono_min_param, monocot_branch_length_to_keep)

if nargin < 6
    split_mode = "auto";
    shave = 40;
    initial_blur_factor=20;
    branch_threshold=125;
    endpoint_threshold=175;
    sd_threshold=75;
    trim_threshold=0.01;
    mono_min_param=30;
    monocot_branch_length_to_keep=200;
  end





% Define image dimensions
pixel_area_um2 = pixel_length_um*pixel_length_um;
pixel_area_mm2 = pixel_area_um2/1000000;

image_area_mm2 = (y_pixels*x_pixels)*(pixel_area_mm2);

X = imread(filename);
% Convert to grayscale
X_gray = rgb2gray(X);
sd = std2(X_gray);

if (strcmp(split_mode, "auto")) || (strcmp(split_mode, "split"))


    % Split image
    upperLeftX = X_gray(1:(y_pixels/2), 1:(x_pixels/2));
    upperRightX = X_gray(1:(y_pixels/2), (x_pixels/2):x_pixels);
    lowerLeftX = X_gray((y_pixels/2):y_pixels, 1:(x_pixels/2));
    lowerRightX = X_gray((y_pixels/2):y_pixels, (x_pixels/2):x_pixels);

    images{1} = upperLeftX;
    images{2} = upperRightX;
    images{3} = lowerLeftX;
    images{4} = lowerRightX;

    % Define blur factor based on image sd - well stained images will have
    % larger variance hence larger sd


    for i = 1:4
        % blur based on sd
        sd = std2(images{i});
        m = mean2(images{i});
        bf_scale1 = (40*(sd/100));
        bf_scale2 = (100/(m*2));
        bf = (bf_scale1*bf_scale2) + 5;

        bf_split{i} = bf;


        % Improve contrast
        X_con = adapthisteq(images{i});
        % Blur
        X_blur = blur(X_con, bf);
        % Binarize
        X_bin = imbinarize(X_blur);
        % Invert
        X_inv = imcomplement(X_bin);
        % Skeletonize image
        X_skel = bwskel(X_inv);

        % Create skeleton mask
        B = bwmorph(X_skel, 'branchpoints');
        E = bwmorph(X_skel, 'endpoints');
        [y,x] = find(E); %Now the numel(x) == number of endpoints
        B_loc = find(B); %numel(B_loc) == number of branchpoints

        % Initial mask
        Dmask = false(size(X_skel));

        % Trim the mask
        for k = 1:numel(x)
            %BWDISTGEODESIC starts at a given point, and then calculates the
            %distance as you travel along the path.
            D = bwdistgeodesic(X_skel,x(k),y(k));
            %Start at a endpoint, start walking and find all pixels that are closer
            %than the nearest branchpoint.
            distanceToBranchPt = min(D(B_loc)); %D(B_loc) gives the lengths of each branchpoint
            if isempty(distanceToBranchPt) % if there are no branch points do nothing
            else
                if distanceToBranchPt > trim_factor %If loose branches are larger than our universal threshold we do nothing with them
                else
                    Dmask(D < distanceToBranchPt) =true; %Else we remove them
                end
            end
        end
        % Redefine our skeletonised image
        skelD = X_skel - Dmask;

        % Initial estimate of VD - will be final estimate if not trimmed based
        % on mode
        vd = (sum(skelD(:))*(pixel_length_um))/image_area_mm2;
        % Calculate the ratio of branchpoints:Initial_VD
        ratio = (numel(B_loc)/vd);

        trim = num2str(trim_factor);

        % If this ratio is too high (>0.01), then trimming will occur again,
        % based on the mode length of loose branch ends instead of the
        % "trim_factor". If we did mode from the begining then we would risk
        % removing data from those with low VD.

        if ratio > 0.01
            trim = "mode";
            % Here we calculate the length of each "loose branch end" i.e. each
            % line that potentially could be trimmed
            loose_size = {}; %Empty cell array for each "loose length"
            for k = 1:numel(x)
                D = bwdistgeodesic(X_skel,x(k),y(k));
                distanceToBranchPt = min(D(B_loc));
                loose_size = [loose_size, distanceToBranchPt]; %#ok<*AGROW> (surpresses warning)
            end
            % Remove infinate values, otherwise we cannot calculate mode
            loose_size( cellfun( @(C) any(isnan(C) | isinf(C)), loose_size ) ) = [];

            % Threshold is based on the mode length of each loose strand - based on
            % fact that the most common size will be the smaller threads because
            % these caputre noise that is ~similar all over each image
            T = (mode(cell2mat(loose_size))+20); % we add 20 pixels just to make sure all are captured

            % We then create a mask containing our skeletonised image; this will be
            % what is taken away from our final image after it is adjusted for
            % loose branches
            Dmask = false(size(X_skel));

            % Adjusting the mask by removing strands smaller than our threshold
            for k = 1:numel(x)
                D = bwdistgeodesic(X_skel,x(k),y(k));
                distanceToBranchPt = min(D(B_loc));
                if isempty(distanceToBranchPt)
                else
                    if distanceToBranchPt > T % T = threshold  % Do nothing
                    else
                        Dmask(D < distanceToBranchPt) =true;
                    end
                end
            end

            %Skeleton after trimming
            skelD = X_skel - Dmask;

            % Re-calculate vein_density
            vd = (sum(skelD(:))*(pixel_length_um))/image_area_mm2;
        end
        skel{i} = skelD;
        vein{i} = vd;
    end

    bf = mean(cell2mat(bf_split));

    %%%% Now we need to "stitch" the skeletons together

    % Before doing so we need to convert, from white to black, the pixels at
    % the edge of some of the skeletons so that we don't capture the same bit
    % twice

    % The edges to convert - lets go clockwise
    % upperLeft (skel{1}) - remove pixels along right edge that are more
    % vertical
    % upperRight (skel{2}) - remove pixels along bottom edge that are more
    % horizontal
    % lowerRight (skel{4}) - remove pixels along left edge that are more
    % vertical
    % lowerLeft (skel{3}) - remove pixels along upper edge that are more
    % horizontal

    % Define start of new edge
    UL_edge = (x_pixels/2)-shave;
    UR_edge = (y_pixels/2)-shave;
    LR_edge = shave;
    LL_edge = shave;

    % Crop out the edge we wish to alter
    UL_crop = skel{1}( :,UL_edge:(x_pixels/2));
    UR_crop = skel{2}(UR_edge:(y_pixels/2), :);
    LR_crop = skel{4}( :,1:LR_edge);
    LL_crop = skel{3}(1:LL_edge, :);

    % Define matrix patterns that we are interested in not removing from the
    % final edge

    % Firstly the masks to use for filtering the vertical crops ones (can be used for UL and LR) - these are masks that
    % will be kept, hence when we filter the crops we want to keep lines that
    % are horizontal

    sev1 = [0 0 0 0
        0 0 0 0
        1 1 1 1
        0 0 0 0
        0 0 0 0];

    sev2 = [0 0 0 0
        0 0 0 0
        1 1 1 0
        0 0 0 1
        0 0 0 0];

    sev3 = [0 0 0 0
        0 0 0 1
        1 1 1 0
        0 0 0 0
        0 0 0 0];

    sev4 = [0 0 0 1
        0 0 1 0
        1 1 0 0
        0 0 0 0
        0 0 0 0];

    sev5 = [0 0 0 0
        0 0 0 0
        1 1 0 0
        0 0 1 0
        0 0 0 1];

    sev6 = [0 0 0 1
        0 1 1 0
        1 0 0 0
        0 0 0 0
        0 0 0 0];

    sev7 = [0 0 0 0
        0 0 0 0
        1 0 0 0
        0 1 1 0
        0 0 0 1];

    sev8 = [0 0 0 0
        0 0 0 0
        1 0 0 0
        1 1 1 0
        0 0 0 1];

    sev9 = [0 0 0 1
        1 1 1 0
        1 0 0 0
        0 0 0 0
        0 0 0 0];

    sev10 = [0 0 1 1
        1 1 1 0
        1 0 0 0
        0 0 0 0
        0 0 0 0];

    sev11 = [0 0 0 0
        0 0 0 0
        1 0 0 0
        1 1 1 0
        0 0 1 1];

    % reflections of these

    sev12 = [0 0 0 0
        0 0 0 0
        0 1 1 1
        1 0 0 0
        0 0 0 0];

    sev13 = [0 0 0 0
        1 0 0 0
        1 1 1 1
        0 0 0 0
        0 0 0 0];

    sev14 = [1 0 0 0
        0 1 0 0
        0 0 1 1
        0 0 0 0
        0 0 0 0];

    sev15 = [0 0 0 0
        0 0 0 0
        0 0 1 1
        0 1 0 0
        1 0 0 0];

    sev16 = [1 0 0 0
        0 1 1 0
        0 0 0 1
        0 0 0 0
        0 0 0 0];

    sev17 = [0 0 0 0
        0 0 0 0
        0 0 0 1
        0 1 1 0
        1 0 0 0];

    sev18 = [0 0 0 0
        0 0 0 0
        0 0 0 1
        0 1 1 1
        1 0 0 0];

    sev19 = [1 0 0 0
        0 1 1 1
        0 0 0 1
        0 0 0 0
        0 0 0 0];

    sev20 = [1 1 0 0
        0 1 1 1
        0 0 0 1
        0 0 0 0
        0 0 0 0];

    sev21 = [0 0 0 0
        0 0 0 0
        0 0 0 1
        0 1 1 1
        1 1 0 0];

    % Filter for these patterns
    % UL
    UL_1 = bwhitmiss(UL_crop, sev1);
    UL_2 = bwhitmiss(UL_crop, sev2);
    UL_3 = bwhitmiss(UL_crop, sev3);
    UL_4 = bwhitmiss(UL_crop, sev4);
    UL_5 = bwhitmiss(UL_crop, sev5);
    UL_6 = bwhitmiss(UL_crop, sev6);
    UL_7 = bwhitmiss(UL_crop, sev7);
    UL_8 = bwhitmiss(UL_crop, sev8);
    UL_9 = bwhitmiss(UL_crop, sev9);
    UL_10 = bwhitmiss(UL_crop, sev10);
    UL_11 = bwhitmiss(UL_crop, sev11);
    UL_12 = bwhitmiss(UL_crop, sev12);
    UL_13 = bwhitmiss(UL_crop, sev13);
    UL_14 = bwhitmiss(UL_crop, sev14);
    UL_15 = bwhitmiss(UL_crop, sev15);
    UL_16 = bwhitmiss(UL_crop, sev16);
    UL_17 = bwhitmiss(UL_crop, sev17);
    UL_18 = bwhitmiss(UL_crop, sev18);
    UL_19 = bwhitmiss(UL_crop, sev19);
    UL_20 = bwhitmiss(UL_crop, sev20);
    UL_21 = bwhitmiss(UL_crop, sev21);

    UL_hm = UL_1 + UL_2 + UL_3 + UL_4 + UL_5 + UL_6 + UL_7 + UL_8 + UL_9 + UL_10 + UL_11 + UL_12 + UL_13 + UL_14 + UL_15 + UL_16 + UL_17 + UL_18 + UL_19 + UL_20 + UL_21;

    % LR
    LR_1 = bwhitmiss(LR_crop, sev1);
    LR_2 = bwhitmiss(LR_crop, sev2);
    LR_3 = bwhitmiss(LR_crop, sev3);
    LR_4 = bwhitmiss(LR_crop, sev4);
    LR_5 = bwhitmiss(LR_crop, sev5);
    LR_6 = bwhitmiss(LR_crop, sev6);
    LR_7 = bwhitmiss(LR_crop, sev7);
    LR_8 = bwhitmiss(LR_crop, sev8);
    LR_9 = bwhitmiss(LR_crop, sev9);
    LR_10 = bwhitmiss(LR_crop, sev10);
    LR_11 = bwhitmiss(LR_crop, sev11);
    LR_12 = bwhitmiss(LR_crop, sev12);
    LR_13 = bwhitmiss(LR_crop, sev13);
    LR_14 = bwhitmiss(LR_crop, sev14);
    LR_15 = bwhitmiss(LR_crop, sev15);
    LR_16 = bwhitmiss(LR_crop, sev16);
    LR_17 = bwhitmiss(LR_crop, sev17);
    LR_18 = bwhitmiss(LR_crop, sev18);
    LR_19 = bwhitmiss(LR_crop, sev19);
    LR_20 = bwhitmiss(LR_crop, sev20);
    LR_21 = bwhitmiss(LR_crop, sev21);

    LR_hm = LR_1 + LR_2 + LR_3 + LR_4 + LR_5 + LR_6 + LR_7 + LR_8 + LR_9 + LR_10 + LR_11 + LR_12 + LR_13 + LR_14 + LR_15 + LR_16 + LR_17 + LR_18 + LR_19 + LR_20 + LR_21;

    %%%%%% Now do this for filtering the horizontal crops - means we want to
    %%%%%% keep verticalish lines
    seh1 = rot90(sev1);
    seh2 = rot90(sev2);
    seh3 = rot90(sev3);
    seh4 = rot90(sev4);
    seh5 = rot90(sev5);
    seh6 = rot90(sev7);
    seh7 = rot90(sev8);
    seh8 = rot90(sev9);
    seh9 = rot90(sev9);
    seh10 = rot90(sev10);
    seh11 = rot90(sev11);
    seh12 = rot90(sev12);
    seh13 = rot90(sev13);
    seh14 = rot90(sev14);
    seh15 = rot90(sev15);
    seh16 = rot90(sev16);
    seh17 = rot90(sev17);
    seh18 = rot90(sev18);
    seh19 = rot90(sev19);
    seh20 = rot90(sev20);
    seh21 = rot90(sev21);

    % UR

    UR_1 = bwhitmiss(UR_crop, seh1);
    UR_2 = bwhitmiss(UR_crop, seh2);
    UR_3 = bwhitmiss(UR_crop, seh3);
    UR_4 = bwhitmiss(UR_crop, seh4);
    UR_5 = bwhitmiss(UR_crop, seh5);
    UR_6 = bwhitmiss(UR_crop, seh6);
    UR_7 = bwhitmiss(UR_crop, seh7);
    UR_8 = bwhitmiss(UR_crop, seh8);
    UR_9 = bwhitmiss(UR_crop, seh9);
    UR_10 = bwhitmiss(UR_crop, seh10);
    UR_11 = bwhitmiss(UR_crop, seh11);
    UR_12 = bwhitmiss(UR_crop, seh12);
    UR_13 = bwhitmiss(UR_crop, seh13);
    UR_14 = bwhitmiss(UR_crop, seh14);
    UR_15 = bwhitmiss(UR_crop, seh15);
    UR_16 = bwhitmiss(UR_crop, seh16);
    UR_17 = bwhitmiss(UR_crop, seh17);
    UR_18 = bwhitmiss(UR_crop, seh18);
    UR_19 = bwhitmiss(UR_crop, seh19);
    UR_20 = bwhitmiss(UR_crop, seh20);
    UR_21 = bwhitmiss(UR_crop, seh21);

    UR_hm = UR_1 + UR_2 + UR_3 + UR_4 + UR_5 + UR_6 + UR_7 + UR_8 + UR_9 + UR_10 + UR_11 + UR_12 + UR_13 + UR_14 + UR_15 + UR_16 + UR_17 + UR_18 + UR_19 + UR_20 + UR_21;

    % LL
    LL_1 = bwhitmiss(LL_crop, seh1);
    LL_2 = bwhitmiss(LL_crop, seh2);
    LL_3 = bwhitmiss(LL_crop, seh3);
    LL_4 = bwhitmiss(LL_crop, seh4);
    LL_5 = bwhitmiss(LL_crop, seh5);
    LL_6 = bwhitmiss(LL_crop, seh6);
    LL_7 = bwhitmiss(LL_crop, seh7);
    LL_8 = bwhitmiss(LL_crop, seh8);
    LL_9 = bwhitmiss(LL_crop, seh9);
    LL_10 = bwhitmiss(LL_crop, seh10);
    LL_11 = bwhitmiss(LL_crop, seh11);
    LL_12 = bwhitmiss(LL_crop, seh12);
    LL_13 = bwhitmiss(LL_crop, seh13);
    LL_14 = bwhitmiss(LL_crop, seh14);
    LL_15 = bwhitmiss(LL_crop, seh15);
    LL_16 = bwhitmiss(LL_crop, seh16);
    LL_17 = bwhitmiss(LL_crop, seh17);
    LL_18 = bwhitmiss(LL_crop, seh18);
    LL_19 = bwhitmiss(LL_crop, seh19);
    LL_20 = bwhitmiss(LL_crop, seh20);
    LL_21 = bwhitmiss(LL_crop, seh21);

    LL_hm = LL_1 + LL_2 + LL_3 + LL_4 + LL_5 + LL_6 + LL_7 + LL_8 + LL_9 + LL_10 + LL_11 + LL_12 + LL_13 + LL_14 + LL_15 + LL_16 + LL_17 + LL_18 + LL_19 + LL_20 + LL_21;

    %%%%%%%%% Need to stitch the cropped images to back in to re-make overall
    %%%%%%%%% images

    UL_crop2 = skel{1}( :,1:UL_edge);
    UR_crop2 = skel{2}(1:UR_edge, :);
    LR_crop2 = skel{4}( :,LR_edge:size(skel{4},2));
    LL_crop2 = skel{3}(LL_edge:size(skel{3},1), :);

    % Now stitch back together
    UL = [UL_crop2 UL_hm];
    UR = [UR_crop2
        UR_hm];
    LR = [LR_hm LR_crop2];
    LL = [LL_hm
        LL_crop2];

    % Make sure all are same dimensions
    UL = UL(1:(y_pixels/2),1:(x_pixels/2));
    UR = UR(1:(y_pixels/2),1:(x_pixels/2));
    LR = LR(1:(y_pixels/2),1:(x_pixels/2));
    LL = LL(1:(y_pixels/2),1:(x_pixels/2));
    % Stitch all together
    X_skel = [UL UR
        LL LR];

    X_skel = X_skel(1:(y_pixels), 1:(x_pixels));

    % Calculate final VD
    vein_density = vein{1} + vein{2} + vein{3} + vein{4};

    split = "split";


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% If the image is poorly stained with a low vein density (hence benefits
    %%% from the splitting) then it will have a high numel(n) (number of
    %%% endpoints), therefore this will be used as the cutoff to decide if an
    %%% image needs further processing.

    B = bwmorph(X_skel, 'branchpoints');
    E = bwmorph(X_skel, 'endpoints');
    B_loc = find(B); %numel(B_loc) == number of branchpoints
    [y,x] = find(E); %Now the numel(x) == number of endpoints




    if ((strcmp(split_mode, "auto")) && (numel(B_loc) < branch_threshold) && (numel(x) < endpoint_threshold)) || ((strcmp(split_mode, "auto")) && (sd > sd_threshold))
        split = "notsplit";
        trim = num2str(trim_factor);
        % Improve contrast
        X_con = adapthisteq(X_gray);
        % Blur
        bf = initial_blur_factor;
        X_blur = blur(X_con, bf);
        % Binarize
        X_bin = imbinarize(X_blur);
        % Invert
        X_inv = imcomplement(X_bin);
        % Skeletonize image
        X_skel = bwskel(X_inv);

        % Create skeleton mask
        B = bwmorph(X_skel, 'branchpoints');
        E = bwmorph(X_skel, 'endpoints');
        [y,x] = find(E); %Now the numel(x) == number of endpoints
        B_loc = find(B); %numel(B_loc) == number of branchpoints

        % Initial mask
        Dmask = false(size(X_skel));

        % Trim the mask
        for k = 1:numel(x)
            %BWDISTGEODESIC starts at a given point, and then calculates the
            %distance as you travel along the path.
            D = bwdistgeodesic(X_skel,x(k),y(k));
            %Start at a endpoint, start walking and find all pixels that are closer
            %than the nearest branchpoint.
            distanceToBranchPt = min(D(B_loc)); %D(B_loc) gives the lengths of each branchpoint
            if isempty(distanceToBranchPt) % if there are no branch points do nothing
            else
                if distanceToBranchPt > trim_factor %If loose branches are larger than our universal threshold we do nothing with them
                else
                    Dmask(D < distanceToBranchPt) =true; %Else we remove them
                end
            end
        end
        % Redefine our skeletonised image
        skelD = X_skel - Dmask;

        % Initial estimate of VD - will be final estimate if not trimmed based
        % on mode
        vd = (sum(skelD(:))*(pixel_length_um))/image_area_mm2;

        vein_density = vd;


        % Calculate the ratio of branchpoints:Initial_VD
        ratio2 = (numel(B_loc)/vd);

        % If this ratio is too high (>trim_threshold), then trimming will occur again,
        % based on the mode length of loose branch ends instead of the
        % "trim_factor". If we did mode from the begining then we would risk
        % removing data from those with low VD.

        if ratio2 > trim_threshold
            trim = "mode";
            % Here we calculate the length of each "loose branch end" i.e. each
            % line that potentially could be trimmed
            loose_size = {}; %Empty cell array for each "loose length"
            for k = 1:numel(x)
                D = bwdistgeodesic(X_skel,x(k),y(k));
                distanceToBranchPt = min(D(B_loc));
                loose_size = [loose_size, distanceToBranchPt]; %#ok<*AGROW> (surpresses warning)
            end
            % Remove infinate values, otherwise we cannot calculate mode
            loose_size( cellfun( @(C) any(isnan(C) | isinf(C)), loose_size ) ) = [];

            % Threshold is based on the mode length of each loose strand - based on
            % fact that the most common size will be the smaller threads because
            % these caputre noise that is ~similar all over each image
            T = (mode(cell2mat(loose_size))+20); % we add 20 pixels just to make sure all are captured

            % We then create a mask containing our skeletonised image; this will be
            % what is taken away from our final image after it is adjusted for
            % loose branches
            Dmask = false(size(X_skel));

            % Adjusting the mask by removing strands smaller than our threshold
            for k = 1:numel(x)
                D = bwdistgeodesic(X_skel,x(k),y(k));
                distanceToBranchPt = min(D(B_loc));
                if isempty(distanceToBranchPt)
                else
                    if distanceToBranchPt > T % T = threshold  % Do nothing
                    else
                        Dmask(D < distanceToBranchPt) =true;
                    end
                end
            end

            %Skeleton after trimming
            skelD = X_skel - Dmask;

            % Re-calculate vein_density
            vd = (sum(skelD(:))*(pixel_length_um))/image_area_mm2;
            vein_density = vd;

            Processed = strcat('Processed_', split, '_tf-', trim, '_bf-', num2str(round(bf)), '_', filename);
            X_skel = skelD;

            figure('Name', Processed)
            imshow(labeloverlay(X_gray,X_skel,'Transparency',0,'Colormap',[ ; 1 0 0]))

        else
            Processed = strcat('Processed_', split, '_tf-', trim, '_bf-', num2str(round(bf)), '_', filename);
            X_skel = skelD;

            figure('Name', Processed)
            imshow(labeloverlay(X_gray,X_skel,'Transparency',0,'Colormap',[ ; 1 0 0]))
        end

    else
        Processed = strcat('Processed_', split, '_tf-', trim, '_bf-', num2str(round(bf)), '_', filename);

        figure('Name', Processed)
        imshow(labeloverlay(X_gray,X_skel,'Transparency',0))



    end





else
    if (strcmp(split_mode, "whole"))
        trim = num2str(trim_factor);
        split = "notsplit";
        % blur based on sd
        % sd = std2(images{i});
        % Improve contrast
        X_con = adapthisteq(X_gray);
        % Blur
        bf = initial_blur_factor;
        X_blur = blur(X_con, bf);
        % Binarize
        X_bin = imbinarize(X_blur);
        % Invert
        X_inv = imcomplement(X_bin);
        % Skeletonize image
        X_skel = bwskel(X_inv);

        % Create skeleton mask
        B = bwmorph(X_skel, 'branchpoints');
        E = bwmorph(X_skel, 'endpoints');
        [y,x] = find(E); %Now the numel(x) == number of endpoints
        B_loc = find(B); %numel(B_loc) == number of branchpoints

        % Initial mask
        Dmask = false(size(X_skel));

        % Trim the mask
        for k = 1:numel(x)
            %BWDISTGEODESIC starts at a given point, and then calculates the
            %distance as you travel along the path.
            D = bwdistgeodesic(X_skel,x(k),y(k));
            %Start at a endpoint, start walking and find all pixels that are closer
            %than the nearest branchpoint.
            distanceToBranchPt = min(D(B_loc)); %D(B_loc) gives the lengths of each branchpoint
            if isempty(distanceToBranchPt) % if there are no branch points do nothing
            else
                if distanceToBranchPt > trim_factor %If loose branches are larger than our universal threshold we do nothing with them
                else
                    Dmask(D < distanceToBranchPt) =true; %Else we remove them
                end
            end
        end
        % Redefine our skeletonised image
        skelD = X_skel - Dmask;

        % Initial estimate of VD - will be final estimate if not trimmed based
        % on mode
        vd = (sum(skelD(:))*(pixel_length_um))/image_area_mm2;
        % Calculate the ratio of branchpoints:Initial_VD
        ratio = (numel(B_loc)/vd);


        % If this ratio is too high (>trim_threshold), then trimming will occur again,
        % based on the mode length of loose branch ends instead of the
        % "trim_factor". If we did mode from the begining then we would risk
        % removing data from those with low VD.

        if ratio > trim_threshold
            trim = "mode";
            % Here we calculate the length of each "loose branch end" i.e. each
            % line that potentially could be trimmed
            loose_size = {}; %Empty cell array for each "loose length"
            for k = 1:numel(x)
                D = bwdistgeodesic(X_skel,x(k),y(k));
                distanceToBranchPt = min(D(B_loc));
                loose_size = [loose_size, distanceToBranchPt]; %#ok<*AGROW> (surpresses warning)
            end
            % Remove infinate values, otherwise we cannot calculate mode
            loose_size( cellfun( @(C) any(isnan(C) | isinf(C)), loose_size ) ) = [];

            % Threshold is based on the mode length of each loose strand - based on
            % fact that the most common size will be the smaller threads because
            % these caputre noise that is ~similar all over each image
            T = (mode(cell2mat(loose_size))+20); % we add 20 pixels just to make sure all are captured

            % We then create a mask containing our skeletonised image; this will be
            % what is taken away from our final image after it is adjusted for
            % loose branches
            Dmask = false(size(X_skel));

            % Adjusting the mask by removing strands smaller than our threshold
            for k = 1:numel(x)
                D = bwdistgeodesic(X_skel,x(k),y(k));
                distanceToBranchPt = min(D(B_loc));
                if isempty(distanceToBranchPt)
                else
                    if distanceToBranchPt > T % T = threshold  % Do nothing
                    else
                        Dmask(D < distanceToBranchPt) =true;
                    end
                end
            end

            %Skeleton after trimming
            skelD = X_skel - Dmask;

            % Re-calculate vein_density
            vd = (sum(skelD(:))*(pixel_length_um))/image_area_mm2;
        end

        vein_density = vd;
        X_skel = skelD;

        Processed = strcat('Processed_', split, '_tf-', trim, '_bf-', num2str(round(bf)), '_', filename);

        figure('Name', Processed)
        imshow(labeloverlay(X_gray,X_skel,'Transparency',0, 'Colormap',[ ; 1 0 0]))


    else
        if (strcmp(split_mode, "monocot"))
            trim = num2str(trim_factor);
            split = "notsplit";
            % blur based on sd
            % sd = std2(images{i});
            bf = initial_blur_factor;
            % Improve contrast
            X_con = adapthisteq(X_gray);
            % Blur
            X_blur = blur(X_con, bf);
            % Binarize
            X_bin = imbinarize(X_blur);
            % Invert
            X_inv = imcomplement(X_bin);
            % Skeletonize image
            X_skel = bwskel(X_inv);

            % Create skeleton mask
            B = bwmorph(X_skel, 'branchpoints');
            E = bwmorph(X_skel, 'endpoints');
            [y,x] = find(E); %Now the numel(x) == number of endpoints
            B_loc = find(B); %numel(B_loc) == number of branchpoints

            % Initial mask
            Dmask = false(size(X_skel));

            % Trim the mask
            for k = 1:numel(x)
                %BWDISTGEODESIC starts at a given point, and then calculates the
                %distance as you travel along the path.
                D = bwdistgeodesic(X_skel,x(k),y(k));
                %Start at a endpoint, start walking and find all pixels that are closer
                %than the nearest branchpoint.
                distanceToBranchPt = min(D(B_loc)); %D(B_loc) gives the lengths of each branchpoint
                if isempty(distanceToBranchPt) % if there are no branch points do nothing
                else
                    if distanceToBranchPt > trim_factor %If loose branches are larger than our universal threshold we do nothing with them
                    else
                        Dmask(D < distanceToBranchPt) =true; %Else we remove them
                    end
                end
            end
            % Redefine our skeletonised image
            skelD = X_skel - Dmask;
            skelD = logical(skelD);

            % Because commissural veins complicate analysis, we want to
            % remove them
            B = bwmorph(skelD, 'branchpoints');
            [y,x] = find(B); % Locations of branchpoints
            B_loc = find(B); %numel(B_loc) == number of branchpoints

            % Initial mask
            Dmask = false(size(skelD));


            m1 = [0 1 0
                1 1 1
                0 0 0];

            m2 = [0 1 0
                1 1 0
                0 1 0];

            m3 = [0 1 0
                0 1 1
                0 1 0];

            m4 = [0 0 0
                1 1 1
                0 1 0];

            m5 = [1 0 1
                0 1 0
                0 1 0];

            m6 = [0 0 1
                1 1 0
                0 0 1];

            m7 = [0 1 0
                0 1 0
                1 0 1];

            m8 = [1 0 0
                0 1 1
                1 0 0];

            m9 = [1 0 1
                0 1 0
                1 0 0];

            m10 = [1 0 1
                0 1 0
                0 0 1];

            m11 = [0 0 1
                0 1 0
                1 0 1];

            m12 = [1 0 1
                0 1 0
                0 0 1];

            m13 = [1 0 0
                0 1 0
                1 0 1];

            m14 = [0 0 1
                0 1 0
                1 0 1];

            m15 = [0 1 0
                1 1 0
                0 0 1];

            m16 = [0 1 0
                0 1 1
                1 0 0];

            m17 = [0 0 1
                1 1 0
                0 1 0];

            m18 = [1 0 0
                0 1 1
                0 1 0];

            m = cell(18, 1);
            m{1} = m1;
            m{2} = m2;
            m{3} = m3;
            m{4} = m4;
            m{5} = m5;
            m{6} = m6;
            m{7} = m7;
            m{8} = m8;
            m{9} = m9;
            m{10} = m10;
            m{11} = m11;
            m{12} = m12;
            m{13} = m13;
            m{14} = m14;
            m{15} = m15;
            m{16} = m16;
            m{17} = m17;
            m{18} = m18;

            cor = [];


            for k = 1:numel(x)
                %BWDISTGEODESIC starts at a given point, and then calculates the
                %distance as you travel along the path.
                D = bwdistgeodesic(skelD,x(k),y(k)); % starting at a given branch point
                lengthOfBranchPt = (D(B_loc));
                if isempty(lengthOfBranchPt) % If no branchpoints do nothing
                else
                    lengthOfBranchPt(lengthOfBranchPt==0) = [];
                    lengthOfBranchPt = min(lengthOfBranchPt);
                    if lengthOfBranchPt > monocot_branch_length_to_keep %If branches are large do nothing
                    else
                        if lengthOfBranchPt < mono_min_param % If branches are too small means it is linking two branches so get rid
                        else
                            Dmask(D < lengthOfBranchPt) = true; %else we remove them

                            temp = (D < lengthOfBranchPt);
                            E2 = bwmorph(temp, 'endpoints');
                            [ye,xe] = find(E2); % co-ordinates of endpoints

                            B2 = bwmorph(temp, 'branchpoints');
                            B_loc2 = find(B2);

                            for j = 1:numel(xe)
                                if (ye(j) == 1 || ye(j) == 2 || ye(j) == 3 || xe(j) == 1 || xe(j) == 2 || xe(j) == 3 || ye(j) == (y_pixels) || ye(j) == (y_pixels-1) || ye(j) == (y_pixels - 2) || xe(j) == (x_pixels) || xe(j) == (x_pixels-1) || xe(j) == (x_pixels - 2)) % Do nothing as this will make it out of bounds
                                else
                                    bm = skelD(ye(j)-3:(ye(j)+3) ,xe(j)-3: (xe(j)+3)); % Extract matrix around the point of interest
                                    for i = 1:length(m)
                                        comp = normxcorr2(m{i},bm); % Compare each of your defined matricies with the overall one
                                        cor(i) = max(max(comp));
                                    end
                                    if max(cor) < 0.999999 % i.e. if there is not a branchpoint in the original skeleton that corresponds to an end point in the removed segment

                                        D2 = bwdistgeodesic(temp,xe(j),ye(j)); % starting at a given end point in the removed bit
                                        lengthOfBranchPt2 = (D2(B_loc2));
                                        if isempty(lengthOfBranchPt2) % makes the code work; do nothing
                                        else
                                            Dmask(D2 < lengthOfBranchPt2) = false; % Remove these pixels from final mask
                                        end
                                    else % Else do nothing
                                    end
                                end
                            end
                        end
                    end
                end
            end

            skelD = skelD - Dmask;
            skelD = logical(skelD);



            % Initial estimate of VD - will be final estimate if not trimmed based
            % on mode
            vd = (sum(skelD(:))*(pixel_length_um))/image_area_mm2;
            % Calculate the ratio of branchpoints:Initial_VD
            ratio = (numel(B_loc)/vd);


            % If this ratio is too high (>trim_threshold), then trimming will occur again,
            % based on the mode length of loose branch ends instead of the
            % "trim_factor". If we did mode from the begining then we would risk
            % removing data from those with low VD.

            if ratio > trim_threshold
                trim = "mode";
                % Here we calculate the length of each "loose branch end" i.e. each
                % line that potentially could be trimmed
                loose_size = {}; %Empty cell array for each "loose length"
                for k = 1:numel(x)
                    D = bwdistgeodesic(X_skel,x(k),y(k));
                    distanceToBranchPt = min(D(B_loc));
                    loose_size = [loose_size, distanceToBranchPt]; %#ok<*AGROW> (surpresses warning)
                end
                % Remove infinate values, otherwise we cannot calculate mode
                loose_size( cellfun( @(C) any(isnan(C) | isinf(C)), loose_size ) ) = [];

                % Threshold is based on the mode length of each loose strand - based on
                % fact that the most common size will be the smaller threads because
                % these caputre noise that is ~similar all over each image
                T = (mode(cell2mat(loose_size))+20); % we add 20 pixels just to make sure all are captured

                % We then create a mask containing our skeletonised image; this will be
                % what is taken away from our final image after it is adjusted for
                % loose branches
                Dmask = false(size(X_skel));

                % Adjusting the mask by removing strands smaller than our threshold
                for k = 1:numel(x)
                    D = bwdistgeodesic(X_skel,x(k),y(k));
                    distanceToBranchPt = min(D(B_loc));
                    if isempty(distanceToBranchPt)
                    else
                        if distanceToBranchPt > T % T = threshold  % Do nothing
                        else
                            Dmask(D < distanceToBranchPt) =true;
                        end
                    end
                end

                % Redefine our skeletonised image
                skelD = X_skel - Dmask;
                skelD = logical(skelD);

                % Because commissural veins complicate analysis, we want to
                % remove them
                B = bwmorph(skelD, 'branchpoints');
                [y,x] = find(B); % Locations of branchpoints
                B_loc = find(B); %numel(B_loc) == number of branchpoints

                % Initial mask
                Dmask = false(size(skelD));


                m1 = [0 1 0
                    1 1 1
                    0 0 0];

                m2 = [0 1 0
                    1 1 0
                    0 1 0];

                m3 = [0 1 0
                    0 1 1
                    0 1 0];

                m4 = [0 0 0
                    1 1 1
                    0 1 0];

                m5 = [1 0 1
                    0 1 0
                    0 1 0];

                m6 = [0 0 1
                    1 1 0
                    0 0 1];

                m7 = [0 1 0
                    0 1 0
                    1 0 1];

                m8 = [1 0 0
                    0 1 1
                    1 0 0];

                m9 = [1 0 1
                    0 1 0
                    1 0 0];

                m10 = [1 0 1
                    0 1 0
                    0 0 1];

                m11 = [0 0 1
                    0 1 0
                    1 0 1];

                m12 = [1 0 1
                    0 1 0
                    0 0 1];

                m13 = [1 0 0
                    0 1 0
                    1 0 1];

                m14 = [0 0 1
                    0 1 0
                    1 0 1];

                m15 = [0 1 0
                    1 1 0
                    0 0 1];

                m16 = [0 1 0
                    0 1 1
                    1 0 0];

                m17 = [0 0 1
                    1 1 0
                    0 1 0];

                m18 = [1 0 0
                    0 1 1
                    0 1 0];

                m = cell(18, 1);
                m{1} = m1;
                m{2} = m2;
                m{3} = m3;
                m{4} = m4;
                m{5} = m5;
                m{6} = m6;
                m{7} = m7;
                m{8} = m8;
                m{9} = m9;
                m{10} = m10;
                m{11} = m11;
                m{12} = m12;
                m{13} = m13;
                m{14} = m14;
                m{15} = m15;
                m{16} = m16;
                m{17} = m17;
                m{18} = m18;

                cor = [];


                for k = 1:numel(x)
                    %BWDISTGEODESIC starts at a given point, and then calculates the
                    %distance as you travel along the path.
                    D = bwdistgeodesic(skelD,x(k),y(k)); % starting at a given branch point
                    lengthOfBranchPt = (D(B_loc));
                    if isempty(lengthOfBranchPt) % If no branchpoints do nothing
                    else
                        lengthOfBranchPt(lengthOfBranchPt==0) = [];
                        lengthOfBranchPt = min(lengthOfBranchPt);
                        if lengthOfBranchPt > monocot_branch_length_to_keep %If branches are large do nothing
                        else
                            if lengthOfBranchPt < mono_min_param % If bramches are too small means it is linking two branches so get rid
                            else
                                Dmask(D < lengthOfBranchPt) = true; %else we remove them

                                temp = (D < lengthOfBranchPt);
                                E2 = bwmorph(temp, 'endpoints');
                                [ye,xe] = find(E2); % co-ordinates of endpoints

                                B2 = bwmorph(temp, 'branchpoints');
                                B_loc2 = find(B2);

                                for j = 1:numel(xe)
                                    if (ye(j) == 1 || ye(j) == 2 || ye(j) == 3 || xe(j) == 1 || xe(j) == 2 || xe(j) == 3 || ye(j) == (y_pixels) || ye(j) == (y_pixels-1) || ye(j) == (y_pixels - 2) || xe(j) == (x_pixels) || xe(j) == (x_pixels-1) || xe(j) == (x_pixels - 2)) % Do nothing as this will make it out of bounds
                                    else
                                        bm = skelD(ye(j)-3:(ye(j)+3) ,xe(j)-3: (xe(j)+3)); % Extract matrix around the point of interest
                                        for i = 1:length(m)
                                            comp = normxcorr2(m{i},bm); % Compare each of your defined matricies with the overall one
                                            cor(i) = max(max(comp));
                                        end
                                        if max(cor) < 0.999999 % i.e. if there is not a branchpoint in the original skeleton that corresponds to an end point in the removed segment

                                            D2 = bwdistgeodesic(temp,xe(j),ye(j)); % starting at a given end point in the removed bit
                                            lengthOfBranchPt2 = (D2(B_loc2));
                                            if isempty(lengthOfBranchPt2) % makes the code work; do nothing
                                            else
                                                Dmask(D2 < lengthOfBranchPt2) = false; % Remove these pixels from final mask
                                            end
                                        else % Else do nothing
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

                % Re-calculate vein_density
                vd = (sum(skelD(:))*(pixel_length_um))/image_area_mm2;
            end

            vein_density = vd;
            X_skel = skelD;

            Processed = strcat('Processed_', split, '_tf-', trim, '_bf-', num2str(round(bf)), '_', filename);

            figure('Name', Processed)
            imshow(labeloverlay(X_gray,X_skel,'Transparency',0, 'Colormap',[ ; 1 0 0]))

        end


    end
end
