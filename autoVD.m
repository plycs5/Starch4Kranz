function [vein_density, ratio_Unitrim, grayscale_mean, grayscale_SD] = autoVD(filename, universal_trim_factor, blur_factor, loose_VD_R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read in image and convert to grayscale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = imread(filename);
X_gray = rgb2gray(X);

m = mean2(X_gray);

grayscale_mean = m;
grayscale_SD = std2(X_gray);

m = mean2(X_gray);
X_blur = X_gray;


if m <110
     X_blur = blur(X_gray, blur_factor);
 else 
     if m >=110
         X_blur = blur(X_gray, (2*blur_factor));
     end
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Binarise and skelatonise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Binarize
X_blur_bin = imbinarize(X_blur);

% %Displaying binarized image (commented out)
% subplot 221
% imagesc(X_gray)
% subplot 222
% imagesc(X_gray_bin)
% subplot 223
% imagesc(X_blur)
% subplot 224
% imagesc(X_blur_bin)
% colormap 'gray'

% Invert binarized image
X_blur_bin_inv = imcomplement(X_blur_bin);

% Skeletonize image
X_blur_bin_inv_skel = bwskel(X_blur_bin_inv);

X_blur_bin_inv_skel_unprocessed = X_blur_bin_inv_skel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Perform universal trim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate parameters of branchpoints
B = bwmorph(X_blur_bin_inv_skel, 'branchpoints');
E = bwmorph(X_blur_bin_inv_skel, 'endpoints');
[y,x] = find(E); %Now the numel(x) == number of endpoints 
B_loc = find(B); %numel(B_loc) == number of branchpoints


%We then create a mask containing our skelotised image, this will be
%what is taken away from our final image after it is adjusted for
%smaller strands
Dmask = false(size(X_blur_bin_inv_skel));

for k = 1:numel(x)
    %BWDISTGEODESIC starts at a given point, and then calculates the
    %distance as you travel along the path.
    D = bwdistgeodesic(X_blur_bin_inv_skel,x(k),y(k));
    %Start at a endpoint, start walking and find all pixels that are closer
    %than the nearest branchpoint.
    distanceToBranchPt = min(D(B_loc)); %D(B_loc) gives the lengths of each branchpoint
    if distanceToBranchPt > universal_trim_factor %If loose branches are larger than our universal threshold we do nothing with them
    else
        Dmask(D < distanceToBranchPt) =true; %Else we remove them
    end
end

% Redefine our skeletonised image
skelD = X_blur_bin_inv_skel - Dmask; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate estimate of VD and ratio against (number of loose
%branches)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%
%This will be our VD if ratio < loose_VD_R; i.e. with universal trim
%applied only
vein_density = (sum(skelD(:))*(14/20))/1.541; %The transformation here is specific to my image dimensions and gives the output in um/mm^2
%%%%%%%%%%%%

branchpoint_to_VD_ratio_uniTrim = (numel(B_loc)/vein_density);
ratio_Unitrim = branchpoint_to_VD_ratio_uniTrim;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-trim based on loose_branch_length mode if loose_VD_R1 is not satisfied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ratio_Unitrim > loose_VD_R
    %Here we calculate the length of each "loose branch point" i.e. each
    %line that potentially could be trimmed
    loose_size = {}; %Empty cell array for each "loose length"
    for k = 1:numel(x)
        %BWDISTGEODESIC starts at a given point, and then calculates the
        %distance as you travel along the path.
        D = bwdistgeodesic(X_blur_bin_inv_skel,x(k),y(k));
        %Start at a endpoint, start walking and find all pixels that are closer
        %than the nearest branchpoint.
        distanceToBranchPt = min(D(B_loc)); %D(B_loc) gives the lengths of each branchpoint
        loose_size = [loose_size, distanceToBranchPt];
    end
    %Remove infinate values - do not let us calculate mode
    loose_size( cellfun( @(C) any(isnan(C) | isinf(C)), loose_size ) ) = [];
    
    %Threshold is based on the mode length of each loose strand - based on
    %fact that the most common size will be the smaller threads because
    %these caputre noise that is ~similar all over each image
    T = (mode(cell2mat(loose_size))+20);
    
    
    %We then create a mask containing our skelotised image, this will be
    %what is taken away from our final image after it is adjusted for
    %smaller strands
    Dmask = false(size(X_blur_bin_inv_skel));
    
    %Adjusting the mask by removing strands smaller than our threshold
    for k = 1:numel(x)
        D = bwdistgeodesic(X_blur_bin_inv_skel,x(k),y(k));
        distanceToBranchPt = min(D(B_loc));
        if distanceToBranchPt > T %T = threshold size
        else
            Dmask(D < distanceToBranchPt) =true;
        end
    end
    
    %Skeleton after trimming
    skelD = X_blur_bin_inv_skel - Dmask;
    
    
    %Calculate estimate of VD
    
    %%%%%%%%%%%%
    %This will be the VD if, after Re-Trimming (based on mode LBE length), the ratio < 0.01 i.e. a trim was required but no additional blur
    vein_density = (sum(skelD(:))*(14/20))/1.541;
    %%%%%%%%%%%%
    
     
end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display output - for now this will display both the trimmed and non
%trimmed image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Processed = strcat('Processed_', filename);

Unprocessed = strcat('Unprocessed_', filename);

figure('Name', Processed)
imshow(labeloverlay(X_gray,skelD,'Transparency',0))

%figure('Name', Unprocessed)
%imshow(labeloverlay(X_gray,X_blur_bin_inv_skel_unprocessed,'Transparency',0))



