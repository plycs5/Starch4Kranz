%Running autoVD function
%Run this in the same directory as images
close all
clear all


images = dir('*tif');

%User defined inputs

UTF = 90; %Universal trim Factor
BTF = 15; %Initial blur factor
LVR = 0.01; %Ratio of loose branches to VD used to accept image or see if further trimming is required



vd = {}; % pre-assign cell array for outputsratio_untrim
R_Unitrim = {};
R_Trim = {};
Gmean = {};
GSD = {};


image_name = {}; % Creates an empty cell array

for file = images'
    
    [vein_density, ratio_Unitrim, grayscale_mean, grayscale_SD] = autoVD(file.name, UTF, BTF, LVR);
    
    vd = [vd, vein_density];
    R_Unitrim = [R_Unitrim, ratio_Unitrim];
    Gmean = [Gmean, grayscale_mean];
    GSD = [GSD, grayscale_SD];
    
    name = file.name;%gives the names as a char
    name = {file.name};%converts the names to a cell
    image_name = [image_name, name]; % adds the output to the empty cell array
    
end


results = [image_name; vd; R_Unitrim; Gmean; GSD];
results = transpose(results);

header = {'Image_Name', 'VD', 'R_Unitrim', 'R_Trim', 'R_Reblur_Unitrim', 'R_Reblur_Trim', "Mean", "SD"};

results = [header; results];
