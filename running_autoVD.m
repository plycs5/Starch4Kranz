% Running vein_density250822

% Run this in the same directory as images
%close all
clear all


images = dir('*.tif');

% pre-assign cell array for image nemes and results
image_name = {}; 

vd = {}; 

for file = images'
    
    [vein_density] = vein_density161122(file.name, "monocot", 90, 0.7, 1536, 2048);
    
    name = file.name; % gives the names as a char
    name = {file.name}; % converts the names to a cell
    image_name = [image_name, name]; % adds the output to the empty cell array
    
    vd = [vd, vein_density];

end


results = [image_name; vd];
results = transpose(results);

header = {'Image_Name', 'VD_um_mm2'};

results = [header; results];
