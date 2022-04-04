# autoVD
MATLAB script for extracting vein density measurements from starch-stained leaf samples


## autoVD.m
This script contains functions to grayscale (if necessary) image, blur, binarise, and skeletonise the image. In addition it will calculate the number of Loose Branch Ends (LBEs) and calculate the ratio of this to estimated vein density. Based on this value trimming of LBEs will be done according to a Universal Trim Factor (UTF; defined by user - default of 90) or the mode length of LBEs.

## running_autoVD.m
Contains all that is needed to run autoVD.m script. User must define UTF (default 90), the ratio (default 0.01) and blur factor (default 15).

## blur.m
Blurs grayscale image - this code is copied from here:
https://uk.mathworks.com/matlabcentral/answers/472573-write-a-function-called-blur-that-blurs-the-input-image?s_tid=mwa_osa_a
I have asked how to reference this author but yet have recieved no reply.
