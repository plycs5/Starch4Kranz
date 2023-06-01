# -*- coding: utf-8 -*-
"""
For analysing starch-stained images automically
"""

# Import libraries
import cv2
import numpy as np
import matplotlib.pyplot as plt
import skimage
import scipy.ndimage as ndi
from scipy.spatial import distance
from plantcv import plantcv as pcv
import statistics

# Deine additional functions
###### Round numbers
def round_to_odd(number):
    rounded_number = round(number)
    if rounded_number % 2 == 0:
        return rounded_number + 1
    else:
        return rounded_number
    
######
###### Finding connected pixels
def find_connected_pixels(skeleton_img, start_coord, end_coord, exp_elucidian_dist):
    height, width = skeleton_img.shape
    visited = np.zeros((height, width), dtype=bool)
    connected_pixels = []  # Initialize connected_pixels outside dfs

    def swapPositions(lst, pos1, pos2):
        lst[pos1], lst[pos2] = lst[pos2], lst[pos1]
        return lst

    def dfs(coord):
        visited[coord] = True

        if coord == end_coord:
            return True

        neighbors = [
            (coord[0] - 1, coord[1]),  # above
            (coord[0] + 1, coord[1]),  # below
            (coord[0], coord[1] - 1),  # left
            (coord[0], coord[1] + 1),  # right
            (coord[0] - 1, coord[1] - 1),  # Upper left diagonal
            (coord[0] - 1, coord[1] + 1),  # Upper Right diagonal
            (coord[0] + 1, coord[1] - 1),  # Lower left diagonal
            (coord[0] + 1, coord[1] + 1),  # Lower right diagonal
        ]

        for swap in range(8):
            neighbors = swapPositions(neighbors, swap, 0)
            for neighbor_coord in neighbors:
                if (
                    0 <= neighbor_coord[0] < height
                    and 0 <= neighbor_coord[1] < width
                    and skeleton_img[neighbor_coord] == 255
                    and not visited[neighbor_coord]
                ):
                    connected_pixels.append(neighbor_coord)  # Append to connected_pixels
                    if len(connected_pixels) < (exp_elucidian_dist + 75):
                        if dfs(neighbor_coord):  # Recursively call dfs
                            return True
                    connected_pixels.pop()  # Remove the last appended coordinate if it doesn't satisfy the condition

        return False

    dfs(start_coord)

    return connected_pixels

###################
### Start of function

def Starch4Kranz(filename, trim_factor, pixel_length_um, x_pixels, y_pixels, split_mode = "auto", shave = 40, initial_blur_factor=20, branch_threshold=125, endpoint_threshold=175, sd_threshold=75, trim_threshold=0.01, mono_min_param=30, monocot_branch_length_to_keep=200):

    
    # Define image dimensions
    pixel_area_um2 = pixel_length_um*pixel_length_um;
    pixel_area_mm2 = pixel_area_um2/1000000;
    
    image_area_mm2 = (y_pixels*x_pixels)*(pixel_area_mm2);
    
    # Read in image
    X = cv2.imread(filename)
    
    # Convert to grayscale
    X_gray = cv2.cvtColor(X, cv2.COLOR_BGR2GRAY)
    
    # Compute standard deviation of the grayscale image
    sd = X_gray.std()
    
    # Start script for when "auto" and "split" are selected:
    
    if (split_mode == "auto") or (split_mode == "split"):
        
        # Split image
        upperLeftX = X_gray[0:(y_pixels//2), 0:(x_pixels//2)]
        upperRightX = X_gray[0:(y_pixels//2), (x_pixels//2):x_pixels]
        lowerLeftX = X_gray[(y_pixels//2):y_pixels, 0:(x_pixels//2)]
        lowerRightX = X_gray[(y_pixels//2):y_pixels, (x_pixels//2):x_pixels]
    
        images = [upperLeftX, upperRightX, lowerLeftX, lowerRightX]
        # Make empty arrays of 4 to fill the final results of each segment
        bf_split = [0,0,0,0]
        skel = [0,0,0,0]
        vein = [0,0,0,0]
        
        #Define blur factor based on image sd - well stained images will have larger variance hence larger sd
        
        for i in range(4):
            # blur based on sd
            sd = np.std(images[i])
            m = np.mean(images[i])
            bf_scale1 = (40 * (sd / 100))
            bf_scale2 = (100 / (m * 2))
            bf = round_to_odd((bf_scale1 * bf_scale2) + 5)
    
            bf_split[i] = bf 
            
            # Improve contrast
            X_con = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8)).apply(images[i])
            # Blur
            X_blur = cv2.blur(X_con, ksize=(bf, bf))
            
            # Binarize
            _, X_bin = cv2.threshold(X_blur, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
            # Invert
            X_inv = cv2.bitwise_not(X_bin)
            # Make so binary true and false not binary 255 and 0
            X_inv = X_inv.astype(dtype=bool)
            # Skeletonize image
            X_skel = skimage.morphology.skeletonize(X_inv)
            # Convert back to unit8
            X_skel = X_skel.astype(np.uint8)
            X_skel*=255
                    
            # Search branchpoints
            B = pcv.morphology.find_branch_pts(X_skel)
            B = B.astype(bool)
            
            # Search endpoints
            
            E = pcv.morphology.find_tips(X_skel)
            E = E.astype(bool)            
            
            # Initial skeleton trim
            T = trim_factor
            skelD, segmented_img, segment_objects = pcv.morphology.prune(skel_img=X_skel, size=T)
           
            # Calculate initial VD
            skel_tot = skelD > 20 # Returns bool
            vd = (skel_tot.sum()*(pixel_length_um))/image_area_mm2 
            
            # Calculate the ratio of branchpoints:Initial_VD
            #Based on branchpoints before trimming as this gives a better view of how many superfluous ends there are
            ratio = (B.sum()/vd);
            
                
            # If this ratio is too high (>trim_threshold), then trimming will occur again,
            # based on the mode length of loose branch ends instead of the
            # "trim_factor". If we did mode from the begining then we would risk
            # removing data from those with low VD.
            if (ratio > trim_threshold):
                       
                skelD2, seg_img2, seg_obj2 = pcv.morphology.prune(skel_img=X_skel, size = 0)
                branch_ends, core_trace = pcv.morphology.segment_sort(skel_img=skelD2, objects=seg_obj2)
                # Get lengths of branch ends
                labeled_img  = pcv.morphology.segment_path_length(segmented_img=seg_img2, objects=branch_ends, label="default")
                branch_lengths = pcv.outputs.observations['default']['segment_path_length']['value']
            
                # Set trimming threshold - we add a number of pixels proportionate to the image length so we make sure to capture all of them, i.e. we take the mode length + a bit extra as the trimming factot
                T = statistics.mode(branch_lengths) + (x_pixels/100)
                # Trim
                skelD, segmented_img, segment_objects = pcv.morphology.prune(skel_img=X_skel, size=T)
                # Re-calculate VD
                skel_tot = skelD > 20 # Returns bool
                vd = (skel_tot.sum()*(pixel_length_um))/image_area_mm2 
    
            # Return the values for the ith segment
            skel[i] = skelD
            vein[i] = vd
        
        # Reccord the average bf for all images
        bf = statistics.mean(bf_split)
            
        #### Now we need to "stitch" the skeletons together
    
        # Before doing so we need to convert, from white to black, the pixels at
        # the edge of some of the skeletons so that we don't capture the same bit
        # twice
    
        # The edges to convert - lets go clockwise
        # upperLeft (skel[0]) - remove pixels along right edge that are more vertical
        # upperRight (skel1]) - remove pixels along bottom edge that are more horizontal
        # lowerLeft (skel[2]) - remove pixels along upper edge that are more horizontal
        # lowerRight (skel[3]) - remove pixels along  left edge that are more vertical
        
        # Crop out the edge we wish to alter
        UL_crop = skel[0][:,-shave:]
        UR_crop = skel[1][-shave:,:]
        LL_crop = skel[2][0:shave,:]
        LR_crop = skel[3][:,0:shave]
        
        # Define matrix patterns that we are interested in not removing from the
        # final edge
    
        # Firstly the masks to use for filtering the vertical crops ones (can be used for UL and LR) - these are masks that
        # will be kept, hence when we filter the crops we want to keep lines that
        # are horizontal
    
        # Define matrix patterns that we are interested in not removing from the
        # final edge
    
        # Firstly the masks to use for filtering the vertical crops ones (can be used for UL and LR) - these are masks that
        # will be kept, hence when we filter the crops we want to keep lines that
        # are horizontal
    
        sev1 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [1, 1, 1, 1],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev2 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [1, 1, 1, 0],
            [0, 0, 0, 1],
            [0, 0, 0, 0]])
    
        sev3 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 1],
            [1, 1, 1, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev4 = np.array([[0, 0, 0, 1],
            [0, 0, 1, 0],
            [1, 1, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev5 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]])
    
        sev6 = np.array([[0, 0, 0, 1],
            [0, 1, 1, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev7 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [1, 0, 0, 0],
            [0, 1, 1, 0],
            [0, 0, 0, 1]])
    
        sev8 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 1, 1, 0],
            [0, 0, 0, 1]])
    
        sev9 = np.array([[0, 0, 0, 1],
            [1, 1, 1, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev10 = np.array([[0, 0, 1, 1],
            [1, 1, 1, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev11 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 1, 1, 0],
            [0, 0, 1, 1]])
    
        # reflections of these
    
        sev12 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 1, 1, 1],
            [1, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev13 = np.array([[0, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 1, 1, 1],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev14 = np.array([[1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev15 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 1, 1],
            [0, 1, 0, 0],
            [1, 0, 0, 0]])
    
        sev16 = np.array([[1, 0, 0, 0],
            [0, 1, 1, 0],
            [0, 0, 0, 1],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev17 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 1, 1, 0],
            [1, 0, 0, 0]])
    
        sev18 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 1, 1, 1],
            [1, 0, 0, 0]])
    
        sev19 = np.array([[1, 0, 0, 0],
            [0, 1, 1, 1],
            [0, 0, 0, 1],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev20 = np.array([[1, 1, 0, 0],
            [0, 1, 1, 1],
            [0, 0, 0, 1],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])
    
        sev21 = np.array([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 1, 1, 1],
            [1, 1, 0, 0]])
        
        sev_list = [sev1, sev2, sev3, sev4, sev5, sev6, sev7, sev8, sev9, sev10, sev11, sev12, sev13, sev14, sev15, sev16, sev17, sev18, sev19, sev20, sev21]
        
        # If our crop images matches the sevs (UL and LR) or the sehs (defined in script) then remove matching pixels
        
        for i,sevs in enumerate(sev_list):
            # i is the index
            # sev is sev_list[i]
            seh = np.rot90(sevs) # UR, and LL are based on the seh (rotated 90 of the sevs)
            UL = ndi.binary_hit_or_miss(UL_crop,sevs)
            LR = ndi.binary_hit_or_miss(LR_crop,sevs)
            UR = ndi.binary_hit_or_miss(UR_crop, seh)
            LL = ndi.binary_hit_or_miss(LL_crop, seh)
            
            if i == 0:
                UL_hm, LR_hm, UR_hm, LL_hm = UL, LR, UR, LL
            else:
                UL_hm += UL
                LR_hm += LR
                UR_hm += UR
                LL_hm += LL
        
       
        ######## Need to stitch the cropped images to back in to re-make overall images
        
        ## Extract other part of image (side not cropped)
        
        UL_crop2 = skel[0][:, 0:-shave]
        UR_crop2 = skel[1][0:-shave, :]
        LL_crop2 = skel[2][shave:, :]
        LR_crop2 = skel[3][:, shave:]
     
        UL = np.concatenate((UL_crop2, UL_hm.astype(np.uint8)*255), axis=1)
        UR = np.concatenate((UR_crop2, UR_hm.astype(np.uint8)*255), axis=0)
        LR = np.concatenate((LR_hm.astype(np.uint8)*255, LR_crop2), axis=1)
        LL = np.concatenate((LL_hm.astype(np.uint8)*255, LL_crop2), axis=0)
    
        X_skel = np.vstack((np.hstack((UL, UR)), np.hstack((LL, LR))))
        
        # Calculate final vd
        
        vein_density = sum(vein)
        
        # State it was split
        split = "split"
    
        ######################################################################
        ### If the image is poorly stained with a low vein density (hence benefits
        ### from the splitting) then it will have a high numel(n) (number of
        ### endpoints), therefore this will be used as the cutoff to decide if an
        ### image needs further processing.
    
       
    
        if ((split_mode == "auto") and (B.sum() < branch_threshold) and (E.sum() < endpoint_threshold)) or ((split_mode == "auto") and (sd > sd_threshold)):
        
            # State it is not split
            split = "whole"
            # Improve contrast
            X_con = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8)).apply(X_gray)
            # Blur
            bf = initial_blur_factor
            X_blur = cv2.blur(X_con, ksize=(bf, bf))
            
            # Binarize
            _, X_bin = cv2.threshold(X_blur, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
            # Invert
            X_inv = cv2.bitwise_not(X_bin)
            # Make so binary true and false not binary 255 and 0
            X_inv = X_inv.astype(dtype=bool)
            # Skeletonize image
            X_skel = skimage.morphology.skeletonize(X_inv)
            # Convert back to unit8
            X_skel = X_skel.astype(np.uint8)
            X_skel*=255
             
            # Trimming
            # Search branchpoints
            B = pcv.morphology.find_branch_pts(X_skel)
            B = B.astype(bool)
            
            # Search endpoints
            
            E = pcv.morphology.find_tips(X_skel)
            E = E.astype(bool)     
            
            
            # Initial skeleton trim
            T = trim_factor
            skelD, segmented_img, segment_objects = pcv.morphology.prune(skel_img=X_skel, size=T)
           
            # Calculate initial VD
            skel_tot = skelD > 20 # Returns bool
            vd = (skel_tot.sum()*(pixel_length_um))/image_area_mm2 
            
            # Calculate the ratio of branchpoints:Initial_VD
            #Based on branchpoints before trimming as this gives a better view of how many superfluous ends there are
            ratio = (B.sum()/vd);
            
            # If this ratio is too high (>trim_threshold), then trimming will occur again,
            # based on the mode length of loose branch ends instead of the
            # "trim_factor". If we did mode from the begining then we would risk
            # removing data from those with low VD.
            
            if (ratio > trim_threshold):
            
                skelD2, seg_img2, seg_obj2 = pcv.morphology.prune(skel_img=X_skel, size = 0)
                branch_ends, core_trace = pcv.morphology.segment_sort(skel_img=skelD2, objects=seg_obj2)
                # Get lengths of branch ends
                labeled_img  = pcv.morphology.segment_path_length(segmented_img=seg_img2, objects=branch_ends, label="default")
                branch_lengths = pcv.outputs.observations['default']['segment_path_length']['value']
            
                # Set trimming threshold - we add a number of pixels proportionate to the image length so we make sure to capture all of them, i.e. we take the mode length + a bit extra as the trimming factot
                T = statistics.mode(branch_lengths) + (x_pixels/100)
                # Trim
                skelD, segmented_img, segment_objects = pcv.morphology.prune(skel_img=X_skel, size=T)
                # Re-calculate VD
                skel_tot = skelD > 20 # Returns bool
                vd = (skel_tot.sum()*(pixel_length_um))/image_area_mm2 
                
                X_skel = skelD
                
                # Define how it has been processed split or not, trim factor, blur factor and 
                Processed = split + '_tf:' + str(round(T)) + '_bf:' + str(round(bf)) + '_FILE:' + filename
                
                # Display image
                pcv.params.line_thickness = 10
                image,_ = pcv.morphology.segment_skeleton(skel_img=X_skel, mask = X_gray)
                
                plt.imshow(image)
                plt.title(Processed)
                plt.show
                
                vein_density = vd
                
            else: # no further processing
            # Define how it has been processed split or not, trim factor, blur factor and 
                Processed = split + '_tf:' + str(round(T)) + '_bf:' + str(round(bf)) + '_FILE:' + filename
            
            # Display image
                pcv.params.line_thickness = 10
                image,_ = pcv.morphology.segment_skeleton(skel_img=X_skel, mask = X_gray)
            
                plt.imshow(image)
                plt.title(Processed)
                plt.show
                
                vein_density = vd
            
        else: # If the figure was fine after splitting
            # Define how it has been processed split or not, trim factor, blur factor and 
            Processed = split + '_tf:' + str(round(T)) + '_bf:' + str(round(bf)) + '_FILE:' + filename
        
            # Display image
            pcv.params.line_thickness = 10
            image,_ = pcv.morphology.segment_skeleton(skel_img=X_skel, mask = X_gray)
        
            plt.imshow(image)
            plt.title(Processed)
            plt.show
            vein_density = vd        
    
    # Now we run for if the split mode is selected as "whole"
    elif (split_mode == "whole"):
        # State it is not split
        split = "whole"
        # Improve contrast
        X_con = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8)).apply(X_gray)
        # Blur
        bf = initial_blur_factor
        X_blur = cv2.blur(X_con, ksize=(bf, bf))
        
        # Binarize
        _, X_bin = cv2.threshold(X_blur, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        # Invert
        X_inv = cv2.bitwise_not(X_bin)
        # Make so binary true and false not binary 255 and 0
        X_inv = X_inv.astype(dtype=bool)
        # Skeletonize image
        X_skel = skimage.morphology.skeletonize(X_inv)
        # Convert back to unit8
        X_skel = X_skel.astype(np.uint8)
        X_skel*=255
    
        # Trimming
        
        # Search branchpoints
        B = pcv.morphology.find_branch_pts(X_skel)
        B = B.astype(bool)
            
        # Search endpoints
        
        E = pcv.morphology.find_tips(X_skel)
        E = E.astype(bool)     
    
        # Initial skeleton trim
        T = trim_factor
        skelD, segmented_img, segment_objects = pcv.morphology.prune(skel_img=X_skel, size=T)
        
        # Calculate initial VD
        skel_tot = skelD > 20 # Returns bool
        vd = (skel_tot.sum()*(pixel_length_um))/image_area_mm2 
        
        # Calculate the ratio of branchpoints:Initial_VD
        #Based on branchpoints before trimming as this gives a better view of how many superfluous ends there are
        ratio = (B.sum()/vd);
        
        # If this ratio is too high (e.g. >trim_threshold), then trimming will occur again,
        # based on the mode length of loose branch ends instead of the
        # "trim_factor". If we did mode from the begining then we would risk
        # removing data from those with low VD.
        
        if (ratio > trim_threshold):
        
            skelD2, seg_img2, seg_obj2 = pcv.morphology.prune(skel_img=X_skel, size = 0)
            branch_ends, core_trace = pcv.morphology.segment_sort(skel_img=skelD2, objects=seg_obj2)
            # Get lengths of branch ends
            labeled_img  = pcv.morphology.segment_path_length(segmented_img=seg_img2, objects=branch_ends, label="default")
            branch_lengths = pcv.outputs.observations['default']['segment_path_length']['value']
        
            # Set trimming threshold - we add a number of pixels proportionate to the image length so we make sure to capture all of them, i.e. we take the mode length + a bit extra as the trimming factot
            T = statistics.mode(branch_lengths) + (x_pixels/100)
            # Trim
            skelD, segmented_img, segment_objects = pcv.morphology.prune(skel_img=X_skel, size=T)
            # Re-calculate VD
            skel_tot = skelD > 20 # Returns bool
            vd = (skel_tot.sum()*(pixel_length_um))/image_area_mm2 
            
            X_skel = skelD
            
            # Define how it has been processed split or not, trim factor, blur factor and 
            Processed = split + '_tf:' + str(round(T)) + '_bf:' + str(round(bf)) + '_FILE:' + filename
            
            # Display image
            pcv.params.line_thickness = 10
            image,_ = pcv.morphology.segment_skeleton(skel_img=X_skel, mask = X_gray)
            
            plt.imshow(image)
            plt.title(Processed)
            plt.show
            
            vein_density = vd
        
        else: # If the figure was fine after no initial trim
        # Define how it has been processed split or not, trim factor, blur factor and 
            Processed = split + '_tf:' + str(round(T)) + '_bf:' + str(round(bf)) + '_FILE:' + filename
    
            # Display image
            pcv.params.line_thickness = 10
            image,_ = pcv.morphology.segment_skeleton(skel_img=X_skel, mask = X_gray)
    
            plt.imshow(image)
            plt.title(Processed)
            plt.show
    
            vein_density = vd
    
    # Now for processing if monocot is selected
    else:
        if (split_mode == "monocot"):
            # State it is not split
            split = "monocot"
            # Improve contrast
            X_con = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8)).apply(X_gray)
            # Blur
            bf = initial_blur_factor
            X_blur = cv2.blur(X_con, ksize=(bf, bf))
            
            # Binarize
            _, X_bin = cv2.threshold(X_blur, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
            # Invert
            X_inv = cv2.bitwise_not(X_bin)
            # Make so binary true and false not binary 255 and 0
            X_inv = X_inv.astype(dtype=bool)
            # Skeletonize image
            X_skel = skimage.morphology.skeletonize(X_inv)
            # Convert back to unit8
            X_skel = X_skel.astype(np.uint8)
            X_skel*=255
        
         
    
            
            # Initial skeleton trim
            T = trim_factor
            skelD, segmented_img, segment_objects = pcv.morphology.prune(skel_img=X_skel, size=T)
            
            # Because commissural veins complicate analysis, we need to remove them
            # Find branch points
            branch_points = pcv.morphology.find_branch_pts(skelD)
            
            # Get coordinates of branchpoints
            bp_coords = np.column_stack(np.where(branch_points > 0))
            
            # Distances between branch_points
            bp_dists = distance.cdist(bp_coords, bp_coords, 'euclidean')
            
            n,_ = bp_coords.shape
            
            # Create empty mask for commisural veins
            height, width = skelD.shape
            comi = np.zeros((height, width))
            comi = comi.astype(np.uint8)
            
            for i in range(n):
                try:
                    min_dist = sorted(set(bp_dists[i]))[1]  # This gets the second smallest    value as the smallest will always be 0
                except IndexError:
                    continue  # Skip this iteration if the index is out of range
                        # if this distance is too small, then we do not want to remove it as sometimes small branches are correct and occur due to other artifacts picked up
                if min_dist < mono_min_param:  # do nothing
                    pass
                elif min_dist > monocot_branch_length_to_keep:  # do nothing if large as may indicate an actual vein which we do not want to remove
                    pass
                else:
                    bp_dists[i] = bp_dists[i] + 100  # can be any number, we just need to change this array so it does not match with it in the next step
                    for j in range(n):
                        if np.any(bp_dists[j] == min_dist):  # This searches the other arrays for the same distance value, when it returns true it means the two branch points are connected by the minimum branch length - i.e. they are closest together and need interrogating
                            # Define coordinates of both
                            bp_y, bp_x = bp_coords[i]
                            bp_y2, bp_x2 = bp_coords[j]
                            start_coord = (bp_y, bp_x)
                            end_coord = (bp_y2, bp_x2)
                        # Identify all the coordinates that connect these in the skeleton
                            connected_pixels = find_connected_pixels(skeleton_img=skelD, start_coord=start_coord, end_coord=end_coord, exp_elucidian_dist=(min_dist+75))  # add 75 pixels to account for differences in the connected length and the elucidean distance
                        # make mask of these pixels
                            height, width = skelD.shape
                            remove = np.zeros((height, width))
                            remove = remove.astype(np.uint8)
                            for p in connected_pixels:
                                remove[p] = 255

                                # create mask for commisural veins
                                comi = comi + remove

            
            # Now remove these pixels from the skeleton completely
            comi = comi > 10 # could be any small number
            comi = comi.astype(np.uint8)
            comi = comi*255
            
            # Initial estimate of VD - will be final estimate if not trimmed based on mode
            skelD = skelD - comi
            skel_tot = skelD > 10 # could be any small number
            vd = (skel_tot.sum()*(pixel_length_um))/image_area_mm2 
            
            # Calculate the ratio of branchpoints:Initial_VD
            #Based on branchpoints before trimming as this gives a better view of how many superfluous ends there are
            B = pcv.morphology.find_branch_pts(X_skel)
            ratio = (B.sum()/vd)
            
            if (ratio > trim_threshold):
            
                skelD2, seg_img2, seg_obj2 = pcv.morphology.prune(skel_img=X_skel, size = 0)
                branch_ends, core_trace = pcv.morphology.segment_sort(skel_img=skelD2, objects=seg_obj2)
                # Get lengths of branch ends
                labeled_img  = pcv.morphology.segment_path_length(segmented_img=seg_img2, objects=branch_ends, label="default")
                branch_lengths = pcv.outputs.observations['default']['segment_path_length']['value']
            
                # Set trimming threshold - we add a number of pixels proportionate to the image length so we make sure to capture all of them, i.e. we take the mode length + a bit extra as the trimming factot
                T = statistics.mode(branch_lengths) + (x_pixels/40) # Arbitrarily higher for veins as the mode tends to be much smaller
                # Trim
                skelD, segmented_img, segment_objects = pcv.morphology.prune(skel_img=X_skel, size=T)
                
                # Remove commisural veins
                # Find branch points
                branch_points = pcv.morphology.find_branch_pts(skelD)
                
                # Get coordinates of branchpoints
                bp_coords = np.column_stack(np.where(branch_points > 0))
                
                # Distances between branch_points
                bp_dists = distance.cdist(bp_coords, bp_coords, 'euclidean')
                
                n,_ = bp_coords.shape
                
                # Create empty mask for commisural veins
                comi = np.zeros((height, width))
                comi = comi.astype(np.uint8)
                
                for i in range(n):
                    min_dist = sorted(set(bp_dists[i]))[1]  # This gets the second smallest value as the smallest will always be 0
                    # if this distance is too small, then we do not want to remove it as sometimes small branches are correct and occur due to other artifacts picked up
                    if min_dist < mono_min_param:  # do nothing
                        pass
                    elif min_dist > monocot_branch_length_to_keep:  # do nothing if large as may indicate an actual vein which we do not want to remove
                        pass
                    else:
                        bp_dists[i] = bp_dists[i] + 100  # can be any number, we just need to change this array so it does not match with it in the next step
                        for j in range(n):
                            if np.any(bp_dists[j] == min_dist):  # This searches the other arrays for the same distance value, when it returns true it means the two branch points are connected by the minimum branch length - i.e. they are closest together and need interrogating
                                # Define coordinates of both
                                bp_y, bp_x = bp_coords[i]
                                bp_y2, bp_x2 = bp_coords[j]
                                start_coord = (bp_y, bp_x)
                                end_coord = (bp_y2, bp_x2)
                                # Identify all the coordinates that connect these in the skeleton
                                connected_pixels = find_connected_pixels(skeleton_img=skelD, start_coord=start_coord, end_coord=end_coord, exp_elucidian_dist = (min_dist+75))#add 75 pixels to account for differences in the connected length and the elucidean distance
                                # make mask of these pixels
                                height, width = skelD.shape
                                remove = np.zeros((height, width))
                                remove = remove.astype(np.uint8)
                                for p in connected_pixels:
                                    remove[p] = 255
    
                            # create mask for commisural veins
                                comi = comi + remove
                            
                
                # Now remove these pixels from the skeleton completely
                comi = comi > 10 # could be any small number
                comi = comi.astype(np.uint8)
                comi = comi*255
                
                # Initial estimate of VD - will be final estimate if not trimmed based on mode
                skelD = skelD - comi
                skel_tot = skelD > 10 # could be any small number
                vd = (skel_tot.sum()*(pixel_length_um))/image_area_mm2 
                
                # Read it back into X_skel
                X_skel = skelD
                
                # Define how it has been processed split or not, trim factor, blur factor and 
                Processed = split + '_tf:' + str(round(T)) + '_bf:' + str(round(bf)) + '_FILE:' + filename
                
                # Display image
                pcv.params.line_thickness = 10
                image,_ = pcv.morphology.segment_skeleton(skel_img=X_skel, mask = X_gray)
                
                plt.imshow(image)
                plt.title(Processed)
                plt.show
            
                vein_density = vd
            
            else: # If the figure was fine after no initial trim
                # Define how it has been processed split or not, trim factor, blur factor and 
                X_skel = skelD
                
                Processed = split + '_tf:' + str(round(T)) + '_bf:' + str(round(bf)) + '_FILE:' + filename
     
                # Display image
                pcv.params.line_thickness = 10
                image,_ = pcv.morphology.segment_skeleton(skel_img=X_skel, mask = X_gray)
     
                plt.imshow(image)
                plt.title(Processed)
                plt.show
                
                vein_density = vd
                
    return vein_density, Processed, image
