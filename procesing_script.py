
from aicsimageio import AICSImage
import numpy as np
from scipy import ndimage as ndi

from skimage import (
    color, feature, filters, measure, morphology, segmentation, util
)

from skimage.filters import meijering, sato, frangi, hessian

from skimage import data, filters, measure, morphology
from skimage import data, restoration, util
import os
import re
import pandas as pd
from aicsimageio.readers import CziReader
import psutil
from joblib import Parallel, delayed
from skimage.segmentation import expand_labels, watershed


#PATH TO DIRECTORY WITH IMAGE FILES - EACH DIRECTORY NAMED AS THE MARKER BEING QUANTIFIED SPLIT INTO TWO FOLDERS BASED ON CONDITION
BASE_PATH = ""
#PATH TO OUTPUT DIRECTORIES
BASE_OUT_PATH = ""

VCAD_CHANNEL = 2
INSULIN_CHANNEL = 3
DAPI_CHANNEL = 0
MARKER_CHANNEL = 1


#RUN EACH PATCH IN PARALLEL
def process(i,gene,path):
    raw_im = im_data[:,i,:,:]

    ####   REMOVE NOISE AND SEGMENT VCAD CHANNEL ######

    kernel = restoration.ellipsoid_kernel((12 , 12), 12 )

    
    background_vcad = restoration.rolling_ball(raw_im[VCAD_CHANNEL,:,:],num_threads=10, kernel=kernel )
    new_rawpatch_vcad = raw_im[VCAD_CHANNEL,:,:] - background_vcad

    smooth = filters.gaussian(new_rawpatch_vcad, sigma=1.5)
    thresh_val = filters.threshold_triangle(smooth)
    thresh = smooth > thresh_val
    footprint=[(np.ones((12, 1)), 1), (np.ones((1, 12)), 1)]
    dilate = morphology.binary_dilation(thresh,footprint)


    footprint=[(np.ones((10, 1)), 1), (np.ones((1, 10)), 1)]
    closing = morphology.binary_closing(dilate,footprint=footprint)


    mask = morphology.remove_small_objects(closing, 600)
    all_mask = np.ones(mask.astype(np.uint8).shape)
    
    labels = measure.label(mask)
    
    ####   REMOVE NOISE IN EACH RAW CHANNEL   ######
    
    background_intensity = restoration.rolling_ball(raw_im[MARKER_CHANNEL,:,:],num_threads=10, kernel=kernel )
    new_rawpatch_intensity = raw_im[MARKER_CHANNEL,:,:] - background_intensity

    background_ins = restoration.rolling_ball(raw_im[INSULIN_CHANNEL,:,:],num_threads=10, kernel=kernel )
    new_rawpatch_ins = raw_im[INSULIN_CHANNEL,:,:] - background_ins

    background_dapi = restoration.rolling_ball(raw_im[DAPI_CHANNEL,:,:],num_threads=10, kernel=kernel )
    new_rawpatch_dapi = raw_im[DAPI_CHANNEL,:,:] - background_dapi      
    

    ####   CALCULATE DENOISED INTENSITY   ######

    props = measure.regionprops(labels,intensity_image=new_rawpatch_intensity)
    properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']

    props_vcad = measure.regionprops(labels,intensity_image=new_rawpatch_vcad)
    properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']

    props_ins = measure.regionprops(labels,intensity_image=new_rawpatch_ins)
    properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']

    props_dapi = measure.regionprops(labels,intensity_image=new_rawpatch_dapi)
    properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']


    ####   CALCULATE RAW INTENSITY   ######

    props_ = measure.regionprops(labels,intensity_image=raw_im[MARKER_CHANNEL,:,:])
    properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']

    props_vcad_ = measure.regionprops(labels,intensity_image=raw_im[VCAD_CHANNEL,:,:])
    properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']

    props_ins_ = measure.regionprops(labels,intensity_image=raw_im[INSULIN_CHANNEL,:,:])
    properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']

    props_dapi_ = measure.regionprops(labels,intensity_image=raw_im[DAPI_CHANNEL,:,:])
    properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']



    df = dict()
    
    for index in range(0, labels.max()):
        label_i = props[index].label
        
        mean_int = getattr(props[index],'intensity_mean')
        area = getattr(props[index],'area')
        ecc = getattr(props[index],'eccentricity')
        mean_int_vcad = getattr(props_vcad[index],'intensity_mean')
        mean_int_ins = getattr(props_ins[index],'intensity_mean')
        mean_int_dapi = getattr(props_dapi[index],'intensity_mean')

        mean_int_ = getattr(props_[index],'intensity_mean')
        mean_int_vcad_ = getattr(props_vcad_[index],'intensity_mean')
        mean_int_ins_ = getattr(props_ins_[index],'intensity_mean')
        mean_int_dapi_ = getattr(props_dapi_[index],'intensity_mean')

        df[label_i] = [area,ecc,mean_int,mean_int_vcad,mean_int_ins,mean_int_dapi,mean_int_,mean_int_vcad_,mean_int_ins_,mean_int_dapi_]


    ####   CALCULATE TOTAL IMAGE INTENSITIES DENOISED   ######
    all_mask = np.ones(new_rawpatch_intensity.astype(np.uint8).shape)
    props_full = measure.regionprops_table(
        all_mask.astype(np.uint8),
        intensity_image=new_rawpatch_intensity,
        properties=('label', 'area', 'intensity_mean')
    )
    props_full_vcad = measure.regionprops_table(
        all_mask.astype(np.uint8),
        intensity_image=new_rawpatch_vcad,
        properties=('label', 'area', 'intensity_mean')
    )
    props_full_ins = measure.regionprops_table(
        all_mask.astype(np.uint8),
        intensity_image=new_rawpatch_ins,
        properties=('label', 'area', 'intensity_mean')
    )
    props_full_dapi = measure.regionprops_table(
        all_mask.astype(np.uint8),
        intensity_image=new_rawpatch_dapi,
        properties=('label', 'area', 'intensity_mean')
    )

    ####   CALCULATE TOTAL IMAGE INTENSITIES RAW   ######
    all_mask = np.ones(raw_im[VCAD_CHANNEL,:,:].astype(np.uint8).shape)
    props_full_raw = measure.regionprops_table(
        all_mask.astype(np.uint8),
        intensity_image=raw_im[MARKER_CHANNEL,:,:],
        properties=('label', 'area', 'intensity_mean')
    )
    props_full_vcad_raw = measure.regionprops_table(
        all_mask.astype(np.uint8),
        intensity_image=raw_im[VCAD_CHANNEL,:,:],
        properties=('label', 'area', 'intensity_mean')
    )
    props_full_ins_raw = measure.regionprops_table(
        all_mask.astype(np.uint8),
        intensity_image=raw_im[INSULIN_CHANNEL,:,:],
        properties=('label', 'area', 'intensity_mean')
    )
    props_full_dapi_raw = measure.regionprops_table(
        all_mask.astype(np.uint8),
        intensity_image=raw_im[DAPI_CHANNEL,:,:],
        properties=('label', 'area', 'intensity_mean')
    )


    
    l1 = [[i for val in list(df.keys())],[int(val) for val in list(df.keys())],[val[0] for val in list(df.values())],[np.round(val[1],3) for val in list(df.values())],[np.round(val[2],3) for val in list(df.values())],[np.round(val[3],3) for val in list(df.values())],[np.round(val[4],3) for val in list(df.values())],[np.round(val[5],3) for val in list(df.values())],[np.round(val[6],3) for val in list(df.values())],[np.round(val[7],3) for val in list(df.values())],[np.round(val[8],3) for val in list(df.values())],[np.round(val[9],3) for val in list(df.values())],[np.round(props_full['intensity_mean'],3)[0] for val in list(df.keys())],[np.round(props_full_vcad['intensity_mean'],3)[0] for val in list(df.keys())],[np.round(props_full_ins['intensity_mean'],3)[0] for val in list(df.keys())],[np.round(props_full_dapi['intensity_mean'],3)[0] for val in list(df.keys())],[np.round(props_full_raw['intensity_mean'],3)[0] for val in list(df.keys())],[np.round(props_full_vcad_raw['intensity_mean'],3)[0] for val in list(df.keys())],[np.round(props_full_ins_raw['intensity_mean'],3)[0] for val in list(df.keys())],[np.round(props_full_dapi_raw['intensity_mean'],3)[0] for val in list(df.keys())]]
    df_ = pd.DataFrame(l1)
    df_ = df_.transpose()

    raw_gene = "RAW_" + gene
    full_gene = "FULL_" + gene
    full_gene_raw = "FULL_" + gene + "_RAW"

    
    df_ = df_.rename(columns = {0 : "Patch",1 : "EC Region", 2 : "Area", 3: "Eccentricity", 4: gene , 5: "VCAD",6:"INS", 7: "DAPI",8:raw_gene,9:"RAW_VCAD",10:"RAW_INS",11:"RAW_DAPI",12:full_gene,13:"FULL_VCAD",14:"FULL_INS",15:"FULL_DAPI",16:full_gene_raw,17:"FULL_VCAD_RAW",18:"FULL_INS_RAW",19:"FULL_DAPI_RAW"})
    #WRITE DATAFRAME FOR EACH PATCH INCLUDING PATCH #, EC Region #, Area, Eccentricity, background normalized and raw intensity for each EC Region segmented and each whole patch
    df_.to_csv(f"{path}/{i}.csv",index=False)



for dir_ in os.listdir(BASE_PATH):
    dir_path = os.path.join(BASE_PATH, dir_)
    gene = dir_
    if os.path.isdir(dir_path):
        qdir = 'OUT_' + dir_
        if not os.path.exists(os.path.join(BASE_OUT_PATH,qdir)):
            os.makedirs(os.path.join(BASE_OUT_PATH,qdir))
        for subdir_ in os.listdir(dir_path):
            subdir_path = os.path.join(dir_path,subdir_)
            if os.path.isdir(subdir_path):
                if not os.path.exists(os.path.join(BASE_OUT_PATH,qdir,subdir_)):
                    os.makedirs(os.path.join(BASE_OUT_PATH,qdir,subdir_))
                    for file in os.listdir(subdir_path):
                        truncfile = file.replace(".czi","")

                        if not os.path.exists(os.path.join(BASE_OUT_PATH,qdir,subdir_,truncfile)):
                            os.makedirs(os.path.join(BASE_OUT_PATH,qdir,subdir_,truncfile))
                        outfolder_path = os.path.join(BASE_OUT_PATH,qdir,subdir_,truncfile)

                    
                        file_path = os.path.join(subdir_path,file)

                        #READ IN WHOLE SLIDE IMAGE TILED 
                        reader = CziReader(file_path)
                        tiles = reader.data.shape[1]
                        im_data = reader.data

                        #RUN EACH PATCH IN PARALLEL
                        current_process = psutil.Process()
                        subproc_before = set([p.pid for p in current_process.children(recursive=True)])
                        grouped_data = Parallel(n_jobs=10)(delayed(process)(i,gene,outfolder_path) for i in range(tiles))
                        subproc_after = set([p.pid for p in current_process.children(recursive=True)])
                        for subproc in subproc_after - subproc_before:
                            print('Killing process with pid {}'.format(subproc))
                            psutil.Process(subproc).terminate()
                        
                        #COMBINE EACH PATCH DATAFRAME INTO ONE WHOLE DATAFRAME
                        out_csv_file_name = file.replace(".czi",".csv")
                        out_csv_file = os.path.join(BASE_OUT_PATH,qdir,subdir_,out_csv_file_name)
                        file_extension = '.csv'
                        csv_file_list = []
                        for root, dirs, files in os.walk(outfolder_path):
                            for name in files:
                                if name.endswith(file_extension):
                                    csvfile_path = os.path.join(root, name)
                                    csv_file_list.append(csvfile_path)
                        
                        big_df = pd.concat((pd.read_csv(f) for f in csv_file_list))
                        big_df = big_df.reset_index(drop=True)
                        big_df.to_csv(out_csv_file)