import os
import sys
sys.path.append('/n/home06/akane/code')
import argparse
from pyto_segmenter import PexSegment, MitoSegment
import numpy as np
import pandas as pd
from skimage import io
import pickle
import re



parser = argparse.ArgumentParser(description = 'Measure the amount of \
                                 aggregates in each cell object.')
parser.add_argument('-d', '--img_dir', required = True, 
                    help = 'directory containing images to segment.')
parser.add_argument('-n', '--array_number', required = True,
                    help = 'the job number within the job array.')
parser.add_argument('-a', '--array_length', required = True,
                    help = 'length of the job array.')
parser.add_argument('-t', '--type', required = False,
                    default = 'apc', help = 'type: apc (default, aggregates per cell) or binary.')

args = parser.parse_args()
print(args)
img_dir = args.img_dir
array_n = int(args.array_number)
array_l = int(args.array_length)
expt_type = args.type

def main():
    os.chdir(img_dir+ '/pickles')
    files = [f for f in os.listdir() if '.pickle' in f.lower()]
    os.chdir(img_dir)
    imgs = [f for f in os.listdir() if '.tif' in f.lower()]
    cfp_pickles = [y for y in files if '447' in y]
    yfp_pickles = [y for y in files if '515' in y]
    yfp_imgs = [y for y in imgs if '515' in y]
    rfp_imgs = [y for y in imgs if '594' in y]
    output_frame = pd.DataFrame({'img': [],
                                 'obj_channel': [],
                                 'obj_number': [],
                                 'puncta': [],
                                 'yfp_mean':[],
                                 'yfp_stdev':[],
                                 'rfp_mean':[],
                                 'rfp_stdev':[]
                                    })
    #generates pickles to be analyzed in this array
    pickles = get_pickle_set(img_dir, array_l, array_n) 
    #generates IDs to properly assign cell and foci.
    (pickle_id, pickle_channel) = get_img_ids(pickles,return_channel = True)    
    yfp_ids = get_img_ids(yfp_pickles)
    cfp_ids = get_img_ids(cfp_pickles)
    yfp_im_ids = get_img_ids(yfp_imgs)
    rfp_ids = get_img_ids(rfp_imgs)
    #generates a cell, foci pair of the correct IDs
    for key, value in pickle_id.items():
        os.chdir(img_dir+ '/pickles')
        foci_cts_dict = {}
        vacuole_pkl = cfp_ids[key]
        foci_pkl = yfp_ids[key]
        print(os.getcwd())
        vacfile = open(vacuole_pkl,'rb')
        vacobj = pickle.load(vacfile)
        focifile = open(foci_pkl,'rb')
        fociobj = pickle.load(focifile)
        vacuoles = vacobj.final_cells
        seg_foci = fociobj.peroxisomes
        for i in np.unique(vacuoles):  # gets the unique IDs of the cells
            if i == 0:  # if it's bgrd
                continue
            else:
                overlapping_foci = np.unique(seg_foci[vacuoles == i])  # gets the unique IDs for each focus in the cell, as well as the bgrd
                num_foci = overlapping_foci.shape[0]-1  # shape[0] is the length of the overlapping_foci array; subtract bgrd 
                foci_cts_dict[i] = num_foci  # add a key:value pair to foci_cts_dict where key is cell # and val is num of foci
        volumes_v2 = {}
        foci = {}
        os.chdir(img_dir)
        print('current YFP image file: ' + yfp_im_ids[key])
        yfp_img = io.imread(yfp_im_ids[key])
        yfp_mean = {}
        yfp_stdev = {}
        print('current RFP image file: ' + rfp_ids[key])
        rfp_img = io.imread(rfp_ids[key])
        rfp_mean = {}
        rfp_stdev = {}
        for obj in vacobj.obj_nums:
            print('     current obj number: ' + str(obj))
            volumes_v2[obj] = len(np.flatnonzero(vacuoles == obj))
            foci[obj] = foci_cts_dict[obj]
            yfp_mean[obj] = np.mean(yfp_img[vacobj.final_cells == obj])
            yfp_stdev[obj] = np.std(yfp_img[vacobj.final_cells == obj])
            rfp_mean[obj] = np.mean(rfp_img[vacobj.final_cells == obj])
            rfp_stdev[obj] = np.std(rfp_img[vacobj.final_cells == obj])
        print('')
        currimg_data = pd.DataFrame({'img': pd.Series(data =
                                                          [vacobj.filename]*len(vacobj.obj_nums),
                                                          index = vacobj.obj_nums),
                                         'obj_channel': pd.Series(data =
                                                                  pickle_channel[vacuole_pkl],
                                                              index = vacobj.obj_nums),
                                         'obj_number': pd.Series(data = vacobj.obj_nums,
                                                                 index = vacobj.obj_nums),
                                         'volume': volumes_v2,
                                         'puncta': foci,
                                        })
        output_frame = pd.concat([output_frame, currimg_data])
    print('')
    print('-----------------------------------------------------------------')
    print('-----------------------------------------------------------------')
    print('')
    print('saving data...')
    if not os.path.isdir(img_dir + '/analysis_output'):
        os.mkdir(img_dir + '/analysis_output')
    output_frame.to_csv(img_dir + '/analysis_output/' + str(array_n) +
                        '_analysis_output.csv')


def get_pickle_set(img_dir, array_l, array_n):
    '''Get the subset of pickles for analysis by a given instance.'''
    os.chdir(img_dir+ '/pickles')
    pickle_list = [p for p in os.listdir() if '447' in p]
    pickle_list.sort()
    pickles_per_job = int(len(pickle_list)/array_l)
    split_pickle_list = []
    for i in range(0, len(pickle_list), pickles_per_job):
        split_pickle_list.append(pickle_list[i:i+pickles_per_job])
    n_leftover = len(pickle_list)%array_l
    if n_leftover != 0:
        leftover_pickles = pickle_list[-n_leftover:]
        for x in range(0,len(leftover_pickles)):
            split_pickle_list[x].append(leftover_pickles[x])
    return(split_pickle_list[array_n])

def get_img_ids(img_files, return_channel = False):
    '''Extracts image filenames lacking wavelength identifiers.'''
    channel_re = re.compile('^w\d+[A-Za-z]*[ .]')
    img_ids = []
    channels = []
    for img in img_files:
        print('_______________________________________________________')
        print('     generating image identifier for ' + img)
        img = img.replace('pickled_','')
        img = img.replace('.pickle','')
        img = img.replace('.tif','')
        split_im = img.split('_')
        rm_channel = '_'.join([i for i in split_im if re.search(channel_re, i)
                               == None])
        channel = [i for i in split_im if re.search(channel_re, i)]
        if len(channel) > 1:
            print('WARNING: more than one match to channel ID string!')
        channel = channel[0].split('.')[0]
        channel = channel[-3:]
        print('     image identifier: ' + rm_channel)
        print('     channel: ' + channel)
        img_ids.append(rm_channel)
        channels.append(channel)
    fname_id_dict = dict(zip(img_ids, img_files))
    channel_dict = dict(zip(img_files, channels))
    print('')
    print('done extracting image identifiers.')
    if not return_channel:
        return(fname_id_dict)
    if return_channel:
        return((fname_id_dict, channel_dict))

if __name__ == '__main__':
    main()

