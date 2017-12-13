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
    mch_imgs = [y for y in files if '594' in y]
    gfp_imgs = [y for y in files if '488' in y]
    output_frame = pd.DataFrame({'img': [],
                                 'obj_channel': [],
                                 'obj_number': [],
                                 'puncta': [],
                                    })
    #generates pickles to be analyzed in this array
    pickles = get_pickle_set(img_dir, array_l, array_n) 
    #generates IDs to properly assign cell and foci.
    (pickle_id, pickle_channel) = get_img_ids(pickles,return_channel = True)    
    gfp_ids = get_img_ids(gfp_imgs)
    mch_ids = get_img_ids(mch_imgs)
    #generates a cell, foci pair of the correct IDs
    for key, value in pickle_id.items():
        foci_cts_dict = {}
        cell_pkl = mch_ids[key]
        foci_pkl = gfp_ids[key]
        cellfile = open(cell_pkl,'rb')
        cellobj = pickle.load(cellfile)
        focifile = open(foci_pkl,'rb')
        fociobj = pickle.load(focifile)
        cells = cellobj.final_cells
        seg_foci = fociobj.peroxisomes
        for i in np.unique(cells):  # gets the unique IDs of the cells
            if i == 0:  # if it's bgrd
                continue
            else:
                overlapping_foci = np.unique(seg_foci[cells == i])  # gets the unique IDs for each focus in the cell, as well as the bgrd
                num_foci = overlapping_foci.shape[0]-1  # shape[0] is the length of the overlapping_foci array; subtract bgrd 
                foci_cts_dict[i] = num_foci  # add a key:value pair to foci_cts_dict where key is cell # and val is num of foci
        volumes_v2 = {}
        foci = {}
        for obj in cellobj.obj_nums:
            print('     current obj number: ' + str(obj))
            volumes_v2[obj] = len(np.flatnonzero(cells == obj))
            foci[obj] = foci_cts_dict[obj]  
        print('')
        currimg_data = pd.DataFrame({'img': pd.Series(data =
                                                          [cellobj.filename]*len(cellobj.obj_nums),
                                                          index = cellobj.obj_nums),
                                         'obj_channel': pd.Series(data =
                                                                  pickle_channel[cell_pkl],
                                                              index = cellobj.obj_nums),
                                         'obj_number': pd.Series(data = cellobj.obj_nums,
                                                                 index = cellobj.obj_nums),
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
    pickle_list = [p for p in os.listdir() if '594' in p]
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

