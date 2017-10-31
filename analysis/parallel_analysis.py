import os
import sys
sys.path.append('/n/denic_lab/Users/nweir/python_packages/')
import argparse
from pyto_segmenter import PexSegment, MitoSegment
import numpy as np
import pandas as pd
from skimage import io
import pickle
import re



parser = argparse.ArgumentParser(description = 'Measure fluorescence \
                                 intensity at pexs and mitochondria.')
parser.add_argument('-d', '--img_dir', required = True, 
                    help = 'directory containing images to segment.')
parser.add_argument('-n', '--array_number', required = True,
                    help = 'the job number within the job array.')
parser.add_argument('-a', '--array_length', required = True,
                    help = 'length of the job array.')
parser.add_argument('-t', '--type', required = False,
                    default = 'yfp', help = 'type: yfp (default) or tft.')

args = parser.parse_args()
print(args)
img_dir = args.img_dir
array_n = int(args.array_number)
array_l = int(args.array_length)
expt_type = args.type

def main():
    os.chdir(img_dir)
    files = [f for f in os.listdir() if '.tif' in f.lower()]
    yfp_imgs = [y for y in files if '515' in y]
    gfp_imgs = [y for y in files if '488' in y]
    if expt_type == 'tft':
        mch_imgs = [m for m in files if '594' in m]
        output_frame = pd.DataFrame({'img': [],
                                     'obj_channel':[],
                                     'obj_number': [],
                                     'volume': [],
                                     'yfp': [],
                                     'mcherry': []
                                    })
    if expt_type == 'yfp':
        output_frame = pd.DataFrame({'img': [],
                                     'obj_channel': [],
                                     'obj_number': [],
                                     'volume': [],
                                     'yfp': []
                                    })
    if expt_type == 'stdev':
        output_frame = pd.DataFrame({'img': [],
                                     'obj_channel': [],
                                     'obj_number': [],
                                     'volume': [],
                                     'gfp_mean': [],
                                     'gfp_stdev': []
                                    })
    pickles = get_pickle_set(img_dir, array_l, array_n)
    if expt_type != 'stdev':
        yfp_ids = get_img_ids(yfp_imgs)
        yfp_ids = {v: k for k, v in yfp_ids.items()}
        if expt_type == 'tft':
            mch_ids = get_img_ids(mch_imgs)
            mch_ids = {v: k for k, v in mch_ids.items()}
    for p in pickles:
        os.chdir(img_dir + '/pickles')
        cfile = open(p,'rb')
        cpickle = pickle.load(cfile)
        print('current segmented object image: ' + cpickle.filename)
        (pickle_id, pickle_channel) = get_img_ids([cpickle.filename],
                                                  return_channel = True)
        pickle_id = pickle_id[cpickle.filename]
        pickle_channel = pickle_channel[cpickle.filename]
        if hasattr(cpickle, 'mitochondria'):
            obj_type = 'mito'
        elif hasattr(cpickle, 'peroxisomes'):
            obj_type = 'pex'
        elif hasattr(cpickle, 'final_cells'):
            obj_type = 'cells'
        if expt_type != 'stdev':
            print('current segmented object identifier: ' + pickle_id)
            print('current segmented object channel: ' + pickle_channel)
            os.chdir(img_dir)
            print('current YFP image file: ' + yfp_ids[pickle_id])
            yfp_img = io.imread(yfp_ids[pickle_id])
            if expt_type == 'tft':
                print('current mCherry image file: ' + mch_ids[pickle_id])
                mch_img = io.imread(mch_ids[pickle_id])
            yfp_vals = {}
            volumes_v2 = {}
            if expt_type == 'tft':
                mch_vals = {}
            for obj in cpickle.obj_nums:
                print('     current obj number: ' + str(obj))
                if obj_type == 'mito':
                    yfp_vals[obj] = sum(yfp_img[cpickle.mitochondria == obj])
                    volumes_v2[obj] = len(np.flatnonzero(cpickle.mitochondria == obj))
                    if expt_type == 'tft':
                        mch_vals[obj] = sum(mch_img[cpickle.mitochondria == obj])
                else:
                    yfp_vals[obj] = sum(yfp_img[cpickle.peroxisomes == obj])
                    volumes_v2[obj] = len(np.flatnonzero(cpickle.peroxisomes == obj))
                    if expt_type == 'tft':
                        mch_vals[obj] = sum(mch_img[cpickle.peroxisomes == obj])
                print('     yfp intensity: ' + str(yfp_vals[obj]))
                if expt_type == 'tft':
                    print('     mcherry intensity: ' + str(mch_vals[obj]))
        else:
            print('current segmented object identifier: ' + pickle_id)
            print('current segmented object channel: ' + pickle_channel)
            os.chdir(img_dir)
            print('current GFP image file: ' + gfp_ids[pickle_id])
            gfp_img = io.imread(gfp_ids[pickle_id])
            gfp_mean = {}
            volumes_v2 = {}
            gfp_stdev = {}
            for obj in cpickle.obj_nums:
                print('     current obj number: ' + str(obj))
                if obj_type == 'cells':
                    gfp_mean[obj] = np.mean(gfp_img[cpickle.final_cells == obj])
                    volumes_v2[obj] = len(np.flatnonzero(cpickle.final_cells == obj))
                    gfp_stdev[obj] = np.std(gfp_img[cpickle.final_cells == obj])

            
        print('')
        print('Appending data to output...')
        if expt_type == 'yfp':
            currimg_data = pd.DataFrame({'img': pd.Series(data =
                                                          [cpickle.filename]*len(cpickle.obj_nums),
                                                          index = cpickle.obj_nums),
                                         'obj_channel': pd.Series(data =
                                                                  [pickle_channel]*len(cpickle.obj_nums),
                                                              index = cpickle.obj_nums),
                                         'obj_number': pd.Series(data = cpickle.obj_nums,
                                                                 index = cpickle.obj_nums),
                                         'volume': volumes_v2,
                                         'yfp': yfp_vals,
                                        })
        elif expt_type == 'tft':
            currimg_data = pd.DataFrame({'img': pd.Series(data =
                                                          [cpickle.filename]*len(cpickle.obj_nums),
                                                          index = cpickle.obj_nums),
                                         'obj_channel': pd.Series(data =
                                                                  [pickle_channel]*len(cpickle.obj_nums),
                                                              index = cpickle.obj_nums),
                                         'obj_number': pd.Series(data = cpickle.obj_nums,
                                                                 index = cpickle.obj_nums),
                                         'volume': volumes_v2,
                                         'yfp': yfp_vals,
                                         'mcherry': mch_vals,
                                        })
        elif expt_type == 'stdev':
            currimg_data = pd.DataFrame({'img': pd.Series(data =
                                                          [cpickle.filename]*len(cpickle.obj_nums),
                                                          index = cpickle.obj_nums),
                                         'obj_channel': pd.Series(data =
                                                                  [pickle_channel]*len(cpickle.obj_nums),
                                                              index = cpickle.obj_nums),
                                         'obj_number': pd.Series(data = cpickle.obj_nums,
                                                                 index = cpickle.obj_nums),
                                         'volume': volumes_v2,
                                         'gfp_mean': gfp_mean,
                                         'gfp_stdev': gfp_stdev,
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
    os.chdir(img_dir + '/pickles')
    pickle_list = [p for p in os.listdir() if '.pickle' in p]
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
    fname_id_dict = dict(zip(img_files, img_ids))
    channel_dict = dict(zip(img_files, channels))
    print('')
    print('done extracting image identifiers.')
    if not return_channel:
        return(fname_id_dict)
    if return_channel:
        return((fname_id_dict, channel_dict))

if __name__ == '__main__':
    main()

