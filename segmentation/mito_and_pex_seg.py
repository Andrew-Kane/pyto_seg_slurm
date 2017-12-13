import os
import sys
sys.path.append('/n/home06/akane/code')
import argparse
from pyto_segmenter import PexSegment, MitoSegment
from scipy.ndimage.morphology import binary_dilation, generate_binary_structure
from scipy.ndimage.morphology import iterate_structure

parser = argparse.ArgumentParser(description = 'Segment peroxisomes from \
                                 images and return pickled objects.')
parser.add_argument('-d', '--img_dir', required = True, 
                    help = 'directory containing images to segment.')
parser.add_argument('-ht', '--high_threshold', required = True,
                    help = 'high threshold for canny segmentation.')
parser.add_argument('-lt', '--low_threshold', required = True,
                    help = 'low threshold for canny segmentation.')
parser.add_argument('-n', '--array_number', required = True,
                    help = 'the job number within the job array.')
parser.add_argument('-a', '--array_length', required = True,
                    help = 'length of the job array.')

args = parser.parse_args()
print(args)
img_dir = args.img_dir
high_threshold = int(args.high_threshold)
low_threshold = int(args.low_threshold)
array_n = int(args.array_number)
array_l = int(args.array_length)

def main():
    os.chdir(img_dir)
    flist = os.listdir()
    imgs = [f for f in flist if '.tif' in f.lower()]
    pex_imgs = [im for im in imgs if '594' in im]
    mito_imgs = [im for im in imgs if '447' in im]
    pex_imgs.sort()
    mito_imgs.sort()
    if len(pex_imgs) != len(mito_imgs):
        raise ValueError('Length of pex and mito img sets do not match.')
    ims_per_job = int(len(pex_imgs)/array_l)
    split_pex_list = []
    split_mito_list = []
    for i in range(0, len(pex_imgs), ims_per_job):
        split_pex_list.append(pex_imgs[i:i+ims_per_job])
        split_mito_list.append(mito_imgs[i:i+ims_per_job])
    n_leftover = len(pex_imgs)%array_l
    if n_leftover != 0:
        leftover_pex = pex_imgs[-n_leftover:]
        leftover_mito = mito_imgs[-n_leftover:]
        for x in range(0,len(leftover_pex)):
            split_pex_list[x].append(leftover_pex[x])
            split_mito_list[x].append(leftover_mito[x])
    mito_list = split_mito_list[array_n]
    pex_list = split_pex_list[array_n]
    print('Peroxisome images:')
    print(pex_list)
    for i in range(0,len(pex_list)):
        os.chdir(img_dir)
        print('SEGMENTING ' + pex_list[i] + ' and ' + mito_list[i])
        pex_segmenter = PexSegment.PexSegmenter(pex_list[i], seg_method = 'canny',
                                               high_threshold = high_threshold,
                                               low_threshold = low_threshold)
        pex_obj = pex_segmenter.segment()
        pex_obj.rm_border_objs()
        pex_obj.pickle(output_dir = img_dir + '/pickles', filename =
                       pex_obj.filename[0:pex_obj.filename.index('.tif')] + '_pex.pickle' )
        os.chdir(img_dir)
        mito_segmenter = MitoSegment.MitoSegmenter(mito_list[i], seg_method = 'canny',
                                                   high_threshold = 250,
                                                   low_threshold = 125,
                                                   min_cutoff = 2300)
        mito_obj = mito_segmenter.segment()
        mito_obj.rm_border_objs()
        mito_obj.mitochondria[pex_obj.threshold_img == 1] = 0 # remove overlap
        mito_obj.pickle(output_dir = img_dir + '/pickles', filename = 
                        mito_obj.filename[0:mito_obj.filename.index('.tif')] + '_mito.pickle' )

        del pex_obj
        del mito_obj

if __name__ == '__main__':
    main()
