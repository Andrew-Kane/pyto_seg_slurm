import os
import sys
sys.path.append('/n/home06/akane/code')
import argparse
from pyto_segmenter import CellSegment


parser = argparse.ArgumentParser(description = 'Segment cells from \
                                 images and return pickled objects.')
parser.add_argument('-d', '--img_dir', required = True, 
                    help = 'directory containing images to segment.')
parser.add_argument('-t', '--threshold', required = True,
                    help = 'threshold for segmentation.')
parser.add_argument('-n', '--array_number', required = True,
                    help = 'the job number within the job array.')
parser.add_argument('-a', '--array_length', required = True,
                    help = 'length of the job array.')

args = parser.parse_args()
print(args)
img_dir = args.img_dir
threshold = int(args.threshold)
array_n = int(args.array_number)
array_l = int(args.array_length)

def main():
    os.chdir(img_dir)
    flist = os.listdir()
    imgs = [f for f in flist if '.tif' in f.lower()]
    cell_imgs = [im for im in imgs if '447' in im]
    cell_imgs.sort()
    ims_per_job = int(len(cell_imgs)/array_l)
    split_cell_list = []
    for i in range(0, len(cell_imgs), ims_per_job):
        split_cell_list.append(cell_imgs[i:i+ims_per_job])
    n_leftover = len(cell_imgs)%array_l
    if n_leftover != 0:
        leftover_cell = cell_imgs[-n_leftover:]
        for x in range(0,len(leftover_cell)):
            split_cell_list[x].append(leftover_cell[x])
    cell_list = split_cell_list[array_n]
    print('Cell images:')
    print(cell_list)
    for i in range(0,len(cell_list)):
        os.chdir(img_dir)
        print('SEGMENTING ' + cell_list[i])
        cell_segmenter = CellSegment.CellSegmenter(cell_list[i],
                                               threshold = threshold)
        cell_obj = cell_segmenter.segment()
        cell_obj.pickle(output_dir = img_dir + '/pickles')
        os.chdir(img_dir)
        del cell_obj

if __name__ == '__main__':
    main()
