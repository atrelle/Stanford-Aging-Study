"""
MRIqc directory structure:
sub-<participant_label>[/ses-<session_label>]/<data_type>/

<data_type> = 'func' | 'anat'
"""

from glob import glob
from os import makedirs, symlink
from os.path import split, splitext

def make_links(subs2run):
    for subject in subs2run:
        dst_dir = '/share/awagner/AM/data/MRIqc/data/sub-{}/ses-01/func/'.format(subject)
        makedirs(dst_dir)

        source_files = glob('/share/awagner/AM/data/{}/bold/*.nii.gz'.format(subject))
        source_files.sort()

        for src in source_files:
            # drop path and extension. run split_text twice because of .nii.gz
            src_file_name = splitext(splitext(split(src)[1])[0])[0]

            dst_file_name = 'sub-{}_task-{}_bold.nii.gz'.format(subject, src_file_name)
            dst = dst_dir + dst_file_name

            symlink(src, dst)

    return
