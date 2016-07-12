#!/usr/bin/env python


# GiniClust/src/GiniClust/__main__.py 


# Author: Gregory Giecold
# Maintainer: Qian Zhu
# Affiliations: Harvard University & Princeton University
# Contact information: ggiecold@jimmy.harvard.edu; qzhu@princeton.edu


from __future__ import print_function

from os import chdir, getcwd, makedirs, path
import subprocess
import tarfile

try:
    import wx
except ImportError:
    print('WARNING: GiniClust: Please install wxPython.')

from gooey import Gooey, GooeyParser
import pkg_resources


def extract_file(path2Rscript, target_directory):

    if path2Rscript.endswith('.tar.gz') or path2Rscript.endswith('.tgz'):
        opener, mode = tarfile.open, 'r:gz'
    else: 
        raise ValueError, "\nERROR: GiniClust: failed to extract {0}; no appropriate extractor could be found".format(path2Rscript)
    
    cwd = getcwd()
    chdir(target_directory)
    
    try:
        file = opener(path2Rscript, mode)
        try: 
            file.extractall()
        finally: 
            file.close()
    finally:
        chdir(cwd)


@Gooey(
    program_name='GiniClust')
def main():

    parser = GooeyParser(
        prog='GiniClust', 
        description="Detecting rare cell-types from single-cell "
            "gene expression data", 
        epilog="Contributors: Lan Jiang, "
            "Qian Zhu and Gregory Giecold.\nFor further help or information, "
            "please contact us at lan_jian@hms.harvard.edu,"
            "qzhu@princeton.edu or ggiecold@jimmy.harvard.edu")
        
    parser.add_argument('-t', '--type', nargs='?', type=str,
        choices=['qPCR', 'RNA-seq'], const='qPCR', default='RNA-seq',
        help='Type of your input genomics dataset')
    parser.add_argument('-o', '--output', nargs='?', type=str,
        default=path.join(getcwd(), 'GiniClust_results'),
        help="Specify GiniClust's output directory")
    parser.add_argument('input', type=str, widget='FileChooser',
        help='Name of the file to process')
    
    args = parser.parse_args()
    if not path.exists(args.output):
        makedirs(args.output)
    
    path2Rscript = pkg_resources.resource_filename(__name__, 
        'GiniClust_Main.R.tar.gz')
    extract_file(path2Rscript, args.output)
    
    newpath = path.join(args.output, 'Rfunction')
    if not path.exists(newpath):
        makedirs(newpath)
        
    Rscripts = ['DE_MAST.R.tar.gz', 'DE_MAST_figures.R.tar.gz',
        'DE_t_test.R.tar.gz', 'GiniClust_Clustering.R.tar.gz',
        'GiniClust_Filtering.R.tar.gz', 'GiniClust_Fitting.R.tar.gz',
        'GiniClust_packages.R.tar.gz', 'GiniClust_parameters.R.tar.gz',
        'GiniClust_Preprocess.R.tar.gz', 'GiniClust_tSNE.R.tar.gz']
    for Rscript in Rscripts:
        path2Rscript = pkg_resouces.resource_filename(__name__,
            path.join('Rscripts', Rscript))
        extract_file(path2Rscript, newpath)
    
    cmd = ['Rscript', path.join(args.output, 'GiniClust_Main.R')]
    cmd += ['-f', args.input, '-t', args.type, '-o', args.output]
    subprocess.check_output(cmd, universal_newlines=True)


if __name__ == '__main__':

    main()
    
     
