#!/usr/bin/env python


# GiniClust/src/GiniClust/__main__.py 


# Author: Gregory Giecold
# Maintainer: Qian Zhu
# Affiliations: Harvard University & Princeton University
# Contact information: ggiecold@jimmy.harvard.edu; qzhu@princeton.edu


from __future__ import print_function

from os import getcwd, path
import subprocess

try:
    import wx
except ImportError:
    print('WARNING: GiniClust: Please install wxPython.')

from gooey import Gooey, GooeyParser


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
    
    command = 'Rscript'
    path2Rscript = path.join(getcwd(), 'GiniClust_Main.R')
    args = parser.parse_args()
    
    cmd = [command, path2Rscript]
    cmd += ['-f', args.input, '-t', args.type, '-o', args.output]
    subprocess.check_output(cmd, universal_newlines=True)


if __name__ == '__main__':

    main()
    
     
