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
    program_name='GiniClust',
    advanced=True)
def main():

    parser = GooeyParser(
        prog='', 
        description="Detecting rare cell-types from single-cell "
            "gene expression data", 
        epilog="Contributors: Lan Jiang, "
            "Qian Zhu and Gregory Giecold.\nFor further help or information, "
            "please contact us at lan_jian@hms.harvard.edu,"
            "qzhu@princeton.edu or ggiecold@jimmy.harvard.edu")
                    
    subparsers = parser.add_subparsers(dest='datatype', 
        help="Type of your input genomics dataset")
    
    qPCR_parser = subparsers.add_parser('pcr')
    qPCR_parser.add_argument('Input', type=str, widget='FileChooser',
        help='Select a file to process:')
    qPCR_parser.add_argument('-e', '--epsilon', nargs='?',
        type=float, const=0.25, default=0.25, 
        help='DBSCAN epsilon parameter:')
    qPCR_parser.add_argument('-m', '--minPts', nargs='?',
        type=int, const=5, default=5,
        help='DBSCAN minPts parameter:')
    qPCR_parser.add_argument('-O', '--Output', nargs='?', type=str,
        default=path.join(getcwd(), 'GiniClust_results'),
        help="Specify GiniClust's output directory:")
    
    RNASeq_parser = subparsers.add_parser('rna')
    RNASeq_parser.add_argument('Input', type=str, widget='FileChooser',
        help='Select a file to process:')
    RNASeq_parser.add_argument('-e', '--epsilon', nargs='?',
        type=float, const=0.5, default=0.5, 
        help='DBSCAN epsilon parameter:')
    RNASeq_parser.add_argument('-m', '--minPts', nargs='?',
        type=int, const=3, default=3,
        help='DBSCAN minPts parameter:')
    RNASeq_parser.add_argument('-O', '--Output', nargs='?', type=str,
        default=path.join(getcwd(), 'GiniClust_results'),
        help="Specify GiniClust's output directory:")
    
    command = 'Rscript'
    path2Rscript = path.join(getcwd(), 'GiniClust_Main.R')

    args = parser.parse_args()
    
    if args.datatype == 'pcr':
        datatype_str = 'qPCR'
    else:
        datatype_str = 'RNA-seq'
        
    cmd = [command, path2Rscript]
    cmd += ['-f', args.Input, '-t', datatype_str, '-o', args.Output]
    subprocess.check_output(cmd, universal_newlines=True)


if __name__ == '__main__':

    main()
    
     
