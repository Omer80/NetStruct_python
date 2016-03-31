# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 10:49:02 2016

@author: Omer Tzuk <omertz@post.bgu.ac.il>
"""

import csv
import numpy as np

def buildGDmatrix(csvfile):
    for row in csvfile:
        print row
    
    
    
def main(args):
    with open(args.filename, 'rb') as csvfile:
        spamreader = csv.reader(csvfile)
        buildGDmatrix(spamreader)
    
def add_parser_arguments(parser):
    parser.add_argument('filename', type=str, nargs='?', default='wildass4pops.csv',
                        help='output file')
    return parser


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser = add_parser_arguments(parser)
    args = parser.parse_args()
    main(args)