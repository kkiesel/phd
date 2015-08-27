#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import argparse
import re
import pickle

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    args = parser.parse_args()


    xSecs = {}
    with open( args.filename ) as f:
        for line in f.readlines():
            match = re.match("\s*(\d+)\s*GeV(.*)±(.*)%.*", line )
            if match:
                m_str, xSec_str, uncert_str = match.groups()
                xSecs[ int(m_str) ] = ( float(xSec_str), float( uncert_str ) )


    output = open( args.filename.replace("txt", "pkl" ), 'wb')
    pickle.dump(xSecs, output)
    output.close()




