#!/usr/bin/env python

from pysvg import parser as svgparser
import sys, getopt, argparse, md5, os
from mysvg import utils as svgu

def convert_dir( fill = None, dirname = '.', stroke = None, strokeWidth = None, dest = None):


    files = [ f for f in os.listdir(dirname) if f[-3:] == 'svg']
    print "{0} SVG Files: ".format(len(files))
    print '\n'.join(files)
                                 
    print "Converting fill to: {0}{1}{2}"\
        .format(fill,
                ", stroke to {0}".format(stroke) if stroke else '',
                ", strokeWidth to {0}".format(strokeWidth) if strokeWidth else '') 


    if dest == None:
        dest = ("{0}-fill{1}{2}"\
                .format(fill,
                        "-{0}-stroke".format(stroke) if stroke else '',
                        "-{0}-width".format(strokeWidth) if strokeWidth else '') )
    
    if not os.path.isdir(dest):
        os.mkdir(dest)

    for f in files:
        tree = svgparser.parse(os.path.join(dirname,f))
        elts = svgu.list_elements(tree)
        for e in elts:
            
            if 'set_fill' in dir(e) and fill:
                e.set_fill(fill)
            if 'set_stroke' in dir(e) and stroke:
                e.set_stroke(stroke)
            if 'set_stroke_width' in dir(e) and strokeWidth:
                e.set_stroke_width("{0}px".format(strokeWidth))
        out = os.path.join(dest, f)
        tree.save(out)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('fill',
                        nargs = 1 , type = str,
                        help = "Select a fill.")

    parser.add_argument('stroke',
                        nargs = '?' , type = str, default = None,
                        help = "Select a stroke.")


    parser.add_argument('strokeWidth',
                        nargs = '?' , type = int, default = None,
                        help = "Select a stroke width.")

    args = parser.parse_args()
    args.fill = args.fill[0]

    convert_dir(args.fill, stroke = args.stroke, strokeWidth = args.strokeWidth,
                dirname = '.')
    exit(0)


        
    
        
