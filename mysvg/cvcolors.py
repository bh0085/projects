#!/usr/bin/env python
import getopt, sys
import colors, utils
from pysvg import parser as svgparser

def usage():
    print '''
Usage svg_cvcolors [-h -f -l -c (color) -t (table)]

Convert colors in an svg file using either a single color or a translation table.

OPTIONS:
  -h:   print this help
  -f #:  convert fills (default 1)
  -l #:  convert lines (default 0)
  -c  :  color name or hex.
  -a  :  color all elements


Usage:

Call cvcolors in order to fill uniformly and/or swap a table of colors in an svg file.

The -a flag controls whether elements having no fill/line specified will be colored by the fallback color. Since Adobe illustrator does not assing fill colors to black objects, you may need to set the -a flag in order to swap black out of filled shapes in illustrator. You could of course color the shape white and be fine...


'''
    return 0
  


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hc:t:af:l:", ["help", "color=","table="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    output = None
    verbose = False

    table_path = None
    color = None
    dofill = 1
    dolines = 0
    colorall = False
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-f"):
            dofill = a        
        elif o in ("-l"):
            dolines = a
        elif o in ("-c", "--color"):
            color = a        
        elif o in ("-t", "--table"):
            table_path = a
        elif o in ("-a"):
            colorall = True
        else:
            assert False, "unhandled option"
    if (table_path == None and color == None):
        usage()
        sys.exit()
    
    if table_path == None:
        table = {}
    else:
        fopen = open(table_path)
        table = dict([[e.strip() for e in re.split(re.compile('\s+'),l)]
                       for l in fopen.readlines()]) 
        
    

    ct = colors.color_names
    if ct.has_key(color):
        color = ct[color]
    for k,v in table.iteritems():
        if ct.has_key(k):
            table[ct[k]] = v;
    for k,v in table.iteritems():
        if ct.has_key(v):
            table[k] = ct[v];

            
    
    svg_in = svgparser.parse(sys.stdin)
    
    elts = utils.list_elements(svg_in)
    if dofill:
        for e in elts:
            if not 'get_fill' in dir(e): continue;
            if not colorall and not e.get_fill():
                continue
            
            if table.has_key(e.get_fill()):
                e.set_fill(table[e.get_fill()])
            else:
                e.set_fill(color)
    if dolines:
        for e in elts:
            if not 'get_color' in dir(e): continue;
            if (not colorall ) and ( not e.get_color()):
                continue
            if table.has_key(e.get_color()):
                e.set_color(table[e.get_color()])
            else:
                e.set_color(color)
                                

    encoding ='ISO-8859-1'
    standalone='no'
    
    sys.stdout.write(svg_in.wrap_xml(svg_in.getXML(), encoding, standalone))
    sys.stdout.close()


if __name__ == "__main__":
    main()
