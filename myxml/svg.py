import re
import parser
import numpy as np
from numpy import *
def parse(fname):
    xml = parser.parse( open(fname).read())
    return xml

def get_elt_by_type(data, dtype):
    return  parser.svg_elt(data,dtype)

def get_polys(data, rescale = True):
    elts = get_elt_by_type(data,'polyline')
    all_points = []
    for e in elts:
        points = e.attrib['points']
        all_points.append( array([[float(val) for val in st.split(',')] 
                                  for st in re.split(re.compile('\s+'),points.strip())])
)
    if rescale:
        maxvals, minvals, deltas ,centers = [],[], [], []
        for d in range(shape(all_points[0])[1]):
            maxvals.append(   max([ np.max(p[:,d]) for p in all_points]))
            minvals.append(   min([ np.min(p[:,d]) for p in all_points]))
            deltas.append(  maxvals[-1] - minvals[-1])
            centers.append( (maxvals[-1] + minvals[-1]) / 2)
        dmax = max(deltas)
        for i,p in enumerate(all_points):
            p -= centers[i]
            p /= dmax


    return all_points

