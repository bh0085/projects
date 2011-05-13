import compbio.utils.colors as mycolors
#from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, \
#     AnnotationBbox
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import itertools as it
import matplotlib.transforms as transforms
from matplotlib.transforms import Affine2D
from matplotlib.collections  import LineCollection
import utils as rutils
from matplotlib.patches import Ellipse

#import mpl_toolkits.axisartist.floating_axes as floating_axes

import numpy as np
#import  mpl_toolkits.axisartist.angle_helper as angle_helper
#from matplotlib.projections import PolarAxes
#from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator, \
#     DictFormatter


#from mpl_toolkits.axes_grid.anchored_artists import AnchoredAuxTransformBox
#from mpl_toolkits.axes_grid.anchored_artists import AnchoredDrawingArea


def fignum(num, size = (8,8)):
    plt.close(num)
    return plt.figure(num, size)

    

def grid_rnas(polys,dims = [20,20], colors = None, size = None ):
    '''grid plot a bunch of RNAs according to the color scheme given.

input:
  polys:   vertices for each RNA such as those pulled from the 
           svgs produced by RNAsubopt im compbio.utils.svg.

  colors:  colors for each rna. Can be None or a vector of length
           equal to n rnas or a vector of vectors (allowing coloring 
           of individual vertices)


'''
    f = plt.gcf()
    n = len(polys)

    xdim = floor(sqrt(len(polys)))
    ydim = ceil(len(polys)/xdim)

    ax = f.add_axes([0,0,1,1], 
                       xlim = [-1,xdim],
                       ylim = [-1,ydim])

    for i, p in enumerate(polys):
        show_rna([mod(i,xdim), 
                  floor(i/xdim) ], p/3, 
                 ax = ax, dims = dims,
                 pkw = dict(color =  colors[i] if colors != None else 'black') )
    return f   

def show_rna( emb, vertices , ax = None,dims = [20,20], pkw= {}, **kwargs):
    if ax == None: ax = plt.gca()
    #v0 = ax.plot(*vertices.T*100, lw = 12, color ='white',
    #              transform = transforms.ScaledTranslation(emb[0],emb[1],ax.transData )
    #              )[0]
    v1 = ax.plot(*vertices.T*dims[0],
                  transform = transforms.ScaledTranslation(emb[0],emb[1],ax.transData ),
                  **pkw)[0]

    #f = ax.figure

    #da = DrawingArea(dims[0],dims[1], 0, 0)
    #p =  plt.plot(vertices[:,0]*50, vertices[:,1]*50)[0]
    #da.add_artist(p)

    #xy = emb[0:2]
    
    #bb = ax.bbox.translated(20,20)
    #box2 = AnchoredAuxTransformBox(ax.transData, loc=2)


    #grid_helper = floating_axes.\
    #    GridHelperCurveLinear(Affine2D(),  extremes=(-1, 1, -1,1))
    #rect = 111
    #ax1 = floating_axes.FloatingAxes(f,[0.,0.,.4,.4],\
    #                                     grid_helper = grid_helper,
    #                                 transform = ax.transData)
    #ax1.set_frame_on(False)
    #
    #
    #
    #
    #el = ax1.plot(*vertices.T)[0]
    #ax.add_artist(ax1)
    #
    #
    #
    #box = AnchoredDrawingArea(100, 100, 0, 0,
    #                           loc=3, pad=0.,
    #                          frameon=True)
    #dx = 1
    #dy = 1
    #offset = transforms.ScaledTranslation(dx, dy,f.dpi_scale_trans)
    #shadow_transform = ax.transData + offset
    #verts2 = plt.plot(*vertices.T*array([1,1])[:,newaxis],
    #                   transform = ax.transData 
    #                   )[0]
    #verts1 = plt.plot(*vertices.T*array([100,100])[:,newaxis],
    #                   transform = ax.transData + offset)[0]
    #print verts1.get_transform()
    #box.drawing_area.add_artist(verts1)
    #
    #print
    #print
    #print verts1.get_transform()
    #                  
    #ax.add_artist(box)
    #
    #
    ##ax2.set_transform(transforms.IdentityTransform())
    #
    ##ax.plot(arange(10),0*arange(10))
    ##p.set_transform(t)
    #
    ##a
    #return (box, ax1)
    #
    #
    ##ab = AnnotationBbox(da, xy,
    ##                    xybox=(0,0),
    ##                    xycoords='data',
    ##                    boxcoords=("offset pixels"),
    ##                    box_alignment=(0.,0.5),
    ##                    arrowprops=dict(arrowstyle="->"))
    ##ax.add_artist(p)


        
    
def show_output(outputs, 
		show = 'conservation',
		save = True):
        mvecs = outputs['all_vecs']['all_time']
        tvecs = outputs['all_vecs']['all_mut']
        fvecs = outputs['all_vecs']['fiftyfifty']

	run_id = outputs['run_id']
	structs = outputs['exemplar_structs']
	ref = outputs['reference_seq']
	
	thermo_pairs = outputs['thermo_pairs']
	thermo_inds  = outputs['thermo_ex_inds']

	run_title = outputs['title']
	fam_name = re.compile('RF\d*').search(run_title).group()

	fig = plt.gcf()
	try: fig.clear()
	except Exception, e: print 'wonky 3d bug'
	fig = plt.gcf()
	try: fig.clear()
	except Exception, e: print 'wonky 3d bug'
	fig.canvas.draw()


	exemplar_inds = sorted(list(set(thermo_inds)))
	struct_colors = dict([(exemplar_inds[i], col) 
			       for i, col in enumerate(mycolors.getct(len(exemplar_inds)))]
                              )

	if show == 'embeddings':

	   
	   exemplars = list(set(thermo_inds))
	   pair_embedding =  compute_embedding(thermo_pairs,
	   			       aff_type = 'pairs',
	   			       do_mve = False,
					       ss_multiplier = None)
	   
	   shape_embedding = compute_embedding(thermo_pairs,
	   			       aff_type = 'easy',
	   			       do_mve = False,
	   			       ss_multiplier = None)
	   show_3d = True
	   #shape_embedding[0] is pca
	   rplots.plot_clusters( thermo_inds, {'shape':shape_embedding[0],
	   			     'pairs':pair_embedding[0]}, 
			  plot3d = show_3d,
			  title = 'projection ({0}) '.format(run_id),
			  save = save,
			  colors = struct_colors)

	elif show == 'conservation':
		ax0 = fig.add_subplot('311')
		lstructs =  [project_lstruct(p, len(ref)) for p in structs]
		seismic.seismic([ abs(l) for l in lstructs] , 
				colors = struct_colors.values(),
				ax = ax0)

		myplots.maketitle(ax0, 'Predicted conservation patterns for {0}'.format(fam_name))

		shapes = array([shape(m) for m in mvecs])
		igood = nonzero(greater(shapes[:,1],0))[0]
		clade_colors = mycolors.getct(len(igood))
		mvg, tvg, fvg = [ [vecs[i] for i in igood] for vecs in [mvecs,tvecs,fvecs]]
		cons_types = array([ mvg, tvg, tvg])
		
		for c in cons_types:
			nrm = sum(c.flatten())
			if nrm == 0: nrm = 1
			c /= sum(c.flatten())
		if shape(cons_types)[1] == 0:
			print 'No good vectors!'
			return
		
		mtype_sums = np.sum(np.sum(cons_types,3),0)	
		stype_sums = np.sum(np.sum(cons_types,3),0).T


		ax1 = fig.add_subplot('312')		
		seismic.seismic(stype_sums , 
				colors = struct_colors.values(),
				ax = ax1)

		myplots.maketitle(ax1,'Observed conservation (struct v. clade) patterns for {0}'\
					  .format(fam_name),
				  )

		
		ax2 = fig.add_subplot('313')
		
		seismic.seismic(mtype_sums , 
				ax = ax2, colors = clade_colors, stacked = True,
				label_y = False)

		#myplots.maketitle(ax2, 'Observed conservation (clade v. struct) patterns for {0}'\
		#			   .format(run_title)
		#		   )
		ax2.annotate('Observed conservation (clade v. struct) patterns for {0}'\
				     .format(run_title),
			     [.5,0],xycoords = 'axes fraction', ha = 'center', va = 'top',
			     size = 'x-large')

		if save: fig.savefig(cfg.dataPath('cs874/figs/cons_profiles/{0}.ps'.format(run_title)))
	       				 	

	

	else: raise Exception('show type not implemented: {0}'.format(show))
	

def plot_clusters(inds,
                  embeddings,
                  plot3d = False,
                  title = '',
		  ax_in =None,
		  save = False,
		  colors = None):
        exemplars = list(set(inds))
        if colors == None:
		cluster_colors = dict([(exemplars[i], col) 
                              for i, col in enumerate(mycolors.getct(len(exemplars)))]
                              )
	else:
		cluster_colors = colors

        cols = [cluster_colors[e] for e in inds]
        try: 
		if ax == None: plt.clf()
        except Exception, e: pass
        if ax_in == None: f = plt.gcf()

        for i, k in enumerate(embeddings.keys()):
            embedding = embeddings[k]

	    #if i == 1: raise Exception()
            emb_sig = embedding[:,0:3]
            cluster_vars = array([ var(emb_sig[nonzero(equal(inds, j))[0]])  for j in exemplars])
            indexed_vars = array([ cluster_vars[exemplars.index(j)] for j in inds ])
	    indexed_vars[equal(indexed_vars,0)] = 1

            sizes = 10 *( exp( -1 * ( np.sum((emb_sig - emb_sig[inds,:])**2,1)/indexed_vars)))
            if plot3d:
                if ax_in == None: 
			ax = f.add_subplot('{1}1{0}'.format(i+1, len(embeddings)),projection = '3d')
		else: ax = ax_in
                ax.scatter(array(embedding[:,0],float)
                           ,array(embedding[:,1],float)
                           ,array(embedding[:,2],float), 
                           s = sizes,
                           color = cols)
                ax.set_xticks([])
                ax.set_yticks([])
                for tl in list(it.chain( ax.w_xaxis.get_ticklabels(),
                                    ax.w_yaxis.get_ticklabels(),
                                    ax.w_zaxis.get_ticklabels())): # re-create what autofmt_xdate but with w_xaxis
                    tl.set_visible(False)
                    tl.set_rotation(30)    
            else:
                if ax_in == None: ax = f.add_subplot('{1}1{0}'.format(i+1, len(embeddings)))
		else: ax = ax_in
                ax.scatter(array(embedding[:,0],float)
                           ,array(embedding[:,1],float),
                           s = sizes,
                           color = cols)
                print 'sttring'
            ax.set_title('{0} for subopts in {1}'.format(k, title))
        
        if save: 
		f.savefig(cfg.dataPath('cs874/figs/subopt_embeddings/{0}.ps').format(title))



    


