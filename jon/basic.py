import re
import itertools
from numpy import *
import matplotlib.pyplot as plt


def helloworld():

    print 'hi'

def load(country = '"United States"',
         n_subs = 10,
         do_seismic = False,
         hour_key = 'rpc',
         top_key = 'conversions'):
    '''
takes country, the country to filter rows by. country = 'all' gives no filtering.

takes n_subs.
'''
    report_file = open('/Users/bh0085/Programming/projects/jon/report.csv')
    stat_file = open('/Users/bh0085/Programming/projects/jon/custom_stat.csv')

    #build a list of dictionaries keyed by column name from
    #stat_file.
    l0 = stat_file.readline()
    cols = l0.strip().split(',')
    rows = []
    while 1:
        l = stat_file.readline()
        if l.strip() == '': break

        cre = re.compile('("[^",]*)(,*)([^",]*")')
        l0 = l
        l = re.sub(cre, '\g<1>\g<3>',l)

        rows.append(dict([(cols[i], v) 
                          for i, v in enumerate(l.strip().split(','))]))

    #fix little bugs
    last = rows.pop()
    for r in rows: r['date'] = r['\xef\xbb\xbfdate']
    for r in rows: r['sep_date'] = float(r['date'].split('-')[-1])

    #gathering stats here
    hours = array([e['hour'] for e in rows],float)
    sep_dates = array([e['sep_date'] for e in rows])
    rpcs = [float(e['rpc']) for e in rows]

    #group into (day,subid) bins
    fields = ['rpc', 'conversions', 'clicks', 'payout']
    totals = dict([(k, {}) for k in fields]) 
    tallies = {}
    
    for r in rows:
        if (r['country'] != country) and (country != 'all'): continue
        for f in fields:
            key = (r['affiliate_sub_id1'],
                   r['hour'])
            dict_of_interest = totals[f]
        
            dict_of_interest[key] = dict_of_interest.get(key,0) + float(r[f])
            
        tallies[key] = tallies.get(key,0) + 1

    means = {}
    for f, d in totals.iteritems():
        means[f] = {}
        for k, v in d.iteritems():
            m = v /tallies[k]
            means[f][k] = m

    
    #generate a LIST of mean clickvals, rpc vals
    #using a single set of keys to ensure proper alignment.

    #then zip the product to create input suitable for plotting.
    #[k,2] list ---> [2,k] list
    #exercise for the reader: do the same using numpy.array.T
    plot_data = zip(*[
            [ totals['clicks'][k], totals['rpc'][k] ] 
            for k in means['rpc'].keys()
            ])


    #clicks vs. rpc.
    f = plt.figure(1)
    f.clear()
    ax = f.add_subplot(111)

    ax.scatter(*plot_data)
    ax.annotate('clicks vs. rpc',
                [0,1],
                va = 'top',
                textcoords = 'axes fraction')
    ax.set_xlabel('clicks')
    ax.set_ylabel('rpc')

    subids = list(set([e[0]  for e in totals['clicks'].keys()]))

    #compute hourly breakdown
    hourlies =dict( [(s,dict([(k,zeros(24)) 
                              for k in totals.keys()])) 
                     for s in subids])
    for f,d in totals.iteritems():
        for k,v in d.iteritems():
            hour = float(k[1])
            subid = k[0]
            hourlies[subid][f][hour] += v

    #hourly dict is keyed by subid
    #and has dictionaries for values.
    #each subdictionary catalogs a specific variable over 24 hours.

    


    #I don't have any idea what this stuff does.
    #pick the highest quality subids based on coversions.
    average_conv = [(k,mean(v[top_key])) for k,v in hourlies.iteritems()]
    top = sorted(average_conv, key = lambda x: x[1])[::-1][:n_subs]
    top_keys = [e[0] for e in top]

    tuples = [(k,v[hour_key]) 
              for k,v in hourlies.iteritems()
              if k in top_keys]

    yvals = [t[1] for t in tuples]
    colors = random.rand(len(top),3)

    
    if do_seismic:
     try: 
         import cb.utils.seismic as seismic
         print 'no problem'
         do_seismic = True
     except Exception, e:
         print 'problem (seismic is unavailable)?'
        

    if do_seismic:
       f2 = plt.figure(2)
       ax2 = f2.add_subplot(111)
       seismic.seismic(yvals, colors = colors,
                       labels = top_keys,
                       stacked = False,
                       ax = ax2)

    else:
        f2 = plt.figure(2)
        ax2 = f2.add_subplot(111)
        for i, ys in enumerate(yvals):
            ax2.plot(range(24), ys, color = colors[i],
                     label = top_keys[i])
        ax2.plot(mean(yvals,0),linewidth = 5, color = 'black',
                 label = 'mean for top subids')
        ax2.set_xlabel('hour')
        ax2.set_ylabel(hour_key)
        print '{0} is the best hour (military time)'\
            .format(argmax(mean(yvals,0)))
    ax2.legend()

    

    ax2.annotate('random bullshit.',
                 [0,1],
                 va = 'top',
                 textcoords = 'axes fraction')
    return hourlies, top_keys
