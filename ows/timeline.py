import cb.utils.memo as mem
import cb.utils.plots as myplots
import cb.utils.colors as mycolors
import numpy as np
from numpy import *
import cb.utils.seismic as seismic

hashtags = ['occupywallstreet',
            'ows',
            'occupyboston',
            'owspickuplines']

people = ['occupyinfo']

import tweepy
def last_5(**kwargs):
    dnums= range(1,40)
    delts = [(int(9+floor(d/30)) , 1+int(d % 30)) for d in dnums]  
    days = ['2011-{0}-{1}'.format(*delt) for delt in delts]
    def set_l5(**kwargs):
        days = kwargs.get('days')
        all_results = {}
        for h in hashtags:
            all_results[h] = []
            search = ['#{0} since:{1}'.format(h,d) for d in days]
            for s in search:
                all_results[h].append(tweepy.api.search(s, rpp = 100))
            
        
        return all_results
    return mem.getOrSet(set_l5,**mem.rc(kwargs,
                                        days = days,
                                        name = ','.join(days)[:20]))

def plot_ronn(**kwargs):
    l5 = last_5(**mem.sr(kwargs))
    ys = zeros((len(l5),len(l5.values()[0])))
    for i,h in enumerate(l5.values()):
        for j,d in enumerate(h):
            if len(d) == 0 : continue
            
            cas = [e.created_at for e in d]
            secs = np.sum([(((ca.hour) * 60\
                                 +ca.minute) *60) \
                               + ca.second
                           for ca in cas])
            ys[i,j] = secs

    ys = ys[:,-20:]
    ys = ys - np.min(ys,1)[:,newaxis]
    colors = mycolors.getct(len(ys))
    seismic.seismic(ys, stacked = True, colors = colors)
    

def plot_times(**kwargs):
    vals = last_5(**mem.rc(kwargs))
    
    xs = []
    ys = []
    for v in vals:
        for r in v:

            raise Exception()
            ca = r.created_at
            catime = (((((ca.month * 30) + ca.day) * 24 + ca.hour) * 60\
                          +ca.minute) *60) + ca.second
            ys.append(catime)
            
    xs = range(len(ys))
        
    f = myplots.fignum(3, (4,4))
    ax = f.add_subplot(111)
    ax.scatter(xs,ys)

    return ys

def city_posts(**kwargs):
    def set_cp(**kwargs):
        wc = world_cities()
        lats = [float(w['lat0']) + float(w['lat1'])/60 for w in wc]
        lons = [float(w['lon0']) + float(w['lon1'])/60 for w in wc]
        tags = [','.join([str(lats[i]),str(lons[i]),'500mi']) 
                for i in range(len(lats))]
        res = [tweepy.api.search( '#ows', geocode = t, rpp = 100) for t in tags]
        return {'lats':lats, 'lons':lons, 'posts':res}
    v = mem.getOrSet(set_cp, **mem.rc(kwargs))
    return v
def plot_city_posts(**kwargs):
    cp = city_posts(**mem.sr(kwargs))
    xs = cp['lons']
    ys = cp['lats']
    rs = [len(x) for x in cp['posts']]
        
    f = myplots.fignum(3, (4,4))
    ax = f.add_subplot(111)
    ax.scatter(xs,ys,s=rs)

def world_cities():
    world_cities_text = '''Albany, N.Y.	42	40	73	45	12:00 noon
Albuquerque, N.M.	35	05	106	39	10:00 a.m.
Amarillo, Tex.	35	11	101	50	11:00 a.m.
Anchorage, Alaska	61	13	149	54	8:00 a.m.
Atlanta, Ga.	33	45	84	23	12:00 noon
Austin, Tex.	30	16	97	44	11:00 a.m.
Baker, Ore.	44	47	117	50	9:00 a.m.
Baltimore, Md.	39	18	76	38	12:00 noon
Bangor, Maine	44	48	68	47	12:00 noon
Birmingham, Ala.	33	30	86	50	11:00 a.m.
Bismarck, N.D.	46	48	100	47	11:00 a.m.
Boise, Idaho	43	36	116	13	10:00 a.m.
Boston, Mass.	42	21	71	5	12:00 noon
Buffalo, N.Y.	42	55	78	50	12:00 noon
Calgary, Alba., Can.	51	1	114	1	10:00 a.m.
Carlsbad, N.M.	32	26	104	15	10:00 a.m.
Charleston, S.C.	32	47	79	56	12:00 noon
Charleston, W. Va.	38	21	81	38	12:00 noon
Charlotte, N.C.	35	14	80	50	12:00 noon
Cheyenne, Wyo.	41	9	104	52	10:00 a.m.
Chicago, Ill.	41	50	87	37	11:00 a.m.
Cincinnati, Ohio	39	8	84	30	12:00 noon
Cleveland, Ohio	41	28	81	37	12:00 noon
Columbia, S.C.	34	0	81	2	12:00 noon
Columbus, Ohio	40	0	83	1	12:00 noon
Dallas, Tex.	32	46	96	46	11:00 a.m.
Denver, Colo.	39	45	105	0	10:00 a.m.
Des Moines, Iowa	41	35	93	37	11:00 a.m.
Detroit, Mich.	42	20	83	3	12:00 noon
Dubuque, Iowa	42	31	90	40	11:00 a.m.
Duluth, Minn.	46	49	92	5	11:00 a.m.
Eastport, Maine	44	54	67	0	12:00 noon
Edmonton, Alb., Can.	53	34	113	28	10:00 a.m.
El Centro, Calif.	32	38	115	33	9:00 a.m.
El Paso, Tex.	31	46	106	29	10:00 a.m.
Eugene, Ore.	44	3	123	5	9:00 a.m.
Fargo, N.D.	46	52	96	48	11:00 a.m.
Flagstaff, Ariz.	35	13	111	41	10:00 a.m.
Fort Worth, Tex.	32	43	97	19	11:00 a.m.
Fresno, Calif.	36	44	119	48	9:00 a.m.
Grand Junction, Colo.	39	5	108	33	10:00 a.m.
Grand Rapids, Mich.	42	58	85	40	12:00 noon
Havre, Mont.	48	33	109	43	10:00 a.m.
Helena, Mont.	46	35	112	2	10:00 a.m.
Honolulu, Hawaii	21	18	157	50	7:00 a.m.
Hot Springs, Ark.	34	31	93	3	11:00 a.m.
Houston, Tex.	29	45	95	21	11:00 a.m.
Idaho Falls, Idaho	43	30	112	1	10:00 a.m.
Indianapolis, Ind.	39	46	86	10	12:00 noon
Jackson, Miss.	32	20	90	12	11:00 a.m.
Jacksonville, Fla.	30	22	81	40	12:00 noon
Juneau, Alaska	58	18	134	24	8:00 a.m.
Kansas City, Mo.	39	6	94	35	11:00 a.m.
Key West, Fla.	24	33	81	48	12:00 noon
Kingston, Ont., Can.	44	15	76	30	12:00 noon
Klamath Falls, Ore.	42	10	121	44	9:00 a.m.
Knoxville, Tenn.	35	57	83	56	12:00 noon
Las Vegas, Nev.	36	10	115	12	9:00 a.m.
Lewiston, Idaho	46	24	117	2	9:00 a.m.
Lincoln, Neb.	40	50	96	40	11:00 a.m.
London, Ont., Can.	43	2	81	34	12:00 noon
Long Beach, Calif.	33	46	118	11	9:00 a.m.
Los Angeles, Calif.	34	3	118	15	9:00 a.m.
Louisville, Ky.	38	15	85	46	12:00 noon
Manchester, N.H.	43	0	71	30	12:00 noon
Memphis, Tenn.	35	9	90	3	11:00 a.m.
Miami, Fla.	25	46	80	12	12:00 noon
Milwaukee, Wis.	43	2	87	55	11:00 a.m.
Minneapolis, Minn.	44	59	93	14	11:00 a.m.
Mobile, Ala.	30	42	88	3	11:00 a.m.
Montgomery, Ala.	32	21	86	18	11:00 a.m.
Montpelier, Vt.	44	15	72	32	12:00 noon
Montreal, Que., Can.	45	30	73	35	12:00 noon
Moose Jaw, Sask., Can.	50	37	105	31	11:00 a.m.
Nashville, Tenn.	36	10	86	47	11:00 a.m.
Nelson, B.C., Can.	49	30	117	17	9:00 a.m.
Newark, N.J.	40	44	74	10	12:00 noon
New Haven, Conn.	41	19	72	55	12:00 noon
New Orleans, La.	29	57	90	4	11:00 a.m.
New York, N.Y.	40	47	73	58	12:00 noon
Nome, Alaska	64	25	165	30	8:00 a.m.
Oakland, Calif.	37	48	122	16	9:00 a.m.
Oklahoma City, Okla.	35	26	97	28	11:00 a.m.
Omaha, Neb.	41	15	95	56	11:00 a.m.
Ottawa, Ont., Can.	45	24	75	43	12:00 noon
Philadelphia, Pa.	39	57	75	10	12:00 noon
Phoenix, Ariz.	33	29	112	4	10:00 a.m.
Pierre, S.D.	44	22	100	21	11:00 a.m.
Pittsburgh, Pa.	40	27	79	57	12:00 noon
Portland, Maine	43	40	70	15	12:00 noon
Portland, Ore.	45	31	122	41	9:00 a.m.
Providence, R.I.	41	50	71	24	12:00 noon
Quebec, Que., Can.	46	49	71	11	12:00 noon
Raleigh, N.C.	35	46	78	39	12:00 noon
Reno, Nev.	39	30	119	49	9:00 a.m.
Richfield, Utah	38	46	112	5	10:00 a.m.
Richmond, Va.	37	33	77	29	12:00 noon
Roanoke, Va.	37	17	79	57	12:00 noon
Sacramento, Calif.	38	35	121	30	9:00 a.m.
St. John, N.B., Can.	45	18	66	10	1:00 p.m.
St. Louis, Mo.	38	35	90	12	11:00 a.m.
Salt Lake City, Utah	40	46	111	54	10:00 a.m.
San Antonio, Tex.	29	23	98	33	11:00 a.m.
San Diego, Calif.	32	42	117	10	9:00 a.m.
San Francisco, Calif.	37	47	122	26	9:00 a.m.
San Jose, Calif.	37	20	121	53	9:00 a.m.
San Juan, P.R.	18	30	66	10	1:00 p.m.
Santa Fe, N.M.	35	41	105	57	10:00 a.m.
Savannah, Ga.	32	5	81	5	12:00 noon
Seattle, Wash.	47	37	122	20	9:00 a.m.
Shreveport, La.	32	28	93	42	11:00 a.m.
Sioux Falls, S.D.	43	33	96	44	11:00 a.m.
Sitka, Alaska	57	10	135	15	8:00 a.m.
Spokane, Wash.	47	40	117	26	9:00 a.m.
Springfield, Ill.	39	48	89	38	11:00 a.m.
Springfield, Mass.	42	6	72	34	12:00 noon
Springfield, Mo.	37	13	93	17	11:00 a.m.
Syracuse, N.Y.	43	2	76	8	12:00 noon
Tampa, Fla.	27	57	82	27	12:00 noon
Toledo, Ohio	41	39	83	33	12:00 noon
Toronto, Ont., Can.	43	40	79	24	12:00 noon
Tulsa, Okla.	36	09	95	59	11:00 a.m.
Vancouver, B.C., Can.	49	13	123	06	9:00 a.m.
Victoria, B.C., Can.	48	25	123	21	9:00 a.m.
Virginia Beach, Va.	36	51	75	58	12:00 noon
Washington, D.C.	38	53	77	02	12:00 noon
Wichita, Kan.	37	43	97	17	11:00 a.m.
Wilmington, N.C.	34	14	77	57	12:00 noon
Winnipeg, Man., Can.	49	54	97	7	11:00 a.m.'''
    cols = ['name','lat0','lat1','lon0','lon1','time']
    outs = []
    for c in world_cities_text.splitlines():
        outs.append(dict([(i,e) for i,e in zip(cols,c.split('\t'))]))
    return outs
