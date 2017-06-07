"""
Print out statistics
"""

from collections import OrderedDict
import cksphys.io

def stats():
    d = OrderedDict()
    cache = 1 
    cand = cksphys.io.load_table('cks-cuts',cache=cache)
    stars = cand.groupby('id_starname',as_index=False).first()

    d['n-cand-cks'] = "{}".format( len(cand) )
    d['n-stars-cks'] = "{}".format( len(stars) )

    cand = cand[~cand.isany]
    stars = cand.groupby('id_starname',as_index=False).first()

    d['n-cand-cks-pass'] = "{}".format( len(cand))
    d['n-stars-cks-pass'] = "{}".format( len(stars))

    for k, v in d.iteritems():
        print r"{{{}}}{{{}}}".format(k,v)



