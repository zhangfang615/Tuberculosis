__author__ = 'Fang'

import sys
try:
    import pysam
except Exception as e:
    print >> sys.stderr, str(e)
    sys.exit(-1)
