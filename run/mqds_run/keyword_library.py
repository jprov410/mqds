"""
Module containing keyword library and functions for
processing the input file
"""
modulename='keyword_library'

"""
keyword library
"""
keywords = dict(method='pldm', calculation='redmat', bath='harmonic', nbstep=1000, basis='site', ntraj=10000, nlit=20,
                dump=1, nstate=2, initstate=1, initstatet=1, nslice=1, zpe='0.5', window='0.5', runtime='1000.0',
                temperature='77.0', tdelay1='500.0', tdelay2='0.0', tdelay3='500.0', nstep1=500, nstep2=0, nstep3=500,
                branch1=1, branch2=1, branch3=1)


def isfloat(input):
    try:
        float(input)
        return True
    except ValueError:
        return False

def default(input):
    return keywords[input]
