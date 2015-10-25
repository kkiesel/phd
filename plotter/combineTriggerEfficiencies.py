#!/usr/bin/evn python2
# -*- coding: utf-8 -*-


from math import sqrt


def getInfoFromStr( str ):
    import re
    out = []
    for subStr in str.split(" "):
        subStr = subStr.strip('%')
        try: out.append(float(subStr))
        except: pass
    return out

    #m = re.match( "ε = (.*)% + (.*)% - (.*)%", str )
    #return m.groups()


str1 = "ε = 98.6% + 0.8% - 1.3%"
str2 = "ε = 95.4% + 0.4% - 0.5%"

e1, e1up, e1dn = getInfoFromStr(str1)
e2, e2up, e2dn = getInfoFromStr(str2)


e = e1*e2/100.
eup = sqrt( (e1up*e2)**2 + (e1*e2up)**2 )/100.
edn = sqrt( (e1dn*e2)**2 + (e1*e2dn)**2 )/100.

print "{:.1f}^{{+{:.1f}}}_{{-{:.1f}}}".format(e, eup,edn)
print e, eup, edn



