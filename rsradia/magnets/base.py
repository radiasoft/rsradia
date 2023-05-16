### base.py
### November 2022
### Provides functions & variables basic to rsradia.magnets

from radia import ObjCntStuf, UtiDel

# Define minimum & maximum trangle mesh sizes
TRI_MIN_SIZE = 0.5
TRI_MAX_SIZE = 1.0

# Returns a recursive list of constituent object IDs for Radia object 'radia_id'
def ObjCntStufRec(radia_id):
    sub_ids = []
    for obj in ObjCntStuf(radia_id):
        if len(ObjCntStuf(obj)) > 0:
            sub_ids.extend(ObjCntStufRec(obj))
        else:
            sub_ids.append(obj)
    return sub_ids

# Deletes Radia object 'radia_id' & all constituent objects
def UtiDelStuf(radia_id):
    for sub_id in ObjCntStufRec(radia_id):
        UtiDel(sub_id)
    UtiDel(radia_id)