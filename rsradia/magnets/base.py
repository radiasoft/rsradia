
import radia as rad

# Define minimum & maximum trangle mesh sizes
TRI_MIN_SIZE = 0.5
TRI_MAX_SIZE = 1.0

# Define a basic iron material
IRONMAT = rad.MatSatIsoFrm([20000, 2], [0.1, 2], [0.1, 2])

RADMATERS = {
    "NdFeB": "Neodymium Iron Boron (linear anisotropic)",
    "SmCo5": "Samarium Cobalt (SmCo5, linear anisotropic)",
    "Sm2Co17": "Samarium Cobalt (Sm2Co17,linear anisotropic)",
    "Ferrite": "Ferrite (linear anisotropic)",
    "Xc06": "Industeel (nonlinnear isotropic)",
    "Steel37": "Steel37 (nonlinnear isotropic)",
    "Steel42": "Steel42 (nonlinnear isotropic)",
    "AFK502": "Aperam 502 (nonlinnear isotropic)",
    "AFK1": "Aperam 1 (nonlinnear isotropic)"
}

# Returns a recursive list of constituent object IDs for Radia object 'radia_id'
def ObjCntStufRec(radia_id):
    sub_ids = []
    for obj in rad.ObjCntStuf(radia_id):
        if len(rad.ObjCntStuf(obj)) > 0:
            sub_ids.extend(ObjCntStufRec(obj))
        else:
            sub_ids.append(obj)
    return sub_ids

# Deletes Radia object 'radia_id' & all constituent objects
def UtiDelStuf(radia_id):
    for sub_id in ObjCntStufRec(radia_id):
        rad.UtiDel(sub_id)
    rad.UtiDel(radia_id)
