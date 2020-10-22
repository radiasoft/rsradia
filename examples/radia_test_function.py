import radia as rad
from math import *


def HybridUndCenPart(_gap, _gap_ofst, _nper, _air, _lp, _ch_p, _np, _np_tip, _mp, _cp, _lm, _ch_m_xz, _ch_m_yz,
                     _ch_m_yz_r, _nm, _mm, _cm, _use_ex_sym=False):
    zer = [0, 0, 0]
    grp = rad.ObjCnt([])

    y = _lp[1] / 4
    initM = [0, -1, 0]

    pole = rad.ObjFullMag([_lp[0] / 4, y, -_lp[2] / 2 - _gap / 2 - _ch_p], [_lp[0] / 2, _lp[1] / 2, _lp[2]], zer,
                          [_np[0], int(_np[1] / 2 + 0.5), _np[2]], grp, _mp, _cp)

    if (_ch_p > 0.):  # Pole Tip
        poleTip = rad.ObjThckPgn(_lp[0] / 4, _lp[0] / 2,
                                 [[y - _lp[1] / 4, -_gap / 2 - _ch_p],
                                  [y - _lp[1] / 4, -_gap / 2],
                                  [y + _lp[1] / 4 - _ch_p, -_gap / 2],
                                  [y + _lp[1] / 4, -_gap / 2 - _ch_p]], zer)
        rad.ObjDivMag(poleTip, [_np_tip[0], int(_np_tip[1] / 2 + 0.5), _np_tip[2]]);
        rad.MatApl(poleTip, _mp);
        rad.ObjDrwAtr(poleTip, _cp);
        rad.ObjAddToCnt(grp, [poleTip])

    y += _lp[1] / 4 + _air + _lm[1] / 2

    for i in range(_nper):
        magnet = rad.ObjThckPgn(_lm[0] / 4, _lm[0] / 2,
                                [[y + _lm[1] / 2 - _ch_m_yz_r * _ch_m_yz, -_gap / 2 - _gap_ofst],
                                 [y + _lm[1] / 2, -_gap / 2 - _gap_ofst - _ch_m_yz],
                                 [y + _lm[1] / 2, -_gap / 2 - _gap_ofst - _lm[2] + _ch_m_yz],
                                 [y + _lm[1] / 2 - _ch_m_yz_r * _ch_m_yz, -_gap / 2 - _gap_ofst - _lm[2]],
                                 [y - _lm[1] / 2 + _ch_m_yz_r * _ch_m_yz, -_gap / 2 - _gap_ofst - _lm[2]],
                                 [y - _lm[1] / 2, -_gap / 2 - _gap_ofst - _lm[2] + _ch_m_yz],
                                 [y - _lm[1] / 2, -_gap / 2 - _gap_ofst - _ch_m_yz],
                                 [y - _lm[1] / 2 + _ch_m_yz_r * _ch_m_yz, -_gap / 2 - _gap_ofst]], initM)
        # Cuting Magnet Corners
        magnet = rad.ObjCutMag(magnet, [_lm[0] / 2 - _ch_m_xz, 0, -_gap / 2 - _gap_ofst], [1, 0, 1])[0]
        magnet = rad.ObjCutMag(magnet, [_lm[0] / 2 - _ch_m_xz, 0, -_gap / 2 - _gap_ofst - _lm[2]], [1, 0, -1])[0]

        rad.ObjDivMag(magnet, _nm);
        rad.MatApl(magnet, _mm);
        rad.ObjDrwAtr(magnet, _cm);
        rad.ObjAddToCnt(grp, [magnet])

        initM[1] *= -1
        y += _lm[1] / 2 + _lp[1] / 2 + _air;

        if (i < _nper - 1):
            pole = rad.ObjFullMag([_lp[0] / 4, y, -_lp[2] / 2 - _gap / 2 - _ch_p], [_lp[0] / 2, _lp[1], _lp[2]], zer,
                                  _np, grp, _mp, _cp)

            if (_ch_p > 0.):  # Pole Tip
                poleTip = rad.ObjThckPgn(_lp[0] / 4, _lp[0] / 2,
                                         [[y - _lp[1] / 2, -_gap / 2 - _ch_p],
                                          [y - _lp[1] / 2 + _ch_p, -_gap / 2],
                                          [y + _lp[1] / 2 - _ch_p, -_gap / 2],
                                          [y + _lp[1] / 2, -_gap / 2 - _ch_p]], zer)
                rad.ObjDivMag(poleTip, _np_tip);
                rad.MatApl(poleTip, _mp);
                rad.ObjDrwAtr(poleTip, _cp);
                rad.ObjAddToCnt(grp, [poleTip])

            y += _lm[1] / 2 + _lp[1] / 2 + _air;

    y -= _lp[1] / 4
    pole = rad.ObjFullMag([_lp[0] / 4, y, -_lp[2] / 2 - _gap / 2 - _ch_p], [_lp[0] / 2, _lp[1] / 2, _lp[2]], zer,
                          [_np[0], int(_np[1] / 2 + 0.5), _np[2]], grp, _mp, _cp)
    if (_ch_p > 0.):  # Pole Tip
        poleTip = rad.ObjThckPgn(_lp[0] / 4, _lp[0] / 2,
                                 [[y - _lp[1] / 4, -_gap / 2 - _ch_p],
                                  [y - _lp[1] / 4 + _ch_p, -_gap / 2],
                                  [y + _lp[1] / 4, -_gap / 2],
                                  [y + _lp[1] / 4, -_gap / 2 - _ch_p]], zer)
        rad.ObjDivMag(poleTip, [_np_tip[0], int(_np_tip[1] / 2 + 0.5), _np_tip[2]]);
        rad.MatApl(poleTip, _mp);
        rad.ObjDrwAtr(poleTip, _cp);
        rad.ObjAddToCnt(grp, [poleTip])

    # Symmetries
    if (_use_ex_sym):  # Some "non-physical" mirroring (applicable for calculation of central field only)
        y += _lp[1] / 4
        rad.TrfZerPerp(grp, [0, y, 0], [0, 1, 0])  # Mirror left-right
        rad.TrfZerPerp(grp, [0, 2 * y, 0], [0, 1, 0])

    #     #"Physical" symmetries (applicable also for calculation of total structure with terminations)
    #     rad.TrfZerPerp(grp, zer, [0,1,0]) #Mirror left-right
    #     #Mirror front-back
    #     rad.TrfZerPerp(grp, zer, [1,0,0])
    #     #Mirror top-bottom
    #     rad.TrfZerPara(grp, zer, [0,0,1])

    return grp