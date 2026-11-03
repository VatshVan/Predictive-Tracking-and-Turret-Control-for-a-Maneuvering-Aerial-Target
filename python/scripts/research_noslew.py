"""No-slew check: the turret starts on the target's initial elevation, so the alternative reference laws
are compared without the acquisition slew.  python scripts/research_noslew.py . 300"""
import json, os, sys
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.join(sys.argv[1], 'scripts')); sys.path.insert(0, sys.argv[1])
import monte_carlo as mc
DEFAULT = dict(mc.tc.CFG)

MODES = {'final': {}, 'exact3': {'ref_mode': 'exact_taylor', 'ref_n': 3}, 'contact1': {'ref_mode': 'contact'},
         'mm': {'ref_mode': 'moment_matched'}}
N = int(sys.argv[2])

def job(args):
    mode, s = args
    mc.tc.CFG.clear(); mc.tc.CFG.update(DEFAULT)
    th0 = np.deg2rad(mc.path_features(f'MC{s}')['th0_deg'])
    m = mc._one((s, 100.0, MODES[mode], th0, 'MC', False))
    return mode, s, m['success'], m['n_hits'], m['n_pulses'], m['point_rms_m'], m['t_first_hit']

if __name__ == '__main__':
    jobs = [(md, s) for md in MODES for s in range(N)]
    with Pool(7) as p:
        rows = p.map(job, jobs, chunksize=4)
    out = {}
    for md in MODES:
        r = [x for x in rows if x[0] == md]
        out[md] = dict(n=len(r), ok=sum(x[2] for x in r), hits=sum(x[3] for x in r),
                       pulses=sum(x[4] for x in r), misses=sum(x[4] - x[3] for x in r),
                       rms_mean=float(np.mean([x[5] for x in r])), rms_median=float(np.median([x[5] for x in r])),
                       tfh_mean=float(np.nanmean([x[6] for x in r])))
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(sys.argv[1], 'results', 'research', 'noslew.json'), 'w'), indent=1)
