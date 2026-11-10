#!/usr/bin/env python3
"""Numerically verify CSE-generated force.ur.h against
the template source formulas. For each (LatSet, random inputs), compute
the reference result from the formula and compare with the generated code."""
import re, math, random, sys

LAT = {
 'D2Q5':  dict(q=5,d=2,cs2=1/3,c=[(0,0),(1,0),(-1,0),(0,1),(0,-1)]),
 'D2Q9':  dict(q=9,d=2,cs2=1/3,c=[(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)]),
 'D3Q7':  dict(q=7,d=3,cs2=1/3,c=[(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]),
 'D3Q15': dict(q=15,d=3,cs2=1/3,c=[(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,1,1),(-1,-1,-1),(1,1,-1),(-1,-1,1),(1,-1,1),(-1,1,-1),(-1,1,1),(1,-1,-1)]),
 'D3Q19': dict(q=19,d=3,cs2=1/3,c=[(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,1,0),(-1,-1,0),(1,0,1),(-1,0,-1),(0,1,1),(0,-1,-1),(1,-1,0),(-1,1,0),(1,0,-1),(-1,0,1),(0,1,-1),(0,-1,1)]),
 'D3Q27': dict(q=27,d=3,cs2=1/3,c=[(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,1,0),(-1,-1,0),(1,0,1),(-1,0,-1),(0,1,1),(0,-1,-1),(1,-1,0),(-1,1,0),(1,0,-1),(-1,0,1),(0,1,-1),(0,-1,1),(1,1,1),(-1,-1,-1),(1,1,-1),(-1,-1,1),(1,-1,1),(-1,1,-1),(-1,1,1),(1,-1,-1)]),
}

def w_k(latname, k):
    lat = LAT[latname]
    c2 = sum(x*x for x in lat['c'][k])
    if lat['d'] == 2:
        return {0: 4.0/9.0, 1: 1.0/9.0}.get(c2, 1.0/36.0)
    else:
        return {0: 2.0/9.0, 1: 1.0/9.0, 2: 1.0/72.0}.get(c2, 1.0/216.0)

def ref_force_pop(u, F, latname):
    """Reference: Fi[k] = w_k * F · ((c_k - u)*InvCs2 + (c_k·u)*InvCs4*c_k)
    Template: Fi[i] = w[i] * F * ((c[i] - u) * InvCs2 + (c[i] * u * InvCs4) * c[i])
    where * between vectors is dot product."""
    lat = LAT[latname]
    q, d, cs2 = lat['q'], lat['d'], lat['cs2']
    InvCs2, InvCs4 = 1.0/cs2, 1.0/(cs2**2)
    Fi = []
    for k in range(q):
        ck = lat['c'][k]
        cu = sum(ck[i]*u[i] for i in range(d))
        s = sum(F[i]*((ck[i] - u[i])*InvCs2 + InvCs4*cu*ck[i]) for i in range(d))
        Fi.append(w_k(latname, k) * s)
    return Fi

def ref_scalar_force_pop(u, F, d_idx, latname):
    """Reference: Fi[k] = w_k * F * ((c_k[d]-u[d])*InvCs2 + (c_k·u)*InvCs4*c_k[d])
    Template: v1 = (c[d]-u[d]) * InvCs2; v2 = (c·u * InvCs4) * c[d]; Fi[i] = w[i] * F * (v1 + v2)"""
    lat = LAT[latname]
    q, d, cs2 = lat['q'], lat['d'], lat['cs2']
    InvCs2, InvCs4 = 1.0/cs2, 1.0/(cs2**2)
    Fi = []
    for k in range(q):
        ck = lat['c'][k]
        cu = sum(ck[i]*u[i] for i in range(d))
        v1 = (ck[d_idx] - u[d_idx]) * InvCs2
        v2 = cu * InvCs4 * ck[d_idx]
        Fi.append(w_k(latname, k) * F * (v1 + v2))
    return Fi

def parse_ur(path, struct_name, latname):
    """Parse a .ur.h file and extract the apply() body for a given struct+LatSet."""
    src = open(path).read()
    if struct_name == 'ScalarForcePopImpl':
        pattern = rf'struct {struct_name}<T, {latname}<T>, (\d+)>\{{(.*?)\n\}};'
        matches = list(re.finditer(pattern, src, re.S))
        results = {}
        for m in matches:
            d_val = int(m.group(1))
            body = m.group(2)
            am = re.search(r'compute\([^)]*\)\{(.*?)\n  \}', body, re.S)
            if not am:
                am = re.search(r'compute\([^)]*\)\{\n(.*?)\n\}', body, re.S)
            if am:
                results[d_val] = am.group(1)
        return results
    else:
        pattern = rf'struct {struct_name}<T, {latname}<T>>\{{(.*?)\n\}};'
        m = re.search(pattern, src, re.S)
        if not m:
            return None
        body = m.group(1)
        am = re.search(r'compute\([^)]*\)\{(.*?)\n  \}', body, re.S)
        if not am:
            am = re.search(r'compute\([^)]*\)\{\n(.*?)\n\}', body, re.S)
        if not am:
            am = re.search(r'apply\([^)]*\)\{(.*?)\n  \}', body, re.S)
        if not am:
            am = re.search(r'apply\([^)]*\)\{\n(.*?)\n\}', body, re.S)
        return am.group(1) if am else None

def run_force_pop(body, u, F, latname, scalar_F=False):
    """Execute a ForcePopImpl::compute body."""
    lat = LAT[latname]
    q, d, cs2 = lat['q'], lat['d'], lat['cs2']
    InvCs2, InvCs4 = 1.0/cs2, 1.0/(cs2**2)
    w = [w_k(latname, k) for k in range(q)]
    cvecs = [list(lat['c'][k]) for k in range(q)]
    Fi = [0.0]*q
    env = {'u': list(u), 'w': w, 'cvec': cvecs,
           'InvCs2': InvCs2, 'InvCs4': InvCs4, 'Fi': Fi, 'q': q, 'd': d,
           'math': math, 'D2Q5': None, 'D2Q9': None, 'D3Q7': None,
           'D3Q15': None, 'D3Q19': None, 'D3Q27': None, 'LatSet': None,
           'T': None}
    if scalar_F:
        env['F'] = float(F)
    else:
        env['F'] = list(F)

    for s in [x.strip() for x in body.split(';') if x.strip()]:
        if s.startswith('using') or s.startswith('if') or not s: continue
        s = re.sub(r'^const\s+', '', s)
        # T var = expr
        dm = re.match(r'T\s+(\w+)\s*=\s*(.+)', s)
        if dm:
            env[dm.group(1)] = _eval(dm.group(2), env)
            continue
        # var = expr (after const T stripped)
        dm2 = re.match(r'(\w+)\s*=\s*(.+)', s)
        if dm2 and dm2.group(1) not in ('Fi',):
            env[dm2.group(1)] = _eval(dm2.group(2), env)
            continue
        # Fi[k] = expr
        fm = re.match(r'Fi\[(\d+)\]\s*=\s*(.+)', s)
        if fm:
            k = int(fm.group(1))
            Fi[k] = _eval(fm.group(2), env)
            continue
    return Fi

def _eval(e, env):
    """Simplified expression evaluator."""
    e = e.strip()
    e = e.replace('LatSet::InvCs2', str(env['InvCs2']))
    e = e.replace('LatSet::InvCs4', str(env['InvCs4']))
    for lat in ['D2Q5', 'D2Q9', 'D3Q7', 'D3Q15', 'D3Q19', 'D3Q27']:
        e = e.replace(f'{lat}<T>::InvCs2', str(env['InvCs2']))
        e = e.replace(f'{lat}<T>::InvCs4', str(env['InvCs4']))
    def rw(m):
        return str(env['w'][int(m.group(1))])
    e = re.sub(r'latset::w<[^>]+>\((\d+)\)', rw, e)
    def rc(m):
        return str(env['cvec'][int(m.group(1))][int(m.group(2))])
    e = re.sub(r'latset::c<[^>]+>\((\d+)\)\[(\d+)\]', rc, e)
    def rcu(m):
        k = int(m.group(1))
        return str(sum(env['cvec'][k][i]*env['u'][i] for i in range(env['d'])))
    e = re.sub(r'latset::c<[^>]+>\((\d+)\)\s*\*\s*u', rcu, e)
    e = re.sub(r'T\{([^}]+)\}', r'(\1)', e)
    # F[0], F[1], etc.
    def rF(m):
        f = env['F']
        return str(f[int(m.group(1))]) if isinstance(f, list) else str(f)
    e = re.sub(r'F\[(\d+)\]', rF, e)
    # F as scalar (only if not a vector/list in env)
    if isinstance(env.get('F'), (int, float)):
        e = re.sub(r'(?<!\w)F(?!\w)', str(env['F']), e)
    # Replace bare variable names from env (non-special, single-letter-prefix names)
    for var, val in env.items():
        if var in ('u', 'F', 'w', 'cvec', 'math', 'Fi', 'q', 'd',
                    'InvCs2', 'InvCs4', 'D2Q5', 'D2Q9', 'D3Q7',
                    'D3Q15', 'D3Q19', 'D3Q27', 'LatSet', 'T'):
            continue
        if isinstance(val, (int, float)):
            e = re.sub(rf'(?<!\w){re.escape(var)}(?!\w)', str(val), e)
    # u[d]
    def ru(m):
        return str(env['u'][int(m.group(1))])
    e = re.sub(r'u\[(\d+)\]', ru, e)
    try:
        return eval(e, {"__builtins__": {}}, {'math': math})
    except:
        return 0.0

def main():
    if len(sys.argv) < 3:
        print("Usage: verify_force.py <installed.ur.h> <generated.ur.h>")
        sys.exit(1)
    gen_path = sys.argv[2]
    fails = 0
    rng = random.Random(424242)

    # Test ForcePopImpl
    print("=== force.ur.h ForcePopImpl ===")
    for latname in LAT:
        lat = LAT[latname]
        ok = True
        for trial in range(3):
            u = [rng.uniform(-0.2, 0.4) for _ in range(lat['d'])]
            F = [rng.uniform(-0.1, 0.1) for _ in range(lat['d'])]
            ref = ref_force_pop(u, F, latname)
            body = parse_ur(gen_path, 'ForcePopImpl', latname)
            if body is None:
                print(f"FAIL {latname}: struct not found"); ok = False; fails += 1; break
            try:
                new = run_force_pop(body, u, F, latname)
                if any(abs(a-b) > 1e-9*max(1, abs(a)) for a, b in zip(ref, new)):
                    print(f"MISMATCH {latname}: ref={ref[:3]}... new={new[:3]}...")
                    ok = False; fails += 1; break
            except Exception as ex:
                print(f"ERR {latname}: {ex}"); ok = False; fails += 1; break
        if ok: print(f"OK   {latname}")

    # Test ScalarForcePopImpl
    print("\n=== force.ur.h ScalarForcePopImpl ===")
    for latname in LAT:
        lat = LAT[latname]
        ok = True
        for trial in range(3):
            u = [rng.uniform(-0.2, 0.4) for _ in range(lat['d'])]
            F = rng.uniform(-0.1, 0.1)
            results = parse_ur(gen_path, 'ScalarForcePopImpl', latname)
            if not results:
                print(f"FAIL {latname}: struct not found"); ok = False; fails += 1; break
            for d_idx, body in results.items():
                ref = ref_scalar_force_pop(u, F, d_idx, latname)
                try:
                    new = run_force_pop(body, u, F, latname, scalar_F=True)
                    if any(abs(a-b) > 1e-9*max(1, abs(a)) for a, b in zip(ref, new)):
                        for ii, (ra, nb) in enumerate(zip(ref, new)):
                            if abs(ra-nb) > 1e-9*max(1, abs(ra)):
                                print(f"MISMATCH {latname} d={d_idx} k={ii}: ref={ra} new={nb}")
                        ok = False; fails += 1; break
                except Exception as ex:
                    print(f"ERR {latname} d={d_idx}: {ex}"); ok = False; fails += 1; break
            if not ok: break
        if ok: print(f"OK   {latname}")

    print(f"\n{'ALL PASSED' if fails == 0 else 'SOME FAILED'}")
    sys.exit(1 if fails else 0)

main()
