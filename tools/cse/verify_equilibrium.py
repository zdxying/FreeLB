#!/usr/bin/env python3
"""Numerically compare equilibrium.ur.h specializations.
For each (struct, LatSet) pair, evaluate both old and new apply() functions
with random inputs and verify the feq arrays match."""
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
    """Compute weight for direction k."""
    lat = LAT[latname]
    c2 = sum(x*x for x in lat['c'][k])
    if lat['d'] == 2:
        if c2 == 0: return 4.0/9.0
        elif c2 == 1: return 1.0/9.0
        else: return 1.0/36.0
    else:
        if c2 == 0: return 2.0/9.0
        elif c2 == 1: return 1.0/9.0
        elif c2 == 2: return 1.0/72.0
        else: return 1.0/216.0

def eval_second_order(u, rho, latname):
    """Evaluate SecondOrderImpl::apply using the template source formula."""
    lat = LAT[latname]
    q, d, cs2 = lat['q'], lat['d'], lat['cs2']
    InvCs2 = 1.0/cs2
    InvCs4 = InvCs2*InvCs2
    u2 = sum(x*x for x in u[:d])
    feq = [0.0]*q
    for k in range(q):
        c_k = lat['c'][k]
        uc = sum(u[i]*c_k[i] for i in range(d))
        wk = w_k(latname, k)
        feq[k] = wk * rho * (1.0 + InvCs2*uc + uc*uc*0.5*InvCs4 - InvCs2*u2*0.5)
    return feq

def parse_and_run(path, latname):
    """Parse the .ur.h file and evaluate SecondOrderImpl for a given LatSet."""
    src = open(path).read()
    lat = LAT[latname]
    q, d = lat['q'], lat['d']

    # Find the struct for this LatSet
    pattern = rf'struct SecondOrderImpl<CELL<T, {latname}<T>, TypePack>>\{{(.*?)\n\}};'
    m = re.search(pattern, src, re.S)
    if not m:
        raise ValueError(f"struct for {latname} not found")
    body = m.group(1)

    # Extract the apply function body
    am = re.search(r'apply\([^)]*\)\{(.*?)\n  \}', body, re.S)
    if not am:
        am = re.search(r'apply\([^)]*\)\{\n(.*?)\n\}', body, re.S)
    if not am:
        raise ValueError(f"apply function not found for {latname}")
    apply_body = am.group(1)

    return apply_body

def run_apply_body(body, u, rho, latname):
    """Execute an apply() function body and return feq."""
    lat = LAT[latname]
    q, d, cs2 = lat['q'], lat['d'], lat['cs2']
    InvCs2 = 1.0/cs2
    InvCs4 = InvCs2*InvCs2

    # Build weight and c vectors
    w = [w_k(latname, k) for k in range(q)]
    cvecs = [list(lat['c'][k]) for k in range(q)]

    # Local variables
    feq = [0.0]*q
    env = {'u': u, 'rho': rho, 'w': w, 'cvec': cvecs, 'InvCs2': InvCs2, 'InvCs4': InvCs4,
           'feq': feq, 'q': q, 'd': d, 'math': math}

    # Execute statements
    stmts = [s.strip() for s in body.split(';') if s.strip()]
    for s in stmts:
        # Skip using declarations, if statements, etc.
        if s.startswith('using') or s.startswith('if') or s.startswith('cell'):
            continue
        # Remove const/constexpr qualifiers and optional type
        s = re.sub(r'^(constexpr|const)\s+\w+\s+', '', s)
        s = re.sub(r'^(constexpr|const)\s+', '', s)
        # Handle variable declarations: T var = expr  or  var = expr
        dm = re.match(r'(?:T\s+)?(\w+)\s*=\s*(.+)', s)
        if dm:
            nm, expr = dm.group(1), dm.group(2)
            env[nm] = eval_expr(expr, env)
            continue
        # Handle feq[k] = expr
        fm = re.match(r'feq\[(\d+)\]\s*=\s*(.+)', s)
        if fm:
            k, expr = int(fm.group(1)), fm.group(2)
            feq[k] = eval_expr(expr, env)
            continue

    return feq

def eval_expr(expr, env):
    """Evaluate a simple expression."""
    e = expr.strip()
    # Replace LatSet constants
    e = e.replace('LatSet::InvCs2', str(env['InvCs2']))
    e = e.replace('LatSet::InvCs4', str(env['InvCs4']))
    # Replace latset::w<LatSet>(k)
    def replace_w(m):
        k = int(m.group(1))
        return str(env['w'][k])
    e = re.sub(r'latset::w<LatSet>\((\d+)\)', replace_w, e)
    # Replace latset::c<LatSet>(k)[d]
    def replace_c(m):
        k, d = int(m.group(1)), int(m.group(2))
        return str(env['cvec'][k][d])
    e = re.sub(r'latset::c<LatSet>\((\d+)\)\[(\d+)\]', replace_c, e)
    # Replace T{val}
    e = re.sub(r'T\{([^}]+)\}', r'(\1)', e)
    # Replace u.getnorm2()
    u = env.get('u', [])
    d = env.get('d', len(u))
    u2 = sum(x*x for x in u[:d])
    e = re.sub(r'u\.getnorm2\(\)', str(u2), e)
    # Replace u[d]
    def replace_u(m):
        idx = int(m.group(1))
        return str(env['u'][idx])
    e = re.sub(r'u\[(\d+)\]', replace_u, e)
    # Replace variable references (skip 'u' to avoid clobbering)
    for k, v in env.items():
        if k == 'u':
            continue
        if isinstance(v, (int, float)):
            e = re.sub(r'\b' + re.escape(k) + r'\b', str(v), e)
    # Evaluate
    try:
        return eval(e, {"__builtins__": {}}, {'math': math})
    except:
        return 0.0

def main(oldpath, newpath):
    fails = 0
    rng = random.Random(424242)

    for latname in LAT:
        lat = LAT[latname]
        q, d = lat['q'], lat['d']

        for trial in range(3):
            u = [rng.uniform(-0.2, 0.4) for _ in range(d)]
            rho = rng.uniform(0.5, 1.5)

            # Compute reference using template formula
            ref_feq = eval_second_order(u, rho, latname)

            # Parse and run old code
            try:
                old_body = parse_and_run(oldpath, latname)
                old_feq = run_apply_body(old_body, u, rho, latname)
            except Exception as ex:
                print(f"ERR old {latname}: {ex}")
                fails += 1
                break

            # Parse and run new code
            try:
                new_body = parse_and_run(newpath, latname)
                new_feq = run_apply_body(new_body, u, rho, latname)
            except Exception as ex:
                print(f"ERR new {latname}: {ex}")
                fails += 1
                break

            # Compare
            if any(abs(a-b) > 1e-9*max(1, abs(a)) for a, b in zip(ref_feq, old_feq)):
                print(f"MISMATCH old {latname} trial {trial}: ref={ref_feq} old={old_feq}")
                fails += 1
                break
            if any(abs(a-b) > 1e-9*max(1, abs(a)) for a, b in zip(ref_feq, new_feq)):
                print(f"MISMATCH new {latname} trial {trial}: ref={ref_feq} new={new_feq}")
                fails += 1
                break

        if fails == 0:
            print(f"OK   {latname}")
        else:
            break

    print(f"\n{'PASS' if fails == 0 else 'FAIL'}")
    sys.exit(1 if fails else 0)

main(sys.argv[1], sys.argv[2])
