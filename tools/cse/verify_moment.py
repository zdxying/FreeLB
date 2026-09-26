#!/usr/bin/env python3
"""Numerically compare formulas between old and new moment.ur.h files."""
import re, math, random, sys, ast
import re as _re

LAT = {
 'D2Q5':  dict(q=5,d=2,cs2=1/3,c=[(0,0),(1,0),(-1,0),(0,1),(0,-1)]),
 'D2Q9':  dict(q=9,d=2,cs2=1/3,c=[(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)]),
 'D3Q7':  dict(q=7,d=3,cs2=1/3,c=[(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]),
 'D3Q15': dict(q=15,d=3,cs2=1/3,c=[(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,1,1),(-1,-1,-1),(1,1,-1),(-1,-1,1),(1,-1,1),(-1,1,-1),(-1,1,1),(1,-1,-1)]),
 'D3Q19': dict(q=19,d=3,cs2=1/3,c=[(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,1,0),(-1,-1,0),(1,0,1),(-1,0,-1),(0,1,1),(0,-1,-1),(1,-1,0),(-1,1,0),(1,0,-1),(-1,0,1),(0,1,-1),(0,-1,1)]),
 'D3Q27': dict(q=27,d=3,cs2=1/3,c=[(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,1,0),(-1,-1,0),(1,0,1),(-1,0,-1),(0,1,1),(0,-1,-1),(1,-1,0),(-1,1,0),(1,0,-1),(-1,0,1),(0,1,-1),(0,-1,1),(1,1,1),(-1,-1,-1),(1,1,-1),(-1,-1,1),(1,-1,1),(-1,1,-1),(-1,1,1),(1,-1,-1)]),
}

def parse(path):
    src = open(path).read()
    fns = {}
    # struct Name<CELL<T, LAT<T>, TypePack>[, Extra]*>{
    for m in re.finditer(r'struct (\w+)<CELL<T, (\w+)<T>, TypePack>((?:, \w+)*)>\{(.*?)\n\};', src, re.S):
        name, lat, extra, body = m.group(1), m.group(2), m.group(3), m.group(4)
        am = re.search(r'apply\((.*?)\)\{(.*?)\n  \}', body, re.S)
        if not am:
            am = re.search(r'apply\((.*?)\)\{\n(.*?)\n\}', body, re.S)
        if not am: continue
        key = (name, lat)
        fns[key] = (am.group(1), am.group(2))
    return fns

def bcast(op, a, b):
    if isinstance(a, list) and isinstance(b, list): return [op(x, y) for x, y in zip(a, b)]
    if isinstance(a, list): return [op(x, b) for x in a]
    if isinstance(b, list): return [op(a, y) for y in b]
    return op(a, b)

def veval(node, env):
    if isinstance(node, ast.Expression): return veval(node.body, env)
    if isinstance(node, ast.Constant): return node.value
    if isinstance(node, ast.Name):
        if node.id == 'math': return math
        return env[node.id]
    if isinstance(node, ast.Subscript):
        return veval(node.value, env)[veval(node.slice, env)]
    if isinstance(node, ast.BinOp):
        import operator
        ops = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv}
        return bcast(ops[type(node.op)], veval(node.left, env), veval(node.right, env))
    if isinstance(node, ast.UnaryOp):
        import operator
        return -veval(node.operand, env)
    if isinstance(node, ast.Call):
        fn = veval(node.func, env)
        args = [veval(a, env) for a in node.args]
        if fn is math.sqrt:
            return [math.sqrt(x) for x in args[0]] if isinstance(args[0], list) else fn(*args)
        return fn(*args)
    raise ValueError(f"unsupported node {node}")

def translate(expr):
    e = expr
    e = _re.sub(r'cell\s*\.\s*template\s+get<[^()]*?>\s*\(\)', ' CELLGET ', e)
    e = _re.sub(r'cell\s*\.\s*get[A-Z]\w*\(\)', ' CELLGET ', e)
    e = e.replace(' ', '')
    e = _re.sub(r'cell\[(\d+)\]', r'c[\1]', e)
    e = e.replace('LatSet::cs2', 'CS2')
    e = e.replace('[ForceScheme::scalardir]', '[SD]')
    e = _re.sub(r'T\{([-\d.e]+)\}', r'\1', e)
    e = e.replace('std::sqrt', 'math.sqrt')
    return e

def run(sig, body, latname, cellvals):
    lat = LAT[latname]
    env = {'c': cellvals, 'CS2': lat['cs2'], 'math': math, 'CELLGET': 0.7, 'OMEGA': 0.9}
    u = [0.0]*lat['d']
    tensor = {}
    env['u_value'] = u
    env['tensor'] = tensor
    def split_params(sig):
        out, cur, ang = [], '', 0
        for ch in sig:
            if ch == '<': ang += 1
            elif ch == '>': ang -= 1
            if ch == ',' and ang == 0: out.append(cur); cur = ''
            else: cur += ch
        if cur.strip(): out.append(cur)
        return out
    for decl in split_params(sig):
        decl = decl.strip()
        nm = decl.replace('>>', '> >').split()[-1].replace('&','').strip()
        rng2 = random.Random(424242)   # identical fills for old/new runs
        isarr = 'std::array' in decl
        if nm == 'cell': continue
        elif isarr:
            n = lat['d']*(lat['d']+1)//2
            env[nm] = [rng2.uniform(-0.5, 1.5) for _ in range(n)] if 'const' in decl else [0.0]*max(n, 20)
        elif 'Vector' in decl:
            env[nm] = u if 'u_value' in nm else [rng2.uniform(-0.2,0.4) for _ in range(lat['d'])]
        else: env[nm] = 0.0
    for stmt in body.split(';'):
        s = stmt.strip()
        if not s or s.startswith('if') or s.startswith('using') or s.startswith('cell.template'): continue
        s = re.sub(r'^const\s+\w+\s+', '', s)
        s = re.sub(r'if\s*constexpr\s*\(WriteToField\)\s*', '', s)
        if s.startswith('cell.template'): continue
        if s.startswith('else'):
            s = s[4:].lstrip()
        dm = re.match(r'(?:const\s+)?(?:Vector<[^=]+?>|T)\s+(\w+)\s*=\s*(.+)$', s, re.S)
        if dm and not s.startswith(('tensor','strain_rate','stress')):
            nm = dm.group(1)
            val = veval(ast.parse(translate(dm.group(2)), mode='eval'), env)
            if isinstance(val, list): env[nm] = val
            else: env[nm] = val
            continue
        if s.startswith('return'):
            env['__ret__'] = veval(ast.parse(translate(s[len('return'):].strip()), mode='eval'), env)
            continue
        am = re.match(r'([\w]+)(?:\[([^\]]+)\])?\s*(\+|-|\*|/)?=\s*(.+)$', s, re.S)
        if not am: continue
        tgt, idx, op, expr = am.group(1), am.group(2), am.group(3), translate(am.group(4))
        val = veval(ast.parse(expr, mode='eval'), env)
        def applyop(cur):
            if op == '+': return cur + val
            if op == '-': return cur - val
            if op == '*': return cur * val
            if op == '/': return cur / val
            return val
        if idx is None:
            env[tgt] = applyop(env.get(tgt, 0.0)) if op else val
        else:
            container = env[tgt]
            i = int(eval(idx, {'SD': 0})) if not idx.isdigit() else int(idx)
            if isinstance(container, list): container[i] = applyop(container[i])
            else: container[(tgt,i)] = applyop(container[(tgt,i)])
    return env

def main(oldpath, newpath):
    old, new = parse(oldpath), parse(newpath)
    common = sorted(set(old) & set(new))
    print(f"old={len(old)} new={len(new)} common={len(common)}")
    fails = 0
    for key in common:
        only_old = [k for k in old if k not in new]
        name, lat = key
        L = LAT[lat]
        ok = True
        for trial in range(3):
            cv = [random.uniform(-0.1, 1.2) for _ in range(L['q'])]
            try:
                ro = run(*old[key], lat, list(cv))
                rn = run(*new[key], lat, list(cv))
            except Exception as ex:
                print(f"ERR {key}: {ex}"); ok=False; fails+=1; break
            def outNames(sig):
                names=[]
                ang=0; cur=''
                for ch in sig:
                    if ch=='<': ang+=1
                    elif ch=='>': ang-=1
                    if ch==',' and ang==0: names.append(cur); cur=''
                    else: cur+=ch
                if cur.strip(): names.append(cur)
                res=[]
                for decl in names:
                    d=decl.strip()
                    if 'const' in d or d.endswith('cell'): continue
                    nm=d.replace('>>','> >').split()[-1].replace('&','').strip()
                    if nm!='cell': res.append((nm, 'std::array' in d or 'Vector' in d))
                return res
            def collect(env, outs):
                out=[]
                for nm,isvec in outs:
                    v=env.get(nm, [] if isvec else None)
                    if isvec: out += list(v)
                    elif isinstance(v,float): out.append(v)
                if '__ret__' in env: out.append(env['__ret__'])
                return out
            vo = collect(ro, outNames(old[key][0]))
            vn = collect(rn, outNames(new[key][0]))
            if len(vo)!=len(vn) or any(abs(a-b)>1e-9*max(1,abs(a)) for a,b in zip(vo,vn)):
                print(f"MISMATCH {key}\n  old={vo}\n  new={vn}"); ok=False; fails+=1; break
        if ok: print(f"OK   {name:16s} {lat}")
    print(f"\nonly in OLD: {sorted(set(old)-set(new))}")
    print(f"only in NEW: {sorted(set(new)-set(old))}")
    sys.exit(1 if fails else 0)

main(sys.argv[1], sys.argv[2])
