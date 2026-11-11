#!/usr/bin/env python3
"""Render SHA-checked simulator frames as SVG observation stills.

These are mathematical-state screenshots, NOT photographs of physical vacuum,
atoms or molecules. The 45-site ring is displayed as a 3-by-15 strip.
"""
from hashlib import sha256
from pathlib import Path
from html import escape
import sys

P=1000081
D=45

def frame_lines(path):
    out={}
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        if not line.startswith('MOLECULAR_FRAME|'):
            continue
        d=dict(part.split('=',1) for part in line.split('|')[1:] if '=' in part)
        values=[int(x) for x in d['values'].split(',')]
        if len(values)!=D or any(x<0 or x>=P for x in values):
            raise ValueError('BAD_FRAME_VALUES')
        digest=sha256(b''.join(v.to_bytes(8,'big') for v in values)).hexdigest()
        if digest!=d['sha256']:
            raise ValueError('FRAME_HASH_MISMATCH')
        out[(d['arm'],int(d['coarse_steps']))]=(d,values)
    return out

def image(d,vals):
    arm=d['arm']; n=int(d['coarse_steps'])
    out=['<svg xmlns="http://www.w3.org/2000/svg" width="1250" height="560" viewBox="0 0 1250 560">',
         '<rect width="1250" height="560" fill="#101522"/>',
         f'<text x="48" y="60" fill="white" font-family="monospace" font-size="28">SIMULATED FRAME: {escape(arm)} / {n} ticks</text>',
         '<text x="48" y="101" fill="#c1cad8" font-family="monospace" font-size="17">45-site periodic ring - coordinate-scaled transport, not a physical camera</text>',
         f'<text x="48" y="138" fill="#8ea9ff" font-family="monospace" font-size="15">Underlying clock count: {d["planck_ticks"]} (Planck units assumed)</text>']
    for idx,v in enumerate(vals):
        col,row=idx%15,idx//15
        x,y=48+col*73,177+row*95
        bg='#131c2c' if v==0 else '#15566a' if v==1 else '#1d98b3'
        border='#364256' if v==0 else '#ffd480' if v==1 else '#6fd7ff'
        out.append(f'<rect x="{x}" y="{y}" width="62" height="68" rx="8" fill="{bg}" stroke="{border}" stroke-width="2"/>')
        out.append(f'<text x="{x+7}" y="{y+18}" font-family="monospace" fill="#a2b4cc" font-size="11">{idx:02d}</text>')
        out.append(f'<text x="{x+31}" y="{y+49}" font-family="monospace" fill="white" text-anchor="middle" font-size="26">{v}</text>')
    out.append('<text x="48" y="512" fill="#c1cad8" font-family="monospace" font-size="14">Cyan = pulse occupancy 2 | Gold = injected occupancy 1 | Dark = empty background</text>')
    out.append(f'<text x="48" y="543" fill="#a2b4cc" font-family="monospace" font-size="13">SHA-256 state: {d["sha256"]}</text>')
    out.append('</svg>')
    return '\n'.join(out)+'\n'

def main():
    if len(sys.argv)!=3:
        raise SystemExit('Usage: render_molecular_shots.py INPUT_HBP OUTPUT_FOLDER')
    frames=frame_lines(sys.argv[1])
    folder=Path(sys.argv[2]);folder.mkdir(parents=True,exist_ok=True)
    chosen=(('ideal_vacuum',0),('isolated_pulse',0),('isolated_pulse',1),
            ('isolated_pulse',9),('isolated_pulse',27),('single_perturbation',9))
    for arm,n in chosen:
        d,values=frames[(arm,n)]
        name=f'shot_{arm}_{n:03}.svg'
        s=image(d,values)
        (folder/name).write_text(s,encoding='utf-8')
        print(f'MOLECULAR_SHOT|file={name}|source_sha256={d["sha256"]}|svg_sha256={sha256(s.encode()).hexdigest()}|status=PASS|json=0')
    print('MOLECULAR_SHOTS|count=6|all_source_hashes_verified=1|status=PASS|json=0')
if __name__=='__main__':
    main()
