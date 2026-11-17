#!/usr/bin/env python3
"""Make every image of a figure PDF a plain, opaque, 8-bit true-color image and remove all transparency.

Why: a journal PDF conversion showed the images of a figure as black boxes. They were stored as
palette ("Indexed") images whose palette contains a zero byte. Images stored as plain RGB came through,
with or without a soft mask. Matplotlib writes a palette image whenever an image has at most 256 colors,
so "no alpha" alone does not prevent black boxes.

What this script changes (text, lines and every other object are copied byte for byte):
  1. Every image that has a palette or a soft mask (/SMask) is decoded, composited onto white with its
     soft mask, and stored as 8-bit true color (Flate with PNG row filters: the form Matplotlib itself
     uses for true-color images).
  2. The soft-mask objects are removed.
  3. Lines drawn with partial opacity (Matplotlib alpha=...): the line color is composited onto white
     and the line becomes opaque. This is exact where the line lies on white.
  4. A stroke opacity of 0 (/CA 0, Matplotlib's "no outline"): invisible strokes are dropped.
  After 3 and 4 every graphics state is opaque. Semi-transparent fills, text or images stop the script:
  those must be fixed in the plotting script.

Usage:  python -I pdf_plain_images.py in.pdf out.pdf
Exit code 0 = written and verified; 2 = something is not supported, nothing written.
"""
import io
import os
import re
import sys
import time
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from PIL import Image
from minipdf import PDF, Stream, Ref, Name, content_tokens, write_pdf, reachable, png_unpredict, check_xref


class Unsupported(Exception):
    pass


def space_info(pdf, cs):
    """-> (kind, ncomp, white, base space, hival, lookup);  kind = 'plain' | 'indexed' | 'other'."""
    cs = pdf.resolve(cs)
    if isinstance(cs, str):
        if cs in ('DeviceRGB', 'RGB'):
            return 'plain', 3, 255, Name('DeviceRGB'), None, None
        if cs in ('DeviceGray', 'G'):
            return 'plain', 1, 255, Name('DeviceGray'), None, None
        if cs in ('DeviceCMYK', 'CMYK'):
            return 'plain', 4, 0, Name('DeviceCMYK'), None, None
    elif isinstance(cs, list) and cs:
        head = pdf.resolve(cs[0])
        if head == 'ICCBased':
            n = pdf.resolve(pdf.resolve(cs[1]).dict.get('N'))
            if n in (1, 3, 4):
                return 'plain', n, (0 if n == 4 else 255), cs, None, None
        elif head == 'CalRGB':
            return 'plain', 3, 255, cs, None, None
        elif head == 'CalGray':
            return 'plain', 1, 255, cs, None, None
        elif head in ('Indexed', 'I'):
            kind, n, white, base, _, _ = space_info(pdf, cs[1])
            if kind == 'plain':
                lut = pdf.resolve(cs[3])
                if isinstance(lut, Stream):
                    lut = pdf.decode(lut)
                return 'indexed', n, white, base, pdf.resolve(cs[2]), lut
    return 'other', None, None, cs, None, None


def unpack(data, W, H, bpc, ncomp=1):
    """Samples of an image as (H, W, ncomp) integers 0 .. 2**bpc - 1."""
    if bpc == 8:
        if len(data) < W * H * ncomp:
            raise Unsupported(f'image data too short ({len(data)} < {W * H * ncomp})')
        return np.frombuffer(data, np.uint8)[:W * H * ncomp].reshape(H, W, ncomp)
    if bpc not in (1, 2, 4):
        raise Unsupported(f'{bpc} bits per component')
    rowlen = (W * ncomp * bpc + 7) // 8
    rows = np.frombuffer(data, np.uint8)[:rowlen * H].reshape(H, rowlen)
    bits = np.unpackbits(rows, axis=1)[:, :W * ncomp * bpc].reshape(H, W * ncomp, bpc)
    v = np.zeros((H, W * ncomp), np.uint8)
    for k in range(bpc):
        v = (v << 1) | bits[:, :, k]
    return v.reshape(H, W, ncomp)


def default_decode(pdf, d, ncomp, indexed, bpc):
    dec = pdf.resolve(d.get('Decode'))
    if dec is None:
        return True
    dec = [float(pdf.resolve(x)) for x in dec]
    want = [0.0, float(2 ** bpc - 1)] if indexed else [0.0, 1.0] * ncomp
    return dec == want


def last_filter(pdf, d):
    f = pdf.resolve(d.get('Filter'))
    f = f if isinstance(f, list) else ([] if f is None else [f])
    return str(pdf.resolve(f[-1])) if f else None


def image_array(pdf, st):
    """Pixels of an image in its base color space: (H, W, ncomp) uint8, plus color-space facts."""
    d = st.dict
    W, H = pdf.resolve(d['Width']), pdf.resolve(d['Height'])
    bpc = pdf.resolve(d.get('BitsPerComponent', 8))
    kind, n, white, base, hival, lut = space_info(pdf, d.get('ColorSpace'))
    if kind == 'other':
        raise Unsupported(f'color space {base!r}')
    if not default_decode(pdf, d, n, kind == 'indexed', bpc):
        raise Unsupported('non-default /Decode array')
    if last_filter(pdf, d) in ('DCTDecode', 'DCT'):
        if kind == 'indexed' or n == 4:
            raise Unsupported('JPEG image with palette or CMYK')
        im = Image.open(io.BytesIO(st.raw)).convert('RGB' if n == 3 else 'L')
        arr = np.asarray(im).reshape(H, W, n)
    else:
        data = pdf.decode(st)
        if kind == 'indexed':
            idx = unpack(data, W, H, bpc)[:, :, 0]
            if len(lut) < (hival + 1) * n:
                raise Unsupported('palette shorter than hival + 1 entries')
            if int(idx.max()) > hival:
                raise Unsupported('pixel index above hival')
            pal = np.frombuffer(lut[:(hival + 1) * n], np.uint8).reshape(hival + 1, n)
            arr = pal[idx]
        else:
            if bpc != 8:
                raise Unsupported(f'{bpc} bits per component in a non-palette image')
            arr = unpack(data, W, H, 8, n)
    return np.ascontiguousarray(arr), kind, n, white, base, hival, lut


def smask_array(pdf, sm, W, H):
    sm = pdf.resolve(sm)
    if not isinstance(sm, Stream):
        raise Unsupported('soft mask is not a stream')
    d = sm.dict
    if 'Matte' in d:
        raise Unsupported('premultiplied soft mask (/Matte)')
    if (pdf.resolve(d['Width']), pdf.resolve(d['Height'])) != (W, H):
        raise Unsupported('soft mask size differs from image size')
    kind, n, _, _, _, _ = space_info(pdf, d.get('ColorSpace'))
    bpc = pdf.resolve(d.get('BitsPerComponent', 8))
    if kind != 'plain' or n != 1 or not default_decode(pdf, d, 1, False, bpc):
        raise Unsupported('soft mask is not a plain gray image')
    if last_filter(pdf, d) in ('DCTDecode', 'DCT'):
        a = np.asarray(Image.open(io.BytesIO(sm.raw)).convert('L')).reshape(H, W)
    else:
        a = unpack(pdf.decode(sm), W, H, bpc)[:, :, 0]
        if bpc < 8:
            a = (a.astype(np.uint16) * 255 // (2 ** bpc - 1)).astype(np.uint8)
    return a


def png_rows(arr):
    """(H, W, 1 or 3) uint8 -> zlib data with PNG row filters, exactly as Matplotlib stores images."""
    im = Image.fromarray(arr[:, :, 0] if arr.shape[2] == 1 else arr)
    buf = io.BytesIO()
    im.save(buf, format='png')
    b, pos, out = buf.getvalue(), 8, bytearray()
    while pos < len(b):
        ln, typ = int.from_bytes(b[pos:pos + 4], 'big'), b[pos + 4:pos + 8]
        if typ == b'IDAT':
            out += b[pos + 8:pos + 8 + ln]
        pos += 12 + ln
    return bytes(out)


PATH_OPS = {'m', 'l', 'c', 'v', 'y', 'h', 're'}
PAINT_OPS = {'S', 's', 'f', 'F', 'f*', 'B', 'B*', 'b', 'b*', 'n'}
STROKE_OPS = {'S', 's', 'B', 'B*', 'b', 'b*'}
FILL_OPS = {'f', 'F', 'f*', 'B', 'B*', 'b', 'b*', 'sh'}
TEXT_OPS = {'Tj', 'TJ', "'", '"'}
DROP_STROKE = {'S': b'n', 's': b'n', 'B': b'f', 'B*': b'f*', 'b': b'h f', 'b*': b'h f*'}
COLOR_OP = {1: 'G', 3: 'RG', 4: 'K'}


def num(x):
    s = ('%.6f' % x).rstrip('0').rstrip('.')
    return s if s not in ('', '-0') else '0'


def scan_content(pdf, data, resources, state, events, where, depth=0):
    """Walk a content stream; record every painting operation done with an opacity below 1.
    state = [stroke opacity CA, fill opacity ca, graphics-state name, stroke color (tuple) or None]."""
    res = pdf.resolve(resources) or {}
    gsd = pdf.resolve(res.get('ExtGState')) or {}
    xob = pdf.resolve(res.get('XObject')) or {}
    stack, operands, starts, tr, path_start = [], [], [], 0, None
    used = events.setdefault('used', set())
    for v, s, e in content_tokens(data):
        if not (isinstance(v, Name) and v.startswith('@')):
            operands.append(v)
            starts.append(s)
            continue
        op = v[1:]
        if op == 'q':
            stack.append((list(state), tr))
        elif op == 'Q':
            if stack:
                st, tr = stack.pop()
                state[:] = st
        elif op == 'gs' and operands:
            name = str(operands[-1])
            used.add(name)
            g = pdf.resolve(gsd.get(name)) or {}
            if 'CA' in g:
                state[0] = float(pdf.resolve(g['CA']))
            if 'ca' in g:
                state[1] = float(pdf.resolve(g['ca']))
            state[2] = name
            sm, bm = pdf.resolve(g.get('SMask')), pdf.resolve(g.get('BM'))
            if sm is not None and sm != 'None':
                events['hard'].append((where, name, 'soft mask in a graphics state'))
            if bm not in (None, 'Normal', 'Compatible'):
                events['hard'].append((where, name, f'blend mode {bm}'))
        elif op in ('G', 'RG', 'K'):
            state[3] = tuple(float(x) for x in operands[-{'G': 1, 'RG': 3, 'K': 4}[op]:])
        elif op == 'CS':
            state[3] = {'DeviceGray': (0.0,), 'DeviceRGB': (0.0, 0.0, 0.0),
                        'DeviceCMYK': (0.0, 0.0, 0.0, 1.0)}.get(str(operands[-1]) if operands else '')
        elif op in ('SC', 'SCN'):
            ok = state[3] is not None and len(operands) == len(state[3]) and \
                all(isinstance(x, (int, float)) for x in operands)
            state[3] = tuple(float(x) for x in operands) if ok else None
        elif op == 'Tr' and operands:
            tr = operands[-1]
        elif op == 'Do' and operands:
            x = pdf.resolve(xob.get(str(operands[-1])))
            if isinstance(x, Stream) and x.dict.get('Subtype') == 'Form':
                if depth < 8:
                    scan_content(pdf, pdf.decode(x), x.dict.get('Resources', resources), list(state),
                                 events, where + '>form ' + str(operands[-1]), depth + 1)
            elif state[1] < 1:
                events['fill'].append((where, state[2], state[1], op))
        elif op in PATH_OPS:
            if path_start is None:
                path_start = starts[0] if starts else s
        else:
            strokes = op in STROKE_OPS or (op in TEXT_OPS and tr in (1, 2, 5, 6))
            fills = op in FILL_OPS or (op in TEXT_OPS and tr in (0, 2, 4, 6))
            if strokes and state[0] < 1:
                events['stroke'].append((where, state[2], state[0], op, s, e,
                                         path_start if op in STROKE_OPS else None, state[3]))
            if fills and state[1] < 1:
                events['fill'].append((where, state[2], state[1], op))
            if op in PAINT_OPS:
                path_start = None
        operands, starts = [], []


def main(src, dst):
    pdf = PDF(src)
    if any(w == 'objstm' for w in pdf.where.values()):
        raise Unsupported('PDF with object streams (PDF 1.5+): not handled by this script')
    pages = pdf.pages()
    print(f'{src}: {len(pdf.buf)} B, {len(pdf.objs)} objects, {len(pages)} page(s)')
    new = {}                                                     # num -> replacement object

    # 1. images ---------------------------------------------------------------------------------
    print('\nImages (obj: size, stored as -> written as):')
    nconv = 0
    for onum in sorted(pdf.objs):
        st = pdf.objs[onum]
        if not (isinstance(st, Stream) and st.dict.get('Subtype') == 'Image'):
            continue
        d = st.dict
        if pdf.resolve(d.get('ImageMask')) is True:
            raise Unsupported(f'obj {onum}: stencil mask image')
        if 'Mask' in d:
            raise Unsupported(f'obj {onum}: /Mask (color-key or stencil mask)')
        kind = space_info(pdf, d.get('ColorSpace'))[0]
        if kind != 'indexed' and 'SMask' not in d:
            continue                                             # already plain and opaque: untouched
        arr, kind, n, white, base, hival, lut = image_array(pdf, st)
        H, W = arr.shape[:2]
        if kind == 'indexed':
            nz = lut[:(hival + 1) * n].count(bytes(1))
            was = f'palette of {hival + 1} colors' + (f' with {nz} zero byte(s)' if nz else '')
        else:
            was = f'{n}-component true color'
        note = ''
        if 'SMask' in d:
            a = smask_array(pdf, d['SMask'], W, H)
            a32 = a.astype(np.uint32)[:, :, None]
            out = ((arr.astype(np.uint32) * a32 + (255 - a32) * white + 127) // 255).astype(np.uint8)
            was += ' + soft mask'
            note = (f'; {int((a < 255).sum())} of {W * H} pixels not opaque ({int((a == 0).sum())} fully '
                    f'transparent), {int((out != arr).any(axis=2).sum())} pixel values changed by compositing onto white')
        else:
            out = arr
        if n == 4:
            raw, parms = zlib.compress(out.tobytes(), 9), None
        else:
            raw = png_rows(out)
            parms = {'Predictor': 10, 'Colors': n, 'Columns': W}
            assert png_unpredict(zlib.decompress(raw), W, n, 8) == out.tobytes(), f'obj {onum}: round trip failed'
        nd = {'Type': Name('XObject'), 'Subtype': Name('Image'), 'Width': W, 'Height': H,
              'ColorSpace': base, 'BitsPerComponent': 8}
        for k in ('Interpolate', 'Intent', 'Name', 'OC', 'Metadata', 'StructParent'):
            if k in d:
                nd[k] = d[k]
        nd['Filter'] = Name('FlateDecode')
        if parms:
            nd['DecodeParms'] = parms
        new[onum] = Stream(nd, raw)
        nconv += 1
        print(f'  obj {onum:3d}: {W}x{H}, {was} -> opaque {n}-component 8-bit'
              f' ({len(st.raw)} -> {len(raw)} B){note}')
    print(f'  {nconv} image(s) rewritten')

    # 2. graphics states ------------------------------------------------------------------------
    events = {'stroke': [], 'fill': [], 'hard': []}
    for i, page in enumerate(pages):
        scan_content(pdf, pdf.page_content(page), page.get('Resources'), [1.0, 1.0, None, (0.0,)],
                     events, f'page {i + 1}')
    print('\nGraphics states in use:', ', '.join(sorted(events['used'])) or 'none')
    bad = list(events['hard']) + list(events['fill'])
    bad += [ev[:4] for ev in events['stroke'] if '>form' in ev[0]]
    bad += [ev[:4] for ev in events['stroke'] if ev[2] > 0 and (ev[6] is None or ev[7] is None)]
    if bad:
        for ev in bad[:20]:
            print('  NOT SUPPORTED:', ev)
        raise Unsupported(f'{len(bad)} painting operation(s) with transparency that this script cannot remove')
    drop = [ev for ev in events['stroke'] if ev[2] == 0]
    part = [ev for ev in events['stroke'] if ev[2] > 0]
    print(f'  invisible strokes (opacity 0) dropped: {len(drop)}')
    kinds = {}
    for ev in part:
        blend = tuple(ev[2] * c + (1 - ev[2]) * (0.0 if len(ev[7]) == 4 else 1.0) for c in ev[7])
        kinds.setdefault((ev[1], ev[2], ev[7], blend), 0)
        kinds[(ev[1], ev[2], ev[7], blend)] += 1
    for (name, alpha, col, blend), cnt in kinds.items():
        print(f'  {cnt} line(s) in state /{name}: opacity {alpha:g}, color {" ".join(num(c) for c in col)}'
              f' -> opaque, color {" ".join(num(c) for c in blend)} (composited onto white)')
    if drop or part:
        if len(pages) != 1 or isinstance(pdf.resolve(pages[0].get('Contents')), list):
            raise Unsupported('transparent strokes in a multi-stream or multi-page file')
        cref = pages[0]['Contents']
        data = pdf.page_content(pages[0])
        edits = []                                               # (offset, end, replacement)
        for where, name, alpha, op, s, e, pstart, col in drop:
            edits.append((s, e, DROP_STROKE[op]))
        for where, name, alpha, op, s, e, pstart, col in part:
            cop = COLOR_OP[len(col)]
            blend = [alpha * c + (1 - alpha) * (0.0 if len(col) == 4 else 1.0) for c in col]
            edits.append((pstart, pstart, (' '.join(num(c) for c in blend) + f' {cop}\n').encode()))
            edits.append((e, e, ('\n' + ' '.join(num(c) for c in col) + f' {cop}').encode()))
        buf = bytearray(data)
        for s, e, rep in sorted(edits, key=lambda t: (-t[0], -t[1])):
            buf[s:e] = rep
        cd = {k: v for k, v in pdf.resolve(cref).dict.items() if k not in ('Filter', 'DecodeParms', 'Length')}
        cd['Filter'] = Name('FlateDecode')
        new[cref.num] = Stream(cd, zlib.compress(bytes(buf), 9))
        print(f'  page content rewritten: {len(data)} -> {len(buf)} B ({len(edits)} edits)')

    # every graphics state with an opacity below 1 is now without effect: make it opaque
    containers, gsobjs, changed, mark = set(), set(), [], set()

    def fix(g, label):
        # all opacities are 1 from here on, which is the default: the entries are removed altogether,
        # so that a plain search for /CA or /ca in the file finds nothing
        hit = False
        for k in ('CA', 'ca'):
            if k in g:
                if float(pdf.resolve(g[k])) < 1:
                    changed.append(f'{label} /{k} {pdf.resolve(g[k])}')
                del g[k]
                hit = True
        return hit

    def do_container(c, top):
        for name, g in c.items():
            if isinstance(g, Ref):
                gsobjs.add(g.num)
            elif isinstance(g, dict) and fix(g, f'/{name}'):
                mark.add(top)

    def find(x, top):
        if isinstance(x, Stream):
            find(x.dict, top)
        elif isinstance(x, dict):
            for k, v in x.items():
                if k == 'ExtGState':
                    if isinstance(v, Ref):
                        containers.add(v.num)
                    elif isinstance(v, dict):
                        do_container(v, top)
                else:
                    find(v, top)
        elif isinstance(x, list):
            for v in x:
                find(v, top)

    for onum, o in pdf.objs.items():
        find(o, onum)
    for onum in sorted(containers):
        if isinstance(pdf.objs.get(onum), dict):
            do_container(pdf.objs[onum], onum)
    for onum in sorted(gsobjs):
        if isinstance(pdf.objs.get(onum), dict) and fix(pdf.objs[onum], f'obj {onum}'):
            mark.add(onum)
    for onum in mark:
        new.setdefault(onum, pdf.objs[onum])
    print('  opacity entries below 1 removed from the graphics states:', '; '.join(changed) or 'none')

    # 3. document information: say what was done --------------------------------------------------
    trailer = pdf.trailer()
    if isinstance(trailer.get('Info'), Ref) and isinstance(pdf.objs.get(trailer['Info'].num), dict):
        info = dict(pdf.objs[trailer['Info'].num])
        prod = info.get('Producer', b'')
        prod = prod if isinstance(prod, bytes) else b''
        info['Producer'] = prod + b'; images rewritten as opaque 8-bit true color, no transparency (pdf_plain_images.py)'
        lt = time.localtime()
        off = lt.tm_gmtoff
        info['ModDate'] = (time.strftime('D:%Y%m%d%H%M%S', lt) +
                           '%s%02d\'%02d\'' % ('+' if off >= 0 else '-', abs(off) // 3600, abs(off) % 3600 // 60)).encode()
        new[trailer['Info'].num] = info

    # 4. write: untouched objects byte for byte, unreferenced objects dropped --------------------
    final = dict(pdf.objs)
    final.update(new)
    keep = reachable(final, trailer)
    dropped = sorted(set(final) - keep)
    raw = {n: pdf.buf[pdf.span[n][0]:pdf.span[n][1]] for n in keep if n not in new}
    ver = pdf.buf[5:8] if pdf.buf[:5] == b'%PDF-' else b'1.4'
    size = write_pdf(dst, {n: final[n] for n in keep if n in new}, trailer, version=ver, raw=raw)
    print(f'\n{dst}: {size} B; {len(raw)} objects copied unchanged, {len([n for n in new if n in keep])} rewritten, '
          f'{len(dropped)} dropped (soft masks and stream-length objects)')

    # 5. verify the result ------------------------------------------------------------------------
    chk = PDF(dst)
    nimg = npal = nsm = 0
    for o in chk.objs.values():
        if isinstance(o, Stream) and o.dict.get('Subtype') == 'Image':
            nimg += 1
            npal += space_info(chk, o.dict.get('ColorSpace'))[0] == 'indexed'
            nsm += 'SMask' in o.dict or 'Mask' in o.dict
            chk.decode(o)                                        # every image stream decodes
    ev2 = {'stroke': [], 'fill': [], 'hard': []}
    for i, page in enumerate(chk.pages()):
        scan_content(chk, chk.page_content(page), page.get('Resources'), [1.0, 1.0, None, (0.0,)], ev2, f'page {i + 1}')
    grep = len(re.findall(rb'/SMask|/CA |/ca ', chk.buf))
    below = [m for m in re.findall(rb'/(?:CA|ca) ([0-9.]+)', chk.buf) if float(m) < 1]
    xprob = check_xref(chk.buf)
    for p in xprob:
        print('  cross-reference table:', p)
    ok = not xprob and npal == 0 and nsm == 0 and not below and grep == 0 and not ev2['stroke'] and not ev2['fill'] and not ev2['hard'] \
        and len(chk.pages()) == len(pages) and len(chk.objs) == len(keep)
    print(f'verify: {len(chk.objs)} objects, {nimg} images, {npal} with palette, {nsm} with mask; '
          f'opacity values below 1: {len(below)}; grep "/SMask|/CA |/ca " matches: {grep}; '
          f'result {"OK" if ok else "FAILED"}')
    return 0 if ok else 1


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    if os.path.exists(sys.argv[2]):
        sys.exit(f'{sys.argv[2]} exists: not overwritten')
    try:
        sys.exit(main(sys.argv[1], sys.argv[2]))
    except Unsupported as e:
        print('NOT WRITTEN:', e)
        sys.exit(2)
