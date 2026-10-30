#!/usr/bin/env python
"""Writes cpp/tests/fixtures/cpython/py_text.txt, the oracle of hicx/py_text.hpp.

Random and exhaustive cases for str.isspace / strip / split / split(sep), strict
UTF-8 decoding and universal-newline reading, as the interpreter does them:

    ~/miniconda3/envs/hicx-ref-3.7.6/bin/python cpp/scripts/gen_pytext_fixture.py

Reference: Python 3.12.7. Byte strings are written in hex, '-' for empty.
Records:
  version <python>
  spaces <n> <code point>*n    every code point c with chr(c).isspace()
  t <hex>                      a str, encoded utf-8 with surrogatepass
  s <hex>                      its strip()
  w <n> <hex>*n                its split()
  p <sep hex> <n> <hex>*n      its split(sep)
  r <hex>                      ''.join(split())
  v <hex> <k>                  bytes; k = UnicodeDecodeError.start, or -1
  e <hex> <hex>                bytes decoded utf-8/surrogateescape and its repr()
  reprall <digest>             FNV-1a over repr(chr(c)) for every code point
                               c except surrogates, each followed by a 0 byte
  n <hex> <count> <hex>*count  a valid UTF-8 text file and readlines() of it
                               opened with newline=None
"""
import os
import random
import sys
import tempfile


MASK64 = (1 << 64) - 1


def hx(b):
    return b.hex() if b else '-'


def enc(s):
    return s.encode('utf-8', 'surrogatepass')


def main():
    rng = random.Random(1030)
    spaces = [c for c in range(0x110000) if chr(c).isspace()]
    out = ['version ' + sys.version.split()[0],
           'spaces %d %s' % (len(spaces), ' '.join(map(str, spaces)))]
    pool = ([chr(c) for c in spaces] * 3
            + list('abcXYZ019_-,.;:\t|') * 2
            + ['é', 'ß', '中', '\U0001F600', '​', '⁠', '﻿',
               '᠎', '\ud800', '\udfff', '\u0085', ' ', '　', '\x1b', '\x1c', '\x1f', '\x7f'])
    seps = [',', '\t', '::', 'ab', ' ', '中', ' ', ' x', 'a', '\r\n']
    for _ in range(4000):
        n = rng.choice([0, 1, 2, 3, 5, 8, 12, 20, 40])
        s = ''.join(rng.choice(pool) for _ in range(n))
        sep = rng.choice(seps)
        out.append('t ' + hx(enc(s)))
        out.append('s ' + hx(enc(s.strip())))
        w = s.split()
        out.append('w %d %s' % (len(w), ' '.join(hx(enc(x)) for x in w)))
        p = s.split(sep)
        out.append('p %s %d %s' % (hx(enc(sep)), len(p), ' '.join(hx(enc(x)) for x in p)))
        out.append('r ' + hx(enc(''.join(w))))
    # invalid UTF-8: random bytes, mutated valid text, truncations
    frag = [b'\x00', b'a', b'\x7f', b'\x80', b'\xbf', b'\xc0', b'\xc1', b'\xc2', b'\xdf',
            b'\xe0', b'\xe1', b'\xec', b'\xed', b'\xee', b'\xef', b'\xf0', b'\xf1', b'\xf3',
            b'\xf4', b'\xf5', b'\xf8', b'\xff', b'\x90', b'\x9f', b'\xa0', b'\x8f', b'\xbe']
    for _ in range(6000):
        n = rng.choice([1, 2, 3, 4, 5, 6, 10, 20])
        data = b''.join(rng.choice(frag) for _ in range(n))
        if rng.random() < 0.3:
            data = enc(''.join(rng.choice(pool) for _ in range(rng.randint(1, 8)))
                       .replace('\ud800', 'x').replace('\udfff', 'y'))
            if data and rng.random() < 0.7:
                data = data[:rng.randint(0, len(data))] + bytes([rng.randrange(256)]) + data[rng.randint(0, len(data)):]
        try:
            data.decode('utf-8')
            k = -1
        except UnicodeDecodeError as e:
            k = e.start
        out.append('v %s %d' % (hx(data), k))
    # exhaustive: every 2-, 3- and 4-byte lead with every continuation prefix
    for b0 in range(0x80, 0x100):
        for b1 in range(0x80, 0x100):
            for tail in (b'', b'\x80', b'\x80\x80', b'\xbf\xbf'):
                data = bytes([b0, b1]) + tail
                try:
                    data.decode('utf-8')
                    k = -1
                except UnicodeDecodeError as e:
                    k = e.start
                out.append('v %s %d' % (hx(data), k))
    # repr: random text with quotes, escapes, every category, invalid bytes
    rpool = pool + list('\'"\\\n\r\t\x00\x7f\x80\x9f') + ['\u00ad', '\u2028', '\u2029', '\ue000',
                                                          '\U000f0000', '\U0010ffff', '\u0378', '\ufffe',
                                                          '\U0001f600', '\u200e', '\u0600']
    for _ in range(3000):
        text = ''.join(rng.choice(rpool) for _ in range(rng.choice([0, 1, 2, 4, 8, 16])))
        text = text.replace('\ud800', 'z').replace('\udfff', 'z')
        data = enc(text)
        if rng.random() < 0.4 and data:
            k = rng.randrange(len(data))
            data = data[:k] + bytes([rng.choice([0x80, 0xc0, 0xe4, 0xff, 0xed, 0xf4])]) + data[k:]
        out.append('e %s %s' % (hx(data), hx(repr(data.decode('utf-8', 'surrogateescape')).encode('utf-8'))))
    h = 0xCBF29CE484222325
    for c in range(0x110000):
        if 0xD800 <= c < 0xE000:
            continue
        for byte in repr(chr(c)).encode('utf-8') + b'\x00':
            h = ((h ^ byte) * 0x100000001B3) & MASK64
    out.append('reprall %d' % h)
    # universal newlines
    nl = ['\n', '\r', '\r\n', '\n\r', '\r\r\n']
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'f.txt')
        for _ in range(2000):
            parts = []
            for _ in range(rng.randint(0, 12)):
                parts.append(rng.choice(['', 'a', 'bc', '中', ' ', 'x y']) + rng.choice(nl + ['']))
            s = ''.join(parts)
            with open(path, 'wb') as f:
                f.write(enc(s))
            with open(path, 'r', encoding='utf-8', newline=None) as f:
                lines = f.readlines()
            out.append('n %s %d %s' % (hx(enc(s)), len(lines), ' '.join(hx(enc(x)) for x in lines)))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
