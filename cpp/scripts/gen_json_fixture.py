#!/usr/bin/env python
"""Writes cpp/tests/fixtures/cpython/json.txt, the oracle of the JSON parser
and writer (coolercpp::json::parse / dumps, hicx::json::parse).

Random JSON documents, random mutations of them (so that most are invalid) and
a list of edge cases go through simplejson.loads, the call cooler's
Cooler.info makes, and its dumps:

    ~/miniconda3/envs/hicx-ref-3.7.6/bin/python cpp/scripts/gen_json_fixture.py

Reference: Python 3.12.7, simplejson 3.19.3 (first line). dumps is called with
allow_nan=True: simplejson 3.19 raises ValueError for NaN and the infinities by
default, which the C++ writer does not reproduce (it writes the literals).
Records, text in hex:
  version <python> <simplejson>
  p <text> X                         loads raises JSONDecodeError
  p <text> <value> <dumps>           loads gives value; dumps(value) is dumps
Values are written as one token sequence:
  N  null    T  true    F  false
  i<digits>  int in the int64 range      u<digits>  int in 2^63..2^64-1
  b<bits>    int beyond 64 bits, as the nearest double (the C++ side keeps a
             double; its dumps differs from the integer's, so those records'
             dumps text is not compared)
  d<bits>    float, IEEE bits in hex
  s<hex>     str (UTF-8 with surrogatepass; '-' for empty)
  [ n v*n    list         { n (s<hex> v)*n   dict in insertion order
"""
import importlib.metadata
import math
import random
import struct
import sys

import simplejson


def hx(s):
    b = s.encode('utf-8', 'surrogatepass') if isinstance(s, str) else s
    return b.hex() if b else '-'


def canon(v):
    if v is None:
        return 'N'
    if v is True:
        return 'T'
    if v is False:
        return 'F'
    if isinstance(v, int):
        if -(1 << 63) <= v < (1 << 63):
            return 'i%d' % v
        if v < (1 << 64) and v >= 0:
            return 'u%d' % v
        try:
            f = float(v)
        except OverflowError:
            f = math.inf if v > 0 else -math.inf
        return 'b%016x' % struct.unpack('<Q', struct.pack('<d', f))[0]
    if isinstance(v, float):
        return 'd%016x' % struct.unpack('<Q', struct.pack('<d', v))[0]
    if isinstance(v, str):
        return 's' + hx(v)
    if isinstance(v, list):
        return '[ %d %s' % (len(v), ' '.join(canon(x) for x in v))
    if isinstance(v, dict):
        return '{ %d %s' % (len(v), ' '.join('s%s %s' % (hx(k), canon(x)) for k, x in v.items()))
    raise TypeError(v)


KEYS = ['a', 'b', 'bin-size', 'metadata', '', 'kéy', 'a b', 'x\ty', '中', '\U0001f600', 'q"uote', 'back\\slash']
STRINGS = ['', 'abc', 'x y', 'tab\there', 'nl\nhere', 'é', '中文', '\U0001f600', '\x00', '\x1f', '\x7f',
           ' ', '\ud800', '\udfff', '𐀀', 'q"uote', 'back\\slash', '/slash', '\b\f\r']


def rand_float(rng):
    r = rng.random()
    if r < 0.15:
        return rng.choice([0.0, -0.0, 1.0, 0.5, 1e16, 1e-5, 1e22, 5e-324, 1.7976931348623157e308, 123.456])
    if r < 0.5:
        return struct.unpack('<d', struct.pack('<Q', rng.getrandbits(64)))[0]
    return rng.uniform(-1000, 1000) * 10.0 ** rng.randint(-12, 12)


def rand_value(rng, depth):
    r = rng.random()
    if depth > 3 or r < 0.5:
        k = rng.randrange(8)
        if k == 0:
            return None
        if k == 1:
            return rng.random() < 0.5
        if k == 2:
            return rng.randint(-10**4, 10**4)
        if k == 3:
            return rng.choice([2**63 - 1, -2**63, 2**63, 2**64 - 1, 2**64, 10**30, -10**30, 0, -1])
        if k == 4:
            f = rand_float(rng)
            return f if math.isfinite(f) else 1.5
        if k in (5, 6):
            return rng.choice(STRINGS) if rng.random() < 0.6 else ''.join(rng.choice('abc d_') for _ in range(rng.randint(0, 8)))
        return rng.choice([float('nan'), float('inf'), float('-inf')])
    if r < 0.75:
        return [rand_value(rng, depth + 1) for _ in range(rng.randint(0, 4))]
    return {rng.choice(KEYS): rand_value(rng, depth + 1) for _ in range(rng.randint(0, 4))}


def render(rng, v):
    """JSON text with random spacing and escape choices."""
    ws = lambda: rng.choice(['', '', ' ', '\n', '\t', ' \r\n'])
    if isinstance(v, list):
        return '[' + ws() + (',' + ws()).join(render(rng, x) for x in v) + ws() + ']'
    if isinstance(v, dict):
        return '{' + ws() + (',' + ws()).join(rstr(rng, k) + ws() + ':' + ws() + render(rng, x)
                                              for k, x in v.items()) + ws() + '}'
    if isinstance(v, str):
        return rstr(rng, v)
    if v is None:
        return 'null'
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, float):
        if math.isnan(v):
            return 'NaN'
        if math.isinf(v):
            return 'Infinity' if v > 0 else '-Infinity'
        return rng.choice([repr(v), repr(v).replace('e', 'E'), simplejson.dumps(v)])
    return str(v)


def rstr(rng, s):
    mode = rng.random()
    if mode < 0.4:
        return simplejson.dumps(s)
    out = ['"']
    for ch in s:
        c = ord(ch)
        if ch == '"':
            out.append('\\"')
        elif ch == '\\':
            out.append('\\\\')
        elif c < 0x20 or rng.random() < 0.3:
            if c < 0x10000:
                out.append('\\u%04x' % c if rng.random() < 0.5 else '\\u%04X' % c)
            else:
                n = c - 0x10000
                out.append('\\u%04x\\u%04x' % (0xD800 | (n >> 10), 0xDC00 | (n & 0x3FF)))
        elif ch == '/' and rng.random() < 0.5:
            out.append('\\/')
        else:
            out.append(ch)
    out.append('"')
    return ''.join(out)


MUT = ['{', '}', '[', ']', ',', ':', '"', '\\', ' ', '\n', '0', '1', '-', '+', '.', 'e', 'E', 'n', 'u', 'l', 't',
       'NaN', 'Infinity', '-Infinity', 'null', 'true', '\\u', '\\ud800', '\\udc00', '\\x', '\x01', '\x7f', 'é',
       '﻿', "'", '/', '//', '/*', '1e', '1.', '.5', '01', '-0', '--1', 'nan', 'inf', 'True', 'None', '\t']


def mutate(rng, text):
    for _ in range(rng.randint(1, 3)):
        i = rng.randint(0, len(text))
        r = rng.random()
        if r < 0.4:
            text = text[:i] + rng.choice(MUT) + text[i:]
        elif r < 0.7 and text:
            j = min(len(text), i + rng.randint(1, 3))
            text = text[:i] + text[j:]
        elif text:
            j = min(len(text), i + 1)
            text = text[:i] + rng.choice(MUT) + text[j:]
    return text


def record(text):
    try:
        v = simplejson.loads(text)
    except simplejson.JSONDecodeError:
        return 'p %s X' % hx(text)
    except RecursionError:
        return None
    return 'p %s %s %s' % (hx(text), canon(v), hx(simplejson.dumps(v, allow_nan=True)))


def main():
    rng = random.Random(1033)
    out = ['version %s %s' % (sys.version.split()[0], importlib.metadata.version('simplejson'))]
    edge = ['', ' ', 'null', ' null ', 'nul', 'nulll', 'true', 'True', 'false', 'NaN', '-NaN', 'Infinity', '-Infinity',
            '+Infinity', 'infinity', '-', '+1', '1', '-1', '-0', '0', '00', '01', '-01', '1.', '1.5', '.5', '-.5',
            '1e5', '1E5', '1e+5', '1e-5', '1e', '1e+', '1.5e', '1.e5', '0.0', '-0.0', '0e0', '1e999', '-1e999', '1e-999',
            '9223372036854775807', '9223372036854775808', '-9223372036854775808', '-9223372036854775809',
            '18446744073709551615', '18446744073709551616', '123456789012345678901234567890',
            '""', '"', '"a', '"a"', '"a" "b"', '"\\"', '"\\x"', '"\\u12"', '"\\u12G4"', '"\\ud800"', '"\\udc00"',
            '"\\ud800\\udc00"', '"\\ud800\\u0041"', '"\\udc00\\ud800"', '"\\ud83d\\ude00"', '"\\uD83D\\uDE00"',
            '"a\nb"', '"a\tb"', '"a\x7fb"', '"a\x00b"', '"\\/"', '"\\b\\f\\n\\r\\t"', "'a'", '[', ']', '[]', '[ ]', '[1]',
            '[1,]', '[,1]', '[1 2]', '[1,2]', '[[]]', '[[],[]]', '{', '}', '{}', '{ }', '{"a"}', '{"a":}', '{"a":1}',
            '{"a":1,}', '{,}', '{"a":1 "b":2}', '{"a":1,"a":2}', '{"a":1,"b":2,"a":3}', '{1:2}', '{"a" : 1 }',
            '{"a":{"b":{"c":[1,2,{"d":null}]}}}', '[1] x', '[1]x', '﻿1', ' 1', '1 ', '1 ',
            '/* c */ 1', '// c\n1', '[1,2,3', '{"a":1', 'nan', 'inf', 'None', 'tru', 'fals', '[NaN,Infinity,-Infinity]',
            '{"NaN":NaN}', '" "', '"\ud800"', '["a","b",]', '2018-11-08T20:03:12.737522', 'cooler-0.7.11',
            'HiCMatrix-16.dev', 'unknown', '1_0', '0x10', '1e1_0', '\t\n\r 1 \t\n\r', '\x0b1', '\x0c1', '1\x00']
    for t in edge:
        r = record(t)
        if r:
            out.append(r)
    for depth in (10, 100, 400):
        r = record('[' * depth + ']' * depth)
        if r:
            out.append(r)
    for _ in range(4000):
        text = render(rng, rand_value(rng, 0))
        out.append(record(text))
        out.append(record(mutate(rng, text)))
    print('\n'.join(x for x in out if x))


if __name__ == '__main__':
    main()
