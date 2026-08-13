# Reference-value generator for TimeFrequencyAnalysis gammatone_filter tests.
#
# Route B (bandwidth at given edge attenuation): runs the ACTUAL pyfilterbank
# implementation (pyfilterbank_gammatone.py fetched verbatim from
# https://raw.githubusercontent.com/SiggiGue/pyfilterbank/master/pyfilterbank/gammatone.py,
# fetched 2026-08-13). Note pyfilterbank's attenuation_half_bandwidth_db is signed
# (default -3); TFA's attenuation_db is its magnitude (3 = 3 dB below peak).
#
# Route A (ERB-matching, Hohmann Eqs. 13-17): pyfilterbank does not implement this route
# (its a_gamma branch is commented out), so values are generated from an independent
# line-by-line transcription of the equations below. AMT's hohmann2002filter implements
# the same equations.
import importlib.util
import math

spec = importlib.util.spec_from_file_location("pfb", "pyfilterbank_gammatone.py")
pfb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pfb)

print("== Route B: pyfilterbank design_filter (attenuation_half_bandwidth_db = -3) ==")
cases_b = [
    dict(sample_rate=16000, order=4, centerfrequency=1000.0, band_width=200.0),
    dict(sample_rate=16000, order=4, centerfrequency=1000.0, band_width=None, band_width_factor=1.0),
    dict(sample_rate=44100, order=4, centerfrequency=4000.0, band_width=456.5),
    dict(sample_rate=16000, order=2, centerfrequency=500.0, band_width=100.0),
]
for c in cases_b:
    b, a = pfb.design_filter(attenuation_half_bandwidth_db=-3, **c)
    coef = -a[1]
    print(f"{c}")
    print(f"  coef = {coef.real!r} + {coef.imag!r}im")
    print(f"  norm = {b[0].real!r}")

print()
print("== Route A: independent transcription of Hohmann (2002) Eqs. 13-17 ==")

def route_a(fs, fc, order=4, bw_factor=1.0):
    L, Q = 24.7, 9.265
    erb_aud = (L + fc / Q) * bw_factor                                   # Eq. 13
    a_g = (math.pi * math.factorial(2 * order - 2) * 2.0 ** (-(2 * order - 2))
           / math.factorial(order - 1) ** 2)                             # Eq. 14
    b_damp = erb_aud / a_g                                               # Eq. 15
    lam = math.exp(-2 * math.pi * b_damp / fs)                           # Eq. 16 (sampled decay)
    beta = 2 * math.pi * fc / fs                                         # Eq. 10
    coef = complex(lam * math.cos(beta), lam * math.sin(beta))           # Eq. 1
    k = 2 * (1 - lam) ** order                                           # normalisation, Sec. 2.2
    return coef, k, b_damp, lam

cases_a = [
    dict(fs=16000, fc=1000.0, order=4, bw_factor=1.0),
    dict(fs=16000, fc=100.0, order=4, bw_factor=1.0),
    dict(fs=16000, fc=4000.0, order=4, bw_factor=1.0),
    dict(fs=16000, fc=1000.0, order=4, bw_factor=0.5),
    dict(fs=44100, fc=2000.0, order=3, bw_factor=2.0),
]
for c in cases_a:
    coef, k, b_damp, lam = route_a(**c)
    print(f"{c}")
    print(f"  coef = {coef.real!r} + {coef.imag!r}im")
    print(f"  norm = {k!r}")
    print(f"  b = {b_damp!r}  lambda = {lam!r}")
