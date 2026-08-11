using TimeFrequencyAnalysis
using Test
using LinearAlgebra
using Random

#Deterministic broadband test signal: two tones in noise, 1 second at 8 kHz
Random.seed!(2026)
const fs = 8000.0
const x = cos.(2π*440.0.*(1:8000)./fs) .+ 0.5.*cos.(2π*1000.0.*(1:8000)./fs) .+ 0.3.*randn(8000)
const signal = waveform(x,fs)
