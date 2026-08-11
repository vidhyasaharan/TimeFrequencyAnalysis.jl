#Correlation analyses on signals with known structure: cross correlation must locate an
#embedded segment, autocorrelation of noise must peak at lag zero with the mean-square
#value, and the short-time functions (acf/nacf) are checked against hand-computable
#fixtures (constant and perfectly periodic signals) across all their call forms.

@testset "xcorr" begin
    #h is a segment cut out of x at hstart, so the sliding dot product must peak exactly
    #where the segment came from (offset by the zero-index position z), and the peak
    #value must equal the energy of the segment (its dot product with itself)
    len = 1000
    hlen = 100
    hstart = 101
    x = randn(Float64,len)
    h = x[hstart:hstart+hlen-1]
    z = 10
    y = xcorr(x,h,z)
    y1 = Vector{Float64}(undef,length(x))
    xcorr!(y1, x, h, z)
    @test typeof(y) == Vector{Float64}
    @test length(y) == length(x) #one correlation value per lag, as many lags as samples
    @test y1 == y #in-place and allocating forms agree
    @test argmax(y) == hstart + z - 1 #peak where the segment was cut out
    @test maximum(y) ≈ sum(abs2,h) #peak value is the segment energy
end

@testset "acorr" begin
    #White noise: the autocorrelation must peak at lag zero (index 1), where the
    #window-length normalisation makes the value the mean square of the signal
    len = 1000
    p = 100
    x = randn(Float64, len)
    rxx = Vector{Float64}(undef,p)
    acorr!(rxx,x)
    rx = acorr(x,p)
    @test typeof(rx) == Vector{Float64}
    @test length(rx) == p #one value per lag 0..p-1
    @test rxx == rx #in-place and allocating forms agree
    @test argmax(rx) == 1
    @test maximum(rx) ≈ sum(abs2,x)/len
end

@testset "acf" begin
    #On a constant (all-ones) signal every product is 1, so acf reduces to the number of
    #overlapping samples: win_size - k for lag k — an exact hand-computable expectation
    afs = 8000.0
    s = waveform(ones(Float64, 800), afs)
    win_size = 100
    @test acf(s, 1, 0; win_size) == win_size
    @test acf(s, 1, 10; win_size) == win_size - 10
    #The vector-lag form must map the scalar form over the lags
    ks = [0, 5, 10]
    @test acf(s, 1, ks; win_size) == [acf(s, 1, k; win_size) for k in ks]
    #The time-based forms must agree with the index-based forms after converting
    #seconds to samples at fs
    t = 0.05; lag = 0.005; win_dur = 0.02
    i = time2nsamples(t, afs); k = time2nsamples(lag, afs); ws = time2nsamples(win_dur, afs)
    @test acf(s, t, lag; win_dur) == acf(s, i, k; win_size = ws)
    @test acf(s, t, [lag, 2lag]; win_dur) == [acf(s, t, lag; win_dur), acf(s, t, 2lag; win_dur)]
end

@testset "nacf" begin
    #A pure cosine with an exact 40-sample period: a window correlated with itself one
    #period later is identical (nacf = 1), and half a period later is negated (nacf = -1)
    afs = 8000.0
    period = 40 #a 200 Hz tone at 8 kHz
    s = waveform(cos.(2pi/period .* (0:3999)), afs)
    win_size = 4period
    @test nacf(s, 1, period; win_size) ≈ 1
    @test nacf(s, 1, period ÷ 2; win_size) ≈ -1
    #The norm scaling (Cauchy-Schwarz) bounds the values to [-1, 1] for any signal
    sr = waveform(randn(Float64, 2000), afs)
    vals = nacf(sr, 1, collect(1:50); win_size = 200)
    @test all(abs.(vals) .<= 1)
    #A positive nconst inflates the denominator, damping the value towards zero
    @test 0 < nacf(s, 1, period; win_size, nconst = 10.0) < nacf(s, 1, period; win_size)
    #The time-based form must agree with the index-based form
    @test nacf(s, 100/afs, period/afs; win_dur = win_size/afs) ≈ nacf(s, 100, period; win_size)
    #Keyword-only batch form: one row per lag in min_lag:max_lag, one column per window
    #position, each entry matching the corresponding single-point call
    cf = nacf(s; min_lag = 1, max_lag = 50, win_size = 200, win_shift = 80)
    @test size(cf,1) == 50
    @test size(cf,2) == number_signal_frames(s, 200 + 50, 80)
    @test cf[period, 1] ≈ nacf(s, 1, period; win_size = 200)
    @test eltype(cf) == Float64
end
