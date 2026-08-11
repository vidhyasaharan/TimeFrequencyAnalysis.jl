#Window generation: each supported window type must match its closed-form definition
#w(n) = α₀ - (1-α₀)cos(2πn/(len-1)) evaluated sample by sample, the default type must be
#Hann, and unrecognised types must warn and fall back to Hann rather than error.

@testset "window" begin
    flen = 5

    #Rectangular window: all ones
    win = window(flen;wtype="rect")
    @test typeof(win)<:Vector{Float64}
    @test length(win) == flen
    for i=1:flen
        @test win[i] == 1.0
    end

    #Hann window: α₀ = 0.5 in the closed form above
    win = window(flen;wtype="hanning")
    @test typeof(win)<:Vector{Float64}
    @test length(win) == flen
    α₀ = 0.5
    α₁ = 1 - α₀
    N = flen-1
    wfun(n) = α₀ - α₁*cos(2π*n/N)
    for i=1:flen
        @test win[i] ≈ wfun(i-1)
    end

    #Hamming window: α₀ = 0.54
    win = window(flen;wtype="hamming")
    @test typeof(win)<:Vector{Float64}
    @test length(win) == flen
    α₀ = 0.54
    α₁ = 1 - α₀
    wfun1(n) = α₀ - α₁*cos(2π*n/N)
    for i=1:flen
        @test win[i] ≈ wfun1(i-1)
    end

    #Default window type is Hann
    win1 = window(flen)
    win2 = window(flen;wtype="hanning")
    @test win1 == win2

    #Unrecognised window types warn and fall back to Hann
    win3 = @test_logs (:warn, r"not recognised") window(flen; wtype="bogus")
    @test win3 == win2
end
