#Construction and field integrity of the container types (signal, framed_signal,
#spectrum, timefreq): inputs are stored/converted as documented, defaults are applied,
#and derived quantities (frame counts, frame centre times) are computed correctly.

@testset "signal" begin
    #Multichannel (matrix) input is not supported as such: the constructor must keep only
    #the first channel and warn, storing samples as a vector and fs in the signal precision
    xm = rand(Float64,3,10000)
    s = @test_logs (:warn, r"first channel") signal(xm,1000) #multichannel input takes the first channel with a warning
    @test typeof(s.x) == Vector{Float64}
    @test length(s.x) == 10000
    @test s.fs == 1000.0
    @test typeof(s.fs) == Float64

    #A vector that is already of the right element type must be stored without copying
    xv = rand(Float64,100)
    s = signal(xv,50)
    @test s isa signal{Float64}
    @test s.x === xv #vector input of the right type is stored without copying
end


@testset "framed_signal" begin
    #The array+fs and signal constructors must produce identical framing geometry
    #(default frame duration and shift), with positive integer geometry fields
    frames = framed_signal(x,fs) #from a plain array
    frames_alt = framed_signal(sig) #from the signal object
    @test frames == frames_alt
    @test frames isa framed_signal{Float64}
    @test typeof(frames.signal.x) <: Vector{Float64}
    @test length(frames.signal.x) > 0
    @test typeof(frames.frame_length) <: Int
    @test typeof(frames.frame_shift) <: Int
    @test typeof(frames.num_signal_frames) <: Int
    @test typeof(frames.num_frames) <: Int
    @test frames.frame_length > 0
    @test frames.frame_shift > 0
    @test frames.num_signal_frames > 0
    @test frames.num_frames > 0
    #num_frames includes trailing partial (padded) frames, so it can only exceed the
    #number of frames fully contained in the signal, which itself cannot be fewer than
    #the sliding-window count (length - frame_length)/frame_shift
    @test frames.num_frames >= frames.num_signal_frames
    @test frames.num_signal_frames >= (length(frames.signal.x)-frames.frame_length)/frames.frame_shift
end


@testset "spectrum" begin
    #The spectrum container must store the analysed signal, complex components and the
    #frequency grid (converted to the signal precision), with an empty default title
    xs = rand(Float64,220)
    c = rand(Complex{Float64},111)
    f = convert.(Float64,collect(1:length(c)))
    sp = spectrum(signal(xs,220),c,f)
    @test sp.signal isa signal{Float64}
    @test typeof(sp.components) == Vector{Complex{Float64}}
    @test typeof(sp.frqs) == Vector{Float64}
    @test typeof(sp.title) <: AbstractString
    @test sp.title == "" #default title
end


@testset "timefreq" begin
    #Direct construction from a signal stores all fields verbatim and records that no
    #framing was involved (frames === nothing)
    sw = signal(rand(Float64,1000),1000)
    comps = rand(Float64,5,4)
    frqs = collect(1.0:5.0)
    t = collect(0.1:0.1:0.4)
    tf = timefreq(sw, comps, frqs, t, "test") #construct directly from a signal (no frames)
    @test tf.signal == sw
    @test tf.frames === nothing
    @test tf.components == comps
    @test tf.frqs == frqs
    @test tf.time == t
    @test tf.title == "test"

    #Frame-based construction computes frame centre times: one time per full frame,
    #starting at half a frame length and advancing by the frame shift (in seconds)
    frames = framed_signal(x,fs,0.02,0.01)
    comps2 = rand(Float64, 7, frames.num_signal_frames)
    tf2 = timefreq(frames, comps2, collect(1.0:7.0))
    @test tf2.frames === frames
    @test length(tf2.time) == frames.num_signal_frames
    @test tf2.time[1] ≈ frames.frame_length/(2*fs) #first frame centre
    @test tf2.time[2] - tf2.time[1] ≈ frames.frame_shift/fs #spacing equals the frame shift
end


#modfreq holds a modulation spectrum: rows are carrier channels, columns modulation rates. It is
#deliberately not anchored to a signal (see the docstring), so its precision comes from the
#components rather than from a signal field, and its invariants are checked at construction.
@testset "modfreq" begin
    ncar, nmod = 4, 7
    C = reshape(collect(1.0:(ncar*nmod)), ncar, nmod)
    fcs = [80.0, 160.0, 320.0, 640.0]
    mfr = collect(range(0.0, 15.0; length = nmod))
    bws = 2 .* erb.(fcs)

    mf = modfreq(C, fcs, mfr, bws, 312.5, 10.0, 11, "Modulation Spectrum")
    @test mf isa modfreq{Float64,Float64}
    @test size(mf.components) == (length(mf.fcs), length(mf.mod_frqs))
    @test mf.fcs == fcs && mf.mod_frqs == mfr && mf.bws == bws
    @test mf.fs_env == 312.5 && mf.frame_dur == 10.0 && mf.nframes == 11
    @test mf.title == "Modulation Spectrum"

    #The short form leaves the title empty, as timefreq's does
    @test modfreq(C, fcs, mfr, bws, 312.5, 10.0, 11).title === nothing

    #Every field is converted to the precision carried by the components, whatever it was given as
    mf32 = modfreq(Float32.(C), fcs, mfr, bws, 312.5, 10.0, 11)
    @test mf32 isa modfreq{Float32,Float32}
    @test eltype(mf32.fcs) == Float32
    @test eltype(mf32.mod_frqs) == Float32
    @test eltype(mf32.bws) == Float32
    @test typeof(mf32.fs_env) == Float32
    @test typeof(mf32.frame_dur) == Float32
    @test mf32.nframes isa Int

    #A complex component matrix keeps T real: real(S) is T for both real and complex components
    mfc = modfreq(ComplexF64.(C), fcs, mfr, bws, 312.5, 10.0, 11)
    @test mfc isa modfreq{Float64,ComplexF64}
    @test eltype(mfc.fcs) == Float64

    #Axes must match the component matrix, and the orientation is not interchangeable: a
    #transposed matrix is rejected rather than silently accepted
    @test_throws ErrorException modfreq(permutedims(C), fcs, mfr, bws, 312.5, 10.0, 11)
    @test_throws ErrorException modfreq(C, fcs[1:3], mfr, bws, 312.5, 10.0, 11)
    @test_throws ErrorException modfreq(C, fcs, mfr[1:5], bws, 312.5, 10.0, 11)
    @test_throws ErrorException modfreq(C, fcs, mfr, bws[1:2], 312.5, 10.0, 11)   #one bw per carrier
    @test_throws ErrorException modfreq(C, fcs, mfr, bws, 0, 10.0, 11)            #fs_env positive
    @test_throws ErrorException modfreq(C, fcs, mfr, bws, 312.5, -1, 11)          #frame_dur positive
    @test_throws ErrorException modfreq(C, fcs, mfr, bws, 312.5, 10.0, 0)         #at least one frame
end
