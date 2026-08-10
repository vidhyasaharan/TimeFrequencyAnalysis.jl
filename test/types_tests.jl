@testset "waveform" begin
    xm = rand(Float,3,10000)
    sig = @test_logs (:warn, r"first channel") waveform(xm,1000) #multichannel input takes the first channel with a warning
    @test typeof(sig.x) == Vector{Float}
    @test length(sig.x) == 10000
    @test sig.fs == 1000.0
    @test typeof(sig.fs) == Float

    xv = rand(Float,100)
    sig = waveform(xv,50)
    @test sig isa waveform{Float}
    @test sig.x === xv #vector input of the right type is stored without copying
end


@testset "framed_signal" begin
    frames = framed_signal(x,fs) #from a plain array
    frames_alt = framed_signal(signal) #from the waveform object
    @test frames == frames_alt
    @test frames isa framed_signal{Float}
    @test typeof(frames.signal.x) <: Vector{Float}
    @test length(frames.signal.x) > 0
    @test typeof(frames.frame_length) <: Int
    @test typeof(frames.frame_shift) <: Int
    @test typeof(frames.num_signal_frames) <: Int
    @test typeof(frames.num_frames) <: Int
    @test frames.frame_length > 0
    @test frames.frame_shift > 0
    @test frames.num_signal_frames > 0
    @test frames.num_frames > 0
    @test frames.num_frames >= frames.num_signal_frames
    @test frames.num_signal_frames >= (length(frames.signal.x)-frames.frame_length)/frames.frame_shift
end


@testset "spectrum" begin
    xs = rand(Float,220)
    c = rand(Complex{Float},111)
    f = convert.(Float,collect(1:length(c)))
    sp = spectrum(waveform(xs,220),c,f)
    @test sp.signal isa waveform{Float}
    @test typeof(sp.components) == Vector{Complex{Float}}
    @test typeof(sp.frqs) == Vector{Float}
    @test typeof(sp.title) <: AbstractString
    @test sp.title == "" #default title
end


@testset "timefreq" begin
    sw = waveform(rand(Float,1000),1000)
    comps = rand(Float,5,4)
    frqs = collect(1.0:5.0)
    t = collect(0.1:0.1:0.4)
    tf = timefreq(sw, comps, frqs, t, "test") #construct directly from a waveform (no frames)
    @test tf.signal == sw
    @test tf.frames === nothing
    @test tf.components == comps
    @test tf.frqs == frqs
    @test tf.time == t
    @test tf.title == "test"

    #Frame-based construction computes frame centre times
    frames = framed_signal(x,fs,0.02,0.01)
    comps2 = rand(Float, 7, frames.num_signal_frames)
    tf2 = timefreq(frames, comps2, collect(1.0:7.0))
    @test tf2.frames === frames
    @test length(tf2.time) == frames.num_signal_frames
    @test tf2.time[1] ≈ frames.frame_length/(2*fs) #first frame centre
    @test tf2.time[2] - tf2.time[1] ≈ frames.frame_shift/fs #spacing equals the frame shift
end
