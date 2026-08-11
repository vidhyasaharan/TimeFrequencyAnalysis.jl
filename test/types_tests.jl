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
