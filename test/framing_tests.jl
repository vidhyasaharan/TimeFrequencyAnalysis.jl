#Framing: splitting signals into (overlapping) frames. Verifies frame extraction against
#hand-computed slices, agreement between the sample-count/duration/waveform call forms,
#frame counting, per-frame energy, and the padding helpers used for trailing frames.

@testset "enframe" begin
    #100 samples with frame length 10 and shift 5 give (100-10)/5 + 1 = 19 full frames;
    #column i of the output must be exactly the i-th length-10 slice of the input
    xs = collect(1.0:100.0)
    y = enframe(xs, 10, 5)
    @test size(y) == (10, number_signal_frames(xs, 10, 5))
    @test size(y,2) == 19
    for i ∈ axes(y,2)
        @test y[:,i] == xs[(i-1)*5 .+ (1:10)] #column i is the i-th frame
    end

    #Duration-based and waveform forms give the same frames (0.1 s / 0.05 s at 100 Hz is
    #the same 10-sample frame and 5-sample shift)
    y2 = enframe(xs, 100.0, 0.1, 0.05)
    @test y2 == y
    y3 = enframe(waveform(xs,100.0), 0.1, 0.05)
    @test y3 == y

    #enframe! derives the frame shift from the preallocated matrix size
    y4 = Matrix{Float64}(undef, 10, 19)
    enframe!(y4, xs)
    @test y4 == y
end


@testset "view/extract_frame" begin
    #view_frame and extract_frame must return the same samples as direct slicing;
    #extract_frame materialises them as a Vector (a copy rather than a view)
    frames = framed_signal(x,fs)
    @test view_frame(frames,1) == x[1:frames.frame_length]
    @test extract_frame(frames,1) == view_frame(frames,1)
    @test typeof(extract_frame(frames,1)) == Vector{Float64}
    @test extract_frame(frames,2) == x[frames.frame_shift .+ (1:frames.frame_length)]

    #enframe agrees with view_frame for every full frame
    y = enframe(x, frames.frame_length, frames.frame_shift)
    for i ∈ 1:frames.num_signal_frames
        @test view_frame(frames,i) == y[:,i]
    end

    #Frames beyond num_signal_frames are zero padded to full length but must still hold
    #the trailing real samples (not be entirely zeros)
    lastframe = view_frame(frames, frames.num_frames)
    @test length(lastframe) == frames.frame_length
    @test maximum(abs.(lastframe)) > 0
end


@testset "number_signal_frames" begin
    #All call forms (array or waveform, sample counts or durations) must agree on the
    #number of frames fully contained in the signal: 19 for 100 samples with frame 10 /
    #shift 5, and a single frame when the frame spans the whole signal
    xs = zeros(Float64, 100)
    @test number_signal_frames(xs, 10, 5) == 19
    @test number_signal_frames(xs, 100, 1) == 1
    @test number_signal_frames(xs, 100.0, 0.1, 0.05) == 19
    @test number_signal_frames(waveform(xs,100.0), 10, 5) == 19
    @test number_signal_frames(waveform(xs,100.0), 0.1, 0.05) == 19
end


@testset "frame_energy" begin
    #10 impulses, one every 10 samples; 20-sample non-overlapping frames hold 2 each, so
    #every frame energy (sum of squared samples) is exactly 2, and normalising by the
    #maximum gives all ones
    it = impulse_train(10, 100)
    frames = framed_signal(it, 100.0, 0.2, 0.2)
    energy = frame_energy(frames)
    @test length(energy) == frames.num_frames
    @test all(energy .== 2.0)
    @test frame_energy(frames; normalised = true) == ones(Float64, frames.num_frames)
end


@testset "padding" begin
    #zero_pad surrounds the signal with n zeros on each side; symmetric_pad reflects the
    #signal about its end points (excluding the end sample itself); unpad_vector strips
    #the padding back off, recovering the original in both cases
    xs = collect(1.0:5.0)
    zp = TimeFrequencyAnalysis.zero_pad(xs, 2)
    @test zp == [0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0]
    @test eltype(zp) == Float64
    sp = TimeFrequencyAnalysis.symmetric_pad(xs, 2)
    @test sp == [3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0]
    @test TimeFrequencyAnalysis.unpad_vector(zp, 2) == xs
    @test TimeFrequencyAnalysis.unpad_vector(sp, 2) == xs
end
