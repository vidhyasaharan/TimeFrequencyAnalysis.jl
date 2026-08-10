@testset "enframe" begin
    xs = collect(1.0:100.0)
    y = enframe(xs, 10, 5)
    @test size(y) == (10, number_signal_frames(xs, 10, 5))
    @test size(y,2) == 19
    for i ∈ axes(y,2)
        @test y[:,i] == xs[(i-1)*5 .+ (1:10)] #column i is the i-th frame
    end

    #Duration-based and waveform forms give the same frames
    y2 = enframe(xs, 100.0, 0.1, 0.05)
    @test y2 == y
    y3 = enframe(waveform(xs,100.0), 0.1, 0.05)
    @test y3 == y

    #enframe! derives the frame shift from the matrix size
    y4 = Matrix{Float}(undef, 10, 19)
    enframe!(y4, xs)
    @test y4 == y
end


@testset "view/extract_frame" begin
    frames = framed_signal(x,fs)
    @test view_frame(frames,1) == x[1:frames.frame_length]
    @test extract_frame(frames,1) == view_frame(frames,1)
    @test typeof(extract_frame(frames,1)) == Vector{Float}
    @test extract_frame(frames,2) == x[frames.frame_shift .+ (1:frames.frame_length)]

    #enframe agrees with view_frame for every full frame
    y = enframe(x, frames.frame_length, frames.frame_shift)
    for i ∈ 1:frames.num_signal_frames
        @test view_frame(frames,i) == y[:,i]
    end

    #Frames beyond num_signal_frames are zero padded to full length
    lastframe = view_frame(frames, frames.num_frames)
    @test length(lastframe) == frames.frame_length
    @test maximum(abs.(lastframe)) > 0
end


@testset "number_signal_frames" begin
    xs = zeros(Float, 100)
    @test number_signal_frames(xs, 10, 5) == 19
    @test number_signal_frames(xs, 100, 1) == 1
    @test number_signal_frames(xs, 100.0, 0.1, 0.05) == 19
    @test number_signal_frames(waveform(xs,100.0), 10, 5) == 19
    @test number_signal_frames(waveform(xs,100.0), 0.1, 0.05) == 19
end


@testset "frame_energy" begin
    #10 impulses, one every 10 samples; 20-sample non-overlapping frames hold 2 each
    it = impulse_train(10, 100)
    frames = framed_signal(it, 100.0, 0.2, 0.2)
    energy = frame_energy(frames)
    @test length(energy) == frames.num_frames
    @test all(energy .== 2.0)
    @test frame_energy(frames; normalised = true) == ones(Float, frames.num_frames)
end


@testset "padding" begin
    xs = collect(1.0:5.0)
    zp = TimeFrequencyAnalysis.zero_pad(xs, 2)
    @test zp == [0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0]
    @test eltype(zp) == Float
    sp = TimeFrequencyAnalysis.symmetric_pad(xs, 2)
    @test sp == [3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0]
    @test TimeFrequencyAnalysis.unpad_vector(zp, 2) == xs
    @test TimeFrequencyAnalysis.unpad_vector(sp, 2) == xs
end
