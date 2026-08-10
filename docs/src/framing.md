```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Framing

A [`framed_signal`](@ref) describes how a signal divides into (typically overlapping)
frames; the functions on this page materialise those frames and compute frame-wise
quantities.

- [`view_frame`](@ref) — create a [view](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-views) of one frame
- [`extract_frame`](@ref) — extract one frame as a new vector
- [`enframe`](@ref) — store all frames of a signal as the columns of a matrix
- [`number_signal_frames`](@ref) — number of frames that fit inside a signal
- [`frame_energy`](@ref) — energy of each frame

## `view_frame`

```@docs
view_frame
```

## `extract_frame`

```@docs
extract_frame
```

## `enframe`

```@docs
enframe
enframe!
```

## `number_signal_frames`

```@docs
number_signal_frames
```

## `frame_energy`

```@docs
frame_energy
```
