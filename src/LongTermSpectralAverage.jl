module LongTermSpectralAverage

using Dates
using DSP
using GLMakie
using ProgressMeter
using SignalAnalysis
using SignalAnalysis: SampledSignal
using Statistics
using ThreadsX
using WAV

export ltsa, ltsa_plot

## ltsa

function ltsa(xs::AbstractVector{String}; n::Int = 1024, noverlap::Int = 512)
    fss = getindex.(wavread.(xs; subrange = 1), 2)
    fs = first(fss)
    any(fss .!= fs) && throw(ArgumentError("WAV files have different sample rates"))
    numtimes = length(xs)
    frequencies = DSP.rfftfreq(n, fs)
    prog = Progress(numtimes)
    P = ThreadsX.map(1:numtimes) do i
        next!(prog)
        x, _ = try 
            wavread(xs[i])
        catch e
            throw("$(e): $(xs[i]) is not readable")
        end
            size(x, 2) > 1 && (x = vec(mean(x; dims = 2))) # array of sensors
        removedc!(x)
            spec = spectrogram(x, n, noverlap; fs = fs)
        mean(pow2db.(power(spec)); dims = 2) |> vec
    end |> a -> hcat(a...)
    P, frequencies
end

function ltsa(xs::AbstractVector{<:VecOrMat{T}}; fs::Real = 1, n::Int = 1024, noverlap::Int = 512) where {T}
    numtimes = length(xs)
    frequencies = DSP.rfftfreq(n, fs)
    prog = Progress(numtimes)
    P = ThreadsX.map(1:numtimes) do i
        next!(prog)
        x = xs[i]
        size(x, 2) > 1 && (x = vec(mean(x; dims = 2))) # array of sensors
        removedc!(x)
        spec = spectrogram(x, n, noverlap; fs = fs)
        mean(pow2db.(power(spec)); dims = 2) |> vec
    end |> a -> hcat(a...)
    P, frequencies
    end

ltsa(ss::AbstractVector{SampledSignal}; kwargs...) = ltsa(samples(ss); fs = framerate(ss), kwargs...) 

## ltsa_plot 

function ltsa_plot(xs::AbstractVector{T}, 
                   starttimes::AbstractVector{DateTime}; 
                   kwargs...) where {T}
    P, frequencies = ltsa(xs; kwargs...)
    ltsa_plot(P, starttimes, frequencies)
end

function ltsa_plot(P::AbstractMatrix{T}, starttimes::AbstractVector{DateTime}, frequencies::AbstractVector) where {T}
    numtimes = size(P, 2)
    unixtimes = datetime2unix.(starttimes) 

    fig = Figure(resolution = (1200, 800))
    ax = Axis(fig[1,1]; ylabel = "Frequency (Hz)", xticklabelrotation = π / 4)
    heatmap!(ax, unixtimes, frequencies, collect(P'))
    δ = numtimes > 20 ? numtimes ÷ 10 : 2
    ax.xticks = (unixtimes[1:δ:end], string.(starttimes[1:δ:end]))
    fig
end

end # module LongTermSpectralAverage
