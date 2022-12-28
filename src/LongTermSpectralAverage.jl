module LongTermSpectralAverage

using Dates
using DSP
using GLMakie
using ProgressMeter
using SignalAnalysis: SampledSignal
using Statistics
using WAV

export ltsa, ltsa_plot

function ltsa(xs::AbstractVector{String}; n::Int = 1024, noverlap::Int = 512)
    numtimes = length(xs)
    _, fs = wavread(first(xs); subrange=1)
    freqs = DSP.rfftfreq(n, fs)
    P = zeros(length(freqs), numtimes)
    prog = Progress(numtimes)
    Threads.@threads for i ∈ 1:numtimes
        try 
            x, _ = wavread(xs[i])
            size(x, 2) > 1 && (x = vec(mean(x; dims = 2))) # array of sensors
            spec = spectrogram(x, n, noverlap; fs = fs)
            P[:,i] = mean(pow2db.(power(spec)); dims = 2) |> vec
        catch
            @warn "$(wavpath) is not readable"
        end
        next!(prog)
    end
    P, freqs
end

function ltsa(xs::AbstractVector{<:VecOrMat{T}}; fs::Real = 1, n::Int = 1024, noverlap::Int = 512) where {T}
    numtimes = length(xs)
    freqs = DSP.rfftfreq(n, fs)
    P = zeros(length(freqs), numtimes)
    prog = Progress(numtimes)
    Threads.@threads for i ∈ 1:numtimes
        x = xs[i]
        size(x, 2) > 1 && (x = vec(mean(x; dims = 2))) # array of sensors
        spec = spectrogram(x, n, noverlap; fs = fs)
        P[:,i] = mean(pow2db.(power(spec)); dims = 2) |> vec
        next!(prog)
    end
    P, freqs
end
ltsa(ss::AbstractVector{SampledSignal}; n::Int = 1024, noverlap::Int = 512) = ltsa(samples(ss); 
                                                                                   fs = framerate(ss), 
                                                                                   n = n, 
                                                                                   noverlap = noverlap)

function ltsa_plot(xs::AbstractVector{T}, 
                   starttimes::AbstractVector{DateTime}; 
                   kwargs...) where {T}
    P, freqs = ltsa(xs; kwargs...)
    ltsa_plot(P, starttimes, freqs)
end

function ltsa_plot(P::AbstractMatrix{T}, starttimes::AbstractVector{DateTime}, freqs::AbstractVector) where {T}
    numtimes = size(P, 2)
    unixtimes = datetime2unix.(starttimes) 

    fig = Figure(resolution = (1200, 800))
    ax = Axis(fig[1,1]; ylabel = "Frequency (Hz)", xticklabelrotation = π / 4)
    heatmap!(ax, unixtimes, freqs, collect(P'))
    step = numtimes ÷ 10
    ax.xticks = (unixtimes[1:step:end], string.(starttimes[1:step:end]))
    fig
end

end # module LongTermSpectralAverage
