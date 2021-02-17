export compute_fft, compute_cc, corrplot

@deprecate compute_fft() rfft()
@deprecate compute_cc() correlate()

function corrplot(C::CorrData)
  @warn "corrplot has been deprecated. Please use: `using Plots;plot(C)` to plot CorrData."
end

"""

  compute_fft(R)

Computes windowed rfft of ambient noise data. Returns FFTData structure.

# Arguments
- `R::RawData`: RawData structure
"""
function compute_fft(R::RawData)
    FFT = rfft(R.x,1)
    return FFTData(R.name, R.id,R.loc, R.fs, R.gain, R.freqmin, R.freqmax,
                   R.cc_len, R.cc_step, R.whitened, R.time_norm, R.resp,
                   R.misc, R.notes, R.t, FFT)
end

"""
    compute_cc(FFT1::FFTData, FFT2::FFTData, maxlag::Float64;
               corr_type::String="cross-correlation")

Cross-correlate ambient noise data in the frequency domain.

Cross-correlation can be done using one of three options:

- Cross-correlation: ``C_{AB}(ω) = u_A(ω) u^∗_B(ω)``
- Coherence: ``C_{AB}(ω) = \frac{u_A(ω) u^∗_B(ω)}{∣ u_A(ω) ∣  ∣ u_B(ω) ∣}``
- Deconvolution: ``C_{AB}(ω) = \frac{u_A(ω) u^∗_B(ω)}{∣ u_B(ω) ∣^2}``

# Arguments
- `FFT1::FFTData`: FFTData object of fft'd ambient noise data.
- `FFT2::FFTData`: FFTData object of fft'd ambient noise data.
- `maxlag::Float64`: Maximum lag time (in seconds) in cross-correlation to save,
                     e.g. `maxlag = 20.` will save lag times = -20.:20. s.
- `corr_type::String`: Type of correlation: `cross-correlation`, `coherence` or
                       `deconv`.
"""
function compute_cc(FFT1::FFTData, FFT2::FFTData, maxlag::Float64;
                    corr_type::String="CC")

    comp = FFT1.name[end] * FFT2.name[end]
    # get intersect of dates; return nothing if no intersect
    inter = intersect(FFT1.t,FFT2.t)
    if length(inter) == 0
        throw(ArgumentError("No common windows for $(FFT1.name)-$(FFT2.name) $(FFT1.id)"))
    end

    N = convert(Int,round(FFT1.cc_len * FFT1.fs)) # number of data points
    ind1 = findall(x -> x ∈ inter, FFT1.t)
    ind2 = findall(x -> x ∈ inter, FFT2.t)
    corr = correlate(@views(FFT1.fft[:,ind1]), @views(FFT2.fft[:,ind2]), N,
                     convert(Int,round(maxlag * FFT1.fs)))
    rotated = false

    return CorrData(FFT1, FFT2, comp, rotated, corr_type,
                    maxlag, inter, corr)

end
