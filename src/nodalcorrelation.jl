export scalablexcorr, nodalxcorr, autocorrelate

function scalablexcorr(A::AbstractArray,Ntau::Int,kthresh::AbstractFloat)
    V,S,U = svd(A)
    kind = findall(S .> maximum(S) * 0.05)
    kmax = maximum(kind)

    Nt,Ns = size(A)

    # check that Ntau < Nt / 2 
    Ntau = min(Ntau,Nt ÷ 2)

    # preallocate view 
    V0 = @view V[Ntau÷2:Nt-Ntau÷2-1,1:kmax]
    U0 = @view U[:,1:kmax]
    if isa(V,AbstractGPUArray)
        VT0 = transpose(V[Ntau÷2:Nt-Ntau÷2-1,1:kmax])
        UT0 = transpose(U[:,1:kmax])
    else
        VT0 = V0'
        UT0 = U0'
    end
    
    Cout = similar(A,Ns,Ns,Ntau)
    
    for itau = -Ntau ÷ 2 + 1 : Ntau ÷ 2 
        W = VT0 * V[itau + Ntau÷2:Nt + itau - Ntau÷2-1,1:kmax]
        Cout[:,:,itau+Ntau÷2] .= U0 * W * UT0
    end
    return Cout
end

function scalablexcorr(N::NodalData,Ntau::Real;kthresh=0.05)
    if isa(Ntau,AbstractFloat)
        Ntau = convert(Int,round(Ntau * N.fs[1]))
    end
    return scalablexcorr(N.data,Ntau,kthresh)
end

function nodalxcorr(N::NodalData,maxlag::Real)
    if isa(maxlag,AbstractFloat)
        maxlag = convert(Int,round(maxlag * N.fs[1]))
    end
    return nodalxcorr(N.data,maxlag)

end

function nodalxcorr(A::AbstractArray,maxlag::Int)
    Nt,Ns = size(A)
    Ncorr = Ns * (Ns -1) ÷ 2 
    Cout = similar(A,maxlag * 2 + 1,Ncorr)
    FFTs = rfft(A,1)
    cstart = 0
    cend = 0
    for ii = 1:Ns-1
            cstart = cend + 1 
            cend = cstart + Ns - ii - 1
            Cout[:,cstart:cend] .= correlate(FFTs[:,ii],FFTs[:,ii+1:end],Nt,maxlag)
    end
    return Cout
end

function autocorrelate(N::NodalData,maxlag::Real)
    if isa(maxlag,AbstractFloat)
        maxlag = convert(Int,round(maxlag * N.fs[1]))
    end
    return autocorrelate(N.data,maxlag)
end

function autocorrelate(A::AbstractArray,maxlag::Int64)
    Nt,Ns = size(A)
    FFT = rfft(A,1)
    corrT = irfft(conj.(FFT) .* FFT,Nt,1)
    # return corr[-maxlag:maxlag]
    t = vcat(0:Int(Nt  / 2)-1, -Int(Nt  / 2):-1)
    ind = findall(abs.(t) .<= maxlag)
    newind = fftshift(ind,1)
    return corrT[newind,:]
end