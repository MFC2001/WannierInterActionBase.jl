function Base.read(io::IO, ::Type{wannier90_hr})
    # @assert atol ≥ 0 "atol should be positive!"
    #first line is comment.
    readline(io)

    norb = parse(Int, readline(io))

    nR = parse(Int, readline(io))
    nline = (nR - 1) ÷ 15
    ndegen = Vector{Int}(undef, nR)
    for i in 1:nline
        ndegen[((i-1)*15+1):(i*15)] .= parse.(Int, split(readline(io)))
    end
    ndegen[(nline*15+1):end] .= parse.(Int, split(readline(io)))

    data = readdlm(io)

    R = Vector{ReducedCoordinates{Int}}(undef, nR)
    nv = size(data, 1)
    is = Vector{Int}(undef, nv)
    js = Vector{Int}(undef, nv)
    iRs = Vector{Int}(undef, nv)
    vals = Vector{ComplexF64}(undef, nv)
    iR = 1
    R[1] = ReducedCoordinates{Int}(data[1, 1], data[1, 2], data[1, 3])
    for ih in Base.OneTo(nv)
        path_val_R = ReducedCoordinates{Int}(data[ih, 1], data[ih, 2], data[ih, 3])
        is[ih] = Int(data[ih, 4])
        js[ih] = Int(data[ih, 5])
        vals[ih] = complex(data[ih, 6], data[ih, 7])
        if path_val_R ≠ R[iR]
            iR += 1
            if iR > nR
                error("R is more than ndegen.")
            end
            R[iR] = path_val_R
        end
        iRs[ih] = iR
    end

    @assert all(i->1 ≤ i ≤ norb, is) "Invalid is detected: i ∈ [1, $n]."
    @assert all(j->1 ≤ j ≤ norb, js) "Invalid is detected: j ∈ [1, $n]."

    if all(iszero ∘ imag, vals)
        vals = real(vals)
    end

    return HRFull(is, js, iRs, vals, norb, R; ndegen, inbound=true, dofilter=true)
end
