import Base.BlasInt

# Equivalent of mjdgges
# S and T are modified in place
function qz!(S::Matrix{Float64}, T::Matrix{Float64})
    function selctg(alphar::Ptr{Float64}, alphai::Ptr{Float64}, beta::Ptr{Float64})
        const criterium::Float64 = 1.0 + 1e-6
        ai = unsafe_ref(alphai, 1)
        ar = unsafe_ref(alphar, 1)
        b = unsafe_ref(beta, 1)

        convert(BlasInt, ar^2 + ai^2 < criterium * b^2)
    end

    selctg2 = cfunction(selctg, BlasInt, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}))
    (n,) = size(S)
    @assert size(S) == (n,n)
    @assert size(T) == (n,n)
    @assert stride(S, 1) == 1
    @assert stride(T, 1) == 1
    
    alphar = Array(Float64, n)
    alphai = Array(Float64, n)
    beta = Array(Float64, n)
    lwork = 16*n+16 # Same heuristic choice than mjdgges
    work = Array(Float64, lwork)
    vsl = Array(Float64, n, n)
    vsr = Array(Float64, n, n)
    bwork = Array(BlasInt, n)
    sdim = Array(BlasInt, 1)
    info = Array(BlasInt, 1)
    
    ccall((:dgges_, Base.liblapack_name), Void, (Ptr{Uint8}, Ptr{Uint8}, Ptr{Uint8}, Ptr{Void},
                                           Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                                           Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
                                           Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                                           Ptr{Float64}, Ptr{BlasInt},
                                           Ptr{Float64}, Ptr{BlasInt},
                                           Ptr{Float64}, Ptr{BlasInt},
                                           Ptr{BlasInt}, Ptr{BlasInt}),
          &'N', &'V', &'S', selctg2, &n, S, &stride(S, 2), T, &stride(T, 2),
          sdim, alphar, alphai, beta, vsl, &stride(vsl, 2),
          vsr, &stride(vsr, 2), work, &lwork, bwork, info)

    @assert info[1] == 0

    eigs = Array(Complex128, n)
    for i = 1:n
        eigs[i] = alphar[i]/beta[i] + (alphai[i] == 0 && beta[i] == 0 ? 0.0 : alphai[i]/beta[i])*im
    end
    
    (vsr, sdim[1], eigs)
end

function decision_rules(m::Model, calib::Dict{Symbol, Float64}, steady_state::Vector{Float64})
    @assert length(steady_state) == m.n_endo
    n_dynvars = m.n_back_mixed + m.n_endo + m.n_fwrd_mixed
    yy = [ steady_state[m.zeta_back_mixed]; steady_state; steady_state[m.zeta_fwrd_mixed] ]
    x = zeros(m.n_exo)
    p = calib2vec(m, calib)
    (r, jacob) = m.dynamic_mf(yy, x, p)
    @assert size(jacob) == (m.n_endo, n_dynvars + m.n_exo)

    A = jacob[:, 1:n_dynvars]

    # Deal with static variables
    if m.n_static > 0
        S = jacob[:, m.n_back_mixed + m.zeta_static]
        (Q,R) = qr(S)
        @assert rank(R) == m.n_static
        A = qTmulQR(qrfact(S), A)
    end

    D = [ A[m.n_static+1:end, m.n_back_mixed + m.zeta_back_mixed] A[m.n_static+1:end, m.n_back_mixed+m.n_endo+1:end];
         eye(m.n_back_mixed, m.n_back_mixed)[m.beta_back, :] zeros(m.n_mixed, m.n_fwrd_mixed) ]

    Atilde0plus = A[m.n_static+1:end, m.n_back_mixed + m.zeta_fwrd_mixed]
    Atilde0plus[:, m.beta_fwrd] = zeros(m.n_dynamic, m.n_mixed)
    E = [ -A[m.n_static+1:end, 1:m.n_back_mixed] -Atilde0plus;
         zeros(m.n_mixed, m.n_back_mixed) eye(m.n_fwrd_mixed, m.n_fwrd_mixed)[m.beta_back, :] ]

    (Z, sdim, eigs) = qz!(E, D)

    # BK order condition
    @assert sdim == m.n_back_mixed

    # Compute DR wrt backward endogenous
    Zp = Z'
    Z22 = Zp[m.n_back_mixed+1:end, m.n_back_mixed+1:end]
    Z21 = Zp[m.n_back_mixed+1:end, 1:m.n_back_mixed]

    gyplus = -Z22\Z21

    Z11p = Z[1:m.n_back_mixed, 1:m.n_back_mixed]
    T11 = D[1:m.n_back_mixed, 1:m.n_back_mixed]
    S11 = E[1:m.n_back_mixed, 1:m.n_back_mixed]
    gyminus = Z11p*(T11\S11)/Z11p

    gy = Array(Float64, m.n_endo, m.n_back_mixed)
    gy[m.zeta_back_mixed, :] = gyminus
    gy[m.zeta_fwrd_mixed, :] = gyplus
    
    if m.n_static > 0
        Abreve0s = A[1:m.n_static, m.n_back_mixed + m.zeta_static]
        Abreve0d = A[1:m.n_static, m.n_back_mixed + m.zeta_dynamic]
        Abreveminus = A[1:m.n_static, 1:m.n_back_mixed]
        Abreveplus = A[1:m.n_static, m.n_back_mixed+m.n_endo+1:end]
        gy[m.zeta_static, :] = -Abreve0s \ (Abreveplus*gyplus*gyminus + Abreve0d*gy[m.zeta_dynamic, :] + Abreveminus)
    end

    # Compute DR wrt shocks
    B = jacob[:, m.n_back_mixed + (1:m.n_endo)]
    B[:, m.zeta_back_mixed] += jacob[:, m.n_back_mixed+m.n_endo+(1:m.n_fwrd_mixed)]*gyplus
    gu = -B\jacob[:,m.n_back_mixed+m.n_endo+m.n_fwrd_mixed+1:end]

    (gy, gu, eigs)
end
