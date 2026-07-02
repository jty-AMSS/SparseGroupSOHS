using LinearAlgebra
using FFTW

# Elements of the cyclic group Z_N, written additively.
struct CyclicGroupElement
    value::BigInt  # Representative in {0, ..., N-1}.
    N::BigInt       # Group order.
# Constructor: reduce value modulo N into {0, ..., N-1}.
    function CyclicGroupElement(value::Integer, N::Number)
        N_big = BigInt(N)
        val_big = BigInt(value) % N_big
        val_big = val_big < 0 ? val_big + N_big : val_big 
        new(val_big, N_big)
    end
end

# display_value
function display_value(x::CyclicGroupElement)
    if x.value > x.N ÷ 2
        return x.value - x.N
    else
        return x.value
    end
end

# add
function Base.:+(x::CyclicGroupElement, y::CyclicGroupElement)
    @assert x.N == y.N "Group elements must belong to the same group."
    CyclicGroupElement(x.value + y.value, x.N)
end

# minus
function Base.:-(x::CyclicGroupElement, y::CyclicGroupElement)
    @assert x.N == y.N "Group elements must belong to the same group."
    CyclicGroupElement(x.value - y.value, x.N)
end

# Involution
function Base.conj(x::CyclicGroupElement)
    CyclicGroupElement(x.N - x.value, x.N)
end
#eqv
function Base.:(==)(a::CyclicGroupElement, b::CyclicGroupElement)
    return a.N == b.N && a.value == b.value
end


#show
function Base.show(io::IO, x::CyclicGroupElement)
    print(io, "($(display_value(x)))")
end



function findGroupkey(dict::Dict,x)
     for (k, v) in dict
       if k==x
        return k
       end
    end
    return false
end

 
# Elements of the group algebra C[Z_N].
struct CyclicGroupAlgebraElement
    coefficients::Dict{CyclicGroupElement,ComplexF64}  
    N::BigInt                                           
end

function CyclicGroupAlgebraElement(coeffs::Dict{<:Integer,<:Number}, N::Integer)
    N_big = BigInt(N)
    group_coeffs = Dict{CyclicGroupElement,ComplexF64}()
    for (k, v) in coeffs
        x = CyclicGroupElement(k, N_big)
        group_coeffs[x] = ComplexF64(v)
    end
    CyclicGroupAlgebraElement(group_coeffs, N_big)
end

#  Constructor for the zero element.
CyclicGroupAlgebraElement(N::Integer) = CyclicGroupAlgebraElement(Dict{CyclicGroupElement,ComplexF64}(), BigInt(N))

# Constructor for the identity element 1*e_0.
function CyclicGroupAlgebraElement(N::Integer, is_identity::Bool)
    if is_identity
        x0 = CyclicGroupElement(0, N)
        CyclicGroupAlgebraElement(Dict(x0 => 1.0 + 0.0im), BigInt(N))
    else
        CyclicGroupAlgebraElement(N)
    end
end

#Convert to a dense coefficient vector. This is intended only when N is manageable.
function to_dense(a::CyclicGroupAlgebraElement; max_size::Integer=big(10)^100)
    N = a.N
    if N > max_size
        error("N is too large ($N) to convert to a dense coefficient vector.")
    end

    arr = zeros(ComplexF64, Int(N))
    for (x, c) in a.coefficients
        arr[Int(x.value)+1] = c  
    end
    arr
end

# Convert a dense coefficient vector back to a group algebra element.
function from_dense(arr::Vector, N; threshold::Float64=1e-10)
    coeffs = Dict{CyclicGroupElement,ComplexF64}()
    N_big = BigInt(N)
    for (i, c) in enumerate(arr)
        if abs(c) > threshold  # Keep only numerically significant coefficients.
            x = CyclicGroupElement(i - 1, N_big)  # Add 1 because Julia arrays are 1-indexed.
            coeffs[x] = c
        end
    end
    CyclicGroupAlgebraElement(coeffs, N_big)
end

# add
function Base.:+(a::CyclicGroupAlgebraElement, b::CyclicGroupAlgebraElement)
    @assert a.N == b.N "Group elements must belong to the same group."
    coeffs = copy(a.coefficients)
    for (x, c) in b.coefficients
        t=x
        if (x=findGroupkey(coeffs,x))!=false
            coeffs[x] += c
            if iszero(coeffs[x])
                delete!(coeffs, x)
            end
        else
            coeffs[t] = c
        end
    end
    CyclicGroupAlgebraElement(coeffs, a.N)
end

# minus
function Base.:-(a::CyclicGroupAlgebraElement, b::CyclicGroupAlgebraElement)
    @assert a.N == b.N "Group elements must belong to the same group."
    coeffs = copy(a.coefficients)
    for (x, c) in b.coefficients
        t=x
        if (x=findGroupkey(coeffs,x))!=false
            coeffs[x] -= c
            if iszero(coeffs[x])
                delete!(coeffs, x)
            end
        else
            coeffs[t] = -c
        end
    end
    CyclicGroupAlgebraElement(coeffs, a.N)
end

# Addition by a scalar,
function Base.:+(a::CyclicGroupAlgebraElement, c::Number)
    coeffs = copy(a.coefficients)
    c_complex = ComplexF64(c)
    x0 = CyclicGroupElement(0, a.N)
    if (x=findGroupkey(coeffs,x0))!=false
        coeffs[x] += c_complex
        if iszero(coeffs[x])
            delete!(coeffs, x)
        end
    else
        coeffs[x0] = c_complex
    end

    CyclicGroupAlgebraElement(coeffs, a.N)
end

Base.:+(c::Number, a::CyclicGroupAlgebraElement) = a + c

# Subtraction by a scalar.
function Base.:-(a::CyclicGroupAlgebraElement, c::Number)
    coeffs = copy(a.coefficients)
    c_complex = ComplexF64(c)
    x0 = CyclicGroupElement(0, a.N)
    t=x0
    if (x0=findGroupkey(coeffs,x0))!=false
        coeffs[x0] -= c_complex
        if iszero(coeffs[x0])
            delete!(coeffs, x0)
        end
    else
        coeffs[t] = -c_complex
    end

    CyclicGroupAlgebraElement(coeffs, a.N)
end

function Base.:-(c::Number, a::CyclicGroupAlgebraElement)
    coeffs = Dict{CyclicGroupElement,ComplexF64}(x => -v for (x, v) in a.coefficients)
    x0 = CyclicGroupElement(0, a.N)
    if (x0=findGroupkey(coeffs,x0))!=false
        coeffs[x0] += ComplexF64(c)
    else
        x0 = CyclicGroupElement(0, a.N)
        coeffs[x0] = ComplexF64(c)
    end

    CyclicGroupAlgebraElement(coeffs, a.N)
end

# Multiplication by convolution.
function Base.:*(a::CyclicGroupAlgebraElement, b::CyclicGroupAlgebraElement)
    @assert a.N == b.N "Group elements must belong to the same group."
    coeffs = Dict{CyclicGroupElement,ComplexF64}()
    N = a.N

    for (x, c1) in a.coefficients
        for (y, c2) in b.coefficients
            z = x + y 
            product = c1 * c2
            if (z=findGroupkey(coeffs,z))!=false
                coeffs[z] += product
                if iszero(coeffs[z])
                    delete!(coeffs, z)
                end
            else
                coeffs[x+y] = product
            end
        end
    end

    CyclicGroupAlgebraElement(coeffs, N)
end


# scalar Multiplication
function Base.:*(a::CyclicGroupAlgebraElement, c::Number)
    coeffs = copy(a.coefficients)
    for (x, t) in coeffs
        coeffs[x] *= c
    end
    CyclicGroupAlgebraElement(coeffs, a.N)

end


function Base.:*(c::Number, a::CyclicGroupAlgebraElement)
    coeffs = copy(a.coefficients)
    for (x, t) in coeffs
        coeffs[x] *= c
    end
    CyclicGroupAlgebraElement(coeffs, a.N)

end




# Involution of a group algebra element.
function Base.conj(a::CyclicGroupAlgebraElement)
    # Conjugate the coefficients and apply the group involution x^* = -x mod N.
    coeffs = Dict{CyclicGroupElement,ComplexF64}(conj(x) => conj(c) for (x, c) in a.coefficients)
    CyclicGroupAlgebraElement(coeffs, a.N)
end


function Base.adjoint(a::CyclicGroupAlgebraElement)
return conj(a)
end


#   Discrete Fourier transform via FFTW.
function fft(a::CyclicGroupAlgebraElement; max_size::Integer=big(10)^100, threshold::Float64=1e-10)
    N = a.N

    # for huge N
    if N > max_size
        coeffs = Dict{CyclicGroupElement,ComplexF64}()
        elements = collect(keys(a.coefficients))

        for y in elements
            val = 0.0im
            for (x, c) in a.coefficients
                exponent = -2π * im * x.value * y.value / N
                val += c * exp(exponent)
            end
            coeffs[y] = val / sqrt(N)  # Normalization
        end
        return CyclicGroupAlgebraElement(coeffs, N)
    else
        # For manageable N, compute the full transform using FFTW.
        dense = to_dense(a, max_size=max_size)
        fft_result = FFTW.fft(dense) / sqrt(N)  # FFTW.fft is unnormalized; apply the unitary normalization manually.
        return from_dense(fft_result, N, threshold=threshold)
    end
end

# Inverse unitary discrete Fourier transform using FFTW.
function ifft(a::CyclicGroupAlgebraElement; max_size::Integer=big(10)^100, threshold::Float64=1e-10)
    N = a.N

    # for huge N
    if N > max_size
        coeffs = Dict{CyclicGroupElement,ComplexF64}()
        elements = collect(keys(a.coefficients))

        for x in elements
            val = 0.0im
            for (y, c) in a.coefficients
                exponent = 2π * im * x.value * y.value / N
                val += c * exp(exponent)
            end
            coeffs[x] = val / sqrt(N)  # Normalization
        end
        return CyclicGroupAlgebraElement(coeffs, N)
    else
        # For small N, compute the full inverse transform using FFTW.
        dense = to_dense(a, max_size=max_size)
        ifft_result = FFTW.ifft(dense) * sqrt(N)  # Normalization
        return from_dense(ifft_result, N, threshold=threshold)
    end
end

# Display and printing
function Base.show(io::IO, a::CyclicGroupAlgebraElement)
    if isempty(a.coefficients)
        print(io, "0 (∈ C[Z_$(a.N)])")
        return
    end

# Sort by display value.
    sorted_terms = sort(collect(a.coefficients), by=term -> display_value(term[1]))

    terms = []
    for (x, c) in sorted_terms
        # Format complex coefficients.
        c_str = if isreal(c) && imag(c) ≈ 0
            "$(real(c))"
        else
            "($c)"
        end
        push!(terms, "$c_str*$x")
    end

    print(io, join(terms, " + "), " (∈ C[Z_$(a.N)])")
end

# Compact printing
function Base.print(io::IO, a::CyclicGroupAlgebraElement)
    if isempty(a.coefficients)
        print(io, "0")
        return
    end

    # Sort by display value.
    sorted_terms = sort(collect(a.coefficients), by=term -> display_value(term[1]))

    terms = []
    for (x, c) in sorted_terms
        c_str = if isreal(c) && imag(c) ≈ 0
            "$(round(real(c), digits=4))"
        else
            "($(round(c, digits=4)))"
        end
        push!(terms, "$c_str*$x")
    end

    print(io, join(terms, " + "))
end
