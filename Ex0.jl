include("ZN.jl")
using DynamicPolynomials
using JuMP
using GLPK
using LinearAlgebra
using SparseArrays, Random, COSMO
Random.seed!(1234)
# Number of random instances for each pair (N, alpha).
global MAXtempi=100
#Absolute and relative tolerance used by COSMO.
global eta=1e-6
# alpha used in the experiments.
global alphalist=[1e-1, 1e-2,1e-4,1e-3]
# # Orders of the cyclic groups Z_N used in Table 3 and Figure 1.
global Nlist= [BigInt(10^3),BigInt(10^4),BigInt(10^5),BigInt(10^6),BigInt(10^7)]

# Generate a random real self-adjoint element supported on {0, ±1, ..., ±d}
# and rescale it so that alpha <= f <= 1 pointwise on Z_N.

function RandomPosFunction(N, d, alpha=0)

    D = Dict(0 => 0.0 * im)
    for k = 1:d
        D[k] =  (randn()) 
    end

    f = CyclicGroupAlgebraElement(D, N)
    f = f + f'

    V = real(to_dense(sqrt(N) * ifft(f)))
    M = maximum(V)
    m = minimum(V)
  f=alpha + (1 - alpha) / (M - m) * (f - m) 
    return f
end



function GenerateDegreeBoundGroupElements(N, d)
    S = Vector{CyclicGroupElement}()
    for k = (-d):1:d
        push!(S, CyclicGroupElement(k, N))
    end
    return S
end


# Build the linear constraints for the SOHS/FSOS feasibility SDP:
# dot(A[z], Q) == C_f(z), with Q positive semidefinite.
function GenerateConstractMAtrix(f, S)
    n = length(S)
    b = Dict()
    A = Dict()
    coeffs = f.coefficients
    ProdMatx = Matrix{BigInt}(undef, n, n)
    for i = 1:n
        for j = 1:n
            x = S[j] - S[i]
            ProdMatx[i, j] = x.value
        end
    end
    I = unique(ProdMatx)

    for i in I
        A[i] = (ProdMatx .== i) * 1.0
        x = CyclicGroupElement(i, f.N)
        if (x = findGroupkey(coeffs, x)) != false
            b[i] = coeffs[x]
        else
            b[i] = 0
        end
    end
    return A, b, I
end



function SolveSDP(A, b, I, S, solver="COSMO")
    if solver == "Mosek"
        model = Model(Mosek.Optimizer)
    end
    if solver == "COSMO"
        

        model = JuMP.Model(COSMO.Optimizer)
        set_optimizer_attribute(model, "eps_abs", eta)
        set_optimizer_attribute(model, "eps_rel", eta)
        set_optimizer_attribute(model, "verbose", false)
    end
    if solver=="SDPT3"
        model = Model(SDPT3.Optimizer)
    end
    set_silent(model)

    @variable(model, Q[1:length(S), 1:length(S)] in PSDCone()); 
    for i in I
        @constraint(model, dot(A[i], Q) == b[i]);

    end
    optimize!(model)

    # output
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal solution found")
    else
        println("Optimal solution not found. Termination status: ", termination_status(model))
    end
    return value.(Q), termination_status(model), model
end



function NumericalExample(N, d, alpha, solver="COSMO")

    f = RandomPosFunction(N, d, alpha);
    l2err = 1;
    k = 0;
    states='a';
    while   !(states==OPTIMAL::TerminationStatusCode) #
        k += 1
        S = GenerateDegreeBoundGroupElements(f.N, k * d)
        A, b, I = GenerateConstractMAtrix(f, S)
        
        Q, states, model = SolveSDP(A, b, I, S, solver)
        Fsos = CyclicGroupAlgebraElement(f.N)
        for i = 1:length(S)
            for j = 1:length(S)
                Fsos = Fsos + CyclicGroupAlgebraElement(Dict(S[j] - S[i] => Q[i, j]), f.N)
            end
        end
        e = Fsos - f;
        l2err = (norm(to_dense(e)))
        println([N, d, k, alpha, l2err])
    end

    return N, d, k, alpha, l2err,states

end

model = 0

global LinesOfTable
LinesOfTable = 0
Tab = []

dlist=[5]
for N in Nlist
        for d in dlist
        for alpha in alphalist
            for tempi = 1:MAXtempi
                global LinesOfTable =LinesOfTable+ 1
                # println([N, d, alpha])
                N, d, k, alpha, l2err, states = NumericalExample(N, d, alpha);
                push!(Tab, [N, d, k, alpha, l2err])
            end
        end
    end
end

Klist = []
for i = 1:length(Tab)
    push!(Klist, Tab[i][3])
end

println("K-1 norm:")
println(norm(Klist.-1))
Tabk=Dict()
for N in Nlist 
    for d in dlist
        for alpha in  alphalist
    Tabk[N,d,round(alpha*1e5)]=0
    end
    end
end
for i=1:length(Tab)
Tabk[Tab[i][1],Tab[i][2],round(Tab[i][4]*1e5)]+=(1+10*Tab[i][3])
end

#Print the result
println("N\t\t\\alpha\t\tmean value of sparsity")
for (keys, val) in Tabk
    N_val = keys[1]
    alpha_key = keys[3]
    k = Int(5 - log10(alpha_key))
    mean_val = val / MAXtempi
    println("$N_val\t\t10^(-$k)\t\t$mean_val")
end
