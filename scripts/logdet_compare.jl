push!(LOAD_PATH, pwd()*"/src")

using CMBFields
using BayesLensSPTpol: class
cls = class();

nside = 32
T = Float64
P = Flat{1,nside}

ϕ = simulate(Cℓ_to_cov(FlatS0FourierDiagCov{T,P}, cls[:ell], cls[:ϕϕ]))

"""
Gets a matrix representation of an operator in the T->T basis
TODO: needs some tweaks to work generally then move into main source
"""
function matrix{T<:Field}(::Type{T}, L::LinearFieldOp)
    hcat(((convert(T,L*T((x=zeros(nside,nside); x[i]=1; x))))[:] for i=1:nside^2)...);
end

L_taylens = FlatS0TaylensOp(ϕ,order=10)
L_lenseflow = LenseFlowOp(ϕ,ode45{1e-6,1e-6})

@show logdet(matrix(FlatS0Map{T,P},L_taylens))
@show logdet(matrix(FlatS0Map{T,P},L_lenseflow))
