
# A DataSet defines a posterior and stores all the 
# non-argument variables needed to compute it
abstract type DataSet end


function (ds::DataSet)(;θ...)
    DataSet(map(fieldvalues(ds)) do v
        (v isa ParamDependentOp) ? v(;θ...) : v
    end...)
end

adapt_structure(to, ds::DataSet) = DataSet(adapt(to, fieldvalues(ds))...)


@kwdef mutable struct DefaultDataSet{F} <: DataSet
    
    d :: F           # data
    Cϕ               # ϕ covariance
    Cf               # unlensed field covariance
    Cf̃ = nothing     # lensed field covariance (not always needed)
    Cn               # noise covariance
    M  = 1           # user mask
    B  = 1           # beam and instrumental transfer functions
    Cn̂ = Cn          # approximate noise covariance, diagonal in same basis as Cf
    M̂  = M           # approximate user mask, diagonal in same basis as Cf
    B̂  = B           # approximate beam and instrumental transfer functions, diagonal in same basis as Cf
    D  = IdentityOp  # mixing matrix for mixed parametrization
    G  = IdentityOp  # reparametrization for ϕ
    L  = alloc_cache(LenseFlow(similar(diag(Cϕ))),d) # a CachedLenseFlow which will be reused for memory
    
    DefaultDataSet(::UndefInitializer) = new()
    
end



function check_hat_operators(ds::DataSet)
    @unpack B̂, M̂, Cn̂, Cf = ds()
    @assert(all([(L isa Scalar) || (L isa typeof(Cf)) || (Cf isa FlatIEBCov && L isa DiagOp{<:FlatIEBFourier}) for L in [B̂,M̂,Cn̂]]),
            "B̂, M̂, Cn̂ should be scalars or the same type as Cf")
end


    
@doc doc"""
    resimulate(ds::DataSet; f=..., ϕ=...)
    
Resimulate the data in a given dataset, potentially at a fixed f and/or ϕ (both
are resimulated if not provided)
"""
function resimulate(ds::DataSet; f=simulate(ds.Cf), ϕ=simulate(ds.Cϕ), n=simulate(ds.Cn), f̃=ds.L(ϕ)*f, return_truths=false)
    @unpack M,P,B = ds
    @set! ds.d = M*P*B*f̃ + n
    return_truths ? @namedtuple(ds,f,ϕ,n,f̃) : ds
end



function get_pol_keys(pol)
    pol = Symbol(pol)
    (S, spectra, F_qu, F_eb, nF) = @match pol begin
        :I  => (S0,  (:TT,),            FlatMap,    FlatFourier,    1)
        :P  => (S2,  (:EE,:BB),         FlatQUMap,  FlatEBFourier,  2)
        :IP => (S02, (:TT,:EE,:BB,:TE), FlatIQUMap, FlatIEBFourier, 3)
        _   => throw(ArgumentError("`pol` should be one of :I, :P, or :IP"))
    end
    @namedtuple(S, spectra, F_qu, F_eb, nF)
end


function load_cmb!(ds::DataSet; θpix, Nside, pol, T, storage, rfid, Cℓ, _...)

    # the biggest ℓ on the 2D fourier grid
    ℓmax = round(Int,ceil(√2*fieldinfo(Flat(θpix=θpix,Nside=Nside)).nyq)+1)
    
    # CMB Cℓs
    if Cℓ == nothing
        Cℓ = camb(r=rfid, ℓmax=ℓmax)
    end
    
    # covariances
    @unpack S, spectra = get_pol_keys(pol)
    P = Flat(Nside=Nside, θpix=θpix, ∂mode=∂mode)
    Cϕ₀   = adapt(storage, Cℓ_to_Cov(P, T, S0, (Cℓ.total.ϕϕ)))
    Cfs   = adapt(storage, Cℓ_to_Cov(P, T, S,  (Cℓ.unlensed_scalar[k] for k in spectra)...))
    Cft   = adapt(storage, Cℓ_to_Cov(P, T, S,  (Cℓ.tensor[k]          for k in spectra)...))
    ds.Cf̃ = adapt(storage, Cℓ_to_Cov(P, T, S,  (Cℓ.total[k]           for k in spectra)...))
    ds.Cf = ParamDependentOp((mem; r=rfid, _...)->(mem .= Cfs + T(r/rfid)*Cft), similar(Cfs))
    ds.Cϕ = ParamDependentOp((mem; Aϕ=1  , _...)->(mem .= T(Aϕ) .* Cϕ₀), similar(Cϕ₀))
    
    ds
    
end

function load_noise!(ds::DataSet; μKarcminT, ℓknee, αknee, ℓmax)
    
    # noise Cℓs (these are non-debeamed, hence beamFWHM=0 below; the beam comes in via the B operator)
    if (Cℓn == nothing)
        Cℓn = noiseCℓs(μKarcminT=μKarcminT, beamFWHM=0, ℓknee=ℓknee, αknee=αknee, ℓmax=ℓmax)
    end
    
    # covariances
    Cn̂  = adapt(storage, Cℓ_to_Cov(Pix_data, T, S,  (Cℓn[k]                for k in ks)...))
    if (Cn == nothing); Cn = Cn̂; end

end


function load_beam!(ds::DataSet; beamFWHM)

    # beam
    if (B == nothing)
        B̂ = B = adapt(storage, Cℓ_to_Cov(Pix, T, S, ((k==:TE ? 0 : 1) * sqrt(beamCℓs(beamFWHM=beamFWHM)) for k=ks)..., units=1))
    end

end


function load_mask!(ds::DataSet; θpix, Nside, pol, T, storage, M, M̂, pixel_mask_kwargs)
    
    @unpack spectra, nF = get_pol_keys(pol)
    
    # data mask
    if (M == nothing)
        M̂ = M = adapt(storage, Cℓ_to_Cov(Pix, T, S, ((k==:TE ? 0 : 1) * bandpass_mask.diag.Wℓ for k in spectra)...; units=1))
        if (pixel_mask_kwargs != nothing)
            M = M * adapt(storage, Diagonal(F{Pix_data}(repeated(T.(make_mask(Nside,θpix; pixel_mask_kwargs...).Ix),nF)...)))
        end
    end
    if diag(M̂) isa BandPass
        M̂ = Diagonal(M̂ * one(diag(Cf)))
    end

end

function default_mixing!(ds::DataSet; D, G)
    
    @unpack Cf, Cϕ, Cn̂ = ds
    
    # D mixing matrix
    if (D == nothing)
        σ²len = T(deg2rad(5/60)^2)
        D = ParamDependentOp(
            function (mem;r=rfid,_...)
                Cfr = Cf(mem,r=r)
                mem .= sqrt(Diagonal(diag(Cfr) .+ σ²len .+ 2*diag(Cn̂)) * pinv(Cfr))
            end,
            similar(Cf())
        )
    end

    if (G == nothing)
        Nϕ = quadratic_estimate(ds,(pol in (:P,:IP) ? :EB : :TT)).Nϕ / Nϕ_fac
        G₀ = @. nan2zero(sqrt(1 + 2/($Cϕ()/Nϕ)))
        G = ParamDependentOp((;Aϕ=1,_...)->(@. nan2zero(sqrt(1 + 2/(($(Cϕ(Aϕ=Aϕ))/Nϕ)))/G₀)))
    end
    ds.G = G
    
end


function seed_for_storage!(seed, storage)
    if (seed != nothing)
        if storage == Array; Random.seed!(seed)
        elseif storage == CuArray; CuArrays.CURAND.seed!(seed)
        else; error("Don't know how to set seed for storage=$storage"); end
    end
end


function load_sim_dataset!(
    
    ds :: DataSet;
    
    # basic configuration
    θpix,
    Nside,
    pol,
    T = Float32,
    storage = Array,
    
    # noise parameters, or set Cℓn or even Cn directly
    μKarcminT = 3,
    ℓknee = 100,
    αknee = 3,
    Cℓn = nothing,
    Cn = nothing,
    
    # beam parameters, or set B directly
    beamFWHM = 0,
    B = nothing,
    
    # mask parameters, or set M directly
    pixel_mask_kwargs = nothing,
    bandpass_mask = LowPass(3000),
    M = nothing, M̂ = nothing,
    
    # theory
    rfid = 0.05,
    Cℓ = nothing,
    
    seed = nothing,
    D = nothing,
    G = nothing,
    Nϕ_fac = 2,
    ϕ=nothing, f=nothing, f̃=nothing, Bf̃=nothing, n=nothing, d=nothing, # can override any of these simulated fields
    L = LenseFlow,
    ∂mode = fourier∂
    )
    
    load_cmb!(ds; Base.@locals()...)
    load_mask!(ds; Base.@locals()...)
    load_noise!(ds; Base.@locals()...)
    load_beam!(ds; Base.@locals()...)
    default_mixing!(ds; Base.@locals()...)
    seed_for_storage!(seed, storage)
    resimulate!(ds)
    ds
    
end
