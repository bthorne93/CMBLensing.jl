# these are in a different file than healpix.jl because of
# https://github.com/timholy/Revise.jl/issues/290


# Healpy doesn't yet have the zbounds option, so call into the Healpix Fortran
# libraries directly. 

# needing to manually dlopen this is probably a bug in Healpix or gfortan, not
# sure which...
using Libdl
@init try
    Libdl.dlopen("libgomp.so.1",Libdl.RTLD_GLOBAL)
catch
    @warn "Failed to load libgomp. Healpix support may not work."
end

@generated function map2alm(maps::Array{T,N}; Nside=hp.npix2nside(size(maps,1)), ℓmax=2Nside, mmax=ℓmax, zbounds=[-1,1]) where {N,T<:Union{Float32,Float64}}
    (spin, Tspin) = if (N==1)
        (), () 
    elseif (N==2)
        (2,), (Ref{Int32},)
    else
        error("maps should be Npix-×-1 or Npix-×-2")
    end
    
    
    fn_name = "__alm_tools_MOD_map2alm_$(N==1 ? "sc" : "spin")_$((T==Float32) ? "s" : "d")"
    quote
        
        # for spin-2, Fortran map2alm needs full 12Nside^2-×-2 map array even
        # with zbounds
        if N==2 && size(maps,1)!=12Nside^2
            fullmaps = Array{T,N}(undef,12Nside^2,2)
            fullmaps[1:size(maps,1),:] .= maps
            maps = fullmaps
        end
        
        aℓms = Array{Complex{T}}(undef, ($N,ℓmax+1,mmax+1))
        aℓms .= NaN
        ccall(
            ($fn_name, "libhealpix.so"), Nothing,
            (Ref{Int32}, Ref{Int32}, Ref{Int32}, $(Tspin...), Ref{T}, Ref{Complex{T}}, Ref{Float64}, Ref{Nothing}),
            Nside, ℓmax, mmax, $(spin...), maps, aℓms, Float64.(zbounds), C_NULL
        )
        aℓms
    end
end

@generated function alm2map!(maps::Array{T,N}, aℓms::Array{Complex{T},3}; Nside=hp.npix2nside(size(maps,1)), ℓmax=(size(aℓms,2)-1), mmax=(size(aℓms,3)-1), zbounds=[-1,1]) where {N,T<:Union{Float32,Float64}}
    (spin, Tspin) = if (N==1)
        (), () 
    elseif (N==2)
        (2,), (Ref{Int32},)
    else
        error("maps should be Npix-×-1 or Npix-×-2")
    end
    fn_name = "__alm_tools_MOD_alm2map_$(N==1 ? "sc_wrapper" : "spin")_$((T==Float32) ? "s" : "d")"
    quote
        
        # for spin-2, Fortran alm2map needs full 12Nside^2-×-2 map array even
        # with zbounds
        if N==2 && size(maps,1)!=12Nside^2
            fullmaps = Array{T,N}(undef,12Nside^2,2)
            fullmaps[1:size(maps,1),:] .= maps
        else
            fullmaps = maps
        end
        
        ccall(
           ($fn_name, "libhealpix.so") , Nothing,
           (Ref{Int32}, Ref{Int32}, Ref{Int32}, $(Tspin...), Ref{Complex{T}}, Ref{T}, Ref{Float64}, Ref{Nothing}),
           Nside, ℓmax, mmax, $(spin...), aℓms, fullmaps, Float64.(zbounds), C_NULL
        )
        
        if N==2 && size(maps,1)!=12Nside^2
            maps .= fullmaps[1:size(maps,1),:]
        end
        maps
    end
end
