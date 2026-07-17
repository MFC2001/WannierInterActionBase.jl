# installation of HDF5.jl

HDF5.jl: Note that it can use local library. One method is use the libhdf5 from anaconda:
```julia
julia> using Preferences, UUIDs

julia> set_preferences!(
    UUID("f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"), # UUID of HDF5.jl
    "libhdf5" => "/conda/lib/libhdf5.so",
    "libhdf5_hl" => "/conda/lib/libhdf5_hl.so", force = true)
```

FFTW.jl: need GLIBC_2.14.
AbstractFFTs.jl: need backend