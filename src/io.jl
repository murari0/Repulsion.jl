using DelimitedFiles
using Dates

"""
    writedata(angles, fmt::String)

Writes data to file when only the file format is specified, using a default filename created from the current time in
yyyymmddHH format
"""
writedata(angles, fmt::String) = 
    writedata(angles, fmt, "grid_"*Dates.format(now(), "yyyymmddHH")) # Default filename from current time

"""
    writedata(angles, fmt::String, name::String)

Writes data to file using the file format specifier and filename provided.

The file format specifier is case-insensitive, and falls back to TSV if an unrecognised file format is requested. CSV and TSV
files are written using [`writedlm`](@ref) from the standard library package `DelimitedFiles`. If the MAT or JLD/JLD2 file formats
are requested, the function attempts to use the `matopen` or `jldopen` functions provided by the MAT or JLD/JLD2 packages. If the
required package has not been loaded by the user, the user is notified and the function falls back to TSV.
"""
function writedata(angles, fmt::String, name::String)
    fmt = lowercase(fmt) # User defined file format specifier
    
    if fmt == "csv"
        open(name*".dat", "w") do file
            write(file, "alphas,betas,gammas,weight\n")
            writedlm(file, angles, ',')
        end
    elseif fmt == "tsv"
        open(name*".dat", "w") do file
            write(file, "alphas\tbetas\tgammas\tweight\n")
            writedlm(file, angles, '\t')
        end
    elseif fmt == "mat"
        alphas = [angles[i][1] for i in eachindex(angles)]
        betas = [angles[i][2] for i in eachindex(angles)]
        gammas = [angles[i][3] for i in eachindex(angles)]
        weights = [angles[i][4] for i in eachindex(angles)]
        try
            matopen(name*".mat", "w") do file
                write(file, "alphas", alphas)
                write(file, "betas", betas)
                write(file, "gammas", gammas)
                write(file, "weights", weights)
            end
        catch e
            if isa(e, UndefVarError)
                @warn "Package MAT.jl is not loaded/installed - defaulting to TSV file"
                writedata(angles, "tsv", name)
            else
                rethrow()
            end
        end                
    elseif fmt == "jld" || fmt == "jld2"
        try
            jldopen(name*"."*fmt, "w") do file
                write(file, "angles", angles)
            end
        catch e
            if isa(e, UndefVarError)
                @warn "Package JLD.jl/JLD2.jl is not loaded/installed - defaulting to TSV file"
                writedata(angles, "tsv", name)
            else
                rethrow()
            end
        end
    else
        @warn "Unrecognised file format specifier - defaulting to TSV file" fmt
        writedata(angles, "tsv", name)
    end
end
