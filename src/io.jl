using DelimitedFiles
using Dates

writedata(angles, fmt::String) = 
    writedata(angles, fmt, Dates.format(now(), "yyyymmddHH")) # Default filename from current time

function writedata(angles, fmt::String, name::String)
    fmt = lowercase(fmt) # User defined file format specifier
    
    if fmt == "csv"
        open(name*".dat", "w") do file
            write(file, "alphas,betas,gammas,weight\n")
            writedlm(file, angles, ',')
        end
    elseif fmt == "tsv"
        open(name*".dat", "w") do file
            write(file, "alphas,betas,gammas,weight\n")
            writedlm(file, angles)
        end
    elseif fmt == "mat"
        # To-do
    elseif fmt == "jld2"
        # To-do
    end
end
