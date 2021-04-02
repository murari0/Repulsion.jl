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
    else if fmt == "tsv"
        open(name*".dat", "w") do file
            write(file, "alphas,betas,gammas,weight\n")
            writedlm(file, angles)
        end
    else if fmt == "mat"
        # To-do
    else if fmt == "jld2"
        # To-do
    end
end
