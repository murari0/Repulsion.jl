using DelimitedFiles
using Dates

function writedata(angles, fileopts::String...)
    fmt = lowercase(first(fileopts)) # User defined file format specifier
    if length(fileopts) == 1
        name = Dates.format(now(), "yyyymmddHH") # Construct default filename from current date/hour string
    else if length(fileopts) == 2
        name = fileopts[2] # Accept user defined file name
    else
        @warn "Too many file options. Discarding " fileopts[3:end] # Reject extra arguments
        name = fileopts[2]
    end

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
