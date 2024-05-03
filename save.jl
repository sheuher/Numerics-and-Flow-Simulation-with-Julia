using DelimitedFiles

saveDir = joinpath(@__DIR__, "dat")
if isdir(saveDir) != true; mkdir(saveDir); end

function save(Phi, PhiName::String, nT::Int64)
    
    varName = PhiName
    savePath = string( rpad(varName, 6, '_'), string(nT, pad=6), ".dat")
    savePath = joinpath(@__DIR__, "dat/$(savePath)")

    open(savePath, "w") do file
        writedlm(file, Phi)
    end

    varName
end

macro Name(arg)
    string(arg)
end