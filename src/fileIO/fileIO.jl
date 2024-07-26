
"""
    load(f::File{format"DAT"}, T=Float64)

Load a .dat file via the FileIO package.
"""
function load(f::File{format"DAT"}, T=Float64)

    return readMultipatch(FileIO.filename(f), T)
end


"""
    load(f::File{format"DAT"}, T=Float64)

Load a step file via the FileIO package.

If the file has no ending the magic byte is used to identify the file format.

Supported endings: .stp, .step, .stpnc, .p21, .210
"""
function load(f::File{format"STEP"}, T=Float64)

    return readStep(FileIO.filename(f), T)
end
