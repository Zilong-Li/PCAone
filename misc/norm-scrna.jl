using CodecZstd
using Statistics

function nm(input::AbstractString)
    # N = number of cells, M = number of genes
    N = zeros(UInt32, 1)
    M = zeros(UInt32, 1)
    open(input) do file
        stream = ZstdDecompressorStream(file)
        read!(stream, N)
        read!(stream, M)
        close(stream)
    end
    return N[], M[]
end

# no.of row of csv
function nrow(;csvfile::AbstractString="")
    nrow = 0
    open(csvfile, "r") do f
        while !eof(f)
            nrow += 1
            readline(f)
        end
    end
    nrow
end

# no. of columns of csv
function ncol(;csvfile::AbstractString="")
    counter = 0
    open(csvfile, "r") do f
        xx = readline(f)
        xx = split(xx, ",")
        length(xx)
    end
end

function csv2bin(;csvfile::AbstractString="", binfile::AbstractString="")
    N = zeros(UInt32, 1)
    M = zeros(UInt32, 1)
    N[] = UInt32(nrow(csvfile=csvfile))
    M[] = UInt32(ncol(csvfile=csvfile))
    counter = 0
    open(binfile, "w") do file
        stream = ZstdCompressorStream(file)
        write(stream, N)
        write(stream, M)
        open(csvfile , "r") do f
            while !eof(f)
                counter += 1
                print("\r", counter)
                xx = readline(f)
                xx = split(xx, ",") # Imported as Character
                # Assume input data is Integer
                x = zeros(UInt32, length(xx))
                for i=1:length(x)
                    x[i] = floor(UInt32, parse(Float64, xx[i]))
                end
                write(stream, x)
            end
        end
        close(stream)
    end
    print("\n")
end

function getMedianLibSize(binfile::AbstractString)
    N, M  = nm(binfile)
    tmpN = zeros(UInt32, 1)
    tmpM = zeros(UInt32, 1)
    x = zeros(UInt32, M)
    libsize = zeros(UInt32, M)

    open(binfile) do file
        stream = ZstdDecompressorStream(file)
        read!(stream, tmpN)
        read!(stream, tmpM)
        for n = 1:N
            read!(stream, x)
            for i in 1:length(x)
                libsize[i] += x[i]
            end
        end
        close(stream)
    end

    return libsize
end


function normalizex(x::Array{UInt32,1}, libsize::Array{UInt32,1}, scale::AbstractString = "log")

    xx = convert(Vector{Float64}, x)
    # Normalization
    # for i in 1:length(xx)
    #     xx[i] = median(libsize) * xx[i] / libsize[i]
    # end
    xx = median(libsize) .* (xx ./ libsize)
    if scale == "log"
        xx = log10.(xx .+ 1)
    end

    # centering
    xx = xx .- mean(xx)

    return xx
end


function bin2csv(binfile::AbstractString, outfile::AbstractString)
    N, M  = nm(binfile)
    tmpN = zeros(UInt32, 1)
    tmpM = zeros(UInt32, 1)
    x = zeros(UInt32, M)
    libsize = getMedianLibSize(binfile)

    open(outfile, "w") do file2
        stream2 = ZstdCompressorStream(file2)
        open(binfile) do file
            stream = ZstdDecompressorStream(file)
            read!(stream, tmpN)
            read!(stream, tmpM)
            for n = 1:N
                read!(stream, x)
                # println(Float64.(x))
                x2 = normalizex(x, libsize)
                out = String[]
                for i in 1:length(x2)
                    push!(out, string(x2[i]))
                end
                write(stream2, join(out, ",") * "\n")
            end
          close(stream)
        end
        close(stream2)
    end
end

