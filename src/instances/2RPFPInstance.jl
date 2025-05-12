struct TwoRPFPInstance
    n::Int
    p::Matrix{Int}
    phat::Matrix{Int}
    Γ1::Int
    Γ2::Int
end

function TwoRPFPInstanceReader(path::String, Γ1_percent::Float64, Γ2_percent::Float64)
    # Read the file
    file = open(path, "r")
    lines = readlines(file)
    close(file)
    # remove empty lines
    lines = filter(line -> !isempty(line), lines)
    n = length(lines)
    p = zeros(Int, n, 2)
    phat = zeros(Int, n, 2)
    for i in 1:n
        #strip the line of leading and trailing whitespace
        lines[i] = string.(strip(lines[i]))
        #split the line by tab
        values = string.(split())
        #parse the values
        p1, p2, phat1, phat2 = parse.(Int, values[1:4])
        #create the matrix
        p[i, :] = [p1, p2]
        phat[i, :] = [phat1, phat2]
    end

    return TwoRPFPInstance(n, p, phat, floor(Int, n * Γ1_percent), floor(Int, n * Γ2_percent))
    
end