struct TwoRPFPInstance
    n::Int
    p::Matrix{Float64}
    phat::Matrix{Float64}
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
    p = zeros(Float64, 2, n)
    phat = zeros(Float64, 2, n)
    for i in 1:n
        #strip the line of leading and trailing whitespace
        lines[i] = string.(strip(lines[i]))
        #split the line by tab
        values = string.(split(lines[i], "\t"))
        #parse the values
        p1, p2, phat1, phat2 = parse.(Float64, values[1:4])
        #create the matrix
        p[:, i] = [p1, p2]
        phat[:, i] = [phat1, phat2]
    end

    return TwoRPFPInstance(n, p, phat, floor(Int, n * Γ1_percent), floor(Int, n * Γ2_percent))
    
end

function TwoRPFPInstance(n::Int, G::Float64)
    p_min = 1
    p_max = 100

    p = rand(p_min:p_max, 2, n)

    phat = [rand(floor(p_i * 2 /G):ceil(p_i * 7 /G)) for p_i in p]

    Γ1 = rand(floor(5*10^-3*G*n):ceil(9*10^-3*G*n))
    Γ2 = rand(floor(5*10^-3*G*n):ceil(9*10^-3*G*n))

    return TwoRPFPInstance(n, p, phat, Γ1, Γ2)
    
end