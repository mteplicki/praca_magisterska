function MultiUnrelatedInstanceReader(path::String)
    # Read the file
    file = open(path, "r")
    lines = readlines(file)
    close(file)
    # Parse the first line
    # first line is n m Γ
    n, m, _ = parse.(Int, string.(split(lines[1])))
    # next are the n lines with format: 0 p0 \t 1 p1 \t 2 p2 ... m pm
    p = zeros(Int, m, n)
    # there are no phat values in the file, we will not generate them now
    #skip net line
    
    for i in 1:n
        #strip the line of leading and trailing whitespace
        line = string.(strip(lines[i+2]))
        #split the line by tab
        values = string.(split(line, "\t"))
        @show values
        # select every second value
        values_splitted = values[2:2:end]
        # convert to Int
        values_splitted_parsed = parse.(Int, values_splitted)
        p[:, i] = values_splitted_parsed
    end
    return n, m, p
end