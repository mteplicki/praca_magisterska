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
    for i in 1:n
        #strip the line of leading and trailing whitespace
        line = string.(strip(lines[i+1]))
        #split the line by tab
        values = string.(split(line, "\t"))
        values_splitted = [parse(Int, string.(split(v))[2]) for v in values]
        p[i, :] = values_splitted
    end
    return n, m, p
end