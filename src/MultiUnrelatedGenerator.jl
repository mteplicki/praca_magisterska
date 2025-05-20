using MasterThesis
"""
  generate_instances(dir::String, k::Int, n::Int, m::Int, G::Float64;
                     p_min=1, p_max=100)

Wygeneruj i zapisz `k` losowych instancji m×n (parametr G) do katalogu `dir`.
Pliki będą nazwane np. `inst_1.bin`, …, `inst_k.bin`.
"""
function generate_instances(dir::AbstractString, k::Int, n::Int, m::Int, G::Float64;
                             p_min::Int=1, p_max::Int=100)
    isdir(dir) || mkpath(dir)
    for i in 1:k
        inst = MultiUnrelatedInstance(n, m, G, p_min, p_max)
        path = joinpath(dir, "inst_$(n)x$(m)_G$(G)_$(i).bin")
        save_instance(path, inst)
    end
end