using MasterThesis

"""
  generate_instances(dir::String, k::Int, n::Int, R::Float64, T::Float64, G::Float64)
  Wygeneruj i zapisz `k` losowych instancji single machine due dates (parametr G) do katalogu `dir`.
  Będzie używać podprocedury `SingleMachineDueDates` z `src/instances/SingleMachineDueDates.jl`.
    Pliki będą nazwane np. `inst_n=010_R=R_T=t_G=G_k=k.bin`, ... `inst_n=100_R=R_T=t_G=G_k=k.bin`.
    n musi mieć co najmniej 3 cyfry.
"""
function generate_instances(dir::AbstractString, k::Int, n::Int, R::Float64, T::Float64, G::Float64)
    isdir(dir) || mkpath(dir)
    for i in 1:k
        inst = SingleMachineDueDates(n, R, T, G)
        path = joinpath(dir, "inst_n=$(lpad(n, 3, '0'))_R=$(R)_T=$(T)_G=$(G)_k=$(i).bin")
        save_instance(path, inst)
    end
end

dir = "./instances/single_tardy_jobs"

for n = [10,20]
    for R = [0.2,0.6]
        for T = [0.2,0.6]
            for G = [10., 100.]
                generate_instances(dir, 3, n, R, T, G)
            end
        end
    end
end

#create file model_type with the content "SingleTardyJobs"
open(dir*"/model_type", "w") do f
    write(f, "SingleTardyJobs")
end


