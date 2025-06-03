using Serialization
using Printf
using FilePathsBase: walkdir, isdir, dirname, joinpath, mkpath
import MasterThesis: TwoRPFPInstanceReader, save_instance

function generate_2RFP(dir_in::String,
                       dir_out::String;
                       Γ1_pct::Float64=0.1,
                       Γ2_pct::Float64=0.1)
    # utwórz katalog wyjściowy jeśli nie istnieje
    isdir(dir_out) || mkpath(dir_out)

    for (root, dirs, files) in walkdir(dir_in)
        for file in files[1:min(5, length(files))]
            # pomijamy ukryte pliki
            startswith(file, ".") && continue

            path_in  = joinpath(root, file)
            # wczytaj instancję 2RPF
            inst = TwoRPFPInstanceReader(path_in, Γ1_pct, Γ2_pct)

            # zbuduj nazwę wyjściową analogicznie do oryginalnej ścieżki
            rel = replace(path_in, dir_in*"/" => "")
            name = replace(rel, r"[\/\.]" => "_") * "_Γ1=$(Γ1_pct)_Γ2=$(Γ2_pct).bin"
            path_out = joinpath(dir_out, name)

            @printf("[%s] -> [%s]\n", path_in, path_out)
            save_instance(path_out, inst)
        end
    end
end

# nowa funkcja: generuje czytelną nazwę wyjściową
function make_output_name(path_in::AbstractString,
                          dir_in::AbstractString;
                          Γ1_pct::Float64,
                          Γ2_pct::Float64)
    rel = relpath(path_in, dir_in)              # ścieżka względna
    parts = splitpath(rel)                      # ["subdir","plik.ext"]
    stem, _ = splitext(parts[end])              # ("plik", ".ext")
    parts[end] = stem
    base = join(parts, "_")                     # "subdir_plik"
    return "$(base)_Γ1=$(Γ1_pct)_Γ2=$(Γ2_pct).bin"
end

# uniwersalna wersja generatora
function generate_2RFP_crossplatform(dir_in::String,
                                     dir_out::String;
                                     Γ1_pct::Float64=0.1,
                                     Γ2_pct::Float64=0.1)
    isdir(dir_out) || mkpath(dir_out)

    for (root, dirs, files) in walkdir(dir_in)
        for file in files[1:min(5, length(files))]
            startswith(file, ".") && continue

            path_in = joinpath(root, file)
            inst    = TwoRPFPInstanceReader(path_in, Γ1_pct, Γ2_pct)

            name     = make_output_name(path_in, dir_in;
                                       Γ1_pct=Γ1_pct, Γ2_pct=Γ2_pct)
            path_out = joinpath(dir_out, name)

            @printf("[%s] -> [%s]\n", path_in, path_out)
            save_instance(path_out, inst)
        end
    end
end