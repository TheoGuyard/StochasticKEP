"""
    kep_pool

Contruct a pool from a `.wmd` and a `.dat` files from PrefLib.
    
# Parameters
* `wmd_file::String` : Absolute path of the `.wmd` file.
* `dat_file::String` : Absolute path of the `.dat` file.
"""
function kep_pool(wmd_file::String, dat_file::String)
    
    wmd_file_name = split(split(wmd_file, '/')[end], '.')[1]
    dat_file_name = split(split(dat_file, '/')[end], '.')[1]

    wmd_file_name == dat_file_name || throw(ArgumentError(".wmd and .dat files do not correspond to the same dataset."))
    isfile(wmd_file) || throw(ArgumentError(".wmd file not found."))
    isfile(dat_file) || throw(ArgumentError(".dat file not found."))

    # Extract the graph structure from the .wmd file using a MetaGraph
    file = readdlm(wmd_file, '\n')
    V = parse(Int, split(file[1],',')[1])
    pool = MetaDiGraph(V, 0)
    for line in file[2:end]
        splitted_line = split(line, ',')
        if length(splitted_line) == 3
            # /!\ Pairs are numbered from 0 in the second part of the file
            source = parse(Int, splitted_line[1]) + 1
            destination = parse(Int, splitted_line[2]) + 1
            weight = parse(Int, splitted_line[3])
            add_edge!(pool, source, destination)
            set_prop!(pool, source, destination, :weight, weight)
        end
    end

    # Extract meta information from the .dat file
    file = readdlm(dat_file, '\n')
    for line in file[2:end]
        splitted_line = split(line, ',')
        pair = parse(Int, splitted_line[1])
        set_prop!(pool, pair, :donor, splitted_line[2])
        set_prop!(pool, pair, :patient, splitted_line[3])
        set_prop!(pool, pair, :wifep, parse(Bool, splitted_line[4]))
        set_prop!(pool, pair, :pra, parse(Float64, splitted_line[5]))
        set_prop!(pool, pair, :altruist, parse(Bool, splitted_line[7]))
    end
    Va = [v for v in vertices(pool) if get_prop(pool, v, :altruist)]
    set_prop!(pool, :nb_altruist, length(Va))

    return pool

end