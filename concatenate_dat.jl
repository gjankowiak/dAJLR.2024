push!(LOAD_PATH, "src")

import Serialization
import Bend

function concat(params_fn_1, params_fn_2; reverse_first::Bool=false, reverse_second::Bool=false, keep_unique::Bool=false, keep_sorted::Bool=false)
    P1 = Serialization.deserialize(params_fn_1)
    P2 = Serialization.deserialize(params_fn_2)

    X_fn_1 = replace(params_fn_1, "_P.dat" => ".dat")
    X_fn_2 = replace(params_fn_2, "_P.dat" => ".dat")

    X1 = Serialization.deserialize(X_fn_1)
    X2 = Serialization.deserialize(X_fn_2)

    if reverse_first
        reverse!(P1)
        reverse!(X1)
    end

    if reverse_second
        reverse!(P2)
        reverse!(X2)
    end

    P = [P1; P2]
    X = [X1; X2]

    if keep_sorted
        idx = sortperm(P, by=x -> x.epsilon)
        P = P[idx]
        X = X[idx]
    end

    if keep_unique
        unique_idx = unique(i -> P[i].epsilon, 1:length(P))

        P = P[unique_idx]
        X = X[unique_idx]
    end

    new_X_fn = replace(params_fn_1, "_P.dat" => "_new.dat")
    new_params_fn = replace(params_fn_1, "_P.dat" => "_new_P.dat")

    r = true
    if isfile(new_X_fn)
        r = Bend.prompt_yes_no("Overwrite "*new_X_fn*" ?")
    end
    if r
        Serialization.serialize(new_X_fn, X)
        Serialization.serialize(new_params_fn, P)
        println("Wrote concatenated results to:")
        println(new_X_fn)
        println(new_params_fn)
    end

end

concat(ARGS[1], ARGS[2]; keep_sorted=true, keep_unique=false)
