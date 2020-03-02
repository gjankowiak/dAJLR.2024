push!(LOAD_PATH, "src")

import Serialization
import Bend

function concat(params_fn_1, params_fn_2)
    P1 = Serialization.deserialize(params_fn_1)
    P2 = Serialization.deserialize(params_fn_2)

    X_fn_1 = replace(params_fn_1, "_P.dat" => ".dat")
    X_fn_2 = replace(params_fn_2, "_P.dat" => ".dat")

    X1 = Serialization.deserialize(X_fn_1)
    X2 = Serialization.deserialize(X_fn_2)

    P = [P1; P2]
    idx = sortperm(P, by=x -> x.epsilon)

    P = P[idx]
    X = [X1; X2]
    X = X[idx]

    unique_idx = unique(i -> P[i].epsilon, 1:length(P))

    P = P[unique_idx]
    X = X[unique_idx]

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

concat(ARGS[1], ARGS[2])
