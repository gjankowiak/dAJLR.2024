import HDF5
import DelimitedFiles
import UnicodePlots

function main()
    f = HDF5.h5open("source_data.hdf5", "r")
    X0 = read(f, "s36X")
    close(f)

    delta_θ_deg = 1 # was 2
    delta_θ = deg2rad(delta_θ_deg)

    delta_θ = 2π / 1000

    N = 720
    t = range(0, 2π, length=N+1)[1:N]

    X0[1:N] .= cos.(3*t)

    X1 = zeros(4*N + 3)
    X1[2N+1:3N] .= X0[N+1:2N]
    X1[3N+1:4N] .= X0[N+1:2N] .+ delta_θ

    X1[end-1:end] .= X0[end-1:end]

    println(UnicodePlots.lineplot(X1[1:2N]))
    println(UnicodePlots.lineplot(X1[1+2N:4N]))

    open("omega0_double_loop.csv", "w") do io
        DelimitedFiles.writedlm(io, X1)
    end
end

main()
