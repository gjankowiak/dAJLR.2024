import DelimitedFiles

const half_N = 720
const N = 2half_N

phase_shift = 2π / 100

output_filename = "omega2_double_circle.csv"

function main()
    t = range(0, 2π, half_N + 1)[1:half_N]

    rhos = zeros(N)

    thetas = cat(t, 2π + phase_shift .+ t; dims=1)

    res = cat(rhos, thetas, zeros(3); dims=1)

    open(output_filename, "w") do io
        DelimitedFiles.writedlm(io, res)
    end
end

main()
