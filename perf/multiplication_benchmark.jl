using Random, IntervalLinearAlgebra, Plots

sizes = Iterators.flatten((2:9, 10:10:90, 100:100:1000)) |> collect

times = zeros(length(sizes), 4)

seed = 42
Random.seed!(seed)

random_interval_matrix(N) = Interval.(randn(N, N) .± abs.(randn(N, N)))

for (j, mult_mode) in enumerate((:fast, :rank1, :slow))

    setmultiplication(mult_mode)
    Random.seed!(seed)
    for (i, N) in enumerate(sizes)
        A = random_interval_matrix(N)
        B = random_interval_matrix(N)
        t = @benchmark $A * $B
        times[i, j] = minimum(t.times)
        @show N
    end

end

# normal multiplication with accurate rounding mode
setmultiplication(:slow)
setrounding(Interval, :accurate)
Random.seed!(seed)
for (i, N) in enumerate(sizes)
    A = random_interval_matrix(N)
    B = random_interval_matrix(N)
    t = @benchmark $A * $B
    times[i, 4] = minimum(t.times)
    @show N
end

# back to default settings
setrounding(Interval, :tight)
setmultiplication(:fast)

# floating point multiplication
times_float = zeros(size(sizes))
for (i, N) in enumerate(sizes)
    A = randn(N, N)
    B = randn(N, N)
    t = @benchmark $A * $B
    times_float[i] = minimum(t.times)
    @show N
end

# plotting
labels = ["fast" "rank1" "slow-tight" "slow-accurate"]
plot(sizes, times/1e6; label=labels, axis=:log, m=:auto, legend=:topleft)
plot!(sizes, times_float/1e6, label="Float64", m=:auto)
xlabel!("size")
ylabel!("time [ms]")
yticks!(10.0 .^(-3:5))

ratios = times ./ times_float
labels = ["fast" "rank1" "slow-tight" "slow-accurate"]
plot(sizes, ratios; label=labels, axis=:log, m=:auto, legend=:topleft)
xlabel!("size")
ylabel!("ratio")
yticks!(10.0 .^(-3:5))

## overestimate benchmarks

radii = 10.0 .^ (-15:15)
sizes = [10, 50, 100, 200]
Nsamples = 10

min_overestimate = zeros(length(radii), length(sizes))
max_overestimate = zeros(length(radii), length(sizes))
mean_overestimate = zeros(length(radii), length(sizes))

# function random_interval_matrix(N, r)
#     Ac = randn(N, N)
#     return Ac .± (r * abs.(Ac))
# end

random_interval_matrix(N, r) = randn(N, N) .± r*abs.(randn(N, N))
Random.seed!(seed)
for (i, N) in enumerate(sizes)
    for (j, r) in enumerate(radii)

        tmp_statistics = zeros(3)
        for _ in 1:Nsamples
            A = random_interval_matrix(N, r)
            B = random_interval_matrix(N, r)

            setmultiplication(:fast)
            Cfast = A * B
            setmultiplication(:slow)
            Cslow = A * B

            ratios = (diam.(Cfast)./diam.(Cslow) .- 1) * 100
            tmp_statistics .+= (minimum(ratios), mean(ratios), maximum(ratios))
        end

        min_overestimate[j, i] = tmp_statistics[1]/Nsamples
        mean_overestimate[j, i] = tmp_statistics[2]/Nsamples
        max_overestimate[j, i] = tmp_statistics[3]/Nsamples
        println("N=$(i)/$(length(sizes)) r=$(j)/$(length(radii))")
    end
end


for (i, N) in enumerate(sizes)
    plot(radii, max_overestimate[:, i], xaxis=:log, m=:o, label="max")
    plot!(radii, min_overestimate[:, i], m=:o, label="min")
    plot!(radii, mean_overestimate[:, i], m=:o, label="mean")
    xlabel!("radius")
    ylabel!("%-overestimate")
    title!("matrix size N=$N")
    savefig("accurate_overestimate_N$(N)_random_radius.pdf")
end
