using BenchmarkTools, CSV, DataFrames
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
set_problem = String.(readdlm("testenomes.txt"))
csv_file = open("persistentinfo5.csv", "w")
#csv_file_benchmark = open("benchlm3_40.csv", "w")
df = DataFrame()
for probname ∈ set_problem
  log_file = open("logpers.txt", "w")
  solved = false
  try
    xk = MinimalLevenbergMarquardt(prob.model, ones(prob.dim), prob.data, prob.dim, 1.0e-3, Matrix{Float64}(I, prob.dim, prob.dim))
    s = MinimalLMPersistent(xk, prob.model, prob.data, prob.dim, prob.nout)
    a = @benchmark MinimalLMPersistent($xk, $prob.model, $prob.data, $prob.dim, $prob.nout)
    row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], median(a.times) / 1e9)])
    df = vcat(df, row)
    CSV.write(csv_file, df)
  catch e
    println("erro: ", e)
    solved = false
    write(log_file, "$probname\n")
  end
  close(log_file)
end
close(csv_file)

#close(csv_file_benchmark)