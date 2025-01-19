using BenchmarkTools, CSV, DataFrames
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
set_problem = String.(readdlm("testenomes.txt"))
csv_file = open("plmtime.csv", "w")
#csv_file_benchmark = open("benchlm3_40.csv", "w")
df = DataFrame()
k = 0
#benchmark_df = DataFrame()
for probname ∈ set_problem
  log_file = open("logpers.txt", "w")
  prob = load_problem(probname)
  solved = false
  try
    xk = MinimalLevenbergMarquardt(prob.model, ones(prob.dim), prob.data, prob.dim, 1.0e-3, Matrix{Float64}(I, prob.dim, prob.dim))
    s = MinimalLMPersistent(xk, prob.model, prob.data, prob.dim, prob.nout)
    a = @benchmark MinimalLMPersistent($xk, $prob.model, $prob.data, $prob.dim, $prob.nout)
    #ndif = norm(prob.solution - s[1])
    k = k + 1
    println(k)
    row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s, median(a.times) / 1e9)])
    df = vcat(df, row)
    #benchmark_row = DataFrame([(probname, prob.npts, prob.nout, minimum(a.times) / 1e9, median(a.times) / 1e9, maximum(a.times) / 1e9)])
    #benchmark_df = vcat(benchmark_df, benchmark_row)
    #df = DataFrame(solution_LOVOCGA = [s], prob_solution = [prob.solution])
    CSV.write(csv_file, df)
    #CSV.write(csv_file_benchmark, benchmark_df)
  catch e
    println("erro: ", e)
    solved = false
    write(log_file, "$probname\n")
  end
  close(log_file)
end
close(csv_file)

#close(csv_file_benchmark)