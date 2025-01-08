using DelimitedFiles, BenchmarkTools, CSV, DataFrames
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
set_problem = String.(readdlm("testenomes.txt"))
csv_file = open("lmlovoinfo.csv", "w")
#csv_file_benchmark = open("benchlm3_40.csv", "w")
df = DataFrame()
k = 0
for probname ∈ set_problem
  log_file = open("logclass.txt", "w")
  prob = load_problem(probname)
  #H = hcat(prob.data, zeros(length(prob.data[:, 1])))
  solved = false
  try
    xk = MinimalLevenbergMarquardt(prob.model, ones(prob.dim), prob.data, prob.dim, 1.0e-3, Matrix{Float64}(I, prob.dim, prob.dim))
    #xk = ones(prob.dim)
    s = MinimalLMLOVO(xk, prob.model, prob.data, prob.dim, prob.nout)
    #a = @benchmark MinimalLMLOVO($xk, $prob.model, $prob.data, $prob.dim, $prob.nout)
    k = k + 1
    println(k)
    row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], s[3])])#, median(a.times) / 1e9)])
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