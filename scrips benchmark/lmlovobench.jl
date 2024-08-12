using DelimitedFiles, BenchmarkTools, CSV, DataFrames

set_problem = String.(readdlm("testenomes.txt"))
csv_file = open("artigolev.csv", "w")
#csv_file_benchmark = open("benchlm3_40.csv", "w")
df = DataFrame()
k = 0
for probname ∈ set_problem
  log_file = open("logclass.txt", "w")
  prob = load_problem(probname)
  #H = hcat(prob.data, zeros(length(prob.data[:, 1])))
  solved = false
  try
    xk = Levenberg(prob.model, ones(prob.dim), prob.data, prob.dim, 1.0e-4)
    #xk = ones(prob.dim)
    s = LMLOVO(xk[1], prob.model, prob.data, prob.dim, prob.nout)
    a = @benchmark LMLOVO($xk[1], $prob.model, $prob.data, $prob.dim, $prob.nout)
    k = k + 1
    println(k)
    row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], s[3], s[4], median(a.times) / 1e9)])
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