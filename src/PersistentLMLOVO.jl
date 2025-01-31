module PersistentLMLOVO

using ForwardDiff, DelimitedFiles, LinearAlgebra

export load_problem, build_problem, visualize, FitProbType, FitOutputType, MinimalLevenbergMarquardt, MinimalLMPersistent, MinimalLMLOVO

import Base.show



struct FitProbType
  name::String
  data::Array{Float64,2}
  npts::Int
  nout::Int
  model::Function
  dim::Int
  cluster::Bool
  noise::Bool
  solution::Array{Float64,1}
  description::String
  sandbox::Vector{Vector{Any}}
  function FitProbType(_name, _data, _npts, _nout, _model, _dim, _cluster, _noise, _solution, _description)
    new(_name, _data, _npts, _nout, _model, _dim, _cluster, _noise, _solution, _description, [["NONE"]])
  end
end

struct FitOutputType
  status::Bool
  solution::Vector{Float64}
  niter::Int
  minimum::Float64
  feval::Int
end

function load_problem(filename::String)
  prob_matrix = readdlm(filename, ':')
  (m, n) = size(prob_matrix)

  if m == 10
    return FitProbType(prob_matrix[1, 2], eval(Meta.parse(prob_matrix[2, 2])), prob_matrix[3, 2], prob_matrix[4, 2], eval(Meta.parse(prob_matrix[5, 2])), prob_matrix[6, 2], prob_matrix[7, 2], prob_matrix[8, 2], eval(Meta.parse(prob_matrix[9, 2])), prob_matrix[10, 2])
  elseif m == 11
    return FitProbType(prob_matrix[1, 2], eval(Meta.parse(prob_matrix[2, 2])), prob_matrix[3, 2], prob_matrix[4, 2], eval(Meta.parse(prob_matrix[5, 2])), prob_matrix[6, 2], prob_matrix[7, 2], prob_matrix[8, 2], eval(Meta.parse(prob_matrix[9, 2])), prob_matrix[10, 2], eval(Meta.parse(prob_matrix[11, 2])))
  else
    error("No type identified!!")
  end
end

function build_problem(probtype::String, limit::Vector{Float64}, params::Vector{Float64})
  if probtype == "quadratic"
    println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
    println("params need to be setup as [coefs of A, b and c,size_prob,npts,nout]")
    nout = Int(params[end])
    npts = Int(params[end-1])
    size_prob = Int(params[end-2])
    dim = length(params[1:end-3])
    n = size_prob
    mA = "["
    sA = [" " for i = 1:n, j = 1:n]
    A = zeros(n, n)
    b = zeros(n)
    c = 0.0
    k = 1
    for i = 1:n
      for j = i:n
        A[i, j] = params[k]
        A[j, i] = A[i, j]
        sA[i, j] = " t[$(k)] "
        sA[j, i] = sA[i, j]
        k += 1
      end
    end
    for i = 1:n
      for j = 1:n
        mA = mA * sA[i, j]
        if j == n
          if i == n
            mA = mA * "]"
          else
            mA = mA * ";"
          end
        end
      end
    end
    mb = "["
    for i = 1:n-1
      mb = mb * " t[$(k)],"
      b[i] = params[k]
      k += 1
    end
    mb = mb * "t[$(k)]]"
    b[n] = params[k]
    k += 1
    mc = " t[$(k)]"
    c = params[k]
    m(x) = x' * A * x + b' * x + c
    #display(A)
    model = "(x,t) -> x'*$(mA)*x + $(mb)'*x + $(mc)"
    r = (limit[2] - limit[1]) / (npts - 1)
    x = rand(limit[1]:r:limit[2], npts)
    for i = 1:n-1
      x = [x rand(limit[1]:r:limit[2], npts)]
    end
    #display(x)
    #global coefs = params[1:end-3]
    y = zeros(npts)
    #
    for i = 1:npts
      y[i] = m(x[i, :])
    end

    k = 1
    iout = []
    while k <= nout
      i = rand([1:npts;])
      if i ∉ iout
        push!(iout, i)
        k = k + 1
      end
    end
    for k = 1:nout
      x[iout[k], :] = x[iout[k], :] + randn(size_prob)
      y[iout[k]] = y[iout[k]] + randn()
    end

    ## A(t) = begin
    #         mA = zeros(size_prob,size_prob)
    #         k = 1
    #         for i=1:size_prob
    #             for j=i:size_prob
    #                 mA[i,j] = t[k]
    #                 mA[j,i] = t[k]
    #                 k += 1
    #             end
    #         end
    #         return mA
    #     end
    # b(t) = begin
    #     vb = zeros(size_prob)
    #     k = Int(((size_prob^2+size_prob)/2)+1)
    #         for i=1:size_prob
    #             vb[i] = t[k]
    #             k += 1
    #         end
    #         return vb
    #     end
    #     c(t) = t[Int(((size_prob^2+size_prob)/2)+size_prob+1)]
    # 
    # model = "(x,t) -> x'*A(t)*x + b(t)'*x + c(t)"


    FileMatrix = ["name :" "Quadratic"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" model; "dim :" dim; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:end-3]]; "description :" "none"]

    open("quadratic_$(params[1])_$(params[end-3])_$(size_prob)_$(dim)_$(npts)_$(nout).csv", "w") do io
      writedlm(io, FileMatrix)
    end
  end
  if probtype == "cubic"
    println("params need to be setup as [vector, npts, nout]")
    p = [params[1], params[2], params[3], params[4]]
    npts = Int(params[5])
    nout = Int(params[6])
    t = range(-5.0, stop=5.0, length=npts)
    x = zeros(npts)
    y = zeros(npts)
    sgn = sign(randn())
    for i = 1:npts
      x[i] = t[i]
      y[i] = p[1] * x[i]^3 + p[2] * x[i]^2 + p[3] * x[i] + p[4] #+ (1.0 + 2 * rand()) * 7.0 * sgn
    end
    k = 1
    iout = []
    while k <= nout
      i = rand([1:npts;])
      if i ∉ iout
        push!(iout, i)
        k = k + 1
      end
    end
    for k = 1:nout
      y[iout[k]] = p[1] * x[iout[k]]^3 + p[2] * x[iout[k]]^2 + p[3] * x[iout[k]] + p[4] + randn() * 200 #rand([0.25*r:0.1*(r); (1 + 0.25) * r])
    end

    FileMatrix = ["name :" "cubic"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> t[1]*x[1]^3 + t[2]*x[1]^2 + t[3]*x[1] + t[4] "; "dim :" 4; "cluster :" "false"; "noise :" "true"; "solution :" [push!(p)]; "description :" "type: cubic model with noise and outliers"]

    open("cubic_$(p[1])_$(p[2])_$(p[3])_$(nout).csv", "w") do io
      writedlm(io, FileMatrix)
    end
  end
  if probtype == "line2d"
    println("params need to be setup as [vector, npts, nout]")
    p = [params[1], params[2]]
    npts = Int(params[3])
    nout = Int(params[4])
    t = range(-15.0, stop=15.0, length=npts)
    x = zeros(npts)
    y = zeros(npts)
    sgn = sign(randn())
    ruid = randn(1, npts)
    for i = 1:npts
      x[i] = t[i]
      y[i] = p[1] * x[i] + p[2] #+ ruid[1,i] #+ (1.0 + 2 * rand()) * 7.0 * sgn
    end
    k = 1
    iout = []
    while k <= nout
      i = rand([1:npts;])
      if i ∉ iout
        push!(iout, i)
        k = k + 1
      end
    end
    for k = 1:nout
      y[iout[k]] = p[1] * x[iout[k]] + p[2] + rand() * 50 #rand([0.25*r:0.1*(r); (1 + 0.25) * r])
    end

    FileMatrix = ["name :" "line2d"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> t[1]*x + t[2]"; "dim :" 2; "cluster :" "false"; "noise :" "false"; "solution :" [push!(p)]; "description :" "type2: line model"]

    open("line2d_$(p[1])_$(p[2])_$(nout).csv", "w") do io
      writedlm(io, FileMatrix)
    end
  end
  if probtype == "sphere2D"
    println("params need to be setup as [center,radious,npts,nout]")
    c = [params[1], params[2]]
    r = params[3]
    npts = Int(params[4])
    x = zeros(npts)
    y = zeros(npts)
    z = zeros(npts)
    ruid = randn(2, npts)
    θ = range(0, stop=2π, length=npts) #Int(ceil(npts/2)))
    #θ2 = range(5*π/4, stop=7*π/4, length= 2*npts)#Int(ceil(npts/2)))
    for k = 1:npts
      x[k] = c[1] + r * cos(θ[k]) #+ ruid[1, k]
      y[k] = c[2] + r * sin(θ[k]) #+ ruid[2, k]
    end
    nout = Int(params[5])
    k = 1
    iout = []
    while k <= nout
      i = rand([1:npts;])
      if i ∉ iout
        push!(iout, i)
        k = k + 1
      end
    end
    #dx = rand() * 0.1 * r # deslocamento aleatório em x
    #dy = rand() * 0.1 * r # deslocamento aleatório em y
    for k = 1:nout
      x[iout[k]] = x[iout[k]] + rand([-(0.5)*r:0.1:(0.5)*r;])
      y[iout[k]] = y[iout[k]] + rand([-(0.5)*r:0.1:(0.5)*r;])   #rand([0.25*r:0.1*(r); (1 + 0.25) * r])
    end
    FileMatrix = ["name :" "sphere2D"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 - t[3]^2"; "dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [push!(c, r)]; "description :" "type3: test sphere2d with noise and outliers"]

    open("sphere2D_$(c[1])_$(c[2])_$(c[3])_$(nout).csv", "w") do io
      writedlm(io, FileMatrix)
    end
  end
  if probtype == "sphere3D"
    println("params need to be setup as [center,radious,npts,nout]")
    c = [params[1], params[2], params[3]]
    r = params[4]
    npts = Int(params[5])
    x = zeros(npts)
    y = zeros(npts)
    z = zeros(npts)
    w = zeros(npts)
    θ = range(0, stop=2π, length=npts)
    φ = range(0, stop=π, length=npts)
    #φ2 = range(5π/6, stop=π, length=npts)
    rd = randn(3, npts)
    #if iseven(npts)==false
    #    l = Int(round(npts/2))
    #   h = Int(ceil(npts/2))
    #else
    #   l = Int(npts/2)
    #  h = Int(npts/2) + 1
    # end
    for k = 1:npts #forma de espiral - ao criar outro forma, se obtem metade dos circulos máximos
      x[k] = c[1] + r * cos(θ[k]) * sin(φ[k]) + rd[1, k]
      y[k] = c[2] + r * sin(θ[k]) * sin(φ[k]) + rd[2, k]
      z[k] = c[3] + r * cos(φ[k]) + rd[3, k]
    end
    nout = Int(params[6])
    k = 1
    iout = []
    while k <= nout
      i = rand([1:npts;])
      if i ∉ iout
        push!(iout, i)
        k = k + 1
      end
    end
    for k = 1:nout
      x[iout[k]] = x[iout[k]] + rand([-(1 + 0.15)*r:0.1:(1+0.15)*r;])
      y[iout[k]] = y[iout[k]] + rand([-(1 + 0.15)*r:0.1:(1+0.15)*r;])
      z[iout[k]] = z[iout[k]] + rand([-(1 + 0.15)*r:0.1:(1+0.15)*r;])
    end
    FileMatrix = ["name :" "sphere3D"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 +(x[3]-t[3])^2 - t[4]^2"; "dim :" 4; "cluster :" "false"; "noise :" "true"; "solution :" [push!(c, r)]; "description :" [[c, c]]]

    open("sphere3D_$(c[1])_$(c[2])_$(c[3])_$(c[4])_$(nout).csv", "w") do io #o que essa linha faz exatamente?
      writedlm(io, FileMatrix)
    end
  end
end

function diferential(model, θ, data, dim)

  cl2(θ) = model(x, θ)

  grad_model!(h, x_, θ) = begin

    global x = x_

    return ForwardDiff.gradient(h, θ)
  end
  (m, n) = size(data)
  J = zeros(m, Int(dim))
  for i = 1:m
    if n == 2
      inputs = data[i, 1]
    else
      inputs = data[i, 1:end-1]
    end
    J[i, :] = grad_model!(cl2, inputs, θ)
  end
  return J
end

function func(x, model, data)
  (m, n) = size(data)
  F = zeros(m)
  for i = 1:m
    if n == 2
      inputs = data[i, 1]
    else
      inputs = data[i, 1:end-1]
    end
    F[i] = model(inputs, x) - data[i, end]
  end
  return F
end

#Para rodar Levenberg insira o modelo, o chute inicial, os dados disponíveis do modelo e a dimensão do modelo ---- > a função func devolve o vetor do modelo aplicado nos pontos data e a função diferential gera a matriz jacobiana do modelo


function MinimalLevenbergMarquardt(model, x, data, dim, ε, Id, λ_min=0.7)
  F = func(x, model, data)
  J = diferential(model, x, data, dim)
  xn = zeros(dim)
  λ = 1.0
  k = 1
  JtF = -(J') * F
  while norm(JtF, 2) > ε && k <= 100
    d = (J' * J + λ * Id) \ (JtF)
    xn .= x .+ d
    Fn = func(xn, model, data)
    if norm(Fn, 2) < norm(F, 2)
      x .= xn
      if λ < λ_min
        λ = λ_min
      else
        λ = λ / 2.0
      end
      F .= Fn
      J = diferential(model, x, data, dim)
      JtF = -(J') * F
    else
      λ = 2.0 * λ
    end
    k = k + 1
  end
  return x
end

function MinimalLMPersistent(xk, model, data, dim, nout, ε=1.0e-4)
  ordres = sort_funcion_res(xk, model, data, nout)
  antres = 0.0
  #k = 1
  #kk = 0
  Id = Matrix{Float64}(I, dim, dim)
  while abs(ordres[2] - antres) > ε
    antres = ordres[2]
    xk = MinimalLevenbergMarquardt(model, xk, ordres[1], dim, ε, Id)
    #kk = kk + xk[2]
    ordres = sort_funcion_res(xk, model, data, nout)
   # k = k + 1
  end
  return xk, ordres[2]
end

function MinimalLMLOVO(xk, model, data, dim, nout, ε=1.0e-4, MAXIT=100)
  ordres = sort_funcion_res(xk, model, data, nout)
  R = func(xk, model, ordres[1])
  J = diferential(model, xk, ordres[1], dim)
  Id = Matrix{Float64}(I, dim, dim)
  λ_min = 0.7
  λ = 1.0
  k = 0
  JtR = -J' * R
  xn = zeros(dim)
  while norm(JtR, 2) > ε && k < MAXIT #&& newdata[2] > 10e-6
    d = (J' * J + λ * Id) \ (JtR)
    xn .= xk .+ d
    Rn = func(xn, model, ordres[1])
    if norm(Rn, 2) < norm(R, 2)
      xk .= xn
      if λ < λ_min
        λ = λ_min
      else
        λ = λ / 2.0
      end
      #R .= Rn 
      ordres = sort_funcion_res(xk, model, data, nout)
      R = func(xk, model, ordres[1])
      J = diferential(model, xk, ordres[1], dim)
      JtR = -J' * R
    else
      λ = 2.0 * λ
    end
    k = k + 1
  end
  return xk, ordres[2]
end


function sort_funcion_res(x, model, data, nout)
  P = data
  (n, m) = size(data)
  v = zeros(n)
  for i = 1:n
    if m == 2
      inputs = data[i, 1]
    else
      inputs = data[i, 1:end-1]
    end
    v[i] = (model(inputs, x) - data[i, end])^2
  end
  indtrust = [1:n;]
  for i = 1:n-nout+1
    for j = i+1:n
      if v[i] > v[j]
        aux = v[j]
        v[j] = v[i]
        v[i] = aux
        aux2 = indtrust[j]
        indtrust[j] = indtrust[i]
        indtrust[i] = aux2
      end
    end
  end
  #    println(indtrust[n-nout+1:n])
  return P[indtrust[1:n-nout], :], sum(v[1:n-nout])
end


function show(io::IO, fout::FitOutputType)

  print(io, "  ▶ Output ◀ \n")
  if Bool(fout.status) == true
    print(io, "  ↳ Status (.status) = Convergent \n")
  else
    print(io, "  ↳ Status (.status) = Divergent \n")
  end
  print(io, "  ↳ Solution (.solution) = $(fout.solution) \n")
  print(io, "  ↳ Number of iterations (.niter) = $(fout.niter) \n")
  print(io, "  ↳ Minimum (.minimum) = $(fout.minimum) \n")
  print(io, "  ↳ Number of function calls (.feval) = $(fout.feval) \n")
end

end