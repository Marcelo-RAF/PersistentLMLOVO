using GLMakie, Colors 


"""Essa função serve para gerar dados próximos a uma Esfera. Primeiro, geramos a esfera (responda de um ajuste) e depois geramos os dados. Isso é só para exemplificar
julia> S = Sphere(Point3f([1.0,2.0,3.0]),10.0)
julia> x,y,z = points_near_sphere(S,20)
"""
function points_near_sphere(S::Sphere{Float32},n::Int)
    θ = range(0, 2π, length=n)
    ϕ = range(0, π, length=n)
    x = zeros(n^2)
    y = zeros(n^2)
    z = zeros(n^2)
    k = 0
    for i=1:n 
        for j=1:n 
            k += 1
            d = randn()
            x[k] = S.center[1] + S.r * sin(ϕ[j]) * cos(θ[i]) + d
            y[k] = S.center[2] + S.r * sin(ϕ[j]) * sin(θ[i]) + d 
            z[k] = S.center[3] + S.r * cos(ϕ[j]) + d
        end
    end
    return x, y, z
end
"""
Essa função faz o plot. 
julia> build_sphere_plot(S,x,y,z)
"""
function build_sphere_plot(S::Sphere{Float32},x::Vector{Float64},y::Vector{Float64},z::Vector{Float64})
    fig = Figure()
    ax = Axis3(fig[1, 1],aspect=(1,1,1))#, title = "Esfera ajustada a um conjunto de pontos")
    scatter!(ax, x, y, z, markersize = 15, strokewidth = 0, color = (:blue,0.7))
    mesh!(ax,S,color = (:red,0.4),transparency = true)
    # esconde os eixos
    hidespines!(ax)
    # esconde a grade
    hidedecorations!(ax)
    fig
end
"""Esta função gera os pontos que estão sobre um plano. Note que os limites de rx e ry desta função devem ser maiores do que a points_near_plane para que o plano ultrapasse os limites dos dados e seja visualizado. Preciso passar os limites de rx e ry para parâmetros"""
function points_in_plane(normal::Vector{Float64},d::Float64,n::Int)
    f(xy) = (-normal[1]*xy[1]-normal[2]*xy[2]-d)/normal[3]
    rx = [-10.0:(20.0/n):10.0;]
    ry = [-10.0:(20.0/n):10.0;]
    xy = [[i,j] for i ∈ rx for j ∈ ry]
    z = f.(xy)
    x = zeros(length(xy))
    y = zeros(length(xy))
    for i=1:length(xy)
        x[i] = xy[i][1]
        y[i] = xy[i][2]
    end
    return x,y,z
end

"""Essa função gera pontos que são os peturbados... serve apenas para testes de gráficos"""
function points_near_plane(normal::Vector{Float64},d::Float64,n::Int)
    f(xy) = (-normal[1]*xy[1]-normal[2]*xy[2]-d)/normal[3]
    rx = [-8.0:(20.0/n):8.0;]
    ry = [-8.0:(20.0/n):8.0;]
    xy = [[i,j] for i ∈ rx for j ∈ ry]
    z = f.(xy)
    x = zeros(length(xy))
    y = zeros(length(xy))
    for i=1:length(xy)
        x[i] = xy[i][1]
        y[i] = xy[i][2]
    end
    return x+0.5*randn(length(x)),y+0.5*randn(length(y)),z+0.5*randn(length(z))
end

"""Funciona assim: 
julia> x,y,z = points_in_plane([1,1,1.0],2.0,40)

julia> xn,yn,zn = points_near_plane([1,1,1.0],2.0,20)

julia> build_plane_plot(x,y,z,xn,yn,zn) # e pronto!

"""
function build_plane_plot(x::Vector{Float64},y::Vector{Float64},z::Vector{Float64},xn::Vector{Float64},yn::Vector{Float64},zn::Vector{Float64}) 
    fig = Figure()
    ax = Axis3(fig[1, 1],aspect=(1,1,1), title = "Plano ajustado a um conjunto de pontos")
    scatter!(ax, xn, yn, zn, markersize = 10, strokewidth = 0, color = (:blue,0.7))
    surface!(ax,x,y,z,  color=fill(RGBA(1.,0.,0.,0.5),100,100), transparency=true)
    # esconde os eixos
    hidespines!(ax)
    # esconde a grade
    hidedecorations!(ax)
    fig
end


function points_in_circle(u::Vector{Float64},v::Vector{Float64},c::Vector{Float64},r::Float64,n_points::Int)
    θ = LinRange(0, 2π, n_points)
    x = [r * cos(θ_i) * u[1] + r * sin(θ_i) * v[1] + c[1] for θ_i in θ]
    y = [r * cos(θ_i) * u[2] + r * sin(θ_i) * v[2] + c[2] for θ_i in θ]
    z = [r * cos(θ_i) * u[3] + r * sin(θ_i) * v[3] + c[3] for θ_i in θ]
    return x,y,z 
end

function points_near_circle(x::Vector{Float64},y::Vector{Float64},z::Vector{Float64})
    v = 0.5*randn(length(x))
    return x+v,y+v,z+v
end 

""" Funciona assim:
julia> x,y,z = points_in_circle([1.0,0,0.0],[0,1,0.0],[0.0,0.0,0.0],6.0,100)
julia> xn,yn,zn = points_near_circle(x,y,z)
julia> build_circle_plot(x,y,z,xn,yn,zn)
"""
function build_circle_plot(x::Vector{Float64},y::Vector{Float64},z::Vector{Float64},xn::Vector{Float64},yn::Vector{Float64},zn::Vector{Float64})
    fig = Figure()  
    ax = Axis3(fig[1, 1],aspect=(1,1,1), title = "Círculo ajustado a um conjunto de pontos")
    limits!(ax, -10, 10, -10, 10, -10, 10)
    scatter!(ax, xn, yn, zn, markersize = 10, strokewidth = 0, color = (:blue,0.7))
    lines!(ax,x,y,z,color = :red, linewidth = 3)
    hidespines!(ax)
    hidedecorations!(ax)
    fig
end

function points_in_line(P::Vector{Float64},u::Vector{Float64},n_points::Int64)
    λ = LinRange(-2, 2, n_points)
    x = [P[1]+λ_i*u[1] for λ_i ∈ λ]    
    y = [P[2]+λ_i*u[2] for λ_i ∈ λ]
    z = [P[3]+λ_i*u[3] for λ_i ∈ λ]
    return x,y,z 
end

function points_near_line(x::Vector{Float64},y::Vector{Float64},z::Vector{Float64})
    v = 0.1*randn(length(x))
    return x+v,y+v,z+v
end 

""" Funciona assim:
julia> x,y,z = points_in_line([1.0,0,0.0],[0,1,0.0],[0.0,0.0,0.0],6.0,100)
julia> xn,yn,zn = points_near_line(x,y,z)
julia> build_line_plot(x,y,z,xn,yn,zn)
"""
function build_line_plot(x::Vector{Float64},y::Vector{Float64},z::Vector{Float64},xn::Vector{Float64},yn::Vector{Float64},zn::Vector{Float64})
    fig = Figure()  
    ax = Axis3(fig[1, 1],aspect=(1,1,1), title = "Reta ajustada a um conjunto de pontos")
    limits!(ax, -3, 3, -3, 3, -3, 3)
    scatter!(ax, xn, yn, zn, markersize = 10, strokewidth = 0, color = (:blue,0.7))
    lines!(ax,x,y,z,color = :red, linewidth = 3)
    hidespines!(ax)
    hidedecorations!(ax)
    fig
end