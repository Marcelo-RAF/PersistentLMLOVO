using GLMakie


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
    ax = Axis3(fig[1, 1],aspect=(1,1,1), title = "Esfera ajustada a um conjunto de pontos")
    scatter!(ax, x, y, z, markersize = 15, strokewidth = 0, color = (:blue,0.7))
    mesh!(ax,S,color = (:red,0.4),transparency = true)
    # esconde os eixos
    hidespines!(ax)
    # esconde a grade
    hidedecorations!(ax)
    fig
end
