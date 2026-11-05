@userplot ColorScatter

@recipe function f(cs::ColorScatter;
                   normal      = nothing,
                   clims       = nothing,
                   colorbar_title = "",
                   markersize  = 8,
                   markerstrokewidth = 0.0,
                   alpha       = 1.0,
                   label       = "")
    values, locations = cs.args
    @assert length(values) == length(locations) "values 和 locations 长度必须一致"

    pts = map(r -> Vec3{Float64}(r), locations)

    # -------- 通用属性 --------
    seriestype          := :scatter
    marker_z            := values
    markersize          := markersize
    markerstrokewidth   := markerstrokewidth
    alpha               := alpha
    label               := label
    colorbar_title      := colorbar_title
    legend              := false

    # 自动颜色范围
    if clims === nothing
        m = max(maximum(abs, values), eps())
        clims := (-m, m)
    else
        clims := clims
    end

    # -------- 分支 --------
    if isnothing(normal)
        # 3D 分支
        xs = getindex.(pts, 1)
        ys = getindex.(pts, 2)
        zs = getindex.(pts, 3)
        return (xs, ys, zs)
    else
        # 投影分支
        frame = plane_frame(normal)
        pu = map(r -> dot(r, frame.e1), pts)
        pv = map(r -> dot(r, frame.e2), pts)
        aspect_ratio := :equal
        # 让坐标轴标签能反映投影方向
        xlabel := "u·r  (u = $(round.(frame.e1, digits=3)))"
        ylabel := "v·r  (v = $(round.(frame.e2, digits=3)))"
        return (pu, pv)
    end
end