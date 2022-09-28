function update_coloring!(p, cs)
    active = cs[:active]
    f = cs[active][:map]
    kwargs = NamedTuple(cs[active])
    final_f = state->f(state; kwargs...)
    p[:colorfunc] = final_f
end

function show_color_settings!(color_settings)
    window = GLMakie.Screen()
    fig = Figure()
    settings_for = color_settings[][:active]
    if settings_for == 1
        make_dialog_color_depth!(window, fig, color_settings)
    elseif settings_for == 2
        make_dialog_color_lineages!(window, fig, color_settings)
    end
    resize_to_layout!(fig)
    display(window, fig)
 
    return nothing
end

function make_dialog_color_depth!(window, fig, color_settings)
    grid = fig[1,1] = GridLayout()
    Label(grid[1,1], text="Depth")
    textbox_roots = Textbox(grid[1,2], validator=Int)
    set_textbox_display!(textbox_roots, string(color_settings[][1][:depth]))

    oac_grid = grid[2,1]
    btn_apply = Button(oac_grid[1,1], label="Apply")

    on(textbox_roots.stored_string) do s
        color_settings[][1][:depth] = parse(Int, s)
    end
    on(btn_apply.clicks) do _
        notify(color_settings)
    end
end
function make_dialog_color_lineages!(window, fig, color_settings)
    grid = fig[1,1] = GridLayout()
    Label(grid[1,1], text="Roots")
    textbox_roots = Textbox(grid[1,2], validator=s->validate_css(Int, s))
    set_textbox_display!(textbox_roots, join(color_settings[][2][:roots], ','))

    oac_grid = grid[2,1]
    btn_apply = Button(oac_grid[1,1], label="Apply")

    on(textbox_roots.stored_string) do s
        color_settings[][2][:roots] = cs_string_to_vec(Int, s)
    end
    on(btn_apply.clicks) do _
        notify(color_settings)
    end
end

set_textbox_display!(tb::Textbox, v::Vector) = set_textbox_display!(tb, join(v,','))
function set_textbox_display!(tb::Textbox, s::AbstractString)
    tb.stored_string[] = s
    tb.displayed_string[] = tb.stored_string[]
end

## -- Inspector -- ##

struct TumorInspector
    fig
    tumor_plot
    slice_plot
    color_settings
    plane_settings
    radius_settings
end

function TumorInspector(state::TumorConfiguration{A}, args...; kwargs...) where A<:RealLattice{Int}
    state_obs = Observable(state)

    color_settings = Observable(Dict())
    color_menu_labels = ["Inner", "Lineages"]
    color_menu_maps = [color_depth, color_lineages]
    for (i,(l,f)) in enumerate(zip(color_menu_labels, color_menu_maps))
        color_settings[][i] = Dict()
        color_settings[][i][:map] = f
        color_settings[][i][:short] = l
        color_settings[][i][:palette] = ColorFunctions.default_palette
    end
    color_settings[][:active] = a = 1
    color_settings[][a][:depth] = 2

    plane_dir = Observable(Vec3f(0,0,1))
    plane_offset = Observable(dot(Lattices.midpointcoord(state_obs[].lattice),plane_dir[]))
    plane_settings = Dict(:dir => plane_dir, :offset => plane_offset)
    plane_origin = @lift $plane_offset*$plane_dir
    plane = @lift Lattices.Plane($plane_origin, $plane_dir)

    fig = Figure(backgroundcolor=:lightgray)
    grid_plots = fig[1,1] = GridLayout()

    grid_controls = fig[1,2] = GridLayout(;tellheight=false)
    tumor_ax = Axis3(grid_plots[1,1], aspect=:data)
    tumor_plot = tumorplot!(tumor_ax, state_obs, args...; kwargs...)
    limits_origin = (0,0,0)
    limits_lengths = Tuple(coord(state.lattice, size(state.lattice)))
    limits!(tumor_ax, Rect(limits_origin..., limits_lengths...))
    
    slice_ax = Axis(grid_plots[2,1], aspect=1)
    slice_plot = tumorplot!(slice_ax, state_obs; plane)
    
    ## Color
    Label(grid_controls[0,1], text="Coloring")

    menu_color = Menu(grid_controls[0,2:3], options=zip(color_menu_labels, eachindex(color_menu_labels)))
    btn_color_settings = Button(grid_controls[0,end+1], label="⋯")

    on(menu_color.selection) do selection
        @show selection
        if !haskey(color_settings[], selection)
            color_settings[][selection] = Dict()
        end
        cs = color_settings[][selection]
        if selection == 1
            if !haskey(cs, :depth)
                cs[:depth] = 2
            end
            cs[:map] = ColorFunctions.color_depth
        elseif selection == 2
            if !haskey(cs, :roots)
                cs[:roots] = [1]
            end
            cs[:map] = ColorFunctions.color_lineages
        end
        color_settings[][:active] = selection
        # update_coloring!(tumor_plot, color_settings)
        notify(color_settings)
    end

    on(btn_color_settings.clicks) do _
        show_color_settings!(color_settings)
    end

    on(color_settings) do cs
        @info "Color settings updated."
        update_coloring!(tumor_plot, cs)
        update_coloring!(slice_plot, cs)
    end

    ## Radius filter controls
    radius_settings = Dict(:r => tumor_plot[:radius], :lesseq => tumor_plot[:radius_eq])
    Label(grid_controls[1,1], text="Filter radius", tellheight=false)
    toggle_radius = Toggle(grid_controls[1,2])
    toggle_radius_eq_le = Toggle(grid_controls[1,3])
    Label(grid_controls[1,4], text="≤/=")
    Label(grid_controls[2,1], text="Radius")
    connect!(tumor_plot[:radius_eq], toggle_radius_eq_le.active)

    text_radius = Textbox(grid_controls[2,3], placeholder="99", validator=Float64)
    btn_radius_plus = Button(grid_controls[2,2], label="+")
    btn_radius_minus = Button(grid_controls[2,4], label="-")
    on(btn_radius_plus.clicks) do _
        r = (tumor_plot[:radius][] += 1)
        text_radius.stored_string[] = string(r)
        text_radius.displayed_string[] = string(r)
    end
    on(btn_radius_minus.clicks) do _
        r = max(0, (tumor_plot[:radius][] - 1))
        tumor_plot[:radius][] = r
        text_radius.stored_string[] = string(r)
        text_radius.displayed_string[] = string(r)
    end
    
    Label(grid_controls[3,1], text="Midpoint")
    text_mp = Textbox(grid_controls[3,2], placeholder="x,y,z", validator=s->validate_css(Float32, s,3))
    btn_reset_mp = Button(grid_controls[3,3], label="Reset")
    connect!(tumor_plot[:filter_radius], toggle_radius.active)
    on(text_radius.stored_string) do s
        @info "radius changed to $s"
        tumor_plot[:radius][] = parse(Float64, s)
    end
    on(text_mp.stored_string) do s
        @info "MP changed to $(xyz_string_to_point3(s))"
        @show typeof(xyz_string_to_point3(s))
        tumor_plot[:midpoint][] = xyz_string_to_point3(s)
    end
    on(btn_reset_mp.clicks) do _
        @info "Reset btn clicked"
        v = Lattices.midpointcoord(tumor_plot[:state][].lattice)
        #tumor_plot[:midpoint][] = v
        text_mp.stored_string[] = join(v,',')
        text_mp.displayed_string[] = text_mp.stored_string[]
    end

    ## Plane filter controls
    Label(grid_controls[4,1], text="Plane direction")
    btn_plane_x = Button(grid_controls[4,2], label="x")
    btn_plane_y = Button(grid_controls[4,3], label="y")
    btn_plane_z = Button(grid_controls[4,4], label="z")
    btns_plane = [btn_plane_x, btn_plane_y, btn_plane_z]
    btn_plane_z.clicks[] = 1
    # Label(grid_controls[5,1], text="Origin")
    Label(grid_controls[5,1], text="Offset")
    text_plane_offset = Textbox(grid_controls[5,2], stored_string=string(plane_offset[]), validator=Float32)
    btn_offset_plus = Button(grid_controls[5,3], label="+")
    btn_offset_minus = Button(grid_controls[5,4], label="-")
    on(btn_offset_plus.clicks) do _
        o = (plane_offset[] += 1)
        text_plane_offset.stored_string[] = string(o)
        text_plane_offset.displayed_string[] = string(o)
    end
    on(btn_offset_minus.clicks) do _
        o = (plane_offset[] -= 1)
        text_plane_offset.stored_string[] = string(o)
        text_plane_offset.displayed_string[] = string(o)
    end
    onany(getproperty.(btns_plane, :clicks)..., text_plane_offset.stored_string) do bx,by,bz,o
        plane_offset[] = parse(Float32, o)
        v = Vec3f(0,0,0)
        dir_label = ""
        bx%2==1 && (v += Vec3f(1,0,0); dir_label*="x" )
        by%2==1 && (v += Vec3f(0,1,0); dir_label*="y")
        bz%2==1 && (v += Vec3f(0,0,1); dir_label*="z")
        plane_dir[] = normalize(v)
        @info "New slice plane $(plane_dir[])"
        # slice_plot.
    end

    return TumorInspector(fig, tumor_plot, slice_plot, color_settings, plane_settings, radius_settings)
end
