# Place-name labels: reading, styling and two layouts for "outer" labels
# (glaciers, ice shelves, settlements, fjords), which are moved off the ice
# and connected to their anchor point by a leader line:
#   :coastal – NASA/SVS style: pushed along the local coast normal to just
#              beyond the ice edge, then slid along the coast to avoid overlaps.
#   :columns – Rignot & Mouginot (2012) style: stacked in left/right margin
#              columns sorted by anchor y.
# All other label types are drawn at their anchor.

struct PlaceLabel
    name::String           # as shown: "name (historical)" where a historical name is given
    lon::Float64
    lat::Float64
    type::String
    tier::Int              # 1 major, 2 detail, 3 gazetteer; 0 = added on request
    source::String
    alt::Vector{String}    # other names (for searching)
end

"""
Read a label CSV (columns name, lat, lon, type, source, tier and optionally alt,
|-separated, and historical). A label with a historical name (Greenland: the
Danish or English name) is shown as "name (historical)"; both names are also
searchable on their own.
"""
function read_labels(path)
    d, h = readdlm(path, ',', Any; header=true, quotes=true)
    col(n, default) = (i = findfirst(==(n), vec(h)); i === nothing ? fill(default, size(d, 1)) : d[:, i])
    str(v) = strip(string(v))
    function label(n, lat, lon, t, src, tier, alt, hist)
        n, hist = str(n), str(hist)
        alt = filter(!isempty, split(str(alt), "|"))
        isempty(hist) && return PlaceLabel(n, Float64(lon), Float64(lat), str(t), Int(tier), str(src), alt)
        return PlaceLabel("$n ($hist)", Float64(lon), Float64(lat), str(t), Int(tier), str(src), unique([n; hist; alt]))
    end
    return [label(r...) for r in zip(col("name", ""), col("lat", 0), col("lon", 0), col("type", ""), col("source", ""),
                                     col("tier", 1), col("alt", ""), col("historical", ""))]
end

"Name normalised for matching: case-folded, accents removed."
normname(s) = Unicode.normalize(String(s); casefold=true, stripmark=true)

matches(l::PlaceLabel, q) = any(occursin(q, normname(n)) for n in vcat(l.name, l.alt))

"Print the labels whose name (or another name) contains `pattern`."
function list_names(labs, pattern="")
    q = normname(pattern)
    sel = sort(filter(l -> matches(l, q), labs); by=l -> (normname(l.name), l.tier))
    println(rpad("name", 40), rpad("type", 10), rpad("tier", 6), rpad("lat", 10), rpad("lon", 10), "source / other names")
    for l in sel
        other = isempty(l.alt) ? "" : " / " * join(l.alt, ", ")
        println(rpad(l.name, 40), rpad(l.type, 10), rpad(l.tier, 6), rpad(round(l.lat, digits=2), 10),
                rpad(round(l.lon, digits=2), 10), l.source, other)
    end
    println(length(sel), " names (tier 1-2: label CSV, 3: gazetteer; add any of them with add=<name>)")
end

"Great-circle distance (km) between two labels."
function distance_km(a::PlaceLabel, b::PlaceLabel)
    φ1, φ2, Δλ = deg2rad(a.lat), deg2rad(b.lat), deg2rad(b.lon - a.lon)
    return 6371.0 * acos(clamp(sin(φ1)*sin(φ2) + cos(φ1)*cos(φ2)*cos(Δλ), -1, 1))
end

"Levenshtein edit distance."
function editdistance(a, b)
    a, b = collect(a), collect(b)
    d = collect(0:length(b))
    for i in eachindex(a)
        prev, d[1] = d[1], i
        for j in eachindex(b)
            prev, d[j+1] = d[j+1], min(d[j+1] + 1, d[j] + 1, prev + (a[i] == b[j] ? 0 : 1))
        end
    end
    return d[end]
end

"""
Labels with the names in `add` (from `labs` or else the gazetteer `gaz`, matched
on any of their names, ignoring case and accents) set to tier 0, so they are
always drawn. Unknown names are reported with the closest matches.
"""
function add_names(labs, gaz, add)
    labs = copy(labs)
    tier0(l) = PlaceLabel(l.name, l.lon, l.lat, l.type, 0, l.source, l.alt)
    for name in add
        q = normname(name)
        exact(l) = any(normname(n) == q for n in vcat(l.name, l.alt))
        k = findfirst(exact, labs)
        if k !== nothing
            labs[k] = tier0(labs[k])
        elseif (g = findfirst(exact, gaz)) !== nothing
            near = [l.name for l in labs if l.type == gaz[g].type && distance_km(l, gaz[g]) < 30]
            isempty(near) || @warn "$(gaz[g].name) may be the same feature as $(join(near, ", ")) of the label CSV"
            push!(labs, tier0(gaz[g]))
        else
            pool = vcat(labs, gaz)
            dist = [minimum(occursin(q, normname(n)) ? 0 : editdistance(q, normname(n)) for n in vcat(l.name, l.alt))
                    for l in pool]
            close = unique(l.name for l in pool[sortperm(dist)[1:min(5, end)]])
            @warn "name not found: $name. Closest: $(join(close, "; "))"
        end
    end
    return labs
end

# Base text style per label type (sizes in figure units, scaled by `scale`)
const LSTYLE = Dict(
    "region"     => (size=20, font=:bold,    color=:gray30,  case=uppercase, outer=false, marker=false),
    "mountains"  => (size=16, font=:italic,  color=:sienna4, case=identity,  outer=false, marker=false),
    "basin"      => (size=15, font=:italic,  color=:gray35,  case=identity,  outer=false, marker=false),
    "dome"       => (size=14, font=:regular, color=:black,   case=identity,  outer=false, marker=true),
    "station"    => (size=13, font=:regular, color=:black,   case=identity,  outer=false, marker=true),
    "icecore"    => (size=13, font=:regular, color=:darkred, case=identity,  outer=false, marker=true),
    "lake"       => (size=13, font=:italic,  color=:navy,    case=identity,  outer=false, marker=true),
    "island"     => (size=13, font=:regular, color=:gray20,  case=identity,  outer=false, marker=false),
    "sea"        => (size=22, font=:italic,  color=:white,   case=identity,  outer=false, marker=false),
    "settlement" => (size=13, font=:regular, color=:black,   case=identity,  outer=true,  marker=true),
    "glacier"    => (size=13, font=:regular, color=:black,   case=identity,  outer=true,  marker=false),
    "iceshelf"   => (size=13, font=:bold,    color=:black,   case=identity,  outer=true,  marker=false),
    "fjord"      => (size=13, font=:italic,  color=:navy,    case=identity,  outer=true,  marker=false),
    "site"       => (size=12, font=:italic,  color=:black,   case=identity,  outer=false, marker=true),
)

"Text with a white halo drawn as a separate layer underneath (keeps glyphs crisp)."
function halotext!(ax, x, y; text, fontsize, halo=true, halowidth=0.25fontsize, halocolor=(:white, 0.8), kw...)
    halo && text!(ax, x, y; text, fontsize, color=(:white, 0.0), strokecolor=halocolor, strokewidth=halowidth, kw...)
    text!(ax, x, y; text, fontsize, kw...)
end

"""
Place area labels (regions, seas, mountains, large ice shelves): each is kept at
its anchor if free, otherwise moved to the first free spot on rings of up to
`maxshift[k]` text heights around it (or the least-overlapping one). `idx` gives
the placement order; `obstacles` are boxes to keep clear.
"""
function place_area!(pos, wh, idx, obstacles; maxshift=fill(3.0, length(pos)))
    placed = Tuple[obstacles...]
    for k in idx
        h = wh[k][2]; p0 = pos[k]
        cands = vcat([(0.0, 0.0)], [(r*h*cos(a), r*h*sin(a)) for r in 0.5:0.5:maxshift[k] for a in (0:11) .* (π/6)])
        best, bestc = (0.0, 0.0), Inf
        for d in cands
            b = boxat(p0 .+ d, wh[k], :center)
            c = sum((overlap_area(b, o) for o in placed); init=0.0)
            c < bestc && ((best, bestc) = (d, c))
            c == 0 && break
        end
        pos[k] = p0 .+ best
        push!(placed, boxat(pos[k], wh[k], :center))
    end
end

"Warn about every pair of overlapping label boxes (and labels hitting `obstacles`)."
function report_overlaps(bx, names, obstacles=Tuple[])
    n = 0
    for a in eachindex(bx), b in a+1:length(bx)
        overlap(bx[a], bx[b]) && (n += 1; @warn "labels overlap: $(names[a]) / $(names[b])")
    end
    for a in eachindex(bx), o in obstacles
        overlap(bx[a], o) && (n += 1; @warn "label overlaps an obstacle: $(names[a])")
    end
    n == 0 && println("labels: no overlaps")
    return n
end

"Point on box b (padded) where a leader line from anchor p should end."
function leader_attach(p, b, pad)
    p[1] > b[2] && return (b[2] + pad, (b[3] + b[4])/2)
    p[1] < b[1] && return (b[1] - pad, (b[3] + b[4])/2)
    return (p[1], p[2] > b[4] ? b[4] + pad : b[3] - pad)
end

"Leader line with a light casing, visible on both dark ocean and white ice."
function leader!(ax, xs, ys; lw, color=(:black, 0.75), casing=(:white, 0.7))
    lines!(ax, xs, ys; color=casing, linewidth=3lw)
    lines!(ax, xs, ys; color, linewidth=lw)
end

"""
Text box (width, height) in km of a text of one or more lines, measured with
Makie's text layout for a theme font (Symbol) or a font file.
"""
function textbox(txt, size, kmpp; font=:regular)
    ft = Makie.to_font(font isa Symbol ? Makie.theme(:fonts)[font][] : font)
    w = maximum(Makie.widths(Makie.text_bb(ln, ft, Float32(size)))[1] for ln in split(txt, '\n'))
    return (w*kmpp, 1.15*size*kmpp*(count(==('\n'), txt) + 1))
end

"Box (xmin, xmax, ymin, ymax) for text at q with horizontal alignment ha."
function boxat(q, wh, ha)
    w, h = wh
    x0 = ha === :left ? q[1] : (ha === :right ? q[1] - w : q[1] - w/2)
    return (x0, x0 + w, q[2] - h/2, q[2] + h/2)
end

"""
Boxes (xmin, xmax, ymin, ymax) of half-width `hw` km every `step` km along a
NaN-separated polyline, as label obstacles for lines.
"""
function polyline_boxes(pts; step=30.0, hw=10.0)
    boxes = Tuple[]
    for k in 1:length(pts)-1
        a, b = pts[k], pts[k+1]
        (any(isnan, a) || any(isnan, b)) && continue
        n = max(1, ceil(Int, hypot((b .- a)...)/step))
        for t in range(0, 1, length=n+1)[1:end-1]
            q = a .+ t .* (b .- a)
            push!(boxes, (q[1] - hw, q[1] + hw, q[2] - hw, q[2] + hw))
        end
    end
    return boxes
end

overlap_area(a, b) = max(0.0, min(a[2], b[2]) - max(a[1], b[1])) * max(0.0, min(a[4], b[4]) - max(a[3], b[3]))
overlap(a, b) = overlap_area(a, b) > 0

"Area of box a lying outside the limits (xmin, xmax, ymin, ymax)."
outside_area(a, lim) = (a[2] - a[1])*(a[4] - a[3]) - overlap_area(a, lim)

halign_for(φ) = cos(φ) > 0.35 ? :left : (cos(φ) < -0.35 ? :right : :center)

"""
Outward coast-normal direction field: minus the gradient of the ice mask
smoothed with Gaussian σ (km), evaluated on a grid of spacing ~`res` km.
Returns a function p -> unit vector (falls back to the radial direction from
`center` where the smoothed gradient vanishes, i.e. deep in the interior).
"""
function coast_normal_field(icemask, x, y; σ=60.0, res=8.0, center=(0.0, 0.0))
    st = max(1, round(Int, res/(x[2] - x[1])))
    xc = x[1:st:end]; yc = y[1:st:end]
    S  = gauss_smooth(Float64.(icemask[1:st:end, 1:st:end]), σ/(xc[2] - xc[1]))
    nx, ny = size(S); d = xc[2] - xc[1]
    return function (p)
        i = clamp(round(Int, (p[1] - xc[1])/d) + 1, 2, nx - 1)
        j = clamp(round(Int, (p[2] - yc[1])/d) + 1, 2, ny - 1)
        gx = S[i+1, j] - S[i-1, j]; gy = S[i, j+1] - S[i, j-1]
        g = hypot(gx, gy)
        if g < 1e-4
            v = p .- center; return v ./ max(hypot(v...), 1e-6)
        end
        return (-gx/g, -gy/g)
    end
end

"""
Distance (km) from `p` along unit vector `v` at which the ice is left for good
(`nclear` consecutive ice-free steps).
"""
function ice_exit_distance(p, v, icemask, x, y; step=4.0, nclear=6, smax=3000.0)
    dx = x[2] - x[1]; dy = y[2] - y[1]
    s = 0.0; sexit = 0.0; clear = 0
    while s < smax
        px, py = p[1] + s*v[1], p[2] + s*v[2]
        i = round(Int, (px - x[1])/dx) + 1; j = round(Int, (py - y[1])/dy) + 1
        (1 <= i <= length(x) && 1 <= j <= length(y)) || break
        if icemask[i, j]
            clear = 0; sexit = s
        else
            clear += 1
            clear >= nclear && break
        end
        s += step
    end
    return sexit
end

"""
Coastal layout. Each outer label is pushed from its anchor along the local
outward normal of its mask (`masks[k]`, e.g. all ice, or grounded ice for
glaciers so they stop at the grounding line) to just beyond the edge (+ offset), then overlaps
(with each other and with the `fixed` boxes of in-place labels) are resolved by
sliding labels along the coast tangent; once a label has slid `tmax` km it is
pushed further out along the normal instead, by at most `maxpush` km (beyond
that it stays put, and the overlap is left to the caller). Labels whose ice exit is more than
`maxlead` km away (e.g. ice streams feeding the big embayments) are demoted to
in-place labels. Returns (positions, halign, outer) for every label.
"""
function layout_coastal(pts, txtwh, outer, masks, x, y; fixed=Tuple[], offset=120.0, maxlead=500.0,
                        niter=4000, step=4.0, tmax=250.0, maxpush=150.0)
    n = length(pts)
    nrm = IdDict(m => coast_normal_field(m, x, y) for m in unique(objectid, masks))
    v = [nrm[masks[k]](pts[k]) for k in 1:n]
    L = [outer[k] ? ice_exit_distance(pts[k], v[k], masks[k], x, y) + offset : 0.0 for k in 1:n]
    outer = [outer[k] && L[k] - offset <= maxlead for k in 1:n]
    Lmax = L .+ maxpush
    tng = [(-vk[2], vk[1]) for vk in v]
    t = zeros(n)
    pos(k) = outer[k] ? (pts[k][1] + L[k]*v[k][1] + t[k]*tng[k][1], pts[k][2] + L[k]*v[k][2] + t[k]*tng[k][2]) : pts[k]
    ha(k)  = outer[k] ? halign_for(atan(v[k][2], v[k][1])) : :center
    away(k, q) = (d = pos(k) .- q; s = d[1]*tng[k][1] + d[2]*tng[k][2]; s >= 0 ? 1.0 : -1.0)
    function nudge!(k, sgn)
        if abs(t[k] + sgn*step) <= tmax
            t[k] += sgn*step
        elseif L[k] + step <= Lmax[k]
            L[k] += step
        end
    end
    ctr(b) = ((b[1] + b[2])/2, (b[3] + b[4])/2)
    idx = findall(outer)
    for _ in 1:niter
        moved = false
        for a in idx, b in idx
            a < b || continue
            ba, bb = boxat(pos(a), txtwh[a], ha(a)), boxat(pos(b), txtwh[b], ha(b))
            overlap(ba, bb) || continue
            # separate the pair in opposite directions along the mean tangent
            tm = tng[a] .+ tng[b]; d = pos(b) .- pos(a)
            sg = d[1]*tm[1] + d[2]*tm[2] >= 0 ? 1.0 : -1.0
            sa = tng[a][1]*tm[1] + tng[a][2]*tm[2] >= 0 ? 1.0 : -1.0
            sb = tng[b][1]*tm[1] + tng[b][2]*tm[2] >= 0 ? 1.0 : -1.0
            nudge!(a, -sg*sa); nudge!(b, sg*sb)
            moved = true
        end
        for a in idx, fb in fixed
            overlap(boxat(pos(a), txtwh[a], ha(a)), fb) || continue
            nudge!(a, away(a, ctr(fb)))
            moved = true
        end
        moved || break
    end
    return [pos(k) for k in 1:n], [ha(k) for k in 1:n], outer
end

"Spread sorted 1-D positions `y` (ascending) to minimum spacing `gap`, staying close to targets."
function spread1d(y, gap)
    p = copy(y)
    for _ in 1:500
        changed = false
        for k in 2:length(p)
            d = p[k] - p[k-1]
            if d < gap
                s = (gap - d)/2
                p[k-1] -= s; p[k] += s; changed = true
            end
        end
        changed || break
    end
    return p
end

"""
Column layout. Outer labels with anchor x < xsplit go to the left column at
x = xleft, the others to the right column at x = xright. With
colalign = :outward text grows away from the map (columns in outer margins);
with :inward it grows towards the map (columns inside the map frame).
"""
function layout_columns(pts, txtwh, outer; xsplit, xleft, xright, gapfac=1.05, colalign=:outward)
    n = length(pts)
    pos = copy(pts); ha = fill(:center, n)
    for s in (-1, 1)
        idx = [k for k in 1:n if outer[k] && (s < 0 ? pts[k][1] < xsplit : pts[k][1] >= xsplit)]
        isempty(idx) && continue
        idx = idx[sortperm([pts[k][2] for k in idx])]
        gap = gapfac*maximum(txtwh[k][2] for k in idx)
        for (k, yk) in zip(idx, spread1d([pts[k][2] for k in idx], gap))
            pos[k] = (s < 0 ? xleft : xright, yk)
            ha[k]  = (s < 0) == (colalign === :outward) ? :right : :left
        end
    end
    return pos, ha, outer
end

"Area (km²) of the connected region of `mask` containing (or within 2 cells of) point p."
function connected_area(mask, x, y, p)
    dx = x[2] - x[1]
    i0 = round(Int, (p[1] - x[1])/dx) + 1; j0 = round(Int, (p[2] - y[1])/dx) + 1
    nx, ny = size(mask)
    seeds = [(i, j) for i in i0-2:i0+2, j in j0-2:j0+2 if 1 <= i <= nx && 1 <= j <= ny && mask[i, j]]
    isempty(seeds) && return 0.0
    seen = falses(nx, ny); stack = [first(seeds)]; seen[first(seeds)...] = true; n = 0
    while !isempty(stack)
        i, j = pop!(stack); n += 1
        for (a, b) in ((i+1, j), (i-1, j), (i, j+1), (i, j-1))
            (1 <= a <= nx && 1 <= b <= ny && mask[a, b] && !seen[a, b]) || continue
            seen[a, b] = true; push!(stack, (a, b))
        end
    end
    return n*dx^2
end

"""
Positions of the labels `keep` (see draw_labels!). Returns the label styles,
texts, font sizes, anchors, box sizes, positions, alignments, the kind flags
(outer, area, demoted) and the final label boxes.
"""
function place_labels(keep; layout, proj, kmpp, scale, obstacles, styles=LSTYLE, icemask=nothing, shelfmask=nothing,
                      groundedmask=nothing, x=nothing, y=nothing, offset=150.0, maxlead=500.0, tmax=250.0,
                      maxpush=150.0, shelf_area_min=60_000.0, xsplit=0.0, xleft=0.0, xright=0.0,
                      colalign=:outward, limits=nothing, maplimits=limits)
    sts  = [styles[l.type] for l in keep]
    txt  = [st.case(l.name) for (l, st) in zip(keep, sts)]
    fs   = [st.size*scale for st in sts]
    pts  = [project(proj, l.lon, l.lat) for l in keep]
    wh   = [textbox(t, f, kmpp; font=st.font) for (t, f, st) in zip(txt, fs, sts)]
    bigshelf = [l.type == "iceshelf" && shelfmask !== nothing &&
                connected_area(shelfmask, x, y, p) >= shelf_area_min for (l, p) in zip(keep, pts)]
    outer0 = [st.outer && !b for (st, b) in zip(sts, bigshelf)]
    # glaciers exit from grounded ice (label near the grounding line), the rest from all ice
    masks = [l.type == "glacier" && groundedmask !== nothing ? groundedmask : icemask for l in keep]
    markerbox(p) = (m = 5*scale*kmpp; (p[1] - m, p[1] + m, p[2] - m, p[2] + m))
    markers = [markerbox(pts[k]) for k in eachindex(keep) if sts[k].marker && !sts[k].outer]

    # area labels first, so they are clear of each other and of the symbols; the
    # least mobile go first (largest first), sea names (which can roam open water) last
    area = [!outer0[k] && !sts[k].marker for k in eachindex(keep)]
    maxshift = [l.type == "sea" ? 5.0 : 3.0 for l in keep]
    aidx = sort(findall(area); by=k -> (maxshift[k], -prod(wh[k])))
    apos = copy(pts)
    place_area!(apos, wh, aidx, vcat(obstacles, markers); maxshift)

    pos, ha, outer = layout === :coastal ?
        layout_coastal(pts, wh, outer0, masks, x, y; offset, maxlead, tmax, maxpush,
                      fixed=vcat(obstacles, markers, [boxat(apos[k], wh[k], :center) for k in aidx])) :
        layout_columns(pts, wh, outer0; xsplit, xleft, xright, colalign)

    # second pass: move area labels that the coastal labels still hit
    pos[aidx] .= apos[aidx]
    place_area!(pos, wh, aidx, vcat(obstacles, markers, [boxat(pos[k], wh[k], ha[k]) for k in eachindex(keep) if outer[k]]);
                maxshift)

    # Obstacles: outer labels at their final spots and fixed in-place labels.
    # Marker labels (stations, domes, demoted glaciers) then pick the first free
    # candidate position around their symbol (right, left, above, below).
    demoted = [outer0[k] && !outer[k] for k in eachindex(keep)]
    # keep outer labels inside the frame, in-place labels inside the map
    if limits !== nothing
        for k in eachindex(keep)
            lim = outer[k] ? limits : maplimits
            b = boxat(pos[k], wh[k], ha[k]); m = 10*kmpp
            sx = max(0.0, lim[1] + m - b[1]) - max(0.0, b[2] - lim[2] + m)
            sy = max(0.0, lim[3] + m - b[3]) - max(0.0, b[4] - lim[4] + m)
            pos[k] = (pos[k][1] + sx, pos[k][2] + sy)
        end
    end
    flex = [!outer[k] && (sts[k].marker || demoted[k]) for k in eachindex(keep)]
    boxes = Tuple[obstacles...]
    for k in eachindex(keep)
        (sts[k].marker || demoted[k]) && push!(boxes, markerbox(pts[k]))    # the symbol itself is an obstacle
        flex[k] && continue
        push!(boxes, boxat(pos[k], wh[k], ha[k]))
    end
    order = sortperm([(keep[k].type == "dome" ? 0 : keep[k].type == "station" ? 1 : 2) for k in eachindex(keep)])
    for k in order
        flex[k] || continue
        p = pts[k]; w, h = wh[k]; g = 6*scale*kmpp
        cands = [((p[1] + g, p[2]), :left), ((p[1] - g, p[2]), :right),
                 ((p[1], p[2] + h/2 + g), :center), ((p[1], p[2] - h/2 - g), :center),
                 ((p[1] + g, p[2] + h), :left), ((p[1] + g, p[2] - h), :left),
                 ((p[1] - g, p[2] + h), :right), ((p[1] - g, p[2] - h), :right)]
        cost(c) = (bb = boxat(c[1], wh[k], c[2]);
                   sum((overlap_area(bb, b) for b in boxes); init=0.0) +
                   (maplimits === nothing ? 0.0 : 10*outside_area(bb, maplimits)))
        c = cands[argmin(cost.(cands))]
        pos[k], ha[k] = c
        push!(boxes, boxat(c[1], wh[k], c[2]))
    end

    return (; sts, txt, fs, pts, wh, pos, ha, outer, area, demoted,
              boxes=[boxat(pos[k], wh[k], ha[k]) for k in eachindex(keep)])
end

"""
Labels to drop so that no two boxes overlap: of each overlapping pair the one
of higher tier, or the later of two labels of tier >= 2, or for two labels of
tier <= 1 the nearest label of a higher tier; also labels of tier >= 2 that hit
an obstacle.
"""
function overlap_losers(bx, keep, obstacles)
    lose = Set{Int}()
    for a in eachindex(bx), b in a+1:length(bx)
        (a in lose || b in lose || !overlap(bx[a], bx[b])) && continue
        ta, tb = keep[a].tier, keep[b].tier
        k = ta > tb ? a : tb > ta ? b : ta >= 2 ? b : 0
        if k == 0
            # two important labels pushed together: make room by dropping the
            # nearest label of a higher tier, if any
            m = ((bx[a][1] + bx[a][2] + bx[b][1] + bx[b][2])/4, (bx[a][3] + bx[a][4] + bx[b][3] + bx[b][4])/4)
            cand = [c for c in eachindex(bx) if keep[c].tier > ta && c ∉ lose]
            dist(c) = hypot((bx[c][1] + bx[c][2])/2 - m[1], (bx[c][3] + bx[c][4])/2 - m[2])
            isempty(cand) || (k = cand[argmin(dist.(cand))])
        end
        k > 0 && push!(lose, k)
    end
    for a in eachindex(bx), o in obstacles
        keep[a].tier >= 2 && overlap(bx[a], o) && push!(lose, a)
    end
    return sort(collect(lose))
end

"""
Draw labels of tier <= maxtier. `layout` is :coastal (needs icemask, x, y) or
:columns (needs xsplit, xleft, xright). `scale` multiplies all font sizes.
`obstacles` are extra boxes (km) to keep clear, e.g. drawn region labels. Ice
shelves whose connected floating area exceeds `shelf_area_min` km² (needs
shelfmask on the x, y grid) are labelled in place on the shelf.

Where labels still overlap after the layout, the less important one (higher
tier; the later one of two tier-2 labels) is dropped and the layout repeated,
so dense label sets show as many names as fit. Overlaps among tier 0-1 labels
are kept and reported.

`styles` maps label types to their style (LSTYLE by default; fonts are theme
names or font files); `ink` and `paper` are the colours of leaders, symbols and
halos. Returns the boxes of the drawn labels.
"""
function draw_labels!(ax, labs, proj; layout, kmpp, scale=1.0, maxtier=1, styles=LSTYLE, types=keys(styles),
                      seacolor=:white, ink=:black, paper=:white, elbow=50.0, obstacles=Tuple[], kw...)
    keep = [l for l in labs if l.tier <= maxtier && l.type in types && haskey(styles, l.type)]
    dropped = String[]
    P = place_labels(keep; layout, proj, kmpp, scale, obstacles, styles, kw...)
    for _ in 1:10
        lose = overlap_losers(P.boxes, keep, obstacles)
        isempty(lose) && break
        append!(dropped, [keep[k].name for k in lose])
        keep = keep[setdiff(eachindex(keep), lose)]
        P = place_labels(keep; layout, proj, kmpp, scale, obstacles, styles, kw...)
    end
    isempty(dropped) || println("labels: no room for ", join(dropped, ", "))
    (; sts, txt, fs, pts, wh, pos, ha, outer, area, demoted) = P
    kind = [outer[k] ? "coastal" : area[k] ? "area" : "symbol" for k in eachindex(keep)]
    report_overlaps(P.boxes, txt .* " (" .* kind .* ")", obstacles)

    for k in eachindex(keep)
        l, st, p, q = keep[k], sts[k], pts[k], pos[k]
        col = l.type == "sea" ? seacolor : st.color
        if outer[k]
            a = leader_attach(p, boxat(q, wh[k], ha[k]), 4*kmpp*scale)
            if layout === :columns
                xe = a[1] + (p[1] > a[1] ? elbow : -elbow)
                leader!(ax, [p[1], xe, a[1]], [p[2], a[2], a[2]]; lw=0.6*scale, color=(ink, 0.75), casing=(paper, 0.7))
            else
                leader!(ax, [p[1], a[1]], [p[2], a[2]]; lw=0.6*scale, color=(ink, 0.75), casing=(paper, 0.7))
            end
            scatter!(ax, [p]; markersize=3.5*scale, color=ink)
        elseif st.marker
            mk = l.type == "icecore" ? :diamond : l.type == "dome" ? :utriangle : l.type == "pole" ? :star5 : :circle
            scatter!(ax, [p]; markersize=8*scale, marker=mk, color=col, strokecolor=paper, strokewidth=0.8*scale)
        elseif demoted[k]
            scatter!(ax, [p]; markersize=4*scale, color=ink)
        end
        halotext!(ax, q[1], q[2]; text=txt[k], fontsize=fs[k], font=st.font, color=col,
                  align=(ha[k], :center), halo=(l.type != "sea"), halocolor=(paper, 0.8))
    end
    return P.boxes
end
