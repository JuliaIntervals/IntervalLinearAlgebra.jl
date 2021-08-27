using Luxor, IntervalLinearAlgebra, LazySets

A = [2..4 -2..1;-1..2 2..4]
b = [-2..2, -2..2]

polytopes = solve(A, b, LinearOettliPrager())

# convert h polytopes to Luxor objects
function luxify(P::HPolytope)

    # convert to vertex representation
    V = convert(VPolygon, P)
    out = vertices_list(P)

    out_luxor = [Luxor.Point(Tuple(p)) for p in out]

    return out_luxor
end

luxor_polytopes = luxify.(polytopes)

# draw logo
function logo(polytopes, fname)
    Drawing(500, 500, fname)
    origin()
    # transparent sethue("grey30")
    # squircle(O, 248, 248, :fill, rt=0.1)
    cols = [Luxor.julia_green, Luxor.julia_blue, Luxor.julia_purple, Luxor.julia_red]
    sizes = [1, 2.5, 4]
    Luxor.scale(1, -1)
    for i in 1:4
        sethue("black")
        poly(polytopes[i]*60, :stroke, close=true)

        sethue(cols[i])
        poly(polytopes[i]*60, :fill, close=true)
        #setline(rescale(sizes[i], 1, 4, 18, 18))
    end

    finish()
    preview()
end


function lockup(logo, fname)
    Drawing(1000, 300, fname)
    # background("grey90")     # otherwise transparent
    origin()
    logosvg = readsvg(logo) # this requires Luxor#master
    panes = Table([300], [300, 700])
    # place SVG
    @layer begin
        Luxor.translate(panes[1])
        Luxor.scale(0.5)
        placeimage(logosvg, O, centered=true)
    end
    # text
    @layer begin
        sethue("black")#sethue("rebeccapurple")
        setline(0.75)
        Luxor.translate(boxmiddleleft(BoundingBox(box(panes, 2))))
        fontsize(50) #80
        fontface("JuliaMono Bold")
        titlex = -40
        Luxor.text("IntervalLinearAlgebra.jl", Point(titlex, 0))
        @layer begin
            sethue("grey60")
            textoutlines("IntervalLinearAlgebra.jl", Point(titlex, 0), :stroke)
        end
        fontsize(40) #45
        fontface("JuliaMono")
        subx = -38
        Luxor.text("Linear algebra done rigorously", Point(subx, 45))
        @layer begin
            sethue("grey60")
            textoutlines("Linear algebra done rigorously", Point(subx, 45), :stroke)
        end
    end
    finish()
    preview()
end

logo(luxor_polytopes, "logo.svg")
lockup("logo.svg","logo-text.svg")
