@everywhere screenx = 7680
@everywhere screeny = 4320
@everywhere screenx = 1920
@everywhere screeny = 1080
#=@everywhere screenx = 800=#
#=@everywhere screeny = 600=#
#=@everywhere screenx = 400=#
#=@everywhere screeny = 300=#
#=@everywhere screenx = 200=#
#=@everywhere screeny = 150=#

@everywhere xcenter = (-0.0)
@everywhere ycenter = (-0.0)

@everywhere xzoom    = ( 2.79914)
@everywhere yzoom    = ( 3.2539551)

@everywhere minx    = (xcenter + xzoom)
@everywhere maxx    = (xcenter - xzoom)
@everywhere miny    = (ycenter + yzoom)
@everywhere maxy    = (ycenter - yzoom)

@everywhere stepx   = (screenx/(maxx-minx))
@everywhere stepy   = (screeny/(maxy-miny))

@everywhere iters   = 750

@everywhere aa      = 1

@everywhere bailout = 66.6

@everywhere function m(z, cord)
    for i = 1:iters

        z = sin(z) * (cord)

        #=z = cos(z) * (cord)=#

        #=z = sin(z) * (cord) + (z^2.0)=#

        if abs(z) > bailout
            w = (i, z)
            return w
        end

    end

    return (0, z)
end

@everywhere function julia_line(c, cord)
    e = Array((Any), screenx)

    for x in 1:screenx
        e[x] = julia(Complex(minx + x*(maxx-minx)/screenx, imag(c)), cord)
    end

    return e
end

@everywhere function julia(c, cord)
    dx = ((2.125 * xzoom) / screenx) / (aa*2 + 1)
    dy = ((2.125 * yzoom) / screeny) / (aa*2 + 1)

    bm = Array((Any), (aa*2 + 1, aa*2 + 1))

    for i in 1:(aa*2 + 1), j in 1:(aa*2 + 1)
        d = complex(dx * (i-aa), dy * (j-aa))
        bm[i, j] = m(c + d, cord)
    end

    it = sum([ bm[i, j][1] for i in 1:(aa*2 + 1), j in 1:(aa*2 + 1) ]) / ((aa*2 + 1)^2)
    zz = sum([ bm[i, j][2] for i in 1:(aa*2 + 1), j in 1:(aa*2 + 1) ]) / ((aa*2 + 1)^2)

    return (it, zz)
end

function load_pal(name)
    file_pal = open(name, "r")

    dump = readline(file_pal)
    dump = readline(file_pal)
    dump = readline(file_pal)
    dump = readline(file_pal)

    pal = Array((Int, Int, Int), 255)

    pivot :: Int = 0

    r, g, b = 0, 0, 0

    while eof(file_pal) == false && pivot < 255
        pivot += 1

        r = int(readline(file_pal))
        g = int(readline(file_pal))
        b = int(readline(file_pal))

        pal[pivot] = (r, g, b)
    end

    return pal
end

function get_color(pal, index)
    if index == 0
        index += 1
    end

    return pal[index]
end

function ppm_write(img)
    out = open("out.ppm", "w")

    write(out, "P6\n")
    x, y = size(img)
    write(out, "$x $y 255\n")

    for j = 1:y, i = 1:x
        p = img[i,j]

        write(out, uint8(p[1]))
        write(out, uint8(p[2]))
        write(out, uint8(p[3]))
    end
end

function initArray_perLine()
    cmap = Array(Complex, (screeny))

    for y in 1:screeny
        cmap[y] = Complex(0.0, miny + y*(maxy-miny)/screeny)
    end

    return cmap
end

function normalizer(bitmap_z)
    bitmap_t = Array(Float64, (screenx, screeny))


    for i in 1:screenx, j in 1:screeny
        if bitmap_z[i, j][1] == 0
            bitmap_t[i, j] = 0.0
        else
            mag = abs(bitmap_z[i, j][2])
            bitmap_t[i, j] = (bitmap_z[i, j][1] + 1 - log( abs(log(mag))) / log(4))
        end
    end


    nn = int( maximum(bitmap_t) )

    histogram = Array(Int, nn)
    fill!(histogram, 0)

    x, y = size(bitmap_t)

    bitmap = bitmap_t #[ int(bitmap_t[i, j] / max) for i in 1:x, j in 1:y ]

    return bitmap, [1]

    for i = 1:x, j = 1:y
        if bitmap[i, j] > 0 && bitmap[i, j] < nn
            histogram[ int(floor(bitmap[i, j])) ] += 1
        end
    end

    total_hist = sum(histogram)

    p   = Array(Float64, (screenx, screeny))
    fda = zeros(Float64, nn)

    for i in 1:nn
        p[i] = (histogram[i] * (nn)) / (x/y) :: Float64
    end

    for i in 1:nn
        for j in 1:i
            fda[i] += int(p[j])
        end
    end

    return bitmap, fda
end

function colorizer(bitmap, fda, pal)
    x, y = size(bitmap)
    bitmap_color = Array((Int, Int, Int), (x,y))

    max  = maximum(fda) / 255
    maxb = maximum(bitmap) / 255

    for i in 1:x, j in 1:y
        if bitmap[i, j] <= 0
            bitmap[i, j] = 1
        end

        if bitmap[i, j] > 255
            bitmap[i, j] = 255
        end

        #=println(fda[ bitmap[i, j] ], " ", bitmap[i, j])=#

        bitmap_color[i, j] = get_color(pal, int(bitmap[i, j] / maxb))

        #=bitmap_color[i, j] = get_color(pal, int( fda[ int(bitmap[i, j]) ] / max))=#

    end

    return bitmap_color
end

function main()
    println("\nStarting")

    tic()

    #=println("Running at ", sum( pmap(peakflops, [ 2000 for i = 1:nworkers()])) / 10^9, " GFlops.")=#

    timert :: Float64 = toq()
    timer  :: Float64 = timert
    println("Measurament took:   \t", timert)

    tic()
    pal = load_pal("pals/friendly.ppm")

    timert = toq()
    timer += timert

    println("Pallete loading took:\t", timert)

    ################################################

    tic()

    cord = Complex( 1.2859512
                  , 0.87555525
                ) * pi/2.0

    cmap = initArray_perLine()

    timert = toq()
    timer += timert

    println("Preallocing took:  \t", timert)

    ################################################

    tic()

    bitmap_zz = pmap(x -> julia_line(x, cord), cmap, err_retry=true, err_stop=false)

    bitmap_z = reshape([ bitmap_zz[i][j] for j in 1:screenx, i in 1:screeny ], (screenx, screeny))

    #=bitmap_z = reshape(bitmap_zz, (screenx, screeny))=#

    timert = toq()
    timer += timert

    println("Iterating took:  \t", timert)

    ################################################

    tic()

    bitmap, fda = normalizer(bitmap_z)

    timert = toq()
    timer += timert

    println("Normalization took:\t", timert)

    ################################################

    tic()

    bitmap_color = colorizer(bitmap, fda, pal)

    timert = toq()
    timer += timert

    println("Colouring took:  \t", timert)

    ################################################

    tic()
    ppm_write(bitmap_color)

    timert = toq()
    timer += timert

    println("Salving took:    \t", timert)

    ################################################

    println("Total time was:  \t", timer)

    println("Done\n\n")
end

main()

#=@profile main()=#
#=Profile.print()=#
