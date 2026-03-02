load "params.gp"
set terminal png size 500,500 #set the dimensions of the image (the frame)
set xlabel "Position (x)" #name the x axis
set ylabel "Position (z)" #name tha y axis
set xrange [0:xmax] #set the boundaries of the x values
set yrange [0:ymax] #set the range as the same as x for the circle to not need to be distorted


do for [i=1:nframes] {
    set output sprintf("frames/frames_%03d.png", i)
    set title sprintf("Time: %.3f s", i*dt*npoints/nframes)
    plot_cmd = "plot "

    do for [j=1:npar] {
        if (j > 1) { 
            plot_cmd = plot_cmd . ", " 
        }
        plot_cmd = plot_cmd . sprintf("'plot2D.dat' every ::%d::%d using %d:%d:(%g) with circles fs solid notitle", i-1, i-1, j+npar, j, r[j])
    }
    eval(plot_cmd)
    percent = i*100.0/nframes
    int_5=int(nframes*0.05)
    if (i%int_5==0){
        print sprintf("\rProgress: %3.0f%% (%d/%d)", percent, i, nframes)
    }
}
set output