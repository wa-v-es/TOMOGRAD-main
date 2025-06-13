#!/bin/bash

##

REGION="-171/-133/56/71"      # Define your region: min_lon/min_lat/max_lon/max_lat
PROJ="m0.1i"               # Define the projection: here it's a mercator with width of 0.1 inch
#
gmt set FONT_ANNOT_PRIMARY 6 FONT_HEADING 5 FONT_TITLE 9 MAP_FRAME_TYPE plain
gmt makecpt -T0/1 -G0/1 -I -CgrayC > gg.cpt
for FILE in gradient_*_no_nan.txt; do

    # Extract azimuth from the filename
    AZIMUTH=$(echo "$FILE" | sed -E 's/gradient_([0-9]+)_no_nan\.txt/\1/')
    awk '{ $3 = ($3 < 0) ? -$3 : $3; print }' $FILE > "${FILE%.txt}_posit.txt"
    # Create a grid from the XYZ data
    gmt blockmean "${FILE%.txt}_posit.txt" -R$REGION -I0.1 > "${FILE%.txt}_blockmean.txt" # Create a block mean grid
    gmt surface "${FILE%.txt}_blockmean.txt" -R$REGION -I0.1 -G"${FILE%.txt}.grd"    # Interpolate the grid

    # Make a plot for each file
    gmt begin "plot_${AZIMUTH}" png    # Start a new plot in PDF format with azimuth in the name
        gmt grdimage "${FILE%.txt}.grd" -R$REGION -J$PROJ -Cgg.cpt  # Use a colormap (like viridis)
        # gmt colorbar -DJBC+o0.1i  # Add a color bar
        gmt coast -Wthin   # Add coastlines
        gmt basemap -Bf2g2 -BWSne+t"Gradient at ${AZIMUTH}"    # Add axes and a title
    gmt end     # Finish the plot and display it
    # exit
done
magick -delay 80 -quiet -quality 50 -loop 0 *.png out.gif
