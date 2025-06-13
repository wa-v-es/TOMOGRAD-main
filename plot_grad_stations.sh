#! /bin/bash
#Script for running TOMOGRAD
#TOMOGRAD written by Chungpeng Zhao of ASU
#This script written by Dan Frost of UC Berkeley
#8.6.2021

#RULES (copied out of src/tomo_gradient_beta02_run_center.f90)
#========================================================================================================
# Calculate the tomography gradient for each layer
# Requirement:
#            1. tomography models are stored in different layer files
#            2. each layer is for one depth range
#            3. first coloumn is latitude ,second coloumn is the longitude,then dvs
#              3.1. first line is the minimum depth range and maximum depth of this layer
#            4. model is stored along latitude value line (different longitude)
#            5. latitude  and longitude order does not matter, as long as the first point
# is either minimum(increasing) or maximu (decreasing) value. otherwise  my location
# searching algorithm would fail.
# Algorithm:
#           1. reorganize the location to 2-D array grids
#           2. for each grid location, calculate great circle sample points for every 1 degree
# distances for every 10 degree azimuth direction.
#           3. search for grid locations surrounding these great circle sample points
#           4. calculate an average dvs using these surrounding grid locations for this sample point
#           5. so in azimuth direction, you have many sample locations, with different dvs and distance
# from the grid location
#           6. Using linear_square to fit a line to get the great for this azimuth directions
#           7. Get the maximum or minimum ( whichever their absolute value is bigger) of these gradients
# for different azimuth directions.
#           8. Then this is the gradient for this grid location, and store the direction of this gradient
#
#OUTPUT lat, lon, gradient, gradient_az, av_standard_deviation (of all directions used)
#=========================================================================================================

#Directories
home_dir="/Users/keyser/Research/TOMOGRAD-main/"
cd $home_dir

##
## this takes the mod of gradient.
# awk '{ $3 = ($3 < 0) ? -$3 : $3; print }' GRADIENT.$model_file > GRADIENT.posit_2.5.$model_file
#<55/71.5>	<-165/-135.5>
# exit
model_file=AlaskaMohoOpt_forGrad_.2_cut.XYZ  #XYZ of tomography data for input (lon, lat, dvs at 1x1 degree spacing)

#-----------
#Plot output
#-----------
interp="-I.1"
range="-R-171/-133/55/71" #-171/-133/56/71
range_cal="-R-170/-135/50/75" #-171/-133/56/71

frame="-Jm0.1i"
# frame="-JS-150/62/10c"

border="-B10g10/5g5WSen"
coast_params="-Dh -W.001c  -A100"
xoff0="-X3c"; yoff0="-Y3c"
yoff1="-Y12c"
scalevar="-D5c/-1c/5.5c/0.25ch"

# gmt set FONT_ANNOT_PRIMARY 6 FONT_HEADING 5 FONT_TITLE 9 MAP_FRAME_TYPE plain

# gmt makecpt -T0/1 -I -CgrayC > gg.cpt #-G0/1

gmt makecpt -Cpolar -T-3/3/0.1 -D -I -Z > tomo.cpt
gmt makecpt -Cocean -T0/0.5/0.05 -D -I > grad.cpt

gmt makecpt -T-1/1/.1 -I -Ccork > gg.cpt #-G0/1

# Set GMT defaults
gmt gmtset FONT_TITLE 12
gmt gmtset MAP_TITLE_OFFSET 0.0c
gmt gmtset FONT_ANNOT_PRIMARY 12
gmt gmtset FONT_ANNOT_SECONDARY 12
gmt gmtset FONT_LABEL 12
gmt gmtset MAP_LABEL_OFFSET 0.1c
###
# gradients in same lat long as input, which was created from Miller Moresi dataset (get_xyx_vals.py)
#had 0.2 deg resolution
#  grad --> blockmean --> surface
# cut grd to remove edges. The grad is calculated for 2.5deg circle, so lots of 'edge' effect
#Mask to remove grad values on ocean
# normalising the masked grd using max Z in masked grd
# grad_.25_mask.grd is copied to moho_ak folder to be plotted using py script

###
#Construct grids
for azimuth in $(seq 0 10 10); do


  # azimuth=60
  outps="GRADIENT_"$azimuth"_new.ps"

  awk 'NR>1 && $4=='$azimuth' {print $2, $1, $3}' GRADIENT.$model_file | gmt blockmean $interp $range_cal | gmt surface $interp $range_cal -GGRADIENT.$azimuth.grd
  # exit
  gmt grdcut GRADIENT.$azimuth.grd -Gcut_GRADIENT.$azimuth.grd -R-159.5/-140/56.4/70
  # gmt grdinfo cut_GRADIENT.grd
  ##
  gmt grdlandmask -R-159.5/-140/56.4/70 -Dc -A10 $interp -NNan/1 -Gmask.grd
  # gmt grdinfo mask.grd
  # exit
  gmt grdmath cut_GRADIENT.$azimuth.grd mask.grd MUL = grad_mask.$azimuth.grd
  gmt grdinfo grad_mask.$azimuth.grd
  # Extract the maximum absolute value (positive or negative)
  v_max=$(gmt grdinfo -C grad_mask.$azimuth.grd | awk '{print ($7 > -$6) ? $7 : -$6}')
  echo $vmax
  gmt grdmath grad_mask.$azimuth.grd $v_max DIV = grad_mask.$azimuth.norm.grd
  # v_max=$(gmt grdinfo -T grad_mask.$azimuth.grd | sed 's/.*\///')
  # gmt grdmath grad_mask.$azimuth.grd 0 SUB $v_max 0 SUB DIV = grad_mask.$azimuth.norm.grd

  # gmt grdinfo normalized.grd
  # awk 'NR>1 {print $2, $1, $3}' $model_file | gmt blockmean $interp $range | gmt surface $interp $range -GTOMO.grd

  #Plot gradient
  gmt grdimage grad_mask.$azimuth.norm.grd $range $frame $border -Cgg.cpt $xoff0 $yoff0 -P -K > $outps
  gmt pscoast $range $frame $coast_params -O -K >> $outps
  gmt psscale $scalevar -Cgg.cpt -B0.2:"Moho gradient $azimuth (% dV/0.1deg)": -O -K >> $outps

  awk '{print $8, $7}' STA_DISTANCE_LOC_gridnumber128.txt | gmt psxy $range $frame $border -St0.3c -W.5p,darkorange -O -P -K >> $outps

  awk '{print $8, $7}' STA_DISTANCE_LOC_gridnumber97.txt | gmt psxy $range $frame $border -St0.3c -W.5p,maroon -O -P >> $outps
  # awk '{print $8, $7}' STA_DISTANCE_LOC_gridnumber112.txt | gmt psxy $range $frame $border -St0.3c -W.5p -O -P >> $outps

  #
  #Plot tomography
  # gmt grdimage TOMO.grd $range $frame $border:."$model_file - gradient over $gradient_radius degrees at azimuth=$plot_azimuth": -Ctomo.cpt $yoff1 -P -O -K >> $outps
  # gmt pscoast $range $frame $coast_params -O -K >> $outps
  # gmt psscale $scalevar -Ctomo.cpt -B1:"Velocity (% dV)": -O >> $outps
  gmt ps2raster -A -Tj -E720 -P -Z -Vq $outps
  # open GRADIENT_$azimuth.jpg
  # exit
  # mv -f *$azimuth*.txt 2.5deg_grad_grds/
  rm -f *GRAD*.grd
  mv -f *$azimuth*.grd 2.5deg_grad_grds/
  rm -f mask.grd

  mv -f *$azimuth*.jpg 2.5deg_grad_grds/
  echo "Azimuth:$azimuth done"
  # exit
done
# rm
# rm temp.grd

# gv $outps &
