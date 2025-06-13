#! /bin/bash

awk '{ print $2,$1,$3} ' AlaskaMohoOpt_forGrad_cut.XYZ > AlaskaMohoOpt_forGrad_cut_format.XYZ
gmt xyz2grd AlaskaMohoOpt_forGrad_cut_format.XYZ -GAlaskaMohoOpt.grd -R-164.9/-134.9/54.9/71.6 -I0.1 -V
gmt gmtinfo AlaskaMohoOpt.grd

for azimuth in $(seq 0 10 360); do
    gmt grdgradient AlaskaMohoOpt.grd -A${azimuth} -N -Ggradient_${azimuth}.grd
    gmt grd2xyz gradient_${azimuth}.grd > gradient_${azimuth}.txt
    awk '!/NaN/' gradient_${azimuth}.txt > gradient_${azimuth}_no_nan.txt
done
# gmt grdgradient AlaskaMohoOpt.grd -A50 -GAlaskaMohoOpt_grad.grd -V
# gmt grd2xyz AlaskaMohoOpt_grad.grd > gradient_data.txt -V
