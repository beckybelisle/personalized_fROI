# example call of FreeSurfer mri_vol2surf

mri_vol2surf --mov ${zstatVolume} \
	     --regheader ${subj} \
	     --surf-fwhm 0 \
	     --hemi ${hemi} \
	     --trgsubject ${subj} \
	     --projdist 0 \
	     --interp nearest \
	     --o ${outDir}/${hemi}.zstat3.dist0.interpNN.gii \
	     --out_type gii
