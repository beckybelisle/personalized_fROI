import os
import nibabel as nib
import numpy as np
from scipy import stats, interpolate
import matplotlib.pyplot as plt
import pandas as pd
from kneed import KneeLocator

def personalized_fROI(sig_map, effectSize_map, parcel, region_name):
    # mask fMRI results with parcel
    sig_parcel = sig_map[parcel == 1]
        
    # percentile-ize
    sig_pct = np.round(np.vectorize(lambda x: stats.percentileofscore(sig_parcel, x))(sig_parcel), 0)
    effectSize_parcel = effectSize_map[parcel == 1]

    # avg effect size at each z-stat %ile
    avg_data = []

    # loop through %ile values, getting avg effect size
    for i in range(101):
        if len(effectSize_parcel[sig_pct == i]) == 0:
            print("no values at this pctile")
            avg_effect = 0
        else:
            avg_effect = np.mean(effectSize_parcel[zstat_pct == i])
        new_row = {'pctile': i, 'effect': avg_effect}
        avg_data.append(new_row)

    avg_df = pd.DataFrame(avg_data)


    y_data = avg_df["effect"].to_numpy()
    x_data = avg_df["pctile"].to_numpy()

    # fit smooth curve (just for visualization)
    tck, u = interpolate.splprep([x_data, y_data], s=100)
    u_new = np.linspace(0, 100, 1)
    x_new, y_new = interpolate.splev(u_new, tck)

    # find knee (using the non-interpolated data)
    # note: additional parameters exist / can be altered if knee is not found (eg, switching to online = True)
    knee_locator = KneeLocator(x_data, y_data, curve='convex', direction='increasing', interp_method = 'polynomial')
    knee_point = knee_locator.knee

    # plotting: avg effect size against significance %ile + the fitted curve + the knee point
    plt.scatter(x_data, y_data, label='Original Data')

    plt.plot(x_new, y_new, '--', label='Fitted Curve')

    # check if a knee point was found before adding to plot
    if knee_point is not None:
        plt.scatter(knee_point, y_data[np.where(x_data == knee_point)], color='orange', label='Knee Point')
        plt.axvline(x=knee_point, color='orange', linestyle='--', linewidth = 2)
        plt.xlabel('Activation Significance (Percentile)', fontsize = 15)
        plt.ylabel('Average fMRI Response', fontsize = 15)
        plt.title('Parcel %s: Optimal fROI Percentile'%(region_name), fontsize = 15)
        print(knee_point)

    return knee_point, sig_pct # the %ile cutoff in terms of activation significance (for typical approach, would be 90)


# ----EXAMPLE USAGE------
knees_df = pd.DataFrame(columns=['Subject', 'Hemi', 'Region', 'Knee'])


if parcelHemi == "lh":
    region_list = ["3", "4", "5", "6", "7", "8"]  # LH
if parcelHemi == "rh":
    region_list = ["1", "2", "3", "7", "8", "10"] # RH

parcel_version = "probthresh0.15_FWHM8_sigpct_0.5"
for subj in subject_list:
    
    # Z STAT (on surface)
    funcDir = "%s/analysis/%s/func-data/LangLoc/overlays/native"%(mainProjDir, subj)
    zstat = nib.load("%s/%s.zstat3.dist0.interpNN.gii"%(funcDir, parcelHemi)).darrays[0].data


    # COPE (on surface)
    cope = nib.load("%s/%s.cope3.dist0.interpNN.gii"%(funcDir, parcelHemi)).darrays[0].data

    # PARCEL (on surface)
    parcelDir = "%s/analysis/%s/lipkin-parcels_surface/%s"%(mainProjDir, subj, parcel_version)

      
    for parcelNum in region_list:
        print("---starting parcel %s-----"%(parcelNum))

        parcels_all = nib.load("%s/%s.parcels.gii"%(parcelDir, parcelHemi)).darrays[0].data
        parcel = parcels_all.copy()
        parcel[parcel != int(parcelNum)] = 0 # isolate specific parcel
        parcel[parcel > 0] = 1 # binarize

  
        knee_point, zstat_pct = personalized_fROI(zstat, cope, parcel, parcelNum)

        knees_df = pd.concat([knees_df, pd.DataFrame({'Subject': [subj], 'Hemi': [parcelHemi], 'Region': [parcelNum],'Knee': [knee_point]})], ignore_index=True)


        # saving the fROI, using the personalized threshold
        outDir = "%s/analysis/%s/surface-parcel-fROIs/04-10-2025"%(mainProjDir,subj)
        isExist = os.path.exists(outDir)
        if not isExist:
            os.makedirs(outDir)

        fROI_parc = zstat_pct.copy() # percentile values (within parcel only)

        # threshold + binarize the z-stat %ile array
        fROI_parc[fROI_parc < knee_point] = 0
        fROI_parc[fROI_parc >= knee_point] = 1


        # fROI was initially defined relative to parcel
        # need to put back into "context" of entire surface (for this hemi)
        fROI_surf = parcel.copy()
        fROI_surf[fROI_surf == 1] = fROI_parc
        
        # initialize fROI GIFTI as the parcel GIFTI
        fROI = nib.load("%s/%s.parcels.gii"%(parcelDir, parcelHemi))
        fROI.darrays[0].data = fROI_surf

        # save the new fROI GIFTI
        nib.save(fROI, "%s/%s.fROI_%s.gii"%(outDir, parcelHemi, parcelNum))
        print("---wrote the fROI gifti---")
