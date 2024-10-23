
PRO ch_plot_all_iso_pops

;+
; NAME:
;     CH_PLOT_ALL_ISO_POPS
;
; PURPOSE:
;     Creates plots of level populations along isoelectronic sequences for
;     a wide range of sequences. Intended for checking atomic data prior to
;     new CHIANTI releases.
;
; CATEGORY:
;     CHIANTI; data validation.
;
; CALLING SEQUENCE:
;     CH_PLOT_ALL_ISO_POPS
;
; INPUTS:
;     None.
;
; OUTPUTS:
;     Creates two directories in the current working directory called
;     'lower_levels' and 'upper_levels'. A number of png files are created
;     in the former directory showing level populations for the lower
;     levels of the ion (typically those of the ground configuration). The
;     'upper_levels' directory contains plots for higher levels typically
;     belonging to the first excited configuration. All plots are created
;     with the routine ch_plot_iso_pops.
;
; CALLS:
;     CH_GET_VERSION, CH_PLOT_ISO_POPS, ZION2NAME, ZION2ELEMENT
;
; EXAMPLE:
;     IDL> ch_plot_all_iso_pops
;
; MODIFICATION HISTORY:
;     Ver. 1, 15-Oct-2024, Peter Young
;-




outdir='lower_levels'
chck=file_info(outdir)
IF chck.exists EQ 0 THEN file_mkdir,outdir

seq= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
lev1=[2,2,2,1,1,1,1,1,1, 2, 2, 1, 1, 1, 1, 1, 1, 2]
lev2=[9,7,8,5,5,5,5,5,3, 5, 5, 5, 5, 5, 5, 5, 3,13]

n=n_elements(seq)

ver=ch_get_version()


FOR i=0,n-1 DO BEGIN
  IF seq[i] LE 12 THEN z_ref=14 ELSE z_ref=26
 ;
  zion2name,z_ref,z_ref-seq[i]+1,ionname
  z2element,seq[i],element
  nlev=lev2[i]-lev1[i]+1
  p=ch_plot_iso_pops(ionname,indgen(nlev)+lev1[i])
  outfile='plot_'+ver+'_'+element+'_'+trim(lev1[i])+'_'+trim(lev2[i])+'.png'
  outfile=concat_dir(outdir,outfile)
  p.save,outfile,width=1000
  p.close
ENDFOR


;---------------------
outdir='upper_levels'
chck=file_info(outdir)
IF chck.exists EQ 0 THEN file_mkdir,outdir

seq= [4,  5, 6, 7, 8, 9,10,12,13,14,15,16,17,18]
lev1=[6,  6, 6, 6, 6, 4, 6, 6, 6, 6, 6, 6, 4,14]
lev2=[10,10,15,13, 9,10,15,10,12,14,10, 9,31,17]

n=n_elements(seq)

ver=ch_get_version()



FOR i=0,n-1 DO BEGIN
  IF seq[i] LE 12 THEN z_ref=14 ELSE z_ref=26
 ;
  zion2name,z_ref,z_ref-seq[i]+1,ionname
  z2element,seq[i],element
  nlev=lev2[i]-lev1[i]+1
  p=ch_plot_iso_pops(ionname,indgen(nlev)+lev1[i])
  outfile='plot_'+ver+'_'+element+'_'+trim(lev1[i])+'_'+trim(lev2[i])+'.png'
  outfile=concat_dir(outdir,outfile)
  p.save,outfile,width=1000
  p.close
ENDFOR

message,/info,/cont,'Plots have been written to the directories lower_levels and upper_levels.'


END
