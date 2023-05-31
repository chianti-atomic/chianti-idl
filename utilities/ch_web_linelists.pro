
PRO ch_web_linelists 


;+
; NAME:
;     CH_WEB_LINELISTS
;
; PURPOSE:
;     Generates the set of line list tables that are linked to the CHIANTI
;     website. Only the latex files are created; the user needs to run
;     pdflatex to generate the pdf files.
;
;     The line lists are created for the "extended-flare" DEM for a
;     pressure of 10^16 and for photospheric abundances. Intensity limits
;     are applied to reduce the number of lines in the tables. Population
;     lookup tables are used to reduce computation time.
;
; CATEGORY:
;     CHIANTI; line lists.
;
; CALLING SEQUENCE:
;     CH_WEB_LINELISTS
;
; INPUTS:
;      None.
;
; OUTPUTS:
;      Creates the following six latex files in the working directory.
;       ch_line_list_v10.0.1_1_50.tex
;       ch_line_list_v10.0.1_50_150.tex
;       ch_line_list_v10.0.1_150_912.tex
;       ch_line_list_v10.0.1_912_2000.tex
;       ch_line_list_v10.0.1_2000_10000.tex
;       ch_line_list_v10.0.1_10000_600000.tex
;
; EXAMPLE:
;      IDL> ch_web_linelists
;
; MODIFICATION HISTORY:
;      Ver.1, 08-Jun-2022, Peter Young
;      Ver.2, 31-May-2023, Peter Young
;        Chnaged abundance file to !abund_file.
;-



ch_ver=ch_get_version()
press=1e16

;
; The DEM is hard-coded to the extended flare DEM.
;
dem_name=concat_dir(!xuvtop,'dem')
dem_name=concat_dir(dem_name,'flare_ext.dem')
chck=file_info(dem_name)
IF chck.exists EQ 0 THEN BEGIN
  message,/info,/cont,'The DEM file was not found. Returning...'
  return
ENDIF 

abund_name=!abund_file
chck=file_info(abund_name)
IF chck.exists EQ 0 THEN BEGIN
  message,/info,/cont,'The abundance file was not found. Returning...'
  return
ENDIF 

w0=1
w1=50
mini=1.16e6/1e3
outfile='ch_line_list_v'+ch_ver+'_'+trim(w0)+'_'+trim(w1)+'.tex'
latex_wvl_dem,w0,w1,outfile=outfile,dem_name=dem_name,abund_name=abund_name,/all, $
              mini=mini,pressure=press,/lookup

w0=50
w1=150
mini=2.28e6/1e3
outfile='ch_line_list_v'+ch_ver+'_'+trim(w0)+'_'+trim(w1)+'.tex'
latex_wvl_dem,w0,w1,outfile=outfile,dem_name=dem_name,abund_name=abund_name,/all, $
              mini=mini,pressure=press,/lookup


w0=150
w1=912
mini=1.27e6/1e3
outfile='ch_line_list_v'+ch_ver+'_'+trim(w0)+'_'+trim(w1)+'.tex'
latex_wvl_dem,w0,w1,outfile=outfile,dem_name=dem_name,abund_name=abund_name,/all, $
              mini=mini,pressure=press,/lookup


w0=912
w1=2000
mini=1.15e6/1e3
outfile='ch_line_list_v'+ch_ver+'_'+trim(w0)+'_'+trim(w1)+'.tex'
latex_wvl_dem,w0,w1,outfile=outfile,dem_name=dem_name,abund_name=abund_name,/all, $
              mini=mini,pressure=press,/lookup


w0=2000
w1=10000
mini=3.29e4/1e3
outfile='ch_line_list_v'+ch_ver+'_'+trim(w0)+'_'+trim(w1)+'.tex'
latex_wvl_dem,w0,w1,outfile=outfile,dem_name=dem_name,abund_name=abund_name,/all, $
              mini=mini,pressure=press,/lookup


w0=10000
w1=600000.
mini=1.0
outfile='ch_line_list_v'+ch_ver+'_'+trim(w0)+'_'+trim(w1)+'.tex'
latex_wvl_dem,w0,w1,outfile=outfile,dem_name=dem_name,abund_name=abund_name,/all, $
              mini=mini,pressure=press,/lookup


END
