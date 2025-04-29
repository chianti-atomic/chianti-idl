
FUNCTION ch_read_atmos, filename, keep_all=keep_all

;+
; NAME:
;     CH_READ_ATMOS
;
; PURPOSE:
;     Reads an atmosphere parameters file that is needed for modeling
;     charge transfer.
;
; CATEGORY:
;     CHIANTI; read.
;
; CALLING SEQUENCE:
;     Result = CH_READ_ATMOS( )
;
; INPUTS:
;     None.
;
; OPTIONAL INPUTS:
;     Filename: The name of the file to read. If not specified, then a
;               widget will appear, asking you to choose from one of the
;               CHIANTI options.
;	
; KEYWORD PARAMETERS:
;     KEEP_ALL: By default, the routine removes data points for which the
;               temperature is not monotonically increasing with height.
;               Setting this keyword will keep all of the data points.
;
; OUTPUTS:
;     A structure with the following tags:
;      .filename  Name of the data file.
;      .temp      Array of temperatures (K).
;      .elec_dens Array of electron densities (cm^-3).
;      .height    Array of heights (km).
;      .pressure  Array of pressures (K cm^-3).
;      .hyd_dens  Array of hydrogen (neutral & ionized) densities.
;      .h_elec    Array of hydrogen density to electron density ratios.
;      .h1_frac   Array of neutral hydrogen ionization fractions.
;      .he1_frac  Array of neutral helium ionization fractions.
;      .he2_frac  Array of neutral helium ionization fractions.
;      .comments  String array giving the file comments.
;      .time_stamp String giving the time the structure was created.
;
;     If the file does not contain data for helium, then atmos_he1 and
;     atmos_he2 will be set to -1.
;
;     If a problem is found, then a value of -1 is returned.
;
; EXAMPLE:
;     IDL> atmos_data=ch_read_atmos()
;
;     IDL> atmos_data=ch_read_atmos('my_atmos_file.dat')
;
; MODIFICATION HISTORY:
;     Ver.1, 23-Apr-2025, Peter Young
;       Code extracted from ch_adv_model_setup.
;-


if n_elements(filename) eq 0 then BEGIN
  search_dir=concat_dir(!xuvtop,'ancillary_data')
  search_dir=concat_dir(search_dir,'advanced_models')
  search_dir=concat_dir(search_dir,'model_atmospheres')
  filename=ch_get_file(path=search_dir, $
                       filter='*.dat', $
                       tit='Select an atmosphere file for the charge transfer calculation') 
ENDIF 

chck=file_info(filename)
IF chck.exists EQ 0 THEN BEGIN
  message,/info,/cont,'The specified file does not exist. Returning...'
  return,-1
ENDIF 


atmos_data=strarr(file_lines(filename))
openr,lf,filename,/get_lun
readf,lf,atmos_data
free_lun,lf

comment_ind=where(atmos_data.startswith('-1'),nc)
if nc gt 0 then nrecs=comment_ind[0] else nrecs=n_elements(atmos_data)

IF nc EQ 0 THEN comments='' ELSE comments=atmos_data[nrecs+1:*]

rdfloat,filename,atmos_temp,atmos_elec,atmos_height,$
        atmos_pressure,atmos_hyd,atmos_h1,atmos_he1,atmos_he2,numline=nrecs,/double


; all atmospheric parameters below the temperature minimum and above the temperature maximum
; have to be removed for temperature interpolation in ionization equilibrium
IF NOT keyword_set(keep_all) THEN BEGIN 
  atmos_min=min(atmos_temp,atmin)
  atmos_max=max(atmos_temp,atmax)
  if atmin ne 0 or atmax ne nrecs-1 then BEGIN
    IF NOT keyword_set(quiet) THEN BEGIN 
      print,'Double-valued temperature array found in atmosphere file,'
      print,'parameters below temperature minimum and above temperature maximum will be removed'
    ENDIF 
    atmos_temp=atmos_temp[atmin:atmax]
    atmos_elec=atmos_elec[atmin:atmax]
    atmos_height=atmos_height[atmin:atmax]
    atmos_pressure=atmos_pressure[atmin:atmax]
    atmos_hyd=atmos_hyd[atmin:atmax]
    atmos_h1=atmos_h1[atmin:atmax]
    if n_elements(atmos_he2) gt 0 then begin
      atmos_he1=atmos_he1[atmin:atmax]
      atmos_he2=atmos_he2[atmin:atmax]
    endif
  ENDIF
ENDIF 

output={ filename: filename, $
         temp: atmos_temp, $
         elec_dens: atmos_elec, $
         height:  atmos_height,$
         pressure: atmos_pressure, $
         hyd_dens: atmos_hyd, $
         h_elec: atmos_hyd/atmos_elec, $
         h1_frac: atmos_h1, $
         comments: comments, $
         time_stamp: systime()}

IF n_elements(atmos_he1) NE 0 THEN output=add_tag(output,atmos_he1,'he1_frac')
IF n_elements(atmos_he2) NE 0 THEN output=add_tag(output,atmos_he2,'he2_frac')

return,output

END
