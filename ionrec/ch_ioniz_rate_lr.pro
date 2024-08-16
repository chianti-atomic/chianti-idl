;+
;
; PROJECT:  CHIANTI
;
;        This program was developed as part of CHIANTI-VIP. 
;        CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;        mantained by Giulio Del Zanna, to develop additional features and
;        provide them to the astrophysics community.
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the
;       University of Cambridge, Goddard Space Flight Center, and University of Michigan.
;
;
; NAME:
;	
;       CH_IONIZ_RATE_LR 
;
;
; PURPOSE:
;
;       Read level-resolved rate coefficient files and interpolate over requested temperature grid	
;
;
; EXPLANATION
;
;       Returns the level-resolved (lr) rate coefficients using interpolation
;       (in log_10) from the values in the file containing the ionization rate coefficients,
;       for either the charge transfer (CT) ionization/recombination or electron impact
;       direct/indirect ionization.
;
;       NOTE: if the temperature grid is outside the temperature range of the rates,
;             estimates are given by using the rates at the first and last temperatures.
;
;             The CT ionization files have the suffix '.ctilvl'
;             The CT recombination files have the suffix '.ctrlvl'
;             The direct impact (DI) ionization files have the suffix '.dilvl'
;             The excitation-autoionization (EA) ionization files have the suffix '.ealvl'
;
;
; CATEGORY
;
;       CHIANTI; atomic rates; ionization; recombination; charge transfer
;
;
; CALLING SEQUENCE:
;
;       IDL> logt=findgen(41)/10.+3.0
;       IDL> t=10.^logt
;       IDL> rt=ch_ioniz_rate_lr(!xuvtop+'/c/c_1/c_1.dilvl',t)
;       IDL> rt=ch_ioniz_rate_lr(!xuvtop+'/c/c_1/c_1.ctilvl',t,/ct,/verbose)
;
;
; INPUTS:
;
;       filename:  the ascii file with the level-resolved rate coefficients.
;
;
; KEYWORDS:
;
;       CT:  read charge transfer (CT) rate coefficients. This must be specified
;            for CT files because without it the routine will still read the file but
;            will return spurious values.
;
;       VERBOSE:  give additional messages
;
;
; OPTIONAL INPUTS:
;
;       temp_in:  the temperatures can be a single value or array.
;
;
; OUTPUTS:
;
;       filename:  file used to obtain rates
;
;       temp:  temperatures for the requested rates, interpolated or otherwise
;
;       rates:  the rate coefficients
;
;       lev_i:  index of initial level before transition
;
;       lev_f:  index of final level after transition
;
;       ion_pot:  the energy difference between initial and final levels of transition
;
;
; OPTIONAL OUTPUTS:
;
;       perturber_z:  atomic number of the perturber involved in CT collision
;
;       perturber_elecs:  number of electrons in perturber before CT collision
;
;
; CALLS:
;
;       INTERPOL
;
;
; PREVIOUS HISTORY:
;
;       NONE
;
;
; WRITTEN:
;         
;       Giulio Del Zanna (GDZ) and Roger Dufresne (RPD)
;       DAMTP, University of Cambridge, 16 Sept 2023 
;
;
; MODIFIED:
;
;       v.2, 17 Oct 2023  GDZ, fixed a typo when setting the rates above the maximum T.
;
;       v.3, GDZ, using a linear interpolation in the log.
;
;       v.4, 29 Feb 2024  RPD, output only in verbose mode warning if interpolation temperature
;               range is outside the range of temperatures in file
;
;       v.5, 30 Apr 2024  RPD, allowed rates to be retrieved from file without interpolation,
;               so that temperature is now a keyword. Fixed typo with interpolating rates
;               when temp_in is a scalar.
;
;
; VERSION:  5
;
;- 

function ch_ioniz_rate_lr,filename,temp_in=temp_in,ct=ct,error=error,verbose=verbose


  error=0
  if n_elements(verbose) eq 0 then verbose=0 

  if not file_exist(filename) then begin 
     print,'% CH_IONIZ_RATE_LR: ERROR, the file '+$
           filename+' was not found ! '
     error=1
     return,-1
  endif 


  ; begin reading rate file
  s=''
  openr,lur,filename,/get_lun

  ; start with the temperatures on the first line
  readf,lur,s

  t=double(str_sep(s,' ',/trim))
  t=t[where(t gt 0.)]           ; IDL fix to remove blanks
  temp_file=t
  logt_file=alog10(t)
  nt_file=n_elements(logt_file)

  
  if n_elements(temp_in) gt 0 then begin

    interp=1
    logt_in=alog10(temp_in)
    nt_in=n_elements(temp_in)

    ; check range:
    if keyword_set(verbose) then $
      if min(logt_in) lt min(logt_file) or max(logt_in) gt max(logt_file) then $
        print,'% CH_IONIZ_RATE_LR: Requested Te outside temperature range in file - will be returning estimated values'

  endif else begin
    
    interp=0
    temp_in=temp_file
  
  endelse
  
  ; now determine the remaining number of lines in the file and set up the various arrays
  nlines=0
  read_file=1
  ; the following works if only one -1 is at the end of the file:
  while read_file eq 1 do begin
    readf,lur,s
    nlines=nlines+1
    if trim(s) eq '' then print,'error:',nlines
    if trim(s) eq -1 then read_file=0
  endwhile 

  nlines=nlines-1             
  if keyword_set(verbose) then print,'% CH_IONIZ_RATE_LR: number of lines in the file:',nlines

  lev_i=intarr(nlines)           ; initial level number 
  lev_f=intarr(nlines)           ;  final level number    
  rates_file=dblarr(nt_file,nlines)
  if interp eq 1 then rates_int=dblarr(nt_in,nlines) ; for interpolated values
  ion_pot=dblarr(nlines)         ; ionization potential, i.e. DE in eV
  
  if keyword_set(ct) then begin 
    perturber_z=intarr(nlines)
    perturber_elecs=intarr(nlines)
  endif
  
  
  ; read the remaining data in the file
  s=''
  point_lun,lur,0
  readf,lur,s                   ; read temperatures
  
  a=0 & b=0 & c=0. & d=0 & e=0. & r=dblarr(nt_file)

  for ii=0,nlines-1 do begin

    readf,lur,s
     
    if keyword_set(ct) then begin
      reads,s,a,b,c,d,e,r

      lev_i[ii]=a
      lev_f[ii]=b
      perturber_z[ii]=c
      perturber_elecs[ii]=d
      ion_pot[ii]=e
    endif else begin
      reads,s,a,b,c,r

      lev_i[ii]=a
      lev_f[ii]=b
      ion_pot[ii]=c
    endelse

    ; read in the rates and interpolate if requested
    if interp eq 0 then rates_file[*,ii]=r $
      
    else begin

      ; avoid zeros for interpolating in log
      bad=where(r le 1.0d-50,nbad)
      if nbad gt 0 then r[bad]=1.0d-50
      
      rates_file[*,ii]=r
      
      ; interpolate rates only within file temperature range, otherwise copy
      ; the first and last values in the file for those outside the range
      ind1=where(logt_in lt min(logt_file),n1)
      ind2=where(logt_in gt max(logt_file),n2)
      ind3=where(logt_in ge min(logt_file) and logt_in le max(logt_file),n3)
      
      if n3 gt 0 then rates_int[ind3,ii]=$
        10.0d0^interpol(alog10(rates_file[*,ii]),logt_file,logt_in[ind3]) 
      if n1 gt 0 then rates_int[ind1,ii]=r[0] ; take first 
      if n2 gt 0 then rates_int[ind2,ii]=r[-1] ; take last

    endelse
    
  endfor

  free_lun,lur

  
  if interp eq 0 then rates_int=rates_file $
  
  ; resolve possible interpolation errors
  else begin  
    bad=where(rates_int lt 0.0,nbad)
    if nbad gt 0 then rates_int[bad]=1.0d-50

    ; check for NaN values 
    bad=where(rates_int ne rates_int,nbad)
    if nbad gt 0 then begin 
       print,'found '+trim(nbad)+' NaN values'
       rates_int[bad]=1.0d-50
    endif
  endelse
  
  
  if keyword_set(ct) then $
    out={filename:filename, temp:temp_in, rates:rates_int, lev_i:lev_i,$
      lev_f:lev_f, perturber_z:perturber_z, perturber_elecs:perturber_elecs, ion_pot:ion_pot} $
  else out={filename:filename,temp:temp_in,rates:rates_int,$
      lev_i:lev_i, lev_f:lev_f, ion_pot:ion_pot}

      
  return,out

end
