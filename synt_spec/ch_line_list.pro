;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;
;
; NAME:
;	CH_LINE_LIST
;
; PURPOSE:
;
;	Create a latex or an ascii file of predicted spectral line intensities and
;       wavelengths corresponding to  selected parameters, as calculated by 
;       CH_SYNTHETIC. Needs as input the line intensity structure calculated by
;       CH_SYNTHETIC (default)  or the SPECTRUM structure output of
;       MAKE_CHIANTI_SPEC.
;
; CALLING SEQUENCE:
;
;       IDL> ch_line_list, transitions, outname, latex=latex, ascii=ascii, $
;       abundfile=abundfile, min_abund=min_abund, $
;       wmin=wmin,wmax=wmax,$
;       SPECTRUM=SPECTRUM, minI=minI,photons=photons,kev=kev, $
;       all=all,no_sort=no_sort, sngl_ion=sngl_ion
;
;
; PROCEDURE:
;
;
; INPUTS:
;
;       The structure created by CH_SYNTHETIC
;
; OPTIONAL INPUTS:
;
;	Wmin:   lower limit of the wavelength/energy range of interest (Angstroms)
;               if kev keyword set, then wmin is in kev	
;	Wmax:   upper limit of the wavelength/energy range of interest (Angstroms)
;               if kev keyword set, then wmax is in kev	
;
;       Mini:   Minimum intensity for line to be included in output
;	
;	SNGL_ION:  specifies a single ion (or a list of ions) to be used instead
;                 of the complete set of ions specified in the structure.
;
;
;       MIN_ABUND:  If set, outputs  only  those elements which 
;                   have an abundance greater than min_abund.  
;
;                   For example, from Allen (1973):
;                   abundance (H)  = 1.
;                   abundance (He) = 0.085
;                   abundance (C)  = 3.3e-4
;                   abundance (O)  = 6.6e-4
;                   abundance (Si) = 3.3e-5
;                   abundance (Fe) = 3.9e-5
;
;
; KEYWORD PARAMETERS:
;
;       LATEX:  Create a latex file (default, exclusive with /ASCII)
;
;       ASCII:  Create an ascii file (exclusive with /LATEX)
;
;	MINI:	Minimum intensity for line to be included in output
;
;       PHOTONS:  units will be in photons rather than ergs
;
;       KEV:  wavelengths will be given in kev rather than Angstroms
;
;       ALL:  if set, then all lines are included.  This means that lines for which
;             only an approximate wavelength is known, denoted by a negative 
;             wavelength value in the .wgfa file, are included.
;             These lines are listed in the file with a * preceding the wavelength.
;
;       NO_SORT:
;             If set, then the lines are *not* sorted in wavelength (or energy).
;
;       SPECTRUM
;
;             If set, IT IS ASSUMED that the input structure is the SPECTRUM
;             structure output of MAKE_CHIANTI_SPEC, where the  line
;             intensities have already been multiplied by the abundance factor!!
;
;
; OUTPUTS:
;
;	A latex (default) or an ascii file with the line list
;
; OPTIONAL OUTPUTS:
;       Max_Int:  The maximum intensity printed in the table.
;
; CALLS: Many SolarSoft routines.
;
;
; COMMON BLOCKS:
;        none.
;
; SIDE EFFECTS:
;
;
; EXAMPLE:
;
;             > ch_line_list, trans,'linelist.tex',/latex, wmin=100.,wmax=200.,/all
;
;
; CATEGORY:
;
;	spectral synthesis.
;
;
; WRITTEN     : 
;       Version 1, Written by: Giulio Del Zanna (GDZ) Oct 31 2001.
;
; MODIFICATION HISTORY:
;
;        V.2, 9-Nov-2001 GDZ. 
;                 Now correctly handles the case when no
;                 abundances are passed to the routine. 
;
;        v.3, 11-Dec-2001, PRY.
;                 Removed calls to get_utc and anytim2cal. Replaced with 
;                 call to systime()
; 
;        v.4, 29-Apr-02, GDZ
;
;                 Fixed a few small bugs, some caused by a change in the
;                 database file format for V4. 
;                 Added only_mini,  file_effarea keywords to be able to use as
;                 input the structure created by MAKE_CHIANTI_SPEC.
;
;        V.5, 22-May-2002, GDZ
;                 generalized directory concatenation to work for
;                 Unix, Windows  and VMS. changed tags.
;                 Changed and added various things, including flabel
;
;        V.6, 12-Aug-02, GDZ
;           Modified the output labeling, and fixed two bugs: 1) when /all was used
;           the keyword /mini was not working. 2) min_abund was not working
;           properly when /spectrum was used.  Reduced size of latex output (was
;           12pt)
;           Changed output in isothermal case (no Tmax given). Better info printed (GDZ)
;
;        V.7, 3-Nov-03  GDZ
;           Modified format e8.2 to e9.2 for Windows compatibility.
;
;        v.8, 18-Jul-2005 GDZ
;           Modified the use of the /kev keyword. Also, now the
;           routine accepts input structure with the units in keV. 
;
;        v.9, 4-Aug-2005 GDZ
;           Corrected a bug introduced in the previous version.
;           Also switched to \documentclass when making the latex file.
;
;        v.10, 8-Jun-2022, Peter Young
;           Added comment to state if population lookup tables were used for
;           the calculation; added MAX_INT= optional output.
;
;
; VERSION     : 10, 8-Jun-2022
;
;-
pro ch_line_list, transitions, outname, latex=latex, ascii=ascii, $
      wmin=wmin,wmax=wmax,$
      SPECTRUM=SPECTRUM, abundfile=abundfile, min_abund=min_abund, $
      minI=minI,photons=photons,kev=kev, $
      all=all,no_sort=no_sort, sngl_ion=sngl_ion, max_int=max_int


;
IF  n_params(0) lt 2 then begin
   print,"  IDL>  ch_line_list,transitions, 'list.tex', /latex [ascii], $ "
   print,'   [abundfile= ,min_abund= ,mini= , wmin=, wmax=,'
   print,'    /SPECTRUM, sngl_ion= ,/all, /no_sort, /photons, /kev]'
   return
ENDIF

IF  keyword_set(latex) THEN BEGIN 
   IF keyword_set(ascii) THEN BEGIN 
      print, 'You can only create either an ascii or a latex file at a time'
      return
   END 
ENDIF ELSE BEGIN 
   IF NOT keyword_set(ascii) THEN BEGIN 
      print, '% CH_LINE_LIST: Assuming a latex output as default...'
      latex = 1
   END 
END 


;
ang2kev=12.39854d
;
if keyword_set(minI) then begin
   min_int=minI
endif else min_int=0.
;

ioneq_name = transitions.ioneq_name

break_file, ioneq_name, disk, dir,ioneq_file , ext
ioneq_file = ioneq_file+ext

ioneq_ref=transitions.ioneq_ref

wvl = transitions.lines(*).wvl
ident_latex = transitions.lines(*).ident_latex
ident_ascii = transitions.lines(*).ident

tmax = transitions.lines(*).tmax
intensity = transitions.lines(*).int 
;snote = transitions.lines(*).snote

;-1 for the theoretical lines, 0 for the observed ones:
flag = transitions.lines(*).flag  

int_units = transitions.int_units

wvl_units = transitions.WVL_UNITS

;check if the line intensities were calulated  with the isothermal approx.

IF tag_exist(transitions, 'DEM') THEN BEGIN 

   isothermal = 0
   dem_ref=transitions.dem_ref
   dem_name =   transitions.dem_name

   break_file, dem_name, disk, dir, dem_file , ext
   dem_file = dem_file+ext

ENDIF ELSE BEGIN 

   isothermal = 1
   iso_t = transitions.logt_isothermal
   iso_em = transitions.logem_isothermal

   iso_t = trim(iso_t)
   iso_em = trim(iso_em)
END 


nlines = n_elements(wvl)

; 'photons/cm2/sr/s'
;  'erg/cm2/sr/s'

IF  keyword_set(photons) THEN  BEGIN
   IF int_units EQ 'photons/cm2/sr/s' THEN BEGIN 
      print,  '% CH_LINE_LIST: Units already in ', int_units
   ENDIF ELSE BEGIN 
      IF int_units EQ 'erg/cm2/sr/s' THEN BEGIN 
         intensity = intensity/1.986e-8*wvl
         int_units = 'photons/cm2/sr/s'
         print, '% CH_LINE_LIST: Units converted to ', int_units
      ENDIF 
   ENDELSE 
ENDIF ELSE BEGIN 
   IF int_units EQ 'photons/cm2/sr/s' THEN BEGIN
;6.626d-27*2.998d+10*1.d+8
      intensity = intensity*1.986e-8/wvl
      int_units = 'erg/cm2/sr/s'
      print, '% CH_LINE_LIST: Units converted to ', int_units
   END 
ENDELSE 

IF n_elements(wmin) EQ 0 THEN wmin = min(wvl) ELSE wmin = float(wmin) > min(wvl)
IF n_elements(wmax) EQ 0 THEN wmax = max(wvl) ELSE wmax = float(wmax) < max(wvl)

IF  keyword_set(kev) THEN  BEGIN 

;
if wvl_units eq 'Angstroms' then begin 
   wvl=ang2kev/wvl
   w_min=ang2kev/wmax
   w_max=ang2kev/wmin
   wmin = w_min
   wmax = w_max
 wvl_units='keV' 
endif  ; else units already in keV

ENDIF 
;


IF NOT keyword_set(SPECTRUM) THEN BEGIN 

;
; multiply for the element abundances:
;
; Choose an abundance file if not passed to the routine
;

   IF n_elements(abundfile) EQ 0 THEN BEGIN 
      dir=concat_dir(!xuvtop, 'abundance')
      abund_name=ch_get_file(path=dir,filter='*.abund',title='Select Abundance File')
   ENDIF ELSE BEGIN 

      IF NOT file_exist(abundfile) THEN BEGIN 
         print, '% CH_LINE_LIST: Error, wrong ABUNDANCE file '
         return
      END 
      abund_name=abundfile
   END

;read abundance file

   read_abund,abund_name,abund,abund_ref

   break_file, ABUND_NAME, disk, dir,abund_file , ext
   abund_file = abund_file+ext

   line_abunds = abund[transitions.lines.iz-1]


; Multiply for the abundance factor:
;-----------------------------------

   intensity = intensity*line_abunds
   max_int=max(intensity)

;Check the min_abund.
;----------------------

   IF n_elements(min_abund) EQ 0  THEN BEGIN 
      index = where(line_abunds GT 0)
      min_abund = min(line_abunds(index))
   ENDIF 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   IF keyword_set(ALL) THEN $
     index = where((wvl GE wmin) AND (wvl LE wmax) AND $
                   (line_abunds GE  min_abund) AND (intensity GE min_int), nlines) ELSE $
     index = where((flag  EQ 0) AND (wvl GE wmin) AND (wvl LE wmax) AND $
                   (line_abunds GE  min_abund) AND (intensity GE min_int) , nlines)

ENDIF  ELSE BEGIN 

   break_file,transitions.ABUND_NAME, disk, dir,abund_file , ext
   abund_file = abund_file+ext
   abund_ref = transitions.abund_ref

;Check the min_abund.
;----------------------

   IF n_elements(min_abund) EQ 0  THEN BEGIN 
      line_abunds = transitions.abund[transitions.lines.iz-1]
      index = where(line_abunds GT 0)
      min_abund = min(line_abunds(index))
   ENDIF 

   IF tag_exist(transitions,'file_effarea') THEN BEGIN  
      file_effarea = transitions.file_effarea
      break_file,file_effarea, disk, dir, fef, ext
      fef = fef+ext
   END

   IF keyword_set(ALL) THEN $
     index = where((wvl GE wmin) AND (wvl LE wmax) AND $
                   (intensity GE min_int), nlines) ELSE $
     index = where((flag  EQ 0) AND (wvl GE wmin) AND (wvl LE wmax) AND $
                   (intensity GE min_int) , nlines)

END  



IF nlines EQ 0 THEN BEGIN 
   print,  '% CH_LINE_LIST: no lines ! '
   return
ENDIF ELSE BEGIN 
   w_min = min(wvl[index])
   w_max = max(wvl[index])
END 


IF n_elements(sngl_ion) GT 0 THEN BEGIN 
   n_match = 0

   flabel = transitions.lines[*].snote 
   FOR i=0,  n_elements(flabel) -1 DO BEGIN 
      spectroscopic2ion,transitions.lines[i].snote, dummy
      flabel[i] = dummy
   END 

   FOR  ilist=0, n_elements(sngl_ion)-1 DO  BEGIN
      index1 = where(flabel EQ sngl_ion[ilist], nn)
      n_match =n_match +nn
   ENDFOR 
   IF n_match EQ 0 THEN BEGIN 
      print,  '% CH_LINE_LIST: no lines for the selected ions ! '
      return
   ENDIF 
ENDIF 


IF NOT keyword_set(no_sort) THEN BEGIN
   index_sort = sort(wvl(index))
   index = index(index_sort)
END



;OPEN THE FILE

openw,luo,outname,/get_lun

;;;;;; HEADER ;;;;;

IF  keyword_set(latex) THEN BEGIN 

   printf,luo,'\documentclass[]{article}'
;
printf,luo, '\usepackage{longtable}'
printf,luo, ' '

   printf,luo,'% \setlength{\textheight}{250mm} %UK'
   printf,luo,'\setlength{\textheight}{220mm}   %US'
   printf,luo,'\setlength{\headheight}{0.0mm}'
   printf,luo,'\setlength{\headsep}{0.0mm}' 
   printf,luo,'\setlength{\topmargin}{-5mm}' 
   printf,luo,'\setlength{\textwidth}{165mm}'
   printf,luo,'\setlength{\oddsidemargin}{-0.5cm}' ;-1.0cm
;
   printf,luo,'\begin{document}'

   printf,luo,'\begin{center}'
   printf,luo,'\large\bf{Predicted XUV Line Intensities}\\'
   printf,luo,' CHIANTI database - Version '+transitions.version+'  \\'
   printf,luo,'\vspace{0.5cm}'

   CASE transitions.model_name OF
      'Constant pressure': BEGIN 
         dummy = 'Calculated with '+transitions.model_name+ '= '+$
           string(transitions.model_pe,'(e9.2)')+ ' (cm$^{-3}$ K) \\' 
      END 
      'Constant density': BEGIN 
         dummy = 'Calculated with '+transitions.model_name+ '= '+$
           string(transitions.model_ne,'(e9.2)')+ ' (cm$^{-3}$) \\'
      END 
      'Function':BEGIN 
         break_file, transitions.model_file, disk, dir, file,ext
         dummy = 'Calculated with the (Te,Ne) file: \verb|'+file+ext+ '| \\'
      END 
   ENDCASE

   printf,luo, dummy


   IF  wvl_units eq 'keV' THEN  $
     printf,luo,'\large\bf{'+strtrim(string(w_min,'(f6.3)'),2)+$
     ' to '+strtrim(string(w_max,'(f6.3)'),2)+' keV}\\' ELSE if wvl_units eq 'Angstroms' then $
     printf,luo,'\large\bf{'+strtrim(string(w_min,'(f6.1)'),2)+$
     ' to '+strtrim(string(w_max,'(f6.1)'),2)+' \AA}\\'

   printf,luo,'\large\bf{Number of lines:  '+strtrim(string(n_elements(index),'(i7)'),2)+'}\\'

   printf,luo,' Minimum intensity = '+string(min_int)+'  \\'
   printf,luo,'\large\bf{Units are: '+int_units+' }\\'

   IF n_elements(fef) GT 0 THEN $
     printf,luo,'\large\bf{Effective area file:} \verb|'+fef+'| \\'


   IF keyword_set(ALL) THEN BEGIN 
      printf,luo,'\large\bf{Lines marked with a * do not have correspondent observed energy levels} \\'
      printf,luo,'\large\bf{and have approximate wavelengths.} \\'
   END

   utc=systime()
   printf,luo,'Calculated: '+utc+ ' \\'

;
;printf,luo,'\today\\'

;
;printf,luo,'{\bf Ions included:}\\'
;printf,luo,elstage
;printf,luo,'\\'
;
;
   printf,luo,'\vspace{0.5cm}'
;
   printf,luo,'{\bf Ionization Fractions file:} \verb|'+ ioneq_file +'| \\'
; get rid of %'s and _'s which kill you in latex
;
   asize=n_elements(ioneq_ref)
   for iref=0,asize-1 do begin
      reflen=strlen(ioneq_ref(iref))
      ipc=strpos(ioneq_ref(iref),'%')
      thisref=strmid(ioneq_ref(iref),ipc+1,reflen-ipc-1)
      thisref=repstr(thisref,'&','\&')
      thisref=repstr(thisref,'_','\_')
      IF trim(thisref) NE '' THEN   printf,luo,thisref+' \\'
   endfor
;
;
   printf,luo,'   \vspace{0.5cm}'
;

   printf,luo,'{\bf Elemental Abundance file:} \verb|'+ abund_file +'| \\'
;
;
; get rid of %'s and _'s which kill you in latex
;
   asize= n_elements(abund_ref) ; size(abund_ref)

;for iref=0,asize(1)-1 do BEGIN
   for iref=0,asize-1 do BEGIN
      reflen=strlen(abund_ref(iref))
      ipc=strpos(abund_ref(iref),'%')
      thisref=strmid(abund_ref(iref),ipc+1,reflen-ipc-1)
      thisref=repstr(thisref,'&','\&')
      thisref=repstr(thisref,'_','\_')
      IF trim(thisref) NE '' THEN printf,luo,thisref+' \\'
   endfor
;
;
   printf,luo,' \vspace{0.2cm}'

   printf,luo,'Minimum abundance = '+string( min_abund)+ '  \\'

   printf,luo,' \vspace{0.5cm}'
;

   IF NOT isothermal THEN BEGIN 

      printf,luo,'{\bf Differential Emission Measure file:} \verb|'+$
        dem_file +'| \\'
;
      asize=n_elements(dem_ref)
      for iref=0,asize-1 do begin
         reflen=strlen(dem_ref(iref))
         ipc=strpos(dem_ref(iref),'%')
         thisref=strmid(dem_ref(iref),ipc+1,reflen-ipc-1)
         thisref=repstr(thisref,'&','\&')
         thisref=repstr(thisref,'_','\_')
         IF trim(thisref) NE '' THEN  printf,luo,thisref+' \\'
      endfor
;
   ENDIF ELSE BEGIN 

      printf,luo,'{\bf calculated with  isothermal approximation: }\\'
      printf,luo, 'log$_{10}$ T='+ arr2str(iso_t, ',',/trim)+'\\'
      printf,luo, 'log$_{10}$ EM='+ arr2str(iso_em, ',',/trim)+'\\'

    END

   IF transitions.lookup EQ 1 THEN BEGIN
     printf,luo,'\vspace{0.5cm}'
     printf,luo,'\textbf{Calculation performed with population lookup tables.}'
   ENDIF 


   printf,luo,'\vspace{0.5cm}'
;
;printf,luo,'\begin{center}'

   printf,luo,'\end{center}'
;
   printf,luo,'\newpage'

   printf,luo,'\begin{longtable}[c]{@{}|lrlcr|@{}}'
   printf,luo,'\caption{\em \normalsize Line List } \label{tab:obs} \\'
   printf,luo,'\hline'

   if  wvl_units eq 'keV' then begin
      printf,luo,' Ion & Energy (KeV)  & Transition & T$_{\rm max}$ & Int\\ \hline'
   endif else if wvl_units eq 'Angstroms' then begin
      printf,luo,' Ion & $\lambda$ (\AA) & Transition & T$_{\rm max}$ & Int\\ \hline'
  end 
;
   printf,luo,'\hline \endfirsthead'
   printf,luo,'\caption[]{\normalsize (continued)}\\ '
   printf,luo,'\hline '
   if wvl_units eq 'keV' then begin
      printf,luo,' Ion & Energy (keV)  & Transition  & T$_{\rm max}$ & Int\\ \hline'
   endif else if wvl_units eq 'Angstroms' then begin
      printf,luo,' Ion & $\lambda$ (\AA) & Transition & T$_{\rm max}$ & Int\\ \hline'
   end
   printf,luo,'\hline \endhead'
;
;

ENDIF ELSE IF keyword_set(ascii) THEN BEGIN 


   printf,luo,'Predicted XUV Line Intensities'
   printf,luo,'CHIANTI database - Version '+transitions.version

   CASE transitions.model_name OF
      'Constant pressure': BEGIN 
         dummy = 'Calculated with '+transitions.model_name+ '= '+$
           string(transitions.model_pe,'(e9.2)')+ ' (cm-3 K) '
      END 
      'Constant density': BEGIN 
         dummy = 'Calculated with '+transitions.model_name+ '= '+$
           string(transitions.model_ne,'(e9.2)')+ ' (cm-3) '
      END
      'Function':BEGIN 
         break_file, transitions.model_file, disk, dir, file,ext
         dummy = 'Calculated with the (Te,Ne) file: '+file+ext
      END 
   ENDCASE

   printf,luo, dummy

   IF  wvl_units eq 'keV' THEN  $
     printf,luo,strtrim(string(w_min,'(f6.3)'),2)+' to '+strtrim(string(w_max,'(f6.3)'),2)+' keV' ELSE $
     if wvl_units eq 'Angstroms' then $
printf,luo,strtrim(string(w_min,'(f6.1)'),2)+' to '+strtrim(string(w_max,'(f6.1)'),2)+' Angstroms'

   printf,luo,'Number of lines:  '+strtrim(string(n_elements(index),'(i7)'),2)

   printf,luo,'Minimum intensity = '+string(min_int)
   printf,luo,'Units are:  '+int_units

   IF n_elements(fef) GT 0 THEN $
     printf,luo,'Effective area file: '+fef

   IF keyword_set(ALL) THEN BEGIN 
      printf,luo,'Lines marked with a * do not have correspondent observed energy levels '
      printf,luo,'and have approximate wavelengths.'
   END 

   get_utc,utc
   utc = anytim2cal(utc)
   printf,luo,'Calculated: '+utc

;
;printf,luo,'{\bf Ions included:}\\'
;printf,luo,elstage

   printf,luo,''

   printf,luo,'Ionization Fractions file: '+ ioneq_file

; get rid of %'s and _'s
   asize=n_elements(ioneq_ref)
   for iref=0,asize-1 do begin
      reflen=strlen(ioneq_ref(iref))
      ipc=strpos(ioneq_ref(iref),'%')
      thisref=strmid(ioneq_ref(iref),ipc+1,reflen-ipc-1)
      thisref=repstr(thisref,'&','\&')
      thisref=repstr(thisref,'_','\_')
      IF trim(thisref) NE '' THEN printf,luo,thisref
   endfor

   printf,luo,''

   printf,luo,' Elemental Abundance file: '+ abund_file 

; get rid of %'s and _'s

   asize= n_elements(abund_ref) ; size(abund_ref)
   for iref=0,asize-1 do BEGIN
      reflen=strlen(abund_ref(iref))
      ipc=strpos(abund_ref(iref),'%')
      thisref=strmid(abund_ref(iref),ipc+1,reflen-ipc-1)
      thisref=repstr(thisref,'&','\&')
      thisref=repstr(thisref,'_','\_')
      IF trim(thisref) NE '' THEN  printf,luo,thisref
   endfor

   printf,luo,'Minimum abundance = '+string(min_abund)


   printf,luo,''

   IF NOT isothermal THEN BEGIN 

      printf,luo,'Differential Emission Measure file: '+dem_file

      asize=n_elements(dem_ref)
      for iref=0,asize-1 do begin
         reflen=strlen(dem_ref(iref))
         ipc=strpos(dem_ref(iref),'%')
         thisref=strmid(dem_ref(iref),ipc+1,reflen-ipc-1)
         thisref=repstr(thisref,'&','\&')
         thisref=repstr(thisref,'_','\_')
         IF trim(thisref) NE '' THEN printf,luo,thisref
      endfor

   ENDIF ELSE BEGIN 

      printf,luo,'calculated with  isothermal approximation:'
      printf,luo, 'log_10 T='+ arr2str(iso_t, ',',/trim)
      printf,luo, 'log_10 EM='+ arr2str(iso_em, ',',/trim)

    END

   IF transitions.lookup EQ 1 THEN BEGIN
     printf,luo,''
     printf,luo,'Calculation performed with population lookup tables.'
   ENDIF 

   printf,luo,''

   
   if  wvl_units eq 'keV' then begin
      printf,luo,'Energy (keV)   Intensity  Ion        Tmax  Transition '
   endif else if wvl_units eq 'Angstroms' then begin
      printf,luo,'Wavelength (A) Intensity  Ion        Tmax  Transition '
  END 

   printf,luo,''

END                             ;ascii header 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


for j=0L,nlines-1 do begin

   i = index(j)

;check that we want this ion:
;---------------------------

   IF n_elements(sngl_ion) GT 0 THEN $
     index1 = where(sngl_ion EQ flabel[i], nn) ELSE $
     nn=1

   IF nn EQ 1 THEN BEGIN 

      if wvl_units eq 'keV' then begin
         wvls=strtrim(string(wvl(i),'(f14.6)'),2)
      endif else if wvl_units eq 'Angstroms' then begin
         wvls=strtrim(string(wvl(i),'(f12.4)'),2)
      END

;if the line does not have observed energy level, then the 
; wavelength is not accurate, and a * is added to the ion designation.
;-------------------------------------------------------

      ion = transitions.lines(i).snote 
      if flag(i)  EQ -1  then ion=ion+' *'

      IF   isothermal THEN  item1='-' ELSE   item1=string(tmax(i),'(f4.1)')

;
      if keyword_set(photons) then begin
         item2=string(intensity(i),'(e10.2)')
      endif else begin
         item2=string(intensity(i),'(e9.2)')
      endelse
;

      IF  keyword_set(latex) THEN $ 
        printf,luo, ion+' & '+wvls+' & ' +ident_latex(i)+' & '+item1+' & ' +item2+ ' \\' ELSE $
        IF  keyword_set(ascii) THEN $
        printf,luo, strpad(wvls,12,/after)+strpad(item2, 10,/after)+strpad(ion,14,/after)+$
        strpad(item1, 7,/after) +  trim(ident_ascii(i))
;

   ENDIF 
ENDFOR

IF  keyword_set(latex) THEN BEGIN 
   printf,luo,'\hline'
   printf,luo,'\end{longtable}'
   printf,luo,'\end{document}'
ENDIF 

free_lun,luo

print, ' The file '+outname+' has been created ! '

IF  keyword_set(latex) THEN print, ' Now latex three times the file '

END 
