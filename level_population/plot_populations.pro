;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. 
;
;
; NAME:
;	PLOT_POPULATIONS
;
; PURPOSE:
;	plot the population of a number of the lowest levels as a function of 
;       electron density for a specific temperature
;
; CATEGORY:
;
;	science.
;
; CALLING SEQUENCE:
;
;       PLOT_POPULATIONS,Ion,T,Nlevels
;
;
; INPUTS:
;
;       Ion:  CHIANTI style name for the ion, i.e., 'c_6' for C VI
;       T:  electron temperature (K) at which the populations are calculated.
;
;
; OPTIONAL INPUT:
;
;       NLEVELS: the maximum number of levels displayed. If not
;                passed, all the levels are displayed.
;
;	LEVELS:  a string specifying the levels to be displayed. E.g.
;                levels='1-5, 8,9,10' or levels='1-10' or levels='1,2,3,7,25' 
;
;       DENSITIES: the array  of electron densities at which the
;                  populations are calculated. 
;
;       Ioneq_File: The name of an ionization equilibrium file, which
;                   is used to populate the level-resolved
;                   ionization/recombination data structure. If not
;                   specified then the file !ioneq_file is used.
;
;       Radtemp: If photon excitation is included (by defining RPHOT),
;                then this input specifies the blackbody radiation
;                temperature in K. If not specified, then it is set to
;                6000 K.
;       Rphot:   Distance from the centre of the star in stellar radius units.
;                That is, RPHOT=1 corresponds to the star's
;                surface. If RPHOT is not specified, then photon
;                excitation will be switched off when pop_solver is
;                called.
;       Path:    This directly specifies the path where the
;                ion's data files are stored. If not set, then
;                the files are taken from the user's CHIANTI
;                distribution. 
;
;       Abund_File: The name of a CHIANTI format element abundance
;                   file. This would be used by pop_solver to compute
;                   the proton-to-electron ratio.
;
;       title:   the title for the plot.
;
;       chars:   the charsize for the plot, default=1.6
;
;       char_label the charsize of the labels, default=1.6
;
;
; KEYWORDS:
;
;       NOIONREC: If set, then level-resolved ionization/recombination
;                rates are not read for the ion, even if they exist.
;
;       NO_RREC: If set, do not read the level-resolved radiative
;                recombination (RR) files.
;       NOPROT:  If set, then proton rates are not read for the ion,
;                even if they exist.
;
;
; OUTPUTS:
;       A plot with the populations, and the option to save it. 
;
; OPTIONAL OUTPUTS:
;                OUTFILE: if specified, an ascii file with the
;                populations is written.
;
; COMMON BLOCKS:
;
;	common elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
;       common wgfa, wvl,gf,a_value
;       common upsilon,t_type,c_ups,splups
;
; CALLS
;       conv_levs
;
;
; EXAMPLE:
;
; to plot populations of the 5 ground configuration levels of Fe XIII
; and store these values in a file 'Fe_XIII.lis' for a temperature of 1.5 MK
;
;             > plot_populations,'fe_13',1.5e+6,5,outfile='Fe_XIII.lis'
;
; MODIFICATION HISTORY:
;
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;       November 1997:  Ken Dere, added psfile keyword
;       September 1999:  Ken Dere, for Version 3, 
;       14-Jul-2000     Peter Young, now calls pop_solver
;	Version 6, 21-Dec-2000, William Thompson
;		Modified for better cross-platform capability.
;
;       V. 7, 21-May-2002, Giulio Del Zanna (GDZ), modified the plotting bit
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;       V 8, 15-July-2002, GDZ 
;         Added keyword  not_interactive
;
;       V 9, 4-Oct-2003, GDZ
;                  modified the input to POP_SOLVER, now it is dealt with an
;                  input structure.
;
;       V 10, 4-May-2005, Enrico Landi (EL)
;                  Modified in order to include ionization and recombination
;                  data in the input to POP_SOLVER, now it allows choice of .ioneq
;                  file needed to include ionization and recombination.
;
;       V 11, 12-Aug-2005, GDZ 
;                  The number of levels is now optional. Also, a check that the
;                  number of levels must be less than the levels in the file is
;                  now enforced.
;
;       V 12, 12-Jun-2009, Enrico Landi
;                  Changed the definition of the temperature array for ion fractions
;                  in the IONREC variable, now taken directly from the output of
;                  READ_IONEQ.PRO
;
;       v 13, 1-May-2014, GDZ 
;                  replaced the splups with the v.8 scups files.
;
;       v 14, 21-Jul-2016, Peter Young
;                  fixed bug when loading splstr structure; added
;                  input ioneq_file=; routine now uses !ioneq_file by
;                  default (no longer need to use widget to select file).
;
;       v.15, 20 Nov 2018, Giulio Del Zanna
;
;             **** removed commons and modified completely the routine.
;
;       v.16,  14-Dec-2018 GDZ, added  keywords and the function conv_levs
;
; VERSION     :  16
;
;-

FUNCTION conv_levs, levels

level_numbers = -1

str = str_sep(levels, ',', /trim)

FOR i=0, n_elements(str)-1 DO BEGIN 

   dummy = fix(str_sep(str[i], '-', /trim))


   IF  n_elements(dummy) EQ 1 THEN level_numbers =[level_numbers, dummy[0]] ELSE $
     IF   n_elements(dummy) EQ 2 THEN BEGIN 

      level_numbers =[level_numbers,dummy[0]+indgen(dummy[1]-dummy[0]+1)]

   ENDIF ELSE BEGIN 
      print, 'Error in the levels parameter.... EXIT'
      return, -1
   END 

ENDFOR 

level_numbers =level_numbers[1:*]

return, level_numbers

END   


pro plot_populations, gname,t,nlevels,levels=levels,densities=densities,$
                      outfile=outfile, $
                      noprot=noprot,radtemp=radtemp,rphot=rphot, $
                      PATH=path, abund_file=abund_file, title=title, $
                      ioneq_file=ioneq_file, chars=chars,char_label=char_label,$
                      noionrec=noionrec, no_rrec=no_rrec


;
  if n_params(0) lt 2 then begin
     print,'   IDL> plot_populations,ion,temp [,n_levels, psfile= ,outfile= ,'
     print,'                radtemp= , rphot=  ]'
     print," i.e. > plot_populations,'c_5',1.e+6,4"
     return
  endif



;  This loads up the ion's atomic data  

  input=ch_setup_ion(gname,rphot=rphot,radtemp=radtemp,noprot=noprot, $
                     ioneq_file=ioneq_file,abund_file=abund_file,path=path, $
                     quiet=quiet, noionrec=noionrec, no_rrec=no_rrec )


; either the nlevels or levels should be defined.

  IF n_elements(nlevels) EQ 0 and n_elements(levels) eq 0 THEN begin  
     nlevels = n_elements(input.ecm)
     level_numbers = indgen(nlevels)+1 
  endif else if n_elements(nlevels) EQ 0 and n_elements(levels) gt 0 THEN begin 
     level_numbers = conv_levs(levels)
     if level_numbers[0] eq -1 then message,'Error - returning'

  endif else if n_elements(nlevels) gt 0 and n_elements(levels) eq 0 THEN begin
     nlevels = nlevels < n_elements(input.ecm)
     level_numbers = indgen(nlevels)+1 
  endif else begin 
     print, '% ERROR, either the maximum number of levels or the levels need to be defined !'
     return
  end 


  term=input.elvlcstr.data.conf+' '+input.elvlcstr.data.level

  term_tex= input.elvlcstr.data.conf_latex+' '+input.elvlcstr.data.level_latex 

; if integer
  nl=n_elements(input.elvlcstr.data)

  s_strings=strarr(nl)
  conf=strarr(nl)

  for ii=0,nl-1 do begin 
     if input.elvlcstr.data[ii].mult eq fix(input.elvlcstr.data[ii].mult) then $
        s_strings[ii]= trim(input.elvlcstr.data[ii].mult) else $
           s_strings[ii]= trim(2*input.elvlcstr.data[ii].mult)+'/2'

     conf[ii]=convert_config(input.elvlcstr.data[ii].conf, /idl)

  endfor

  term_idl = conf+ $
             ' !E'+trim(s_strings)+'!N'+trim(input.elvlcstr.data.L_SYM)+$
             '!I'+trim(input.elvlcstr.data.j_str)+'!N'
  

  print,' ' 
  print,' Levels: '
  FOR  i=0,n_elements(level_numbers)-1 DO $
     print,level_numbers[i],strpad(term(level_numbers[i]-1),30,/after),format='(i5,4x,a30)'


; define the densities

  IF n_elements(densities) eq 0 then densities=10.^(indgen(10)+4)

  nd=n_elements(densities)

;  calculate level populations
;------------------------------

  pop_solver, input, t , densities ,pop, radfunc=radfunc,verbose=verbose
                                ;

  pop = REFORM(pop)
  pops = TRANSPOSE(pop[*,0:nlevels-1])

;
  print,' '
  print,'Density       Populations'
  print,' '

  FOR i=0,nd-1 DO begin 

     print,' '
     print, 'Density: '+trim(densities[i]) 

     FOR j=0, nlevels-1 DO $
        print, string(trim(level_numbers[j]), pops[level_numbers[j]-1,i], $
                      format='(" ",i4,": ",500e10.3)')

  endfor 

;
;; print,' '
;; print,'Density       Log10 Populations'
;; print,' '
;; FOR i=0,nd-1 DO print,format='(e10.3,"   ",500f10.4)',densities[i],alog10(pops[*,i])

;; print, ''
;; print,' Levels: '
;; for i=0,nlevels-1 do print,i+1,strpad(term(i),30,/after),format='(i5,4x,a30)'

  ion2spectroscopic, gname,species

  IF n_elements(title) EQ 0 THEN title=species+' T='+trim(t)+' K '

  xtitle='Electron Density (cm!u-3!n)'
  ytitle='Fractional Population'

  if n_elements(chars) eq 0 then chars=1.6
  if n_elements(char_label) eq 0 then char_label=1.6


  window, 0
  circle_sym,/fill

  dx1 = min(densities)*0.1
  dx2 = max(densities)*0.25

  x_min=min(densities)-2*dx1
  x_max=max(densities)+2*dx2

  y_min=1.e-5
  y_max=1.2

  pr=''
begin_post:

  IF pr EQ 'y' THEN BEGIN 
     thick = 3.
     lthick = 3.
  ENDIF  ELSE BEGIN 
     thick =1.
     lthick = 1.
  END 

  !p.multi=0

  rand=randomu(seed,nd)
  isrt=sort(rand)               ; to try to fix problems with labeling curves

  plot_oo,densities,pops(0,*),psym=-8,$
          xrange=[x_min,x_max], yrange=[y_min,y_max],$
          xstyle=1,ystyle=1, title=title, $
          xtitle=xtitle,ytitle=ytitle, /nodata, $
          thick=thick,xthick=lthick,ythick=lthick, CHARthick=thick,chars=chars


;xyouts, 0.1, 0.9 ,/norm, ' Levels: '

  dy = (y_max-y_min)*0.05
  dyb = dy*0.5*0.1

  int =  pops[level_numbers-1, *]

  left_y = 0
  right_y = 0
  
  num_plot=0
  
  FOR  i=0 ,n_elements(level_numbers)-1 DO  BEGIN 

     oplot, densities, int[i, *], thick=lthick

     indx = where(densities LE x_max AND $
                  densities GE x_min AND $
                  int[i,*] GE y_min AND int[i,*] LE y_max, nn)


     IF nn GE 2 THEN BEGIN 

        left_x =  densities[indx[0]] ; -1.2*dx1
        right_x = densities[indx[nn-1]]


        IF i NE 0 THEN BEGIN 
           dummy = where( abs(left_y- int[i, indx[0]]) LT dyb , nov)
           IF nov GT 0 THEN xlt = left_x - nov*dx1*0.3 ELSE xlt = left_x
        ENDIF ELSE  xlt = left_x

        left_y =[left_y, int[i, indx[0]]]

        XYOUTS, xlt , int[i, indx[0]], trim(level_numbers[i]), $
                charsize=char_label, CHARthick=thick

        IF i NE 0 THEN BEGIN 
           dummy = where( abs(right_y- int[i, indx[nn-1]]) LT dyb, nov)
           IF nov GT 0 THEN xrt =right_x +nov*dx2*0.7 ELSE  xrt=right_x
        ENDIF ELSE  xrt=right_x

        right_y =[right_y, int[i, indx[nn-1]] ]

        XYOUTS, xrt, int[i, indx[nn-1]], ' '+trim(level_numbers[i]),$
                charsize=char_label, CHARthick=thick

        xyouts, 0.2, 0.8-0.04* num_plot,/norm,$
                string(level_numbers[i])+'  '+$
                strpad(term_idl(level_numbers[i]-1),30,/after), CHARthick=thick,chars=char_label
        
        num_plot=num_plot+1
        
     ENDIF  ELSE IF nn EQ 1 THEN BEGIN 

        xc =  densities[indx[0]]

        XYOUTS, xc,  int[i, indx[0]], trim(level_numbers[i]), $
                charsize=char_label, CHARthick=thick

     END 


  ENDFOR

  print2d_plot, pr=pr , x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
                go_to_line=go_to_line, out_name=out_name , /ask_name

  if go_to_line eq 'y' then goto,begin_post

;
; write a file 

  if n_elements(outfile) gt 0 then begin

     openw,luo,outfile,/get_lun
;

     printf,luo,' Temperature=',t
     printf,luo, ' '
;
     for i=0,nd-1 do begin
        printf,luo,densities(i),pops(*,i)
     endfor
;
     printf,luo, ' '
;
     for i=0,nd-1 do begin
        print,densities(i),alog10(pops(*,i))
        printf,luo,densities(i),alog10(pops(*,i))
     endfor
;
     free_lun,luo
;
  endif
;
end
