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
;       T:  electron temperature (K)
;
;
; OPTIONAL INPUT:
;       Nlevels:  the number of levels for which populations are plotted
;                 starts from level 1 (the ground level)
;       Ioneq_file: The name of an ionization equilibrium file. If not
;                 set, then the CHIANTI default file (!ioneq_file) is used.
;	
; KEYWORD PARAMETERS:
;
;	OUTFILE:   the (optional) name of the output file where the listing 
;                  is produced
;	PSFILE:    the (optional) name of the output file where a postscript 
;                  plot produced
;
;	RADTEMP	Specify background radiation temperature (default: 6000 K)
;
;	RPHOT   Distance from the centre of the star in stellar radius units.
;               I.e., RPHOT=1 corresponds to the star's surface. (Default is
;               infinity, i.e., no photoexcitation.)
;
;       NOPROT  If set, then the default setting will be NOT to use 
;               proton rates. This can be changed within the routine.
;
;               NOT_INTERACTIVE Avoid interactive use.
;
; OUTPUTS:
;
;
; COMMON BLOCKS:
;
;	common elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
;       common wgfa, wvl,gf,a_value
;       common upsilon,t_type,c_ups,splups
;
; CALLS
;
;       POP_SOLVER, ION2SPECTROSCOPIC, ION2FILENAME, READ_IP,
;       CONVERTNAME, READ_ELVLC, READ_WGFA2, READ_SCUPS
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
; VERSION     :  14, 21-July-2016
;
;-
pro plot_populations, gname,t,nlevels,outfile=outfile,psfile=psfile, $
                     noprot=noprot,radtemp=radtemp,rphot=rphot, $
                     not_interactive=not_interactive, ioneq_file=ioneq_file
;
;
common elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
common wgfa, wvl,gf,a_value
common upsilon,splstr
COMMON proton, pstr, pe_ratio
COMMON radiative, radt, dilute

;
if n_params(0) lt 2 then begin
  print,'   IDL> plot_populations,ion,temp [,n_levels, psfile= ,outfile= ,'
  print,'                radtemp= , rphot= , /noprot, ioneq_file=, /not_interactive ]'
  print," i.e. > plot_populations,'c_5',1.e+6,4"
  return
endif
;
xtitle='Electron Density (cm!u-3!n)'
ytitle='Fractional Population'

IF n_elements(ioneq_file) EQ 0 THEN BEGIN
  ioneq_file=!ioneq_file
ENDIF ELSE BEGIN
  chck=file_search(ioneq_file,count=count)
  IF count EQ 0 THEN BEGIN
    print,'% PLOT_POPULATIONS: the specified ioneq file was not found. Returning...'
    return
  ENDIF 
ENDELSE 

;
ion2spectroscopic, gname ,species
title='CHIANTI:  '+species
;
;   get data
;
ion2filename,gname,filename
;
ecname=filename+'.elvlc'
wname=filename+'.wgfa'

outname=filename+'.pop'
pname=filename+'.psplups'
IF keyword_set(noprot) THEN BEGIN
  pstr=-1
ENDIF ELSE BEGIN
  read_splups, pname, pstr, pref, /prot
ENDELSE

IF n_elements(rphot) EQ 0 THEN dilute=0d0 ELSE dilute=r2w(rphot)
IF n_elements(radtemp) EQ 0 THEN radt=6d3 ELSE radt=double(radtemp)

;
read_elvlc,ecname,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,ecmthryd,eref
mult=2.*jj+1.
g=where(ecm eq 0.)
if(max(g) gt 0) then begin
   ecm(g)=ecmth(g)
endif
;
read_wgfa2,wname,wlvl1,wlvl2,wvl1,gf1,a_value1,wref
  ntrans=n_elements(wlvl1)
  nlvls=max([wlvl1,wlvl2])
  wvl=fltarr(nlvls,nlvls)
  gf=fltarr(nlvls,nlvls)
  a_value=fltarr(nlvls,nlvls)
  for itrans=0,ntrans-1 do begin
    wl1=wlvl1(itrans)
    wl2=wlvl2(itrans)
    wvl(wl1-1,wl2-1)=wvl1(itrans)
    gf(wl1-1,wl2-1)=gf1(itrans)
    a_value(wl1-1,wl2-1)=a_value(wl1-1,wl2-1)+a_value1(itrans)  ; specifically for 2 trans types
  endfor
  
; read v.8 format files:
  uname=filename+'.scups'
  read_scups,uname,sp
 splref=sp.info.comments
  splstr=sp   ;.data     ; fixed, PRY, 21-Jul-2016

IF n_elements(nlevels) EQ 0 THEN $
nlevels = n_elements(ecm) ELSE nlevels = nlevels < n_elements(ecm)

nxne=16
xne=fltarr(nxne)
pops=fltarr(nlevels,nxne)
;
;
print,' '
for i=0,nlevels-1 do print,i+1,strpad(term(i),30,/after),format='(i5,4x,a30)'
xne = 10.^FINDGEN(nxne)

pe_ratio=proton_dens(alog10(t))

read_ionrec,filename,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status

IF status lt 0 THEN ionrec = {rec:0., ci:0., temp:0., lev_up_rec:0., lev_up_ci:0., status:-1, ioneq:0., temp_ioneq:0.}
IF status gt 0 THEN BEGIN
   convertname,gname,iz,ion
   ;; dir=concat_dir(!xuvtop,'ioneq')
   ;;    ioneq_name=ch_get_file(path=dir,filter='*.ioneq', $
   ;;                           title='Select Ionization Equilibrium File')
   ;; ff = findfile(ioneq_name)
   ;; IF  ff(0)  NE ''  THEN BEGIN 
     read_ioneq,ioneq_file,ioneq_logt,ioneq,ioneq_ref 
     t_ioneq=ioneq_logt                                        ; Defines temperature array for ion abundances in ionrec
   ;; ENDIF ELSE BEGIN
   ;;    err_msg = '%CH_SYNTHETIC: Error,  no ioneq file found ! -- EXIT'
   ;;    print,err_msg
   ;;    return
   ;; END
   IF ion gt 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-2:ion))
   IF ion eq 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-1:ion))  ; No coll.ionization to neutral ions!
   ionrec = {rec:rec_rate, ci:ci_rate, temp:temp_ionrec, lev_up_rec:luprec, $
   lev_up_ci:lupci, status:status, ioneq:ioneq_ionrec, temp_ioneq:t_ioneq}
ENDIF

;  calculate level populations
;------------------------------

;create input structure for pop_solver:

input = {gname:gname, jj:jj, ecm:ecm,ecmth:ecmth, $
 wvl:wvl, a_value:a_value, splstr:splstr, ionrec:ionrec}

input =CREATE_STRUCT(input, 'pe_ratio',pe_ratio,  'prot_struc', pstr, $
                  'dilute', dilute, 'radtemp', radt ) 

pop_solver, input, t,xne,pop

pop = REFORM(pop)
pops = TRANSPOSE(pop[*,0:nlevels-1])

;
print,' '
print,'Density       Populations'
print,' '
FOR i=0,nxne-1 DO print,format='(e10.3,"   ",500e10.3)',xne[i],pops[*,i]

;
print,' '
print,'Density       Log10 Populations'
print,' '
FOR i=0,nxne-1 DO print,format='(e10.3,"   ",500f10.4)',xne[i],alog10(pops[*,i])

print, ''
print,' Levels: '
for i=0,nlevels-1 do print,i+1,strpad(term(i),30,/after),format='(i5,4x,a30)'

window, 0
circle_sym,/fill

x_min=min(xne)
x_max=max(xne)
y_min=1.e-6
y_max=1.1

pr=''
begin_post:

!p.multi=0

; IF keyword_set(log) THEN $

rand=randomu(seed,nxne)
isrt=sort(rand)                 ; to try to fix problems with labeling curves


plot_oo,xne,pops(0,*),psym=-8,$
  xrange=[x_min,x_max], yrange=[y_min,y_max],$
  xstyle=1,ystyle=1,$
  title=title, $
      xtitle=xtitle,ytitle=ytitle, /nodata

xyouts, 0.1, 0.9 ,/norm, ' Levels: '


;plot_oo,xne,pops(0,*),yr=[1.e-6,1.],title=title, $
;      xtitle=xtitle,ytitle=ytitle
;for i=1,nlevels-1 do oplot,xne,pops(i,*)

FOR  ilist=0, nlevels-1 DO  BEGIN 

oplot,xne, pops(ilist,*) ;, psym=-8

 xyouts, 0.1, 0.85-0.03*ilist,/norm,$
   string(ilist+1)+'  '+strpad(term(ilist),30,/after)

FOR iden=1, nxne-2 DO BEGIN 

;  iden=isrt(ilist<(nxne-1))

   xyouts,xne(iden),pops(ilist,iden),string(ilist+1,'(i3)'),chars=1.2,$
     align=0.5,noclip=0
   xyouts,xne(1),pops(ilist,1),string(ilist+1,'(i3)'),chars=1.2 ,$
     align=0.5,noclip=0
   xyouts,xne(nxne-2),pops(ilist,nxne-2),string(ilist+1,'(i3)'),$
     chars=1.2,align=0.5,noclip=0

ENDFOR

ENDFOR


IF NOT keyword_set(not_interactive) THEN BEGIN 

print2d_plot, pr=pr , x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
  go_to_line=go_to_line, out_name=out_name , /ask_name
if go_to_line eq 'y' then goto,begin_post
END 

if keyword_set(psfile) then begin
   dname = !d.name
   set_plot,'ps'
   device,filename=psfile
   landscape
   !p.font=0
   device,/times,/isolatin1,font_size=14
   thick=3
   plot_oo,xne,pops(0,*),yr=[1.e-6,1.],title=title, $
      xthick=thick,ythick=thick,thick=thick, $
      xtitle=xtitle,ytitle=ytitle
   for i=1,nlevels-1 do oplot,xne,pops(i,*),thick=thick
   device,/close
   !p.font=-1
   set_plot, dname
endif


;
;
if keyword_set(outfile) then begin
   openw,luo,outfile,/get_lun
;
   printf,luo,outname
   printf,luo,' Temperature=',t
   printf,luo, ' '
;
   for i=0,nxne-1 do begin
       printf,luo,xne(i),pops(*,i)
   endfor
;
   printf,luo, ' '
;
   for i=0,nxne-1 do begin
   print,xne(i),alog10(pops(*,i))
      printf,luo,xne(i),alog10(pops(*,i))
   endfor
;
   free_lun,luo
;
endif
;
end
