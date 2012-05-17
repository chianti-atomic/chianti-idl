;+
; Project     : SOHO - CDS     
;                   
; Name        : DEM_FIT
;               
; Purpose     : Calculates the Differential Emission Measure DEM(T) using
;		a set of values in common with other routines.
;               
; Category    : diagnostic analysis
;
; Explanation : This routine (called by CHIANTI_DEM.PRO)   performs
;		 a series of iterations in order to
;		find the DEM that minimize the chi^2. The values in common 
;		with CHIANTI_DEM.PRO are used.
;		As a least squares fit to a non-linear
;       	function, see pages 237-239, Bevington, Data Reduction and Error
;       	Analysis for the Physical Sciences.
;               
; Use         : IDL> dem_fit,y,chisqr 
;		It has to be noted that a general use of this routine 
;		is limited.
;    
; Inputs      : the values stored in common. 
;                              
; Opt. Inputs : None
;               
; Outputs     : log_dem_mesh,y ,chisqr
;		i.e. the DEM mesh values, the intensity values as resulted
;		from this model DEM, and the chi^2.
; 
;               
; Opt. Outputs: None
;               
; Keywords    : 
;		FLAMBDA: the initial value of the parameter flambda.
;
;		SCALE:   the initial value of the parameter scale, that 
;			controls the steps of the iteration.
;
;		N_ITER:
;		optional.It is the number of iterations of the fitting routine.
;		If not set, a default value of 20 is assumed. 
;		Changing this value alone might not affect the fit, since 
;		also the value of DCHISQ_MIN is checked during the fit.
;
;		DCHISQ_M:
;		optional. If not set, a default value of DCHISQ_MIN=1.e-5 
;		is assumed. For each iteration, the CHISQ and it's variation 
;		are calculated. As long as the iteration achieves an
;		improvement in CHISQ greater than  DCHISQ_MIN , another 
;		iteration will be performed.
;
;		FAILED: If the fit fails,the routine returns and flags FAILED=1
;
;		QUIET:
;		optional. Set to avoid various messages
;
; Calls       : DEM_DERIV, DEM_FUNCTN, DEM_CHISQR
;
; Common      : 
;               obs,  obs_int,obs_sig,n_obs
;               dem,  d_dem_temp,dem_temp,log_dem_temp,log_t_mesh,log_dem_mesh 
;		contr,ch_tot_contr
;               
; Restrictions: Not always the fit is successful.
;               
; Side effects: None known.
;               
;               
; Prev. Hist. :
;       Written by Ken Dere (NRL) as part of the CHIANTI package 
;       in collaboration with Brunella Monsignori Fossi, Enrico Landi
;       (Arcetri Observatory, Florence), Helen Mason and Peter Young
;       (DAMTP, Cambridge Univ.). Incorporated into the CDS software.  
;
; Written     : 
;       Giulio Del Zanna (GDZ), 
;	UCLAN  (University of Central Lancashire, UK)  5 November 1997
;
; HISTORY:
;
;       Ver 1, GDZ  5-Nov-97
;
;       Ver 2,  EL 6-Apr-05
;               Renamed the variable "deriv" to "deriv1" to avoid conflicts
;               with an IDL routine with the same name.
;       
; Version     : 2.0  6 April 2005
;
;-        
;

;; Auxiliary routine
;;------------------------------------------------------
;
function dem_functn,dum
  common obs, obs_int,obs_sig,n_obs
  common dem, d_dem_temp,dem_temp,log_dem_temp,log_t_mesh,log_dem_mesh
  common contr, ch_tot_contr
;
;
dem=spline(log_t_mesh,log_dem_mesh,log_dem_temp)
dem=10.^dem
dlnt=alog(10.^d_dem_temp)
out=fltarr(n_obs)
for iobs=0,n_obs-1 do begin
   out(iobs)=total(ch_tot_contr(*,iobs)*dem*dem_temp)*dlnt
endfor
return,out
end
;

;; Auxiliary routine
;  _______________________________________
;
function dem_deriv,dum
  common obs, obs_int,obs_sig,n_obs
  common dem, d_dem_temp,dem_temp,log_dem_temp,log_t_mesh,log_dem_mesh
  common contr, ch_tot_contr
;
;
nterms=n_elements(log_dem_mesh)
deriv1=fltarr(n_obs,nterms)
;
dem_spl_sav=log_dem_mesh
delta_a=0.01*log_dem_mesh>1.e-6   ; the > part necessary when restricting a ge 0.
yfit=dem_functn(1.)
;
for j=0,nterms-1 do begin
;
log_dem_mesh(j)=dem_spl_sav(j)+delta_a(j)
deriv1(0,j)=(dem_functn(1.)-yfit)/delta_a(j)
log_dem_mesh(j)=dem_spl_sav(j)
endfor
;
return,deriv1
;
end
;
; ______________________________________
;
function dem_chisqr,dum
  common obs, obs_int,obs_sig,n_obs
  common dem, d_dem_temp,dem_temp,log_dem_temp,log_t_mesh,log_dem_mesh
  common contr, ch_tot_contr
;
dy=obs_int-dem_functn(1.)
chi=total(dy^2/obs_sig^2)
;
return,chi
end

;  _______________________________________
;
pro dem_fit,y,chisqr, flambda=flambda,scale=scale,$
	niter=niter,dchisq_min=dchisq_min,$
	failed=failed,quiet=quiet


  common obs, obs_int,obs_sig,n_obs
  common dem, d_dem_temp,dem_temp,log_dem_temp,log_t_mesh,log_dem_mesh
  common contr, ch_tot_contr


on_error,2


dem_temp_min=min(log_dem_temp)
dem_temp_max=max(log_dem_temp)

if n_elements(flambda) eq 0 then flambda=10.
if n_elements(scale) eq 0 then scale=1.
if n_elements(niter) eq 0 then niter=20
if n_elements(dchisq_min) eq 0 then dchisq_min=1.e-5
if n_elements(failed) eq 0 then failed=0


;
;---------------------------------------------
;        plot the initial  DEM
;----------------------------------------------
;
device,window_state=ws
if(ws(1) ne 1) then  window,1,ysize=425,xsize=525 else wset,1


;THIS is required: if you don't have DEM values at the extremes of the 
;range, the spline produces values that then create problems...

        log_dem=spline(log_t_mesh,log_dem_mesh,log_dem_temp)

	missing_high=where(log_dem_temp  gt max(log_t_mesh), nc )

;replace with the last value....
;--------------------------------
	if nc gt 0 then log_dem(missing_high)=$
		 log_dem_mesh (n_elements(log_t_mesh)-1)

	missing_low=where(log_dem_temp lt min(log_t_mesh),nc)

;replace with the first value....
;--------------------------------
	if nc gt 0 then log_dem(missing_low)= log_dem_mesh(0)

        dem=10.^log_dem

;
plot_oo,dem_temp,dem,xr=[10.^dem_temp_min,10.^dem_temp_max],yr=[1.e15,1.e30],$
			;/ynozero,yr=[min(dem)/100.,max(dem)*100.],$
	ytitle=' DEM [ cm!S!E-5 !NK!S!E-1!N ]',$
	xtitle=' T [ !eo!nK ]',$
	title='CHIANTI SPECTRAL CODE AND INVERSION TECHNIQUE',chars=1.2



;---------------------------------------------
;        fit DEM
;----------------------------------------------
;
chisqr=dem_chisqr(1.)
chisq1=chisqr
dchisq=1.
;
;print,' dem chisqr = ',chisqr
;
y=dem_functn(1.)
;
;for i=0,n_obs-1 do print,i,obs_id(i),obs_wvl(i),obs_int(i),y(i)

weight=1./ ( obs_sig  > 1)

;
nterms=n_elements(log_t_mesh)   ; # of parameters (mesh points)
nfree=n_obs-nterms     ; Degrees of freedom

 if nfree le 0 then  begin 
	failed=1
	message, 'not enough data points.',/cont
	goto,fin
 endif


iter=0


;Start the iterations:
;---------------------
while (iter le niter) and (dchisq ge dchisq_min) do begin

;
  alpha= fltarr(nterms,nterms)
  beeta=  fltarr(nterms)
  deriv1= fltarr(nterms)
  b=     fltarr(nterms)
  sigmaa=fltarr(nterms)
  array= fltarr(nterms,nterms)

b=log_dem_mesh
deriv1=dem_deriv(1.)
;
dy=obs_int-dem_functn(1.)

;
  for i=0,n_obs-1 do begin

    for j=0,nterms-1 do begin
      beeta(j)=beeta(j)+weight(i)*dy(i)*deriv1(i,j)
      for k=0,j do begin
        alpha(j,k)=alpha(j,k)+weight(i)*deriv1(i,j)*deriv1(i,k)
      endfor
    endfor
  endfor
;
alpha=transpose(alpha)
;
;
  chisq1=dem_chisqr(nfree)
;

start_this_iter:
;****************
;
for j=0,nterms-1 do begin
    for k=0,nterms-1 do begin
      array(j,k)=alpha(j,k)/sqrt(alpha(j,j)*alpha(k,k))
    endfor
    array(j,j)=1.+flambda
endfor
;
  array=invert(array)
;
;
  for j=0,nterms-1 do begin
    log_dem_mesh(j)=b(j)
    for k=0,nterms-1 do begin
      log_dem_mesh(j)=log_dem_mesh(j)+beeta(k)*array(j,k)/sqrt(alpha(j,j)*alpha(k,k))*scale
    endfor
  endfor
;

;calculate the new chisqr:
;-------------------------
  chisqr=dem_chisqr(nfree)

        IF chisqr NE chisqr THEN BEGIN

	   IF NOT keyword_set(quiet) THEN print,"NaN encountered"

 	   IF flambda GT 1.0e10 THEN BEGIN
	    failed=1
 	    message, 'Failed to converge,silly lambda value, done '+trim(iter)+' iterations',/cont
	    goto,fin
	   ENDIF 

	   IF NOT keyword_set(quiet) THEN $
		 print,'increasing flambda to ',10.*flambda

  	    flambda=10.*flambda

	    IF scale gt 0.1 THEN BEGIN
	      IF NOT keyword_set(quiet) THEN $
			print,'reducing the scale to ',scale/2. 
	      scale=scale/2. 
	    ENDIF
	
		goto, start_this_iter

        ENDIF

  if( chisqr lt chisq1 ) then begin 

     for j=0,nterms-1 do begin
  		sigmaa(j)=sqrt(array(j,j)/alpha(j,j))
     endfor

  	if flambda eq 0 then begin 
	failed=1
	IF NOT keyword_set(quiet) THEN $
	message,' flambda=0,done '+trim(iter)+' iterations',/cont
  	goto, fin  
  	endif 

     flambda=flambda/10.

  endif else begin 


;the fit has got worse... check few things
;---------------------------------------

	if flambda eq 0 then begin 
	failed=1
        IF NOT keyword_set(quiet) THEN $
	message,' flambda=0, done '+trim(iter)+' iterations',/cont  
	goto, fin
	endif 

       IF flambda GT 1.0e10 THEN BEGIN
        print, "Too many repeats, silly lambda value."
  	failed=1
	message, 'Failed to converge, done '+trim(iter)+' iterations',/cont
	goto,fin
       ENDIF ELSE BEGIN 

	print,'increasing flambda to ',10.*flambda
	  flambda=10.*flambda

       	IF scale gt 0.1 THEN BEGIN
	 IF NOT keyword_set(quiet) THEN print,'reducing the scale to ',scale/2. 
	  scale=scale/2. 
       	ENDIF
	ENDELSE

  goto, start_this_iter

  endelse


 end_this_iter: 

   chisqr=dem_chisqr(nfree)
   dchisq=(chisq1-chisqr)
   chisq1=chisqr

   iter=iter+1

   IF NOT keyword_set(quiet) THEN $
	print,' iter= ',iter,' chisqr =',chisqr,'  delta chisqr =',dchisq
;
  log_dem=spline(log_t_mesh,log_dem_mesh,log_dem_temp)
   dem=10.^log_dem 
;
   oplot,dem_temp,dem,linestyle=2

;stop

;
endwhile

fin:

;   print out chisqr
;-------------------
chisqr=dem_chisqr(1.)
print,' dem chisqr = ',chisqr
;

y=dem_functn(1.)


; overplot the final DEM:
;------------------------
  log_dem=spline(log_t_mesh,log_dem_mesh,log_dem_temp)
   dem=10.^log_dem 
oplot,dem_temp,dem,thick=3


end
;



