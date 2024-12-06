FUNCTION mpdemfunct,p,pm=pm,n_line=n_line,i_obs=i_obs,i_err=i_err, $
        i_mod=i_mod,spl_logt=spl_logt,n_spl=n_spl,t=t,weights=weights
        

;==========================================================================
;+
;
; arguments
; p             dem values at spline knots
; pm[i,*]       emis*T*d(ln t) array for the ith line
; i_obs         observed intensities (scaled)
; i_err         errors in observed intensity
; i_mod         modeled intensities (scaled)
; spl_t         log T values of spline knots
; n_spl         number of spline knots
; t             points to evaluate DEM on
; weights       relative weighting of each line, default = 1.0
;
; HISTORY:
; progver = 'v2007-May-19' ;--- (M.Weber) Just tweaked a bit to fit into the
;                             organization.
; GDZ - removed abundance factor.
;
;-
;==========================================================================

   dem = 10^(spline(spl_logt[0:n_spl-1],p,t))

   i_mod = (dem##pm)
   chisq = float((i_mod-i_obs)*weights/i_err)
   idx = where(weights eq 0.0,count)
   if(count gt 0) then chisq[idx] = 0.0
   return,chisq


END ;======================================================================


pro mpfit_dem, obs_int, responses, logt_responses, obs_err=obs_err,max_logt=max_logt, $
               spl_logt = spl_logt , spl_logdem=spl_logdem,$
               min_logt=min_logt, dt=dt, out_logt=out_logt, out_logdem=out_logdem, $
                maxiter=maxiter,$
               min_limits=min_limits,  max_limits=max_limits, verbose=verbose,solv_factor=solv_factor,$
               MIN_LOGDEM=MIN_LOGDEM, error=error
;+
; Name        : MPFIT_DEM
;
; Purpose     : Calculates the Differential Emission Measure DEM(T)
;              using MPFIT. 
; 
;         INPUT:
;
;              OBS_INT.  input data (e.g. AIA count rates)  must be
;              fltarr(ni)
;             
;              RESPONSES.  responses is an array fltarr(nt,ni)
;              LOGT_RESPONSES.   log(T[K]) values of the responses fltarr(nt)
;
;         OUTPUT:
;              OUT_LOGT, OUT_LOGDEM: output log T and log DEM values.
;              SPL_LOGDEM: spline log DEM values
;
;         OPTIONAL INPUT: 
;              
;              OBS_ERR: uncertainty in the observed intensities. If
;              not defined, assume 20%.
;
;              WEIGHTS  relative weighting of each line, default= 1.0
;
;              SPL_LOGT: log T values for the spline. They must be at
;              most ni-1 values. If these are not defined, then
;              min_logt, max_logt must be defined.
;
;              SPL_LOGDEM: log DEM initial values for the spline.
; 
;              MIN_LOGT:   minimum temperature (log T) for the DEM
;                            calculation
;              MAX_LOGT:   maximum temperature (log T) for the DEM
;                            calculation
;                        Note: if  min_logt, max_logt are defined, the
;                        spline log T values (ni-1) are equally spaced
;                        between the minimum and maximum values.
;
;              DT:        temperature step (log T) for the DEM
;                            calculation (default=0.1)
;
;
;              MIN_LIMITS,  MAX_LIMITS: minimum and maximum limits for
;              the log DEM values at the spline temperatures.
;              Note: if no limits are given, by default the routine
;              applies a minimum value for the log DEM=17.
;
;              SOLV_FACTOR: the responses are multiplied by
;              10^SOLV_FACTOR before MPFIT is run. By default 21.
;
;              MAXITER: maxium number of iterations
;
;              VERBOSE: if set, prints error messages and plots the
;              DEM with the spline points.
;
; Written     :
;              V. 1.0  21 Nov 2016 Giulio Del Zanna (GDZ). This
;              routine is a modification of 
;              XRT_DEM_ITER_NOWIDGET, written by M. Weber,  which
;              is part of the XRT_DEM  package, and available within
;              the SSW XRT path.
;
;
; Modified    : v.2, added ERROR, the spline DEM output and fixed a bug.
;               v.3,  2 Jun 2017, GDZ. A few minor changes.
;               v.4,  8 Aug 2017, GDZ. Reduced number of default nodes
;               to 6. 
;               v.5, 6-Jul-2018 GDZ.  Added some checking and simplifications. 
;               v.6, 6-Nov-2023 GDZ, bug fix: when in error, return the input
;                    spline DEM values.
;
; Version     : 6
;-

error=0

if n_elements(verbose) eq 0 then verbose=0
if n_elements(solv_factor) eq 0 then solv_factor=21
if N_ELEMENTS(dt) eq 0 then dt=0.1
if n_elements(maxiter) eq 0 then maxiter=1000
if n_elements(obs_err) eq 0 then obs_err=0.2 *obs_int   
if n_elements(min_logdem) eq 0 then min_logdem=17.


n_line = N_ELEMENTS(obs_int)

;weights - the initial weights for each line  
if(N_ELEMENTS(weights) eq 0) then weights = replicate(1.0,n_line)


if N_ELEMENTS(spl_logt) gt 0  then begin 
   if n_elements(spl_logt) gt  n_line-1 then begin 
      error=1
      if verbose then print,'Error ! number of spline nodes too big '
      return
   endif else  begin
; redefine the log T into the fine dt grid :      
      n_spl = N_ELEMENTS(spl_logt) 
      max_logt=max(spl_logt) 
      min_logt=min(spl_logt) 
      nt = round((max_logt - min_logt)/dt + 1)
      new_logt = min_logt + findgen(nt)*dt
      
      if verbose then print,'Redefined min and max log T values '+$
                            'on the basis of the spline points'
      
   endelse    
endif else begin 
                                ;  set splines at equally spaced
                                ;  intervals- set 6 nodes as default,
                                ;  but do not exceed the degrees of freedom.
   n_spl =  6 < (n_line-1)         
   spl_logt = FLTARR(n_spl)
   
   if N_ELEMENTS(max_logt) eq 0 or N_ELEMENTS(min_logt) eq 0 then begin 
      error=1
      if verbose then print,'MPFIT_DEM: Error ! no min/max log T values provided. '
      return
   endif 
   
   nt = round((max_logt - min_logt)/dt + 1)
   new_logt = min_logt + findgen(nt)*dt
   
   spl_logt[0:n_spl-1] = min(new_logt) + findgen(n_spl) * $
                         (max(new_logt)-min(new_logt)) / (n_spl-1)
   
   print, 'log T input values set by default =',spl_logt
   
endelse  


if N_ELEMENTS(spl_logdem) gt 0 then begin 
   if N_ELEMENTS(spl_logdem) ne N_ELEMENTS(spl_logt) then begin 
            error=1
            if verbose then print,$
               'MPFIT_DEM: Error ! number of spline DEM values different than input log T '
      return
   endif 
      spl_logdem=float(spl_logdem)
endif else if N_ELEMENTS(spl_logdem) eq 0 then begin 
  
    spl_logdem = replicate(20.0,n_spl)     
   
   if verbose then print, 'log DEM input values at the nodes set by default =',spl_logdem 
   
  
endif 


; now check the min max limits:
if n_elements(max_limits) eq n_spl then $
   for ispl=0,n_spl-1 do spl_logdem[ispl] = spl_logdem[ispl] < (max_limits[ispl])

if n_elements(min_limits) eq n_spl then $
   for ispl=0,n_spl-1 do spl_logdem[ispl] = spl_logdem[ispl] > (min_limits[ispl] )

input_spl_logdem=spl_logdem
 spl_logdem = spl_logdem -solv_factor

; p array  contains  t*emis*dlnt (when multiplied by dem and summed it gives the intensity)
p = FLTARR(n_line,nt)

emis = fltarr(n_line,nt)

for i=0,n_line-1 do begin
; avoid zeroes:    
   emis[i,*] = interpol(10.^solv_factor* responses[*,i], logt_responses, new_logt) > 1e-20    
   p[i,*] = emis[i,*]*10.^new_logt[*]*alog(10^dt)
endfor


i_mod = FLTARR(n_line)          ; modelled intensities 


;   fa - function arguments for use with mpfit
fa = {pm:p, n_line:n_line, i_obs:obs_int, i_err:obs_err,$
      i_mod:i_mod, t:new_logt, spl_logt:spl_logt, n_spl:n_spl, weights:weights}

;; -- Limit the spline knots in log DEM:

parinfo = replicate({limited:[1,0], limits:[min_logdem-solv_factor, 0.]}, n_spl)

if n_elements(min_limits) eq n_spl then parinfo[*].limits[0]=min_limits-solv_factor

if n_elements(max_limits) eq n_spl then begin 
   parinfo[*].limited[1]=1
   parinfo[*].limits[1]=max_limits-solv_factor
endif

if keyword_set(verbose) then quiet=0 else quiet=1


spl_logdem = mpfit('mpdemfunct', spl_logdem, functargs=fa,    $
                   status=status, errmsg=errmsg, bestnorm=bestnorm, $
                   parinfo=parinfo, MAXITER=maxiter, quiet=quiet)

                                ; check status 

if status le 0 then begin 
; FAILURE:
   if keyword_set(verbose) then PRINT, 'MPFIT: Error ',errmsg
   out_logdem=spl_logdem
   out_logdem[*]=0
   
   error=1
;GDZ return the input spline values:
   spl_logdem=input_spl_logdem

   return
   
endif else begin 
   
   ;; Calculate DEM over the T array  
   out_logdem =  spline(spl_logt, spl_logdem, new_logt)
   out_logt=new_logt
   
                                ; calculate the intensities
   ldem = 10.d^out_logdem 
   i_mod = (ldem##p)
   
   di = (i_mod-obs_int)
   chisq = di^2 * weights / (obs_err)^2
   idx = where(weights eq 0.0,count)
   if(count gt 0) then chisq[idx] = 0.0
   
   chisq = total(chisq)

      out_logdem=solv_factor +out_logdem
   spl_logdem=solv_factor +spl_logdem

   if keyword_set(verbose) then begin 
      
      print, '* Total Chi-squared = '+trim(chisq)
      print, 'Obs int=',  arr2str(obs_int,',')
      print, 'Preditced int=', arr2str(i_mod,',')
      print, 'ratio (predicted/obs)=', arr2str(i_mod/obs_int,',')            
      
;;       window,0
      
;;       plot, chars=1.4, out_logt, out_logdem,psym=10 ,/yno ; xr=[4,7],/xst
;;       oplot,  spl_logt, spl_logdem,psym=5,syms=2
      
;; ; over plot the original DEM:
;;       oplot,  spl_logt,input_spl_logdem, psym=6

;; for ii=0, n_elements(spl_logt)-1 do print, spl_logt[ii], spl_logdem[ii]
      
   endif 


endelse                         ; success 


end  
