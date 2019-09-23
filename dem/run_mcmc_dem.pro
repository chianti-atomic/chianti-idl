;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. See www.chiantidatabase.org
;                   
; Name        : run_mcmc_dem
;     		          
; Purpose     : This is a wrapper routine to call the mcmc_dem
;               function.
;               
; Explanation : This routine is called by CHIANTI_DEM.
;               The mcmc_dem  and suite of routines need to be
;               available in the IDL path. To retrieve them, go to
;               http://hea-www.harvard.edu/PINTofALE/
;
;               The  MCMC_DEM() code runs a
;               Markov-Chain Monte-Carlo algorithm on a set of line
;               fluxes and returns an estimate of the DEM that
;               generates the observed fluxes.
;
; 
;              v.1, 25 Jun 2014,  Giulio Del Zanna (GDZ) 
;              
;
; VERSION  1, 25 Jun 2014
;-
function run_mcmc_dem, wvl,int,emiss,$
                          Z=Z,$
                          logt=logt_grid, demrng=demrng,$
                          fsigma=fsigma,ulim=ulim,$
                          diffem=diffem, nsim=nsim,nbatch=nbatch,nburn=nburn,smoot=smoot,$
                          storpar=storpar,storidx=storidx,$
        simprb=simprb,simdem=simdem,demerr=demerr,simflx=simflx,$
        simprd=simprd,nosrch=nosrch,softlim=softlim,$
        sampenv=sampenv,smooscl=smooscl, _extra=_extra, $      
                       savfil=savfil


  recompile, 'mcmc_dem', status=status
   if status eq 0 then begin 
      print, 'Routine MCMC_DEM not found --- EXIT ! '
      return,-1
      
    endif else begin 
     
      out=mcmc_dem(wvl,int, emiss,$
                          Z=Z,$
                          logt=logt_grid, demrng=demrng,$
                          fsigma=fsigma,ulim=ulim,$
                          diffem=diffem, nsim=nsim,nbatch=nbatch,nburn=nburn,smoot=smoot,$
                          storpar=storpar,storidx=storidx,$
        simprb=simprb,simdem=simdem,demerr=demerr,simflx=simflx,$
        simprd=simprd,nosrch=nosrch,softlim=softlim,$
        sampenv=sampenv,smooscl=smooscl, _extra=_extra, $      
                   savfil=savfil
      

     
     return,out
     
  endelse
    
     
end
