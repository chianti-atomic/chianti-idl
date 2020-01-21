;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. See www.chiantidatabase.org
;                   
; Name        : run_data2dem_reg
;     		          
; Purpose     : This is a wrapper routine to call the data2dem_reg
;               function.
;               
; Explanation : This routine is called by CHIANTI_DEM.
;               The data2dem_reg and suite of routines need to be
;               available in the IDL path. To retrieve them, go to
;               http://www.astro.gla.ac.uk/~iain/demreg/ and add them
;               to the IDL path. 
;               The DATA2DEM_REG recovers the DEM using a GSVD approach,
;               see Hannah & Kontar A&A 539, A146 2012. 
; 
;              v.1, 25 Jun 2014,  Giulio Del Zanna (GDZ)
;
;              v.2 17-Sep-2015, Peter Young
;                 Fixed bug identified by Petros Syntelis.
;              
;
; VERSION  1, 17 Sep 2015
;-

function run_data2dem_reg, log_dem_temp, ch_tot_contr, obs_int, obs_sig,$
                        mint=mint,maxt=maxt,nt=nt,$
                       order=order, guess=guess,reg_tweak=reg_tweak,$
                           channels=channels,debug=debug,gloci=gloci, pos=pos, error=error


  recompile, 'data2dem_reg', status=status
   if status eq 0 then begin 
      print, 'Routine DATA2DEM_REG not found --- EXIT ! '
      
      error=1
      
      return,-1
      
   endif else begin 


out=data2dem_reg(log_dem_temp, ch_tot_contr, obs_int, obs_sig,$
                       mint=mint,maxt=maxt,nt=nt,$
                       order=order, guess=guess,reg_tweak=reg_tweak,$
                       channels=channels,debug=debug,gloci=gloci, pos=pos)

error=0

return, out

endelse
   
   
end
