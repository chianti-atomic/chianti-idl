pro print2d_plot, pr=pr, x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
           go_to_line=go_to_line,out_name=out_name, ask_name=ask_name
;+
; Project     : SOHO - CDS     
;                   
; Name        : print2d_plot
;     		          
; Purpose     : this routine optionally changes the ranges of a 1D plot
;               and print the plot to a landscape postscript file.
;               
; Explanation : 
; 		
;		This routine changes the keywords pr, x_min,x_max,y_min,y_max
;		in order to change the range of a 2-D plot, or to open a
;		SET_PLOT,'ps'   and create a postscript.
;               This is done basically repeating the sequence of commands in the
;               calling routine, within a loop.
;
;
;               The routine that calls print2d_plot should have :
;      
;               window,/free    
;                                      ;;to open a window
;
;               x_min=min(lambda) & x_max=max(lambda)   
;               y_min=min(sp) & y_max=max(sp)   
;               
;                                     ;; to set the initial limits of 
;                                     ;; a one-dimensional plot. 
;               pr=''
;               begin_plot_sp:
;
;
;               plot,lambda,sp,xr=[x_min,x_max],yr=[y_min,y_max],$
;                          xstyle=1,ystyle=1 
;
;                                     ;; to plot (etc....)
;                
;               print2d_plot, pr=pr, x_min=x_min,x_max=x_max,$
;                             y_min=y_min,y_max=y_max,$
;	                      go_to_line=go_to_line,$
;                             out_name=out_name,/ask_name  
;
;               or:
;                 print2d_plot, pr=pr, x_min=x_min,x_max=x_max,$
;                             y_min=y_min,y_max=y_max,$
;	                      go_to_line=go_to_line,$
;                             out_name='spectrum.ps' 
;
;                                    ;;to create the ps file spectrum.ps 
;
;               if go_to_line eq 'y' then goto,begin_plot_sp
;
;                
;
;               If you have an image or a stack of images, the limits cannot be
;               changed and you should omit the min,max limits.
;
;
; Use         : IDL>
;
;
; Examples    : 
;		
;		
;		
;
;    
; Inputs      : 
;		
;               
;               
; Opt. Inputs : 
;
;
;               
; Outputs     : 
;		
;		
;		
;	
; Opt. Outputs:
;		
;	
;		
;
; Keywords    : 
;
;
;
; Calls       : 
;		
;		
;		
;		
;		
; Common      : None.
; 		
;		
;
; Restrictions: None.
;
;               
; Side effects: None known yet.
;               
; Category    : 
;               
; Prev. Hist. :
;
;
;      
; Written     : 
;
;       Giulio Del Zanna (GDZ), 
;	UCLAN  (University of Central Lancashire, UK) 
;
; Modified    : Version 1, GDZ Fri Jan 30 12:21:59 1998
;                v.2 GDZ added  set_plot, 'ps',/copy ,/inter 02-Feb-2000
;		Version 3, 21-Dec-2000, William Thompson
;			Modified for better cross-platform capability.
;
;               Ver 4, 1-May-02, GDZ
;               Modified the setup to go back to display.
;
;
; Version     : Version 4, 1-May-02
;-
;


go_to_line='n'


if ( pr ne strlowcase('y') ) then begin


   IF N_ELEMENTS(x_min) NE 0 AND  N_ELEMENTS(x_max) NE 0 AND $
     N_ELEMENTS(y_min) NE 0 AND  N_ELEMENTS(y_max) NE 0 THEN BEGIN

      answer=''
      read,'do you wish to change the plot limits ? [y/N]  ', answer

   ENDIF ELSE answer='n'

   if(answer eq 'y')then begin

      print,'These are the current limits: '
      print,trim(x_min)+','+trim(x_max)+','+trim(y_min)+','+trim(y_max)

      read,'enter new limits (x_min,x_max,y_min,y_max): ',$
        x_min,x_max,y_min,y_max
      go_to_line='y'
      return
   endif

   read,'do you want a ps file of the last  plot ? [y/N]  ',pr

   if (pr eq strlowcase('y') ) then begin 

      if keyword_set (ask_name) then begin
         out_name=' '
         read,'Type the file name: ',out_name
      endif 

      yes_no, 'Do you want a high resolution ? ', res
      IF res THEN   bits = 8 ELSE bits = 4

      set_plot,'ps', /copy, /interpol

      yesno=''
      read,'Do you want a color output ? [y/N] ',yesno
      if strlowcase(yesno) eq 'y' then color=1 ELSE color=0
      
      yesno=''
      read,'Do you want a portrait ? [y/N] ',yesno
      if strlowcase(yesno) eq 'y' then $

      device, filename=out_name ,color=color, $
        bits=bits, $
        xsize=18,ysize=25,xoffset=1,yoffset=1, /port  ELSE $
        device, filename=out_name ,color=color ,bits=bits, /landscape ;, /encaps

                                ;device,xoffset=1.0,yoffset=10.25,/inches,/landscape
      
      go_to_line='y'


      return

   endif

endif else  begin

   print,' '
   print,'The ps file ',out_name,' has been created'
   print,' '


   device,/close

if not have_windows() then case os_family() of
        'Windows': set_plot,'win'
        'MacOS':   set_plot,'mac'
        else: set_plot,'x'
endcase


   !p.multi=0
;   !p.background=255 
;   !p.color=0

   yesno='' 
   read,'Do you want to replot [y/N] ? ',yesno
   if strlowcase(yesno) eq 'y' then begin 
      go_to_line='y'
      pr=''
   endif else go_to_line='n'

   
endelse

return
end

