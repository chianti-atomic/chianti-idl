;+
; PROJECT     : CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the
;       University of Cambridge,  Goddard Space Flight Center, and University of Michigan. 
;
;                   
; NAME        : CH_GET_FILE
;     		          
; PURPOSE     : to select a file from either a selected directory or the working
;               directory, having an extension.
;
;               
; EXPLANATION : a file in  either a selected directory or the working
;               directory, having an extension can be selected using a
;               widget. Note that both directory and extension have to be
;               supplied. If no file is found, an empty string is returned.                   ; 		
;
; USE         : IDL> name = ch_get_file( '~/', '.pro',  tit=' Select a procedure ')
;
;
; EXAMPLES    : dir= concat_dir(!xuvtop),'dem') 
;		dem_name=ch_get_file(path=dir,filter='*.dem',title='Select DEM File')
;		
;    
; INPUTS      : directory, extension 
;		
;               
; OPT. INPUTS : 
;
;               
; OUTPUTS     : the file name
;	
; OPT. OUTPUTS:
;		
;
; KEYWORDS    : title
;
;
; CALLS       : findfile, break_file
;		
; COMMON      : co
; 		
; RESTRICTIONS:  both directory and extension have to be
;               supplied. 
;
;               
; SIDE EFFECTS: 
;               
; CATEGORY    : 
;               
; PREV. HIST. : extracted from CDS/CHIANTI routines.
;
;
;      
; WRITTEN     : 
;
;       Giulio Del Zanna (GDZ), 
;	DAMTP  (University of Cambridge, UK) 
;
; MODIFIED    : Version 1, GDZ 10-Oct-2000
;               V.2, GDZ, corrected a typo at the end of the file.
;               V.3, GDZ,  generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;               V. 4, 19-July-2002, GDZ
;
;                Added the option to select files also with the standard IDL
;                dialaog_pickfile, and changed a few things...
;
;               V.5, 2-Aug-02, GDZ 
;                 reduced the size of the widget.
;
;               V.6, 12-Aug-02, GDZ
;                 corrected for a bug in the directory output.
;
;               V.7, 3-Nov-03  GDZ
;                 Fixed a bug when using Windows, the returned path was not
;                 correct. 
;
;              v.8 15 May 2024  Giulio Del Zanna
;                  reinstated the group leader.
;
; VERSION     :  8
;
;
;-
;

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;The following is the event processor for the widget that allows you to select
;-----------------------------------------------------------------------------
;any *+EXTENSION file present in the selected directory or in the working directory
;-----------------------------------------------------------------------------

pro ch_get_ab_event,ev

common co, abund_name,abund_info, co_directory, co_extension, filter_info, f
;

widget_control,ev.id,get_uvalue=uvalue
;name = strmid(tag_names(ev,/str),7,1000)

if n_elements(uvalue) eq 0 then uvalue = ''
if n_elements(value)  eq 0 then  value = ''

CASE uvalue  of

'D': widget_control ,/destroy,ev.top 


'ABUND': begin


      abfile = f(ev.value-1)     

      abund_name = ''

;avoid wrong bits:

      CASE   strmid(abfile, 0, 3) OF 

         '---': BEGIN 
            WIDGET_CONTROL,abund_info, set_val=''
         END 

         'SEL':BEGIN 
            abund_name = dialog_pickfile(filter=co_extension, tit='Select file to READ')
            IF abund_name NE '' THEN BEGIN 

; WIDGET_CONTROL, filter_info, set_value=co_extension

               break_file, abund_name, disk, abdir, abfile, co_extension
              co_directory = disk+abdir
              abfile = abfile+co_extension
               WIDGET_CONTROL,abund_info, set_val=abfile
            END 
         END 

         ELSE:BEGIN 
             abund_name = concat_dir(co_directory[ev.value-1], abfile)
;             print,'Using  file: ',abund_name
             widget_control, abund_info, set_val=abfile
         END 
      ENDCASE 
END 

;FILTER: BEGIN 
;WIDGET_CONTROL, filter_info, get_value=co_extension
;END 
 
ENDCASE 

END  


FUNCTION  ch_get_file, path=directory, filter=extension,  title=title,$
            modal=modal, group_leader=group_leader


common co, abund_name,abund_info, co_directory, co_extension, filter_info, f

IF(XRegistered("ch_get_file")) THEN return, 0

IF n_elements(directory) EQ 0 THEN $
  cd, current=directory

IF N_ELEMENTS(extension) EQ 0 THEN extension = '*'

;co_directory = directory
co_extension = extension

abund_name = ''

;IF n_params() LT 2 THEN message, 'you have to specify directory, extension !'


IF n_elements(title) EQ 0 THEN title = 'Select File'

IF (N_ELEMENTS(GROUP) EQ 0)     THEN GROUP = 0

;get the abundance file
;----------------------
;abund_name=pickfile(path=directory,filter=extension,title=title)
;

get_ab_base=widget_base(title=title, /col)
                                ;, xsize=600,ysize=200,, /column,/frame,x_scroll=sz(0),y_scroll=sz(1))

;  the TOP widget is where all the user/numerical/text input takes place
;
;top = widget_base(get_ab_base,/row)

; the top middle widget contains...
;top_middle= widget_base(get_ab_base, /column)

;  information/selection of  data
;
cs1 = widget_base(get_ab_base,/row)

junk = {CW_PDMENU_S,FLAGS:0,NAME:''}

af = findfile(concat_dir(directory, extension))
f1 = strarr(n_elements(af))
f1_dir = strarr(n_elements(af))

FOR i=0, n_elements(af)-1 DO BEGIN 
   break_file, af[i], disk, dir, ff, ext
   f1(i) = ff+ext
   f1_dir[i] = disk+dir
END

f=['---CHIANTI DIRECTORY: '+extension+' ---',f1]
co_directory = ['--', f1_dir]

cd, current=dir
af = findfile(concat_dir(dir, extension))

IF n_elements(af) GT 0 THEN BEGIN 

   f2 = strarr(n_elements(af))
   f2_dir = strarr(n_elements(af))

   FOR i=0, n_elements(af)-1 DO BEGIN 
      break_file, af[i], disk, dir, ff, ext
      f2(i) = ff+ext
      f2_dir[i] = disk+dir
   END

   f=[f, '---WORKING DIRECTORY: '+extension+' ---', f2]
   co_directory = [co_directory, '--', f2_dir]

ENDIF 


f=[f,'SELECT FILES WITH WIDGET', '' ]
co_directory = [co_directory,'--', '--']



;if f2(0) ne '' then f = ['---CHIANTI DIRECTORY: '+extension+' ---',f1,$
;  '---WORKING DIRECTORY: '+extension+' ---', f2, 'SELECT FILES WITH WIDGET', ''] ELSE $
;   f=['---CHIANTI DIRECTORY: '+extension+' ---',f1, 'SELECT FILES WITH WIDGET', '']


desc = [{CW_PDMENU_S,FLAGS:1,NAME:title}]
for i=0,n_elements(f)-1 do begin
   desc = [desc,{CW_PDMENU_S,FLAGS:0,NAME:f(i)}]
endfor

desc(n_elements(desc)-1).flags=2

menu = cw_pdmenu(cs1,desc,uvalue='ABUND',font=font)
;/return_name,

cs2 = widget_base(cs1,/column)

;filter_info = cw_field(cs2,title='Filter:',value=extension,$
;                      /row,xsize=20, UVALUE='FILTER')

;dummy0 = widget_base(cs2, /row)
;dummy = widget_label(dummy0, value='Filter:')
;dummy = widget_text(dummy0,value=extension, xsize=10)

abund_info = cw_field(cs2,title='Current  file:',value='    ',$
                      /row,xsize=40)

done=widget_button(cs2,UVALUE='D',VALUE='OK', xsize=10)

;  make the whole thing happen
;
widget_control,get_ab_base,/realize

xmanager,'ch_get_ab',get_ab_base  , modal=modal, group_leader=group_leader


;--------END OF THE WIDGET  STUFF-----------------------------------------


return, abund_name

END

