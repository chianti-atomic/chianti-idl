
;
; wmb_h5tb_append_records
; 
; Purpose: Adds records to the end of the table. 
;
; Description: wmb_h5tb_append_records adds records to the end of the table 
;              named dset_name attached to the object specified by the 
;              identifier loc_id. The dataset is extended to hold the 
;              new records. 
;              
; Parameters:
; 
; loc_id
;     IN: Identifier of the file or group in which the table is located. 
; dset_name
;     IN: The name of the dataset. 
; nrecords
;     IN: The number of records to insert. 
; databuffer
;     IN: Buffer with data. 
;

pro wmb_h5tb_append_records, loc_id, $
                             dset_name, $
                             nrecords, $
                             databuffer

    compile_opt idl2, strictarrsubs

    ; get the original number of records and fields

    wmb_h5tb_get_table_info, loc_id, dset_name, nfields, nrecords_orig

    ; open the dataset

    did = h5d_open(loc_id, dset_name)
    
    ; append the records
    
    wmb_h5tb_common_append_records, did, $
                                    nrecords, $
                                    nrecords_orig, $
                                    databuffer
    
    ; close

    h5d_close, did
    
end
