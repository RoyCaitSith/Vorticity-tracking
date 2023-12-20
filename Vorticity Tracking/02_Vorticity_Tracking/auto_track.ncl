;编辑软件Notepad++
;具体说明请参考tracking_scheme使用说明
;对于参数不懂的请先阅读文献
;RUNDATIN请用RUNDATIN.vor_global_hilat_1deg_revised.in
;adapt.dat0请用revised版本，使用之前请重命名

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
 
begin 
      
    file_dir            = "/raid63/cffeng/13_vorticity_track_hodges_gaussian" 
    file_tra            = "/home/cffeng/TRACK-1.4.9"
    file_ori            = file_dir+"/original_data"  
    file_smo            = file_dir+"/pf_smooth_data" 
    file_bes            = file_dir+"/best_track"
    file_vor            = file_smo+"/vor_wnp_1p000.nc" 
    file_fil            = file_smo+"/vor_wnp_1p000_npf_2p00_gaus.nc"
    file_sav            = "TESTVOR_BPF_1"
    file_out            = "-o=/home/cffeng "
    trac_mul            = "-f=vortest "
;    trac_stp            = "-n=1,76,1 "
    trac_stp            = "-n=1,64,24 "
    star_tim            = "2016010100"

;tracking_scheme使用说明中的运行tracking scheme部分的第一至第四部可以忽略
;因为我使用的是高斯滤波，不是谱滤波
;    system("rm -rf "+file_dir+"/"+file_sav)
;    system("rm -rf /home/cffeng/"+file_sav)
;    system("rm -rf "+file_tra+"/outdat/*")
;    system("cp "+file_fil+" "+file_tra+"/indat/vertavg_T63filt.nc")

;下面参数不懂的请参照tracking_scheme使用说明
;    system("master -c="+file_sav+" -d=now -e=track.linux -i=vertavg_T63filt.nc "+trac_mul+\
;           "-j=RUN_AT.in "+trac_stp+file_out+"-s=RUNDATIN.vor_global_hilat_T63")
;    system("rm -rf "+file_tra+"/indat/vertavg_T63filt.nc")
;    system("mv /home/cffeng/"+file_sav+" "+file_dir)
;    system("gunzip -c "+file_dir+"/"+file_sav+"/ff_trs_pos.gz > "+file_bes+"/ff_trs_pos_low")
    system(file_tra+"/utils/bin/count "+file_bes+"/ff_trs_pos_low 0 0 5 4 0 "+star_tim+" 6")

;将路径画出来看
;    system("ncl "+file_dir+"/03_draw_track.ncl")

end

