;编辑软件Notepad++
;使用的是水平分辨率为1度的ERA-Interim的相对涡度数据
;做925，900，875和850hPa的垂直平均
;不做时间上的滤波
;进行标准差为2度的高斯平滑
;保存的数据也是水平分辨率为1度的数据
;遇到问题请联系冯城烽，lzfchfthdl@163.com

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
 
begin 
      
    file_dir            = "/raid63/cffeng/14_vorticity_track_hodges_all_year" 
    file_ori            = file_dir+"/original_data"  
    file_smo            = file_dir+"/pf_smooth_data" 
    file_gau            = file_dir+"/gaussian"

    do i_year = 1979,2016
;    do i_year = 2016,2016

        file_vor_1p000 = file_ori+"/vorticity_19000101-19001231_all.nc" 
        file_vor_1p000 = str_sub_str(file_vor_1p000,"1900",tostring(i_year))
        print(file_vor_1p000)
                                            
        a = addfile(file_vor_1p000,"r")  
        vo_temp = flt2dble(short2flt(a->vo)) 
        vor_wnp_1p000 = dim_avg_n_Wrap(vo_temp(time|:,level|2:5,latitude|:,longitude|:),1) 
        time = vor_wnp_1p000&time
        lat = vor_wnp_1p000&latitude
        lon = vor_wnp_1p000&longitude
        n_time = dimsizes(time)
        n_lat = dimsizes(lat)
        n_lon = dimsizes(lon)
    
        vor_wnp_1p000_npf = vor_wnp_1p000
    
        gaus_1p00000_2p00 = fbindirread(file_gau+"/gaus_1p00000_2p00",0,(/21,21/),"double")
        vor_wnp_1p000_npf_2p00_gaus = vor_wnp_1p000_npf
        temp = new((/n_time,3*n_lat-2,3*n_lon/),"double")
        vor_wnp_1p000_npf_180 = vor_wnp_1p000_npf
        vor_wnp_1p000_npf_180(:,:,  0:179) = vor_wnp_1p000_npf(:,:,180:359)
        vor_wnp_1p000_npf_180(:,:,180:359) = vor_wnp_1p000_npf(:,:,  0:179)
        temp(:,0*n_lat+0:1*n_lat-2,0*n_lon+0:1*n_lon-1) = vor_wnp_1p000_npf_180(:,0:n_lat-2:-1,:)
        temp(:,1*n_lat-1:2*n_lat-2,0*n_lon+0:1*n_lon-1) = vor_wnp_1p000_npf
        temp(:,2*n_lat-1:3*n_lat-3,0*n_lon+0:1*n_lon-1) = vor_wnp_1p000_npf_180(:,1:n_lat-1:-1,:)
        temp(:,0*n_lat+0:1*n_lat-2,1*n_lon+0:2*n_lon-1) = vor_wnp_1p000_npf_180(:,0:n_lat-2:-1,:)
        temp(:,1*n_lat-1:2*n_lat-2,1*n_lon+0:2*n_lon-1) = vor_wnp_1p000_npf
        temp(:,2*n_lat-1:3*n_lat-3,1*n_lon+0:2*n_lon-1) = vor_wnp_1p000_npf_180(:,1:n_lat-1:-1,:)
        temp(:,0*n_lat+0:1*n_lat-2,2*n_lon+0:3*n_lon-1) = vor_wnp_1p000_npf_180(:,0:n_lat-2:-1,:)
        temp(:,1*n_lat-1:2*n_lat-2,2*n_lon+0:3*n_lon-1) = vor_wnp_1p000_npf
        temp(:,2*n_lat-1:3*n_lat-3,2*n_lon+0:3*n_lon-1) = vor_wnp_1p000_npf_180(:,1:n_lat-1:-1,:)
    
        do k = 0,n_time-1
            do j = 1*n_lat-1,2*n_lat-2
                do i = 1*n_lon+0,2*n_lon-1
                    vor_wnp_1p000_npf_2p00_gaus(k,j-n_lat+1,i-n_lon) = sum(temp(k,j-10:j+10,i-10:i+10)*gaus_1p00000_2p00)
                end do
            end do
        end do
    
        ncdf_out_fname = "/vor_wnp_1p000_19000101-19001231.nc"
        ncdf_out_fname = str_sub_str(ncdf_out_fname,"1900",tostring(i_year))
        system("rm -rf "+file_smo+ncdf_out_fname) 
        ncdf_out = addfile(file_smo+ncdf_out_fname,"c")
        ncdf_out->$"vo"$ = vor_wnp_1p000
        print(ncdf_out_fname)
    
        ncdf_out_fname = "/vor_wnp_1p000_npf_2p00_gaus_19000101-19001231.nc" 
        ncdf_out_fname = str_sub_str(ncdf_out_fname,"1900",tostring(i_year))
        system("rm -rf "+file_smo+ncdf_out_fname) 
        ncdf_out = addfile(file_smo+ncdf_out_fname,"c")
        ncdf_out->$"vo"$ = vor_wnp_1p000_npf_2p00_gaus(:,::-1,:)
        print(ncdf_out_fname)
    
        delete([/file_vor_1p000,a,vo_temp,vor_wnp_1p000/])
        delete([/time,lat,lon,n_time,n_lat,n_lon/])
        delete([/vor_wnp_1p000_npf/])
        delete([/gaus_1p00000_2p00,vor_wnp_1p000_npf_2p00_gaus,vor_wnp_1p000_npf_180,temp/])
        delete([/ncdf_out_fname,ncdf_out/])

    end do

end

