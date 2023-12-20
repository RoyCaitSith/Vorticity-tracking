;该程序用来对扰动进行筛选

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
    
    file_dir      = "/raid63/cffeng/14_vorticity_track_hodges_all_year"
    file_ori      = file_dir+"/original_data"
    file_tra      = file_dir+"/best_track"
    file_gen      = file_dir+"/JMA_genesis"

    file_lnd_sea  = file_ori+"/land_sea_mask.nc"
    a = addfile(file_lnd_sea,"r")
    lsm = short2flt(a->lsm(0,:,:))
    delete(a)
;    printVarSummary(lsm)
;    print(miao)

    do i_year = 1979,2016
;    do i_year = 2001,2016

        print(i_year)

        file_sea_tem = file_ori+"/sst_19000101-19001231.nc"
        file_tra_cho = file_tra+"/ff_trs_pos_low_1900.new"
        file_sav_cho = file_tra+"/ff_trs_pos_low_select_1900_no_lat_limit.new"
;        file_sav_cho = file_tra+"/ff_trs_pos_low_select_1900.new"
        file_tra_gen = file_gen+"/1900.txt"
        file_sea_tem = str_sub_str(file_sea_tem,"1900",tostring(i_year))
        file_tra_cho = str_sub_str(file_tra_cho,"1900",tostring(i_year))
        file_sav_cho = str_sub_str(file_sav_cho,"1900",tostring(i_year))
        file_tra_gen = str_sub_str(file_tra_gen,"1900",tostring(i_year))

        system("rm -rf "+file_sav_cho)

        a = addfile(file_sea_tem,"r")
        sst = short2flt(a->sst)
        time = a->time
        latitude = a->latitude
        longitude = a->longitude
        n_time = dimsizes(time)
        n_latitude = dimsizes(latitude)
        n_longitude = dimsizes(longitude)
        delete(a)

        nrow = numAsciiRow(file_tra_cho)
        besttrack = asciiread(file_tra_cho,nrow,"string")

        nrow_gen = numAsciiRow(file_tra_gen)
        genesis = asciiread(file_tra_gen,nrow_gen,"string")

        n_track = 0
        count = 4
        do while (count.lt.nrow)
            id = stringtointeger(str_get_field(besttrack(count-1),2," "))
            lines = stringtointeger(str_get_field(besttrack(count),2," "))
            count_next = count+lines+2
            tim = new(lines,"double")
            lon = new(lines,"double")
            lat = new(lines,"double")
            vor = new(lines,"double")
            do j = count+1,count_next-2
                tim(j-count-1) = stringtodouble(str_get_field(besttrack(j),1," "))
                lon(j-count-1) = stringtodouble(str_get_field(besttrack(j),2," "))
                lat(j-count-1) = stringtodouble(str_get_field(besttrack(j),3," "))
                vor(j-count-1) = stringtodouble(str_get_field(besttrack(j),4," "))
            end do
			;是否在JMA_genesis的记录之中
            if_genesis_in = False
            do m = 0,nrow_gen-1
                if((id.eq.stringtointeger(str_get_cols(genesis(m),29,33))))
                    if_genesis_in = True
                end if 
            end do
			;是否在研究区域之中
            if_domain_in = False
            do j = count+1,count_next-2
                if((lon(j-count-1).ge.110).and.(lon(j-count-1).le.180).and.\
                   (lat(j-count-1).ge.  0).and.(lat(j-count-1).le. 40))
                    if_domain_in = True
                end if
            end do
			;t0时刻涡度是否超过2.15度
            if_lat_high = False
            max_ind = ind(vor.eq.max(vor))
            if(lat(max_ind).ge.2.15) 
                if_lat_high = True
            end if
			;t0时刻涡度是否超过1.25*10^-5 s^-1
            if_vor_high = False
;            max_ind = ind(vor.eq.max(vor))
;            abs_vor = vor(max_ind)+2.0*7.2722*sin(3.1415926*lat(max_ind)/180.0)
;            if(abs_vor.ge.2.6)
            if(max(vor).ge.1.25)
                if_vor_high = True
            end if
			;t0时刻前48小时平均SST至少为26.5摄氏度
            if_sst_high = False
            max_vor = 0.0
            do j = count+1,count_next-2
                if(vor(j-count-1).ge.max_vor) 
                    find_tim_max_vor = tim(j-count-1)
                    find_tim_max_vor_j = j-count-1
                end if
            end do  
            do k = 0,n_time-1
                if(cd_calendar(time(k),3).eq.find_tim_max_vor)
                    mean_sst_tim = k 
                end if
            end do
            j_str = max((/doubletointeger(find_tim_max_vor_j-8),0/))
            j_end = find_tim_max_vor_j
            mean_sst_temp = new(j_end-j_str+1,"double")
            do j = j_str,j_end
                mean_sst_temp(j-j_str) = linint2_points(longitude,latitude(::-1),sst(j-j_end+mean_sst_tim,::-1,:),True,lon(j),lat(j),0)
            end do
            n_mean_sst_temp = num(.not.ismissing(mean_sst_temp))
            if(n_mean_sst_temp.ne.0)
                mean_sst = avg(mean_sst_temp)-273.15
                if(mean_sst.ge.26.5)
;                if(mean_sst.ge.26)
                   if_sst_high = True
                end if
            else
                mean_sst = -999.0
            end if
            mean_sst@_FillValue = -999.0
;            if_start_on_sea = False
;            if((lsm(90-round(lat(0),3),mod(round(lon(0),3),360)).eq.0).or.\
;               (lsm(90-round(lat(1),3),mod(round(lon(1),3),360)).eq.0).or.\
;               (lsm(90-round(lat(2),3),mod(round(lon(2),3),360)).eq.0).or.\
;               (lsm(90-round(lat(3),3),mod(round(lon(3),3),360)).eq.0).or.\
;               (max(lat(0:3)).lt.30))
;               if_start_on_sea = True
;            end if
;            flag = (if_domain_in.and.if_lat_high.and.if_vor_high.and.if_sst_high).or.if_genesis_in
            flag = (if_domain_in.and.if_vor_high.and.if_sst_high).or.if_genesis_in
            if(flag)
                hlist1 = [/"TRACK_ID     ",id,"    START_TIME",tim(0)/]  
                hlist2 = [/"POINT_NUM   ",lines/]
                alist = [/tim,lon,lat,vor/]
                write_table(file_sav_cho,"a",hlist1,"%s%2i%s%10.0f")
                write_table(file_sav_cho,"a",hlist2,"%s%3i")
                write_table(file_sav_cho,"a",alist,"%10.0f%16.6f%16.6f%16.6f")
            end if
            count = count_next
            delete([/tim,lon,lat,vor/])
            delete(max_ind)
;            delete(abs_vor)
            delete(mean_sst_temp)
        end do

        delete(sst)
        delete([/time,latitude,longitude/])
        delete([/n_time,n_latitude,n_longitude/])
        delete([/nrow,besttrack/])
        delete([/nrow_gen,genesis/])

    end do

end
