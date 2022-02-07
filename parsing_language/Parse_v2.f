c----------------------------------------------------------------
c
c Converts protein primary sequence to 3-letter code where the
c labels F, D, and P predict protein regions that are folded (F),
c intrinsically disordered (D), or phase-separating intrinsically
c disordered (P), respectively.
c
c STW 12/22/2020
c
c----------------------------------------------------------------

      program Parse_v2
      implicit none

      integer length,num,i,j,npep,jj,
     &        window_size,jjj,middle_position,
     &        num_ala,num_cys,num_asp,num_glu,num_phe,
     &        num_gly,num_his,num_ile,num_lys,num_leu,
     &        num_met,num_asn,num_pro,num_gln,num_arg,
     &        num_ser,num_thr,num_val,num_trp,num_tyr,
     &        num_ala_w,num_cys_w,num_asp_w,num_glu_w,num_phe_w,
     &        num_gly_w,num_his_w,num_ile_w,num_lys_w,num_leu_w,
     &        num_met_w,num_asn_w,num_pro_w,num_gln_w,num_arg_w,
     &        num_ser_w,num_thr_w,num_val_w,num_trp_w,num_tyr_w,
     &        count_regions,count_p,count_w,region,
     &        p_region_start,p_region_end,p_start(10000),
     &        p_end(10000),count_regions_p,count_regions_d,
     &        count_regions_f,d_start(10000),d_end(10000),
     &        f_start(10000),f_end(10000),count_domains,low,
     &        domain(10000)
      character code_inp*11000,code(11000)*3,classification(11000)*3
      real pppiia_h,pppiip_h,pppiig_h,pppiic_h,pppiid_h,pppiie_h,
     &     pppiif_h,pppiih_h,pppiii_h,pppiik_h,pppiil_h,pppiim_h,
     &     pppiin_h,pppiiq_h,pppiir_h,pppiis_h,pppiit_h,pppiiv_h,
     &     pppiiw_h,pppiiy_h,slope,intercept,v_line,net_charge_w,
     &     sum_ppii,v_exponent,fppii,v_flory,rh_w,helix,hydr,
     &     helix_a,helix_c,helix_d,helix_e,helix_f,helix_g,helix_h,
     &     helix_i,helix_k,helix_l,helix_m,helix_n,helix_p,helix_q,
     &     helix_r,helix_s,helix_t,helix_v,helix_w,helix_y,
     &     hydr_a,hydr_c,hydr_d,hydr_e,hydr_f,hydr_g,hydr_h,
     &     hydr_i,hydr_k,hydr_l,hydr_m,hydr_n,hydr_p,hydr_q,
     &     hydr_r,hydr_s,hydr_t,hydr_v,hydr_w,hydr_y,m,b,x,y,
     &     percent_p,percent_cutoff,whole_seq_helix,whole_seq_nu,
     &     rh,net_charge,whole_seq_hydr,dist_norm(10000),
     &     Null_dist,PS_dist,p_dist_sum


c PPII bias measured in peptides by Hilser group; used to calculate Rh and then nu.
c Prot Sci 2013, vol 22, pgs 405-417, in a table in supplementary information

      pppiia_h=0.37
      pppiic_h=0.25
      pppiid_h=0.30
      pppiie_h=0.42
      pppiif_h=0.17
      pppiig_h=0.13
      pppiih_h=0.20
      pppiii_h=0.39
      pppiik_h=0.56
      pppiil_h=0.24
      pppiim_h=0.36
      pppiin_h=0.27
      pppiip_h=1.00
      pppiiq_h=0.53
      pppiir_h=0.38
      pppiis_h=0.24
      pppiit_h=0.32
      pppiiv_h=0.39
      pppiiw_h=0.25
      pppiiy_h=0.25

c normalized frequency frequency of alpha-helix (Nagano, J. Mol. Biol. 75, 401-420 (1973)).

      helix_a=1.29
      helix_c=0.94
      helix_d=1.00
      helix_e=1.54
      helix_f=1.23
      helix_g=0.72
      helix_h=1.29
      helix_i=0.94
      helix_k=1.23
      helix_l=1.23
      helix_m=1.23
      helix_n=0.77
      helix_p=0.70
      helix_q=1.10
      helix_r=0.83
      helix_s=0.78
      helix_t=0.87
      helix_v=0.97
      helix_w=1.06
      helix_y=0.63

c Structure-based interactivity scale used to calculate the hydrophobicity profile of a protein
c from its primary sequence. This intertactivity scale was determined from the residue contact matrix
c of single-domain globular proteins (Bastolla et al., Proteins 58, 22-30 (2005)).

      hydr_a=0.0728
      hydr_c=0.3557
      hydr_d=-0.0552
      hydr_e=-0.0295
      hydr_f=0.4201
      hydr_g=-0.0589
      hydr_h=0.0874
      hydr_i=0.3805
      hydr_k=-0.0053
      hydr_l=0.3819
      hydr_m=0.1613
      hydr_n=-0.0390
      hydr_p=-0.0492
      hydr_q=0.0126
      hydr_r=0.0394
      hydr_s=-0.0282
      hydr_t=0.0239
      hydr_v=0.2947
      hydr_w=0.4114
      hydr_y=0.3113

c Read input protein sequence.

      call get_command_argument(1,code_inp)
      if (len_trim(code_inp) == 0) then
      write(*,*)'no input argument, exiting program'
      stop
      endif

      length=len(code_inp)

c  convert any lower case letter to upper case.
    
      do i=1,length
         num=ichar(code_inp(i:i))
         if (num.ge.97.and.num.le.122) code_inp(i:i) = char(num-32)
      enddo

c  Determine sequence length.

      j=0
      do i=1,length
      if (code_inp(i:i).eq.' ') goto 1
         j=j+1
         code(j)=code_inp(i:i)
1     continue
      enddo 
      npep=j

c Define window size.

      window_size=25

c if protein sequence is less than the window size, stop

      if (npep.lt.window_size) then
      write(*,*)' '
      write(*,*)'input sequence is too short'
      write(*,*)'minimum sequence length is ',window_size
      write(*,*)' '
      stop
      endif

c if protein sequence is too long, stop

      if (npep.gt.10000) then
      write(*,*)' '
      write(*,*)'input sequence is too long'
      write(*,*)'maximum sequence length is 10000'
      write(*,*)' '
      stop
      endif

c calculate nu, helix propensity, and hydrophobicity for the whole sequence

      num_ala=0
      num_cys=0
      num_asp=0
      num_glu=0
      num_phe=0
      num_gly=0
      num_his=0
      num_ile=0
      num_lys=0
      num_leu=0
      num_met=0
      num_asn=0
      num_pro=0
      num_gln=0
      num_arg=0
      num_ser=0
      num_thr=0
      num_val=0
      num_trp=0
      num_tyr=0

      DO J=1,NPEP
         IF (CODE(J).EQ.'A') THEN
            num_ala=num_ala+1
         endif
         IF (CODE(J).EQ.'C') THEN
            num_cys=num_cys+1
         endif
         IF (CODE(J).EQ.'D') THEN
            num_asp=num_asp+1
         endif
         IF (CODE(J).EQ.'E') THEN
            num_glu=num_glu+1
         endif
         IF (CODE(J).EQ.'F') THEN
            num_phe=num_phe+1
         endif
         IF (CODE(J).EQ.'G') THEN
            num_gly=num_gly+1
         endif
         IF (CODE(J).EQ.'H') THEN
            num_his=num_his+1
         endif
         IF (CODE(J).EQ.'I') THEN
            num_ile=num_ile+1
         endif
         IF (CODE(J).EQ.'K') THEN
            num_lys=num_lys+1
         endif
         IF (CODE(J).EQ.'L') THEN
            num_leu=num_leu+1
         endif
         IF (CODE(J).EQ.'M') THEN
            num_met=num_met+1
         endif
         IF (CODE(J).EQ.'N') THEN
            num_asn=num_asn+1
         endif
         IF (CODE(J).EQ.'P') THEN
            num_pro=num_pro+1
         endif
         IF (CODE(J).EQ.'Q') THEN
            num_gln=num_gln+1
         endif
         IF (CODE(J).EQ.'R') THEN
            num_arg=num_arg+1
         endif
         IF (CODE(J).EQ.'S') THEN
            num_ser=num_ser+1
         endif
         IF (CODE(J).EQ.'T') THEN
            num_thr=num_thr+1
         endif
         IF (CODE(J).EQ.'V') THEN
            num_val=num_val+1
         endif
         IF (CODE(J).EQ.'W') THEN
            num_trp=num_trp+1
         endif
         IF (CODE(J).EQ.'Y') THEN
            num_tyr=num_tyr+1
         endif
      enddo

      net_charge=abs(num_asp+num_glu
     &                  -num_lys-num_arg)

      sum_ppii=0.0

      sum_ppii=num_ala*pppiia_h+num_cys*pppiic_h
     &        +num_asp*pppiid_h+num_glu*pppiie_h
     &        +num_phe*pppiif_h+num_gly*pppiig_h
     &        +num_his*pppiih_h+num_ile*pppiii_h
     &        +num_lys*pppiik_h+num_leu*pppiil_h
     &        +num_met*pppiim_h+num_asn*pppiin_h
     &        +num_pro*pppiip_h+num_gln*pppiiq_h
     &        +num_arg*pppiir_h+num_ser*pppiis_h
     &        +num_thr*pppiit_h+num_val*pppiiv_h
     &        +num_trp*pppiiw_h+num_tyr*pppiiy_h

      fppii=sum_ppii/real(npep)
      v_exponent=0.503-0.11*log(1.0-fppii)

      rh=2.16*(real(npep)**(v_exponent))
     &    +0.26*real(net_charge)
     &    -0.29*(real(npep)**(0.5))
      v_flory=log(rh/2.16)/log(real(npep))

      helix=0.0

      helix=num_ala*helix_a+num_cys*helix_c
     &        +num_asp*helix_d+num_glu*helix_e
     &        +num_phe*helix_f+num_gly*helix_g
     &        +num_his*helix_h+num_ile*helix_i
     &        +num_lys*helix_k+num_leu*helix_l
     &        +num_met*helix_m+num_asn*helix_n
     &        +num_pro*helix_p+num_gln*helix_q
     &        +num_arg*helix_r+num_ser*helix_s
     &        +num_thr*helix_t+num_val*helix_v
     &        +num_trp*helix_w+num_tyr*helix_y
      helix=helix/real(npep)

      hydr=0.0

      hydr=num_ala*hydr_a+num_cys*hydr_c
     &        +num_asp*hydr_d+num_glu*hydr_e
     &        +num_phe*hydr_f+num_gly*hydr_g
     &        +num_his*hydr_h+num_ile*hydr_i
     &        +num_lys*hydr_k+num_leu*hydr_l
     &        +num_met*hydr_m+num_asn*hydr_n
     &        +num_pro*hydr_p+num_gln*hydr_q
     &        +num_arg*hydr_r+num_ser*hydr_s
     &        +num_thr*hydr_t+num_val*hydr_v
     &        +num_trp*hydr_w+num_tyr*hydr_y
      hydr=hydr/real(npep)

      whole_seq_hydr=hydr
      whole_seq_helix=helix
      whole_seq_nu=v_flory

      p_dist_sum=0.0

c calculate nu, helix propensity, and hydrophobicity for each window

      DO J=1,NPEP
         if (j.le.(NPEP-window_size+1)) then
         middle_position=j+window_size/2

      num_ala_w=0
      num_cys_w=0
      num_asp_w=0
      num_glu_w=0
      num_phe_w=0
      num_gly_w=0
      num_his_w=0
      num_ile_w=0
      num_lys_w=0
      num_leu_w=0
      num_met_w=0
      num_asn_w=0
      num_pro_w=0
      num_gln_w=0
      num_arg_w=0
      num_ser_w=0
      num_thr_w=0
      num_val_w=0
      num_trp_w=0
      num_tyr_w=0

      DO JJJ=J,J+window_size-1
         IF (CODE(JJJ).EQ.'A') THEN
            num_ala_w=num_ala_w+1
         endif
         IF (CODE(JJJ).EQ.'C') THEN
            num_cys_w=num_cys_w+1
         endif
         IF (CODE(JJJ).EQ.'D') THEN
            num_asp_w=num_asp_w+1
         endif
         IF (CODE(JJJ).EQ.'E') THEN
            num_glu_w=num_glu_w+1
         endif
         IF (CODE(JJJ).EQ.'F') THEN
            num_phe_w=num_phe_w+1
         endif
         IF (CODE(JJJ).EQ.'G') THEN
            num_gly_w=num_gly_w+1
         endif
         IF (CODE(JJJ).EQ.'H') THEN
            num_his_w=num_his_w+1
         endif
         IF (CODE(JJJ).EQ.'I') THEN
            num_ile_w=num_ile_w+1
         endif
         IF (CODE(JJJ).EQ.'K') THEN
            num_lys_w=num_lys_w+1
         endif
         IF (CODE(JJJ).EQ.'L') THEN
            num_leu_w=num_leu_w+1
         endif
         IF (CODE(JJJ).EQ.'M') THEN
            num_met_w=num_met_w+1
         endif
         IF (CODE(JJJ).EQ.'N') THEN
            num_asn_w=num_asn_w+1
         endif
         IF (CODE(JJJ).EQ.'P') THEN
            num_pro_w=num_pro_w+1
         endif
         IF (CODE(JJJ).EQ.'Q') THEN
            num_gln_w=num_gln_w+1
         endif
         IF (CODE(JJJ).EQ.'R') THEN
            num_arg_w=num_arg_w+1
         endif
         IF (CODE(JJJ).EQ.'S') THEN
            num_ser_w=num_ser_w+1
         endif
         IF (CODE(JJJ).EQ.'T') THEN
            num_thr_w=num_thr_w+1
         endif
         IF (CODE(JJJ).EQ.'V') THEN
            num_val_w=num_val_w+1
         endif
         IF (CODE(JJJ).EQ.'W') THEN
            num_trp_w=num_trp_w+1
         endif
         IF (CODE(JJJ).EQ.'Y') THEN
            num_tyr_w=num_tyr_w+1
         endif
      enddo

      net_charge_w=abs(num_asp_w+num_glu_w
     &                  -num_lys_w-num_arg_w)

      sum_ppii=0.0

      sum_ppii=num_ala_w*pppiia_h+num_cys_w*pppiic_h
     &        +num_asp_w*pppiid_h+num_glu_w*pppiie_h
     &        +num_phe_w*pppiif_h+num_gly_w*pppiig_h
     &        +num_his_w*pppiih_h+num_ile_w*pppiii_h
     &        +num_lys_w*pppiik_h+num_leu_w*pppiil_h
     &        +num_met_w*pppiim_h+num_asn_w*pppiin_h
     &        +num_pro_w*pppiip_h+num_gln_w*pppiiq_h
     &        +num_arg_w*pppiir_h+num_ser_w*pppiis_h
     &        +num_thr_w*pppiit_h+num_val_w*pppiiv_h
     &        +num_trp_w*pppiiw_h+num_tyr_w*pppiiy_h

      fppii=sum_ppii/real(window_size)
      v_exponent=0.503-0.11*log(1.0-fppii)

c multiplier of 4 on the sequence, to put nu into the length-independent range

      rh_w=2.16*(real(4*window_size)**(v_exponent))
     &    +0.26*real(4*net_charge_w)
     &    -0.29*(real(4*window_size)**(0.5))
      v_flory=log(rh_w/2.16)/log(real(4*window_size))

      helix=0.0

      helix=num_ala_w*helix_a+num_cys_w*helix_c
     &        +num_asp_w*helix_d+num_glu_w*helix_e
     &        +num_phe_w*helix_f+num_gly_w*helix_g
     &        +num_his_w*helix_h+num_ile_w*helix_i
     &        +num_lys_w*helix_k+num_leu_w*helix_l
     &        +num_met_w*helix_m+num_asn_w*helix_n
     &        +num_pro_w*helix_p+num_gln_w*helix_q
     &        +num_arg_w*helix_r+num_ser_w*helix_s
     &        +num_thr_w*helix_t+num_val_w*helix_v
     &        +num_trp_w*helix_w+num_tyr_w*helix_y
      helix=helix/real(window_size)

      hydr=0.0

      hydr=num_ala_w*hydr_a+num_cys_w*hydr_c
     &        +num_asp_w*hydr_d+num_glu_w*hydr_e
     &        +num_phe_w*hydr_f+num_gly_w*hydr_g
     &        +num_his_w*hydr_h+num_ile_w*hydr_i
     &        +num_lys_w*hydr_k+num_leu_w*hydr_l
     &        +num_met_w*hydr_m+num_asn_w*hydr_n
     &        +num_pro_w*hydr_p+num_gln_w*hydr_q
     &        +num_arg_w*hydr_r+num_ser_w*hydr_s
     &        +num_thr_w*hydr_t+num_val_w*hydr_v
     &        +num_trp_w*hydr_w+num_tyr_w*hydr_y
      hydr=hydr/real(window_size)

c Sector boundaries were defined by the mean and standard deviation in nu, helix
c propensity, and hydrophobicity calculated for each of the testing, null,
c and folded sets. For the folded set, mean hydr is 0.1163289 ± 0.01673017
c
c Windows with hydr value less than 0.08286856 (mean - 2*sd) are classified as disordered (D or P).
c Windows with hydr value greater than or equal to 0.08286856 are classified as folded (F).
c
c dist_norm is the distance from the sector boundary, normalized by the distance to the mean.

      if(hydr.ge.(0.08286856)) then
         classification(middle_position)='F'
         dist_norm(middle_position)=(hydr-0.08286856)/(0.01673017*2.0)
         goto 10
      endif

c For the PS set, mean helix propensity is 0.9596848 ± 0.06516711
c For the PS set, mean nu_model is 0.5416 ± 0.01997657
c For the Null set, mean helix propensity is 1.027303 ± 0.0674386
c For the Null set, mean nu_model is 0.5583985 ± 0.02217093 
c
c Boundary between P and D sectors was defined by the line
c
c y=(-0.317841)*x + 0.865043
c
c determined from the PS and Null sets means and standard deviations.
c
c In a plot of helix propensity (x) versus nu_model (y), a window localized
c to the right of this boundary line is in the D sector. A window localized
c to the left of the boundary is in the P sector.

c distance from the boundary line of the null set means
      m=-1.0/(-0.317841)
      b=0.5583985-m*1.027303
      x=(b-0.865043)/(-0.317841-m)
      y=m*x+b
      Null_dist=sqrt((1.027303-x)*(1.027303-x)+
     &              (0.5583985-y)*(0.5583985-y))

c distance from the boundary line of the testing set means
      m=-1.0/(-0.317841)
      b=0.5416-m*0.9596848
      x=(b-0.865043)/(-0.317841-m)
      y=m*x+b
      PS_dist=sqrt((0.9596848-x)*(0.9596848-x)+
     &              (0.5416-y)*(0.5416-y))

c Here, x and y define the point on the boundary line that makes a
c perpendicular when paired with the point defined by the window 
c of v_flory (i.e., nu_model) and helix (i.e., helix propensity).
c m and b, below, define the equation of this perpendicular line. 
c Because the line is perpendicular to the boundary, its slope will be
c the negative reciprocal of the boundary slope.
      m=-1.0/(-0.317841)
      b=v_flory-m*helix
c x and y define the intersect point of the boundary (y=(-0.317841)*x + 0.865043)
c and perpendicular (y=mx+b). Two equations, two unknowns, so, easily
c solved.
      x=(b-0.865043)/(-0.317841-m)
      y=m*x+b

      if(((v_flory-0.865043)/(-0.317841)).le.(helix)) then
         classification(middle_position)='D'
         dist_norm(middle_position)=
     &      sqrt((helix-x)*(helix-x)+(v_flory-y)*(v_flory-y))/
     &      Null_dist
         goto 10
      else
         classification(middle_position)='P'
         dist_norm(middle_position)=
     &      sqrt((helix-x)*(helix-x)+(v_flory-y)*(v_flory-y))/
     &      PS_dist
            p_dist_sum=p_dist_sum+dist_norm(middle_position)
         goto 10
      endif

10    continue
      endif
      enddo

      do j=1,window_size/2
         classification(j)=classification((window_size/2)+1)
         dist_norm(j)=dist_norm((window_size/2)+1)
      enddo
      do j=npep-(window_size/2)+1,npep
         classification(j)=classification(npep-(window_size/2))
         dist_norm(j)=dist_norm(npep-(window_size/2))
      enddo


      write(*,*)' '
      write(*,'(11000a1)')(classification(j),j=1,npep)

      write(*,*)' '
      write(*,'("Sequence length ",i6)')npep
      write(*,'("Whole sequence nu_model ",f6.3)')whole_seq_nu
      write(*,'("Whole sequence helix prop ",f6.3)')whole_seq_helix
      write(*,'("Whole sequence hydrophobicity  ",f6.3)')whole_seq_hydr
      write(*,'("p_dist_sum  ",f10.3)')p_dist_sum

      open (7,file='residue_level_dist_norm.csv')
      do j=1,npep
      write(7,'(i6,", ",a1,", ",a1,",",f6.3)')j,code(j),
     & classification(j),dist_norm(j)
      enddo
      close(7)

      percent_cutoff=0.90
c find PS (blue) regions 20 residues or longer and labeled P at the percent_cutoff or higher
      i=1
      count_regions=0
20    continue
      count_p=0
      count_w=0
      do j=i,i+19
      region=0
      count_w=count_w+1
      if(classification(j).eq.'P') count_p=count_p+1
      enddo
30    continue
      percent_p=real(count_p)/real(count_w)
      if(percent_p.ge.percent_cutoff) then
         region=1
         p_region_start=i
         p_region_end=j
         if(j.lt.npep) then
         j=j+1
         count_w=count_w+1
         if(classification(j).eq.'P') then
            count_p=count_p+1
         endif
         goto 30
         endif
      endif
      if(region.eq.1) then
         i=j
         count_regions=count_regions+1
         p_start(count_regions)=p_region_start
         p_end(count_regions)=p_region_end
      endif
      if(region.eq.0) i=i+1
      if((i+19).le.npep) goto 20
      count_regions_p=count_regions

c find ID (red) regions 20 residues or longer and labeled D at the percent_cutoff or higher
      i=1
      count_regions=0
40    continue
      count_p=0
      count_w=0
      do j=i,i+19
      region=0
      count_w=count_w+1
      if(classification(j).eq.'D') count_p=count_p+1
      enddo
50    continue
      percent_p=real(count_p)/real(count_w)
      if(percent_p.ge.percent_cutoff) then
         region=1
         p_region_start=i
         p_region_end=j
         if(j.lt.npep) then
         j=j+1
         count_w=count_w+1
         if(classification(j).eq.'D') then
            count_p=count_p+1
         endif
         goto 50
         endif
      endif
      if(region.eq.1) then
         i=j
         count_regions=count_regions+1
         d_start(count_regions)=p_region_start
         d_end(count_regions)=p_region_end
      endif
      if(region.eq.0) i=i+1
      if((i+19).le.npep) goto 40
      count_regions_d=count_regions

c find F (black) regions 20 residues or longer and labeled F at the percent_cutoff or higher
      i=1
      count_regions=0
60    continue
      count_p=0
      count_w=0
      do j=i,i+19
      region=0
      count_w=count_w+1
      if(classification(j).eq.'F') count_p=count_p+1
      enddo
70    continue
      percent_p=real(count_p)/real(count_w)
      if(percent_p.ge.percent_cutoff) then
         region=1
         p_region_start=i
         p_region_end=j
         if(j.lt.npep) then
         j=j+1
         count_w=count_w+1
         if(classification(j).eq.'F') then
            count_p=count_p+1
         endif
         goto 70
         endif
      endif
      if(region.eq.1) then
         i=j
         count_regions=count_regions+1
         f_start(count_regions)=p_region_start
         f_end(count_regions)=p_region_end
      endif
      if(region.eq.0) i=i+1
      if((i+19).le.npep) goto 60
      count_regions_f=count_regions

c find first domain (earliest in sequence) and identify its first residue
      count_domains=count_regions_p+count_regions_d+count_regions_f
      low=1000000
      j=1
      if(j.le.count_domains) then
      do i=1,count_regions_p
         if(p_start(i).lt.low) low=p_start(i)
      enddo
      do i=1,count_regions_d
         if(d_start(i).lt.low) low=d_start(i)
      enddo
      do i=1,count_regions_f
         if(f_start(i).lt.low) low=f_start(i)
      enddo
      domain(j)=low
      endif

c find first residue of each additional domain in successive order
80    j=j+1
      low=1000000
      if(j.le.count_domains) then
      do i=1,count_regions_p
         if(p_start(i).lt.low.and.p_start(i).gt.domain(j-1)) 
     &      low=p_start(i)
      enddo
      do i=1,count_regions_d
         if(d_start(i).lt.low.and.d_start(i).gt.domain(j-1)) 
     &      low=d_start(i)
      enddo
      do i=1,count_regions_f
         if(f_start(i).lt.low.and.f_start(i).gt.domain(j-1)) 
     &      low=f_start(i)
      enddo
      domain(j)=low
      goto 80
      endif

c correct for domain-domain overlap when present
      if(count_domains.gt.1) then
c find overlapping domains and split the overlap, which is percent_cutoff*20/2
      do i=1,count_domains-1
      do j=1,count_regions_p
         if(p_start(j).eq.domain(i)) then
         if(p_end(j).ge.domain(i+1)) then
            p_end(j)=domain(i+1)+int((1.0-percent_cutoff)*20.0/2.0)
            do jj=1,count_regions_d
            if(d_start(jj).eq.domain(i+1)) then
               d_start(jj)=p_end(j)+1
               domain(i+1)=d_start(jj)
            endif
            enddo
            do jj=1,count_regions_f
            if(f_start(jj).eq.domain(i+1)) then
               f_start(jj)=p_end(j)+1
               domain(i+1)=f_start(jj)
            endif
            enddo
         endif
         endif
      enddo
      do j=1,count_regions_d
         if(d_start(j).eq.domain(i)) then
         if(d_end(j).ge.domain(i+1)) then
            d_end(j)=domain(i+1)+int((1.0-percent_cutoff)*20.0/2.0)
            do jj=1,count_regions_p
            if(p_start(jj).eq.domain(i+1)) then
               p_start(jj)=d_end(j)+1
               domain(i+1)=p_start(jj)
            endif
            enddo
            do jj=1,count_regions_f
            if(f_start(jj).eq.domain(i+1)) then
               f_start(jj)=d_end(j)+1
               domain(i+1)=f_start(jj)
            endif
            enddo
         endif
         endif
      enddo
      do j=1,count_regions_f
         if(f_start(j).eq.domain(i)) then
         if(f_end(j).ge.domain(i+1)) then
            f_end(j)=domain(i+1)+int((1.0-percent_cutoff)*20.0/2.0)
            do jj=1,count_regions_d
            if(d_start(jj).eq.domain(i+1)) then
               d_start(jj)=f_end(j)+1
               domain(i+1)=d_start(jj)
            endif
            enddo
            do jj=1,count_regions_p
            if(p_start(jj).eq.domain(i+1)) then
               p_start(jj)=f_end(j)+1
               domain(i+1)=p_start(jj)
            endif
            enddo
         endif
         endif
      enddo
      enddo
      endif

      write(*,*)' '
      write(*,'("Number of PS (blue) regions = ",i6)')count_regions_p
      write(*,'("region  first_residue  last_residue")')
      do i=1,count_regions_p
      write(*,'(i4,7x,i6,7x,i6)') i,p_start(i),p_end(i)
      enddo
      write(*,*)' '
      write(*,'("Number of ID (red) regions = ",i7)')count_regions_d
      write(*,'("region  first_residue  last_residue")')
      do i=1,count_regions_d
      write(*,'(i4,7x,i6,7x,i6)') i,d_start(i),d_end(i)
      enddo
      write(*,*)' '
      write(*,'("Number of F (black) regions = ",i5)')count_regions_f
      write(*,'("region  first_residue  last_residue")')
      do i=1,count_regions_f
      write(*,'(i4,7x,i6,7x,i6)') i,f_start(i),f_end(i)
      enddo

      write(*,*)' '
      do i=1,count_domains
      do j=1,count_regions_p
         if(p_start(j).eq.domain(i)) then
         write(*,'("region ",i4,", PS (blue)",",  first residue ",i5,
     &    ",  last residue ",i5,",  length ",i6)')i,p_start(j),
     &    p_end(j),(p_end(j)-p_start(j)+1)
         endif
      enddo
      do j=1,count_regions_d
         if(d_start(j).eq.domain(i)) then
         write(*,'("region ",i4,", ID  (red)",",  first residue ",i5,
     &    ",  last residue ",i5,",  length ",i6)')i,d_start(j),
     &    d_end(j),(d_end(j)-d_start(j)+1)
         endif
      enddo
      do j=1,count_regions_f
         if(f_start(j).eq.domain(i)) then
         write(*,'("region ",i4,", F (black)",",  first residue ",i5,
     &    ",  last residue ",i5,",  length ",i6)')i,f_start(j),
     &    f_end(j),(f_end(j)-f_start(j)+1)
         endif
      enddo
      enddo

      end