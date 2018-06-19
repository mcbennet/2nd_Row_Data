***,H2
memory,1000,m

! 1 x2,y2,z2
! 2 x    (l=1)
! 3 y    (l=-1)
! 4 xy
! 5 z    (l=0)
! 6 xz
! 7 yz
! 8 xyz

print,orbitals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LONE ATOM GROUNDSTATES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! Si !!!!!!!!!!!!!!!!!

basis={
ecp,Si,10,2,0;
6;
1,1.22418085,4.0;
2,2.05337336,40.72596063;
2,1.7141285,-48.11509746;
2,2.41395005,-37.28006653;
2,2.32084434,38.28006653;
3,1.35299631,4.89672339;
6;
2,1.13070385,-7.68509694;
2,1.16859753,13.98411213;
2,2.36994226,-116498.38332824;
2,2.167341,-9121.48068622;
2,2.44879942,31941.11999828;
2,2.32322104,93679.74429067;
6;
2,1.86811003,41248.64599856;
2,2.10179754,-41245.51022334;
2,1.33467919,-60.37864776;
2,2.29835912,4180.55486914;
2,1.93345601,-142125.41164262;
2,1.99192523,138006.23630568;

include,si_aug-cc-pCV5Z_2s2p2d.basis
}

geometry={Si}

{RHF
  wf,4,7,2
  occ,1,0,1,0,1,0,0,0
  open,1.3,1.5
}
si_hf_gs=energy

PUT,XML,wf_tn-df.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
si_cc_gs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

basis={
ecp,O,2,2,0;
6;
1,8.86932353,6.0;
2,6.05326172,-28.04199287;
2,5.51480979,11.15704031;
2,10.77878397,180.8243251;
2,10.23693413,-179.8243251;
3,7.90462675,53.21594115;
6;
2,7.28893859,-9212.20980516;
2,6.0597119,9226.8656795;
2,10.83143357,58203.26727502;
2,5.75281092,-5120.48607364;
2,10.51155711,-93321.50266843;
2,9.72227746,40239.72318888;
6;
2,7.43321349,10001.55649464;
2,5.85047476,-10012.86801601;
2,5.79011164,8554.95973537;
2,8.08750969,-20342.33136146;
2,8.4322992,11739.44079236;
2,4.71055456,48.9283704;

include,o_aug-cc-pCV5Z.basis
}

geometry={O}

{RHF
  wf,6,6,2
  occ,1,1,1,0,1,0,0,0
  open,1.2,1.5
}
o_hf_gs=energy

{uccsd(t);
   maxit,500;
   core
}
o_cc_gs=energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,Si,10,2,0;
6;
1,1.22418085,4.0;
2,2.05337336,40.72596063;
2,1.7141285,-48.11509746;
2,2.41395005,-37.28006653;
2,2.32084434,38.28006653;
3,1.35299631,4.89672339;
6;
2,1.13070385,-7.68509694;
2,1.16859753,13.98411213;
2,2.36994226,-116498.38332824;
2,2.167341,-9121.48068622;
2,2.44879942,31941.11999828;
2,2.32322104,93679.74429067;
6;
2,1.86811003,41248.64599856;
2,2.10179754,-41245.51022334;
2,1.33467919,-60.37864776;
2,2.29835912,4180.55486914;
2,1.93345601,-142125.41164262;
2,1.99192523,138006.23630568;

include,si_aug-cc-pCV5Z_2s2p2d.basis
}

geometry={Si}

!           +3 +3 +3 +2 +2 +2 +2 +2 +2 +2 +2 +1 +1 +1 +1 +1 +0 +0 +0 +0 +0 +0 -1           
nelectrons=[ 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5]
symm      =[ 1, 3, 4, 1, 3, 1, 7, 1, 4, 2, 8, 3, 7, 1, 4, 8, 1, 1, 8, 5, 2, 7, 8]
ss        =[ 1, 1, 1, 0, 2, 0, 2, 0, 2, 2, 2, 1, 3, 1, 1, 3, 0, 0, 4, 2, 2, 2, 3]
cceval    =[ 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1]

!   +3+3+3+2+2+2+2+2+2+2+2+1+1+1+1+1+0+0+0+0+0+0-1
Ag =[3,2,2,3,3,2,2,2,3,2,2,3,3,3,3,2,3,3,3,3,3,3,3]
B3u=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,1,2]
B2u=[1,2,1,1,2,2,2,1,1,2,1,2,2,2,1,2,2,1,2,2,2,1,2]
B1g=[0,0,1,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,1,1,0]
B1u=[1,1,1,1,1,1,2,2,1,1,2,1,2,1,1,2,1,2,2,2,1,1,2]
B2g=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0]
B3g=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
Au =[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

sing_occ=[1.1,,,,,,,            \ ! +3
,1.3,,,,,,,                     \ ! +3
,1.4,,,,,,,                     \ ! +3
,,,,,,,,                        \ ! +2
,1.1,1.3,,,,,,                  \ ! +2
,,,,,,,,                        \ ! +2
,1.3,1.5,,,,,,                  \ ! +2
,,,,,,,,                        \ ! +2
,1.1,1.4,,,,,,                  \ ! +2
,1.3,1.4,,,,,,                  \ ! +2
,1.5,1.4,,,,,,                  \ ! +2
,1.3,,,,,,,                     \ ! +1
,1.1,1.3,1.5,,,,,               \ ! +1
,1.1,,,,,,,                     \ ! +1
,1.4,,,,,,,                     \ ! +1
,1.3,1.5,1.2,,,,,               \ ! +1
,,,,,,,,                        \ ! +0
,,,,,,,,                        \ ! +0
,1.1,1.3,1.5,1.2,,,,            \ ! +0
,1.1,1.5,,,,,,                  \ ! +0
,1.3,1.4,,,,,,                  \ ! +0
,1.6,1.4,,,,,,                  \ ! +0
,1.2,1.3,1.5,,,,,]                ! -1

DO s=1,#symm

   ! HARTREE-FOCK
   {RHF;maxit,300;
     start,atden
     wf,nelectrons(s),symm(s),ss(s) 
     occ,Ag(s)-2,B3u(s)-1,B2u(s)-1,B1g(s),B1u(s)-1,B2g(s),B3g(s),Au(s)    
     open,sing_occ(8*s-8+1),sing_occ(8*s-8+2),sing_occ(8*s-8+3),sing_occ(8*s-8+4),sing_occ(8*s-8+5),sing_occ(8*s-8+6),sing_occ(8*s-8+7),sing_occ(8*s-8+8);
   }
   hf_energy(s)=energy
   hf_gap(s)=hf_energy(s)-si_hf_gs
   hf_val(s)=hf_energy(s)-si_hf_cs

   IF (cceval(s).EQ.1) THEN

      {uccsd(t)
       core
      }
      cc_energy(s)=energy
      cc_gap(s)=cc_energy(s)-si_cc_gs
      cc_val(s)=cc_energy(s)-si_cc_cs

   END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Si2     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO s=0,6

 y = 1.75+0.125*s
 
 geometry={
   2
  Si Dimer
  Si          0.0000000000        0.0000000000       -y*0.5
  Si          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,8,4,2   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 si2_hf_bd=energy-2.0*si_hf_gs

 {uccsd(t);
    maxit,500;
    core
 }
 si2_cc_bd=energy-2.0*si_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     SiO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,Si,10,2,0;
6;
1,1.22418085,4.0;
2,2.05337336,40.72596063;
2,1.7141285,-48.11509746;
2,2.41395005,-37.28006653;
2,2.32084434,38.28006653;
3,1.35299631,4.89672339;
6;
2,1.13070385,-7.68509694;
2,1.16859753,13.98411213;
2,2.36994226,-116498.38332824;
2,2.167341,-9121.48068622;
2,2.44879942,31941.11999828;
2,2.32322104,93679.74429067;
6;
2,1.86811003,41248.64599856;
2,2.10179754,-41245.51022334;
2,1.33467919,-60.37864776;
2,2.29835912,4180.55486914;
2,1.93345601,-142125.41164262;
2,1.99192523,138006.23630568;

include,si_aug-cc-pCV5Z_2s2p2d.basis

ecp,O,2,2,0;
6;
1,8.86932353,6.0;
2,6.05326172,-28.04199287;
2,5.51480979,11.15704031;
2,10.77878397,180.8243251;
2,10.23693413,-179.8243251;
3,7.90462675,53.21594115;
6;
2,7.28893859,-9212.20980516;
2,6.0597119,9226.8656795;
2,10.83143357,58203.26727502;
2,5.75281092,-5120.48607364;
2,10.51155711,-93321.50266843;
2,9.72227746,40239.72318888;
6;
2,7.43321349,10001.55649464;
2,5.85047476,-10012.86801601;
2,5.79011164,8554.95973537;
2,8.08750969,-20342.33136146;
2,8.4322992,11739.44079236;
2,4.71055456,48.9283704;

include,o_aug-cc-pCV5Z.basis
}

DO s=0,6

 z = 1.15+0.1*s
 
 geometry={
   2
  Si Oxide
  Si          0.0000000000        0.0000000000       -z*0.5
   O          0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,10,1,0   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 sio_hf_bd=energy-si_hf_gs-o_hf_gs

 pop;

 {uccsd(t);
    maxit,500;
    core
 }
 sio_cc_bd=energy-si_cc_gs-o_cc_gs

END DO
