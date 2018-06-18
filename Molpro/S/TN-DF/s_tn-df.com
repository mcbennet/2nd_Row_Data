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

!!!!!!!!!!!! S !!!!!!!!!!!!!!!!!

basis={
ecp,S,10,2,0;
6;
1,2.51977085,6.0;
2,3.22007986,-84.83332404;
2,4.71655238,70.54487302;
2,4.39998291,3581.56671658;
2,4.41784559,-3580.56671658;
3,2.54586294,15.11862509;
6;
2,4.61819246,-231.72652822;
2,2.30938314,244.26248418;
2,2.6507245,-920.53494189;
2,3.3111907,2410.83323256;
2,3.80226712,-2429.46016726;
2,4.46824294,940.1625125;
6;
2,4.48874898,957.88712772;
2,3.37845034,-950.12559451;
2,3.83307173,6481.0599021;
2,1.92699416,-157.23448173;
2,1.98946862,186.18956071;
2,3.99439281,-6509.01396292;

include,/remote/mbennet/projects/ECP/3p/S/analysis/bases/s_aug-cc-pwCV5Z.basis
}

geometry={S}

{RHF
  wf,6,6,2
  occ,1,1,1,0,1,0,0,0
  open,1.5,1.2
}
s_hf_gs=energy

PUT,XML,wf_tn-df.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
s_cc_gs=energy

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

include,/remote/mbennet/projects/ECP/3p/S/analysis/bases/o_aug-cc-pCV5Z.basis
}

geometry={O}

{RHF
  wf,6,7,2
  occ,1,1,1,0,1,0,0,0
  open,1.3,1.5
}
o_hf_gs=energy

{uccsd(t);
   maxit,500;
   core
}
o_cc_gs=energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,S,10,2,0;
6;
1,2.51977085,6.0;
2,3.22007986,-84.83332404;
2,4.71655238,70.54487302;
2,4.39998291,3581.56671658;
2,4.41784559,-3580.56671658;
3,2.54586294,15.11862509;
6;
2,4.61819246,-231.72652822;
2,2.30938314,244.26248418;
2,2.6507245,-920.53494189;
2,3.3111907,2410.83323256;
2,3.80226712,-2429.46016726;
2,4.46824294,940.1625125;
6;
2,4.48874898,957.88712772;
2,3.37845034,-950.12559451;
2,3.83307173,6481.0599021;
2,1.92699416,-157.23448173;
2,1.98946862,186.18956071;
2,3.99439281,-6509.01396292;

include,/remote/mbennet/projects/ECP/3p/S/analysis/bases/s_aug-cc-pwCV5Z.basis
}

geometry={S}

!           +5 +5 +5 +4 +4 +4 +4 +3 +3 +2 +2 +2 +1 +1 +1 +0 +0 +0 -1           
nelectrons=[11,11,11,12,12,12,12,13,13,14,14,14,15,15,15,16,16,16,17]
symm      =[ 1, 3, 4, 1, 3, 1, 7, 2, 4, 1, 4, 8, 2, 8, 5, 1, 5, 2, 2]
ss        =[ 1, 1, 1, 0, 2, 0, 2, 1, 3, 0, 2, 4, 1, 3, 5, 0, 4, 6, 1]
eval      =[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
cceval    =[ 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

!   +5+5+5+4+4+4+4+3+3+2+2+2+1+1+1+0+0+0-1           
Ag =[3,2,2,3,3,2,2,3,3,3,3,3,3,3,3,3,3,3,3] ! x^2, y^2, z^2
B3u=[1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,1,2,2,2] ! x
B2u=[1,2,1,1,2,2,2,1,2,1,2,2,2,2,2,2,2,2,2] ! y
B1g=[0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0] ! xy
B1u=[1,1,1,1,1,1,2,1,1,1,1,2,1,2,2,2,2,2,2] ! z
B2g=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0] ! xz
B3g=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] ! yz
Au =[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] ! xyz

sing_occ=[1.1,,,,,,,            \ ! +5
,1.3,,,,,,,                     \ ! +5
,1.4,,,,,,,                     \ ! +5
,,,,,,,,                        \ ! +4
,1.1,1.3,,,,,,                  \ ! +4
,,,,,,,,                        \ ! +4
,1.3,1.5,,,,,,                  \ ! +4
,1.2,,,,,,,                     \ ! +3
,1.1,1.2,1.3,,,,,               \ ! +3
,,,,,,,,                        \ ! +2
,1.2,1.3,,,,,,                  \ ! +2
,1.1,1.2,1.3,1.5,,,,            \ ! +2
,1.2,,,,,,,                     \ ! +1
,1.2,1.3,1.5,,,,,               \ ! +1
,1.1,1.4,1.2,1.3,1.5,,,         \ ! +1 *
,,,,,,,,                        \ ! +0
,1.4,1.2,1.3,1.5,,,,            \ ! +0 *
,1.1,1.4,1.6,1.2,1.3,1.5,,      \ ! +0 *
,1.2,,,,,,,]                      ! -1 *

DO s=1,#symm

   ! HARTREE-FOCK
   {RHF;maxit,300;
     start,atden
     wf,nelectrons(s)-10,symm(s),ss(s) 
     occ,Ag(s)-2,B3u(s)-1,B2u(s)-1,B1g(s),B1u(s)-1,B2g(s),B3g(s),Au(s)    
     open,sing_occ(8*s-8+1),sing_occ(8*s-8+2),sing_occ(8*s-8+3),sing_occ(8*s-8+4),sing_occ(8*s-8+5),sing_occ(8*s-8+6),sing_occ(8*s-8+7),sing_occ(8*s-8+8);
   }
   hf_energy(s)=energy
   hf_gap(s)=hf_energy(s)-s_hf_gs
   hf_val(s)=hf_energy(s)-s_hf_cs

   IF (cceval(s).EQ.1) THEN

      {uccsd(t)
       core
      }
      cc_energy(s)=energy
      cc_gap(s)=cc_energy(s)-s_cc_gs
      cc_val(s)=cc_energy(s)-s_cc_cs

   END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     S2     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO s=0,6

 y = 1.5+0.1*s
 
 geometry={
   2
   S Dimer
   S          0.0000000000        0.0000000000       -y*0.5
   S          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,12,4,2   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 s2_hf_bd=energy-2.0*s_hf_gs

 {uccsd(t);
    maxit,500;
    core
 }
 s2_cc_bd=energy-2.0*s_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     SO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,S,10,2,0;
6;
1,2.51977085,6.0;
2,3.22007986,-84.83332404;
2,4.71655238,70.54487302;
2,4.39998291,3581.56671658;
2,4.41784559,-3580.56671658;
3,2.54586294,15.11862509;
6;
2,4.61819246,-231.72652822;
2,2.30938314,244.26248418;
2,2.6507245,-920.53494189;
2,3.3111907,2410.83323256;
2,3.80226712,-2429.46016726;
2,4.46824294,940.1625125;
6;
2,4.48874898,957.88712772;
2,3.37845034,-950.12559451;
2,3.83307173,6481.0599021;
2,1.92699416,-157.23448173;
2,1.98946862,186.18956071;
2,3.99439281,-6509.01396292;

include,/remote/mbennet/projects/ECP/3p/S/analysis/bases/s_aug-cc-pwCV5Z.basis

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

include,/remote/mbennet/projects/ECP/3p/S/analysis/bases/o_aug-cc-pCV5Z.basis
}

DO s=0,6

 z = 1.15+0.1*s
 
 geometry={
   2
   S Oxide
   S          0.0000000000        0.0000000000       -z*0.5
   O          0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,12,4,2   ! nelec,symm,spin ( B2g * B3g = B1g ) 
    occ,3,2,2,0
 }
 so_hf_bd=energy-s_hf_gs-o_hf_gs

 pop;

 {uccsd(t);
    maxit,500;
    core
 }
 so_cc_bd=energy-s_cc_gs-o_cc_gs

END DO
