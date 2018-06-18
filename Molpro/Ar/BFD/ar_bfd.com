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

!!!!!!!!!!!! Ar !!!!!!!!!!!!!!!!!

basis={
ecp,Ar,10,2,0;
3;
1,3.09403094,8.0;
3,6.53700323,24.75224749;
2,3.35769859,-20.38446872;
1;
2,3.68203169,30.67006675;
1;
2,3.45735664,20.84338017;

include,/remote/mbennet/projects/ECP/3p/Ar/analysis/bases/ar_aug-cc-pwCV5Z.basis
}

geometry={Ar}

{RHF
  wf,8,1,0
  occ,1,1,1,0,1,0,0,0
}
ar_hf_gs=energy

PUT,XML,wf_bfd.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
ar_cc_gs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

basis={
ECP,H,0,1,0
3;
2,   21.77696655044365 ,     -10.85192405303825 ;
 1,   21.24359508259891 ,       1.00000000000000;
 3,   21.24359508259891 ,      21.24359508259891;
1; 2,                   1.,                      0.;

include,/remote/mbennet/projects/ECP/3p/Ar/analysis/bases/h_aug-cc-pV5Z.basis
}

geometry={H}

{RHF
  wf,1,1,1
  occ,1,0,0,0,0,0,0,0
  open,1.1
}
h_hf_gs=energy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,Ar,10,2,0;
3;
1,3.09403094,8.0;
3,6.53700323,24.75224749;
2,3.35769859,-20.38446872;
1;
2,3.68203169,30.67006675;
1;
2,3.45735664,20.84338017;

include,/remote/mbennet/projects/ECP/3p/Ar/analysis/bases/ar_aug-cc-pwCV5Z.basis
}

geometry={Ar}

!           +7 +6 +6 +5 +5 +4 +4 +4 +3 +3 +3 +2 +2 +2 +2 +1 +1 +1 +1  0  0  0           
nelectrons=[11,12,12,13,13,14,14,14,15,15,15,16,16,16,16,17,17,17,17,18,18,18]
symm      =[ 1, 1, 5, 2, 4, 4, 1, 8, 8, 5, 2, 2, 5, 7, 1, 2, 1, 2, 8, 3, 6, 8]
ss        =[ 1, 0, 2, 1, 3, 2, 0, 4, 3, 5, 1, 6, 4, 2, 0, 1, 3, 5, 7, 2, 4, 6]
eval      =[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
cceval    =[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

!   +7+6+6+5+5+4+4+4+3+3+3+2+2+2+2+1+1+1+1 0 0 0           
Ag =[3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3] ! x^2, y^2, z^2
B3u=[1,1,1,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2] ! x
B2u=[1,1,1,1,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2] ! y
B1g=[0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,1,1,1,1,1,1] ! xy
B1u=[1,1,2,1,1,1,1,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2] ! z
B2g=[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,1,1] ! xz
B3g=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1] ! yz
Au =[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] ! xyz

sing_occ=[1.1,,,,,,,            \ ! +7
,,,,,,,,                        \ ! +6
,1.1,1.5,,,,,,                  \ ! +6
,1.2,,,,,,,                     \ ! +5
,1.1,1.2,1.3,,,,,               \ ! +5
,1.2,1.3,,,,,,                  \ ! +4
,,,,,,,,                        \ ! +4
,1.1,1.2,1.3,1.5,,,,            \ ! +4
,1.2,1.3,1.5,,,,,               \ ! +3
,1.1,1.4,1.2,1.3,1.5,,,         \ ! +3 *
,1.2,,,,,,,                     \ ! +3
,1.1,1.4,1.6,1.2,1.3,1.5,,      \ ! +2 *
,1.4,1.2,1.3,1.5,,,,            \ ! +2 *
,1.3,1.5,,,,,,                  \ ! +2
,,,,,,,,                        \ ! +2
,1.2,,,,,,,                     \ ! +1
,1.4,1.2,1.3,,,,,               \ ! +1 *
,1.4,1.6,1.2,1.3,1.5,,,         \ ! +1 *
,1.1,1.4,1.6,1.2,1.3,1.5,1.7,   \ ! +1 *
,1.4,1.2,,,,,,                  \ !  0 *
,1.4,1.6,1.2,1.3,,,,            \ !  0 *
,1.4,1.6,1.2,1.3,1.5,1.7,,]       !  0 *

DO s=1,#symm

   ! HARTREE-FOCK
   {RHF;maxit,300;
     start,atden
     wf,nelectrons(s)-10,symm(s),ss(s) 
     occ,Ag(s)-2,B3u(s)-1,B2u(s)-1,B1g(s),B1u(s)-1,B2g(s),B3g(s),Au(s)    
     open,sing_occ(8*s-8+1),sing_occ(8*s-8+2),sing_occ(8*s-8+3),sing_occ(8*s-8+4),sing_occ(8*s-8+5),sing_occ(8*s-8+6),sing_occ(8*s-8+7),sing_occ(8*s-8+8);
   }
   hf_energy(s)=energy
   hf_gap(s)=hf_energy(s)-ar_hf_gs
   hf_val(s)=hf_energy(s)-ar_hf_cs

   IF (cceval(s).EQ.1) THEN

      {uccsd(t)
       core
      }
      cc_energy(s)=energy
      cc_gap(s)=cc_energy(s)-ar_cc_gs
      cc_val(s)=cc_energy(s)-ar_cc_cs

   END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     ArH+     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
basis={
ecp,Ar,10,2,0;
3;
1,3.09403094,8.0;
3,6.53700323,24.75224749;
2,3.35769859,-20.38446872;
1;
2,3.68203169,30.67006675;
1;
2,3.45735664,20.84338017;

include,/remote/mbennet/projects/ECP/3p/Ar/analysis/bases/ar_aug-cc-pwCV5Z.basis

ECP,H,0,1,0
3;
2,   21.77696655044365 ,     -10.85192405303825 ;
 1,   21.24359508259891 ,       1.00000000000000;
 3,   21.24359508259891 ,      21.24359508259891;
1; 2,                   1.,                      0.;

include,/remote/mbennet/projects/ECP/3p/Ar/analysis/bases/h_aug-cc-pV5Z.basis
}

DO s=0,7

 y = 0.8+0.1*s

 geometry={
   2
   ArH+
   Ar          0.0000000000        0.0000000000       -y*0.5
    H          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,8,1,0   ! nelec,symm,spin ( B2g * B3g = B1g )
 }
 arh_hf_bd=energy-ar_hf_gs

 {uccsd(t);
    maxit,500;
    core
 }
 arh_cc_bd=energy-ar_cc_gs

END DO

