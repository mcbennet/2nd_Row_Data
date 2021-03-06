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
3;
1,1.80721061,4.0;
3,9.99633089,7.22884246;
2,2.50043232,-13.0672559;
1;
2,2.26686403,21.20531613;
1;
2,2.11659661,15.43693603;

include,si_aug-cc-pCV5Z_2s2p2d.basis
}

geometry={Si}

{RHF
  wf,4,7,2
  occ,1,0,1,0,1,0,0,0
  open,1.3,1.5
}
si_hf_gs=energy

PUT,XML,wf_bfd.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
si_cc_gs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

basis={
ecp,O,2,1,0;
3;
1,9.29793903,6.0;
3,8.86492204,55.787634;
2,8.62925665,-38.81978498;
1;
2,8.71924452,38.41914135;

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
3;
1,1.80721061,4.0;
3,9.99633089,7.22884246;
2,2.50043232,-13.0672559;
1;
2,2.26686403,21.20531613;
1;
2,2.11659661,15.43693603;

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

 z = 1.75+0.125*s
 
 geometry={
   2
  Si Dimer
  Si          0.0000000000        0.0000000000       -z*0.5
  Si          0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,8,4,2   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 si2_hf_bd=energy
 si2_hf_bd=si2_hf_bd-2.0*si_hf_gs

 {uccsd(t);
    maxit,500;
    core
 }
 si2_cc_bd=energy
 si2_cc_bd=si2_cc_bd-2.0*si_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     SiO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,Si,10,2,0;
3;
1,1.80721061,4.0;
3,9.99633089,7.22884246;
2,2.50043232,-13.0672559;
1;
2,2.26686403,21.20531613;
1;
2,2.11659661,15.43693603;

include,si_aug-cc-pCV5Z_2s2p2d.basis

ecp,O,2,1,0;
3;
1,9.29793903,6.0;
3,8.86492204,55.787634;
2,8.62925665,-38.81978498;
1;
2,8.71924452,38.41914135;

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
 sio_hf_bd=energy
 sio_hf_bd=sio_hf_bd-si_hf_gs-o_hf_gs

 pop;

 {uccsd(t);
    maxit,500;
    core
 }
 sio_cc_bd=energy
 sio_cc_bd=sio_cc_bd-si_cc_gs-o_cc_gs

END DO
