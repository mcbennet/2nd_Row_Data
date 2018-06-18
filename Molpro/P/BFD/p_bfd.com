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

!!!!!!!!!!!! P !!!!!!!!!!!!!!!!!

basis={
ecp,P,10,2,0;
3;
1,2.0262281,5.0;
3,9.95970113,10.13114051;
2,2.74841795,-14.94375088;
1;
2,2.60470698,23.6247948;
1;
2,2.549579,18.18547203;

include,/remote/mbennet/projects/ECP/3p/P/analysis/bases/p_aug-cc-pCV5Z_1s1p1d.basis
}

geometry={P}

{RHF
  wf,5,8,3
  occ,1,1,1,0,1,0,0,0
  open,1.3,1.5,1.2
}
p_hf_gs=energy

PUT,XML,wf_bfd.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
p_cc_gs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

basis={
ecp,O,2,1,0;
3;
1,9.29793903,6.0;
3,8.86492204,55.787634;
2,8.62925665,-38.81978498;
1;
2,8.71924452,38.41914135;

include,/remote/mbennet/projects/ECP/3p/P/analysis/bases/o_aug-cc-pCV5Z.basis
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
ecp,P,10,2,0;
3;
1,2.0262281,5.0;
3,9.95970113,10.13114051;
2,2.74841795,-14.94375088;
1;
2,2.60470698,23.6247948;
1;
2,2.549579,18.18547203;

include,/remote/mbennet/projects/ECP/3p/P/analysis/bases/p_aug-cc-pCV5Z_1s1p1d.basis
}

geometry={P}

!           +4 +4 +4 +3 +3 +3 +3 +3 +3 +3 +3 +3 +2 +2 +2 +2 +2 +2 +2 +1 +1 +1 +1 +1 +0 +0 +0 +0 +0 -1           
nelectrons=[ 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6]
symm      =[ 1, 3, 4, 1, 3, 7, 1, 1, 4, 2, 6, 1, 3, 4, 7, 2, 8, 6, 1, 7, 1, 2, 6, 8, 6, 8, 5, 2, 4, 6]
ss        =[ 1, 1, 1, 0, 2, 2, 0, 0, 2, 2, 2, 0, 1, 1, 3, 3, 3, 3, 3, 2, 0, 2, 2, 4, 3, 3, 5, 1, 1, 2]
eval      =[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
cceval    =[ 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1]

!   +4+4+4+3+3+3+3+3+3+3+3+3+2+2+2+2+2+2+2+1+1+1+1+1+0+0+0+0+0-1           
Ag =[3,2,2,3,3,2,2,2,3,2,2,2,3,3,3,3,2,3,2,3,3,3,3,3,3,3,3,3,3,3] ! x^2, y^2, z^2
B3u=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,2,1,1,2,1,1,2,2,1,2] ! x
B2u=[1,2,1,1,2,2,2,1,1,2,1,1,2,1,2,2,2,1,1,2,1,2,1,2,2,2,2,2,2,2] ! y
B1g=[0,0,1,0,0,0,0,0,1,1,1,1,0,1,0,1,0,1,1,0,0,1,1,0,1,1,1,0,1,0] ! xy
B1u=[1,1,1,1,1,2,1,2,1,1,1,1,1,1,2,1,2,1,1,2,1,1,1,2,2,1,2,1,1,2] ! z
B2g=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0] ! xz
B3g=[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,1,0,0,1,0,0,0,0] ! yz
Au =[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] ! xyz

sing_occ=[1.1,,,,,,,     ,1.3,,,,,,,        ,1.4,,,,,,,                              \ ! +4
,,,,,,,,                 ,1.1,1.3,,,,,,     ,1.3,1.5,,,,,,     ,,,,,,,,              \ ! +3
,,,,,,,,                 ,1.1,1.4,,,,,,   ,1.3,1.4,,,,,,  ,1.4,1.7,,,,,,  ,,,,,,,,   \ ! +3
,1.3,,,,,,,              ,1.4,,,,,,,        ,1.1,1.3,1.5,,,,,  ,1.1,1.3,1.4,,,,,     \ ! +2
,1.3,1.5,1.2,,,,,        ,1.1,1.4,1.7,,,,,  ,1.4,1.7,1.6,,,,,                        \ ! +2
,1.3,1.5,,,,,,                  \ ! +1
,,,,,,,,                        \ ! +1
,1.3,1.4,,,,,,                  \ ! +1
,1.4,1.7,,,,,,                  \ ! +1
,1.1,1.2,1.3,1.5,,,,            \ ! +1
,1.3,1.5,1.4,,,,,               \ ! +0
,1.3,1.4,1.7,,,,,               \ ! +0
,1.1,1.4,1.2,1.3,1.5,,,         \ ! +0 *
,1.2,,,,,,,                     \ ! +0
,1.4,,,,,,,                     \ ! +0
,1.2,1.5,,,,,,]                   ! -1 *


DO s=1,#symm

   ! HARTREE-FOCK
   {RHF;maxit,300;
     start,atden
     wf,nelectrons(s),symm(s),ss(s) 
     occ,Ag(s)-2,B3u(s)-1,B2u(s)-1,B1g(s),B1u(s)-1,B2g(s),B3g(s),Au(s)    
     open,sing_occ(8*s-8+1),sing_occ(8*s-8+2),sing_occ(8*s-8+3),sing_occ(8*s-8+4),sing_occ(8*s-8+5),sing_occ(8*s-8+6),sing_occ(8*s-8+7),sing_occ(8*s-8+8);
   }
   hf_energy(s)=energy
   hf_gap(s)=hf_energy(s)-p_hf_gs
   hf_val(s)=hf_energy(s)-p_hf_cs

   IF (cceval(s).EQ.1) THEN

      {uccsd(t)
       core
      }
      cc_energy(s)=energy
      cc_gap(s)=cc_energy(s)-p_cc_gs
      cc_val(s)=cc_energy(s)-p_cc_cs

   END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     P2     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO s=0,6

 y = 1.55+0.1*s
 
 geometry={
   2
   P Dimer
   P          0.0000000000        0.0000000000       -y*0.5
   P          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,10,1,0   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 p2_hf_bd=energy-2.0*p_hf_gs

 {uccsd(t);
    maxit,500;
    core
 }
 p2_cc_bd=energy-2.0*p_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     PO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,P,10,2,0;
3;
1,2.0262281,5.0;
3,9.95970113,10.13114051;
2,2.74841795,-14.94375088;
1;
2,2.60470698,23.6247948;
1;
2,2.549579,18.18547203;

include,/remote/mbennet/projects/ECP/3p/P/analysis/bases/p_aug-cc-pCV5Z_1s1p1d.basis

ecp,O,2,1,0;
3;
1,9.29793903,6.0;
3,8.86492204,55.787634;
2,8.62925665,-38.81978498;
1;
2,8.71924452,38.41914135;
include,/remote/mbennet/projects/ECP/3p/P/analysis/bases/o_aug-cc-pCV5Z.basis
}

DO s=0,6

 z = 1.2+0.1*s
 
 geometry={
   2
   P Oxide
   P          0.0000000000        0.0000000000       -z*0.5
   O          0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,11,2,1   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 po_hf_bd=energy-p_hf_gs-o_hf_gs

 pop;

 {uccsd(t);
    maxit,500;
    core
 }
 po_cc_bd=energy-p_cc_gs-o_cc_gs

END DO
