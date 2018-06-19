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

! Douglas-Kroll
SET,DKROLL=1 
SET,DKHO=10 
SET,DKHP=2 

basis={
include,al_aug-cc-pwCV5Z.basis

include,o_aug-cc-pCV5Z.basis
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LONE ATOM CORE ENERGIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! Al !!!!!!!!!!!!!!!!!

geometry={Al}

{RHF
  wf,10,1,0
  occ,2,1,1,0,1,0,0,0
}
al_hf_cs=energy
al_cc_cs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

geometry={O}

{RHF
  wf,2,1,0
  occ,1,0,0,0,0,0,0,0
}
o_hf_cs=energy
o_cc_cs=energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LONE ATOM GROUNDSTATES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! Al !!!!!!!!!!!!!!!!!

geometry={Al}

{RHF
  wf,13,3,1
  occ,3,1,2,0,1,0,0,0
  open,2.3
}
al_hf_gs=energy

PUT,XML,wf_ae.xml,append;keepsph;

{uccsd(t);
   maxit,500;
}
al_cc_gs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

geometry={O}

{RHF
  wf,8,6,2
  occ,2,1,1,0,1,0,0,0
  open,1.2,1.5
}
o_hf_gs=energy

{uccsd(t);
   maxit,500;
}
o_cc_gs=energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

geometry={Al}

!           +2 +1 +0 +0 -1           
nelectrons=[11,12,13,13,14]
symm      =[ 1, 1, 7, 4, 7]
ss        =[ 1, 0, 3, 1, 2]
cceval    =[ 1, 1, 1, 1, 1]
 
!   +2+1+0+0-1
Ag =[3,3,3,3,3]
B3u=[1,1,1,1,1]
B2u=[1,1,2,1,2]
B1g=[0,0,0,1,0]
B1u=[1,1,2,1,2]
B2g=[0,0,0,0,0]
B3g=[0,0,0,0,0]
Au =[0,0,0,0,0]

sing_occ=[3.1,,,,,,,            \ ! +2
,,,,,,,,                        \ ! +1
,3.1,2.3,2.5,,,,,               \ ! +0
,1.4,,,,,,,                     \ ! +0
,2.3,2.5,,,,,,]                   ! -1


DO s=1,#symm

   ! HARTREE-FOCK
   {RHF;maxit,300;
     start,atden
     wf,nelectrons(s),symm(s),ss(s) 
     occ,Ag(s),B3u(s),B2u(s),B1g(s),B1u(s),B2g(s),B3g(s),Au(s)    
     open,sing_occ(8*s-8+1),sing_occ(8*s-8+2),sing_occ(8*s-8+3),sing_occ(8*s-8+4),sing_occ(8*s-8+5),sing_occ(8*s-8+6),sing_occ(8*s-8+7),sing_occ(8*s-8+8);
   }
   hf_energy(s)=energy
   hf_gap(s)=hf_energy(s)-al_hf_gs
   hf_val(s)=hf_energy(s)-al_hf_cs

   IF (cceval(s).EQ.1) THEN

      {uccsd(t)
      }
      cc_energy(s)=energy
      cc_gap(s)=cc_energy(s)-al_cc_gs
      cc_val(s)=cc_energy(s)-al_cc_cs

   END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Al2     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO s=0,6

 y = 2.1+0.2*s
 
 geometry={
   2
   Al Dimer
   Al          0.0000000000        0.0000000000       -y*0.5
   Al          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,26,2,2   ! nelec,symm,spin ( B2g * B3g = B1g )
    occ,5,2,1,0,4,1,1,0
 }
 al2_hf_bd=energy-2.0*al_hf_gs

 {uccsd(t);
    maxit,500;
 }
 al2_cc_bd=energy-2.0*al_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     MgO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO s=0,6

 z = 1.22+0.1*s
 
 geometry={
   2
   Al Oxide
   Al          0.0000000000        0.0000000000       -z*0.5
    O          0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,21,1,1   ! nelec,symm,spin ( B2g * B3g = B1g ) 
    occ,7,2,2,0
    open,7.1
 }
 alo_hf_bd=energy-al_hf_gs-o_hf_gs

 pop;

 {uccsd(t);
    maxit,500;
 }
 alo_cc_bd=energy-al_cc_gs-o_cc_gs

END DO
