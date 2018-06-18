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
include,/remote/mbennet/projects/ECP/3p/Si/analysis/bases/si_aug-cc-pCV5Z_2s2p2d.basis

include,/remote/mbennet/projects/ECP/3p/Si/analysis/bases/o_aug-cc-pCV5Z.basis
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LONE ATOM CORE ENERGIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! Si !!!!!!!!!!!!!!!!!

geometry={Si}

{RHF
  wf,10,1,0
  occ,2,1,1,0,1,0,0,0
}
si_hf_cs=energy

{uccsd(t);
   maxit,500;
   core
}
si_cc_cs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

geometry={O}

{RHF
  wf,2,1,0
  occ,1,0,0,0,0,0,0,0
}
o_hf_cs=energy

{uccsd(t);
   maxit,500;
   core
}
o_cc_cs=energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LONE ATOM GROUNDSTATES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! Si !!!!!!!!!!!!!!!!!

geometry={Si}

{RHF
  wf,14,7,2
  occ,3,1,2,0,2,0,0,0
  open,2.3,2.5
}
si_hf_gs=energy

PUT,XML,wf_ae.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core,2,1,1,0,1,0,0,0
}
si_uc_gs=energy

{uccsd(t);
   maxit,500;
   core
}
si_cc_gs=energy

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
   core,1,0,0,0,0,0,0,0
}
o_uc_gs=energy

{uccsd(t);
   maxit,500;
   core
}
o_cc_gs=energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

geometry={Si}

!           +3 +3 +3 +2 +2 +2 +2 +2 +2 +2 +2 +1 +1 +1 +1 +1 +0 +0 +0 +0 +0 +0 -1           
nelectrons=[11,11,11,12,12,12,12,12,12,12,12,13,13,13,13,13,14,14,14,14,14,14,15]
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

sing_occ=[3.1,,,,,,,            \ ! +3
,2.3,,,,,,,                     \ ! +3
,1.4,,,,,,,                     \ ! +3
,,,,,,,,                        \ ! +2
,3.1,2.3,,,,,,                  \ ! +2
,,,,,,,,                        \ ! +2
,2.3,2.5,,,,,,                  \ ! +2
,,,,,,,,                        \ ! +2
,3.1,1.4,,,,,,                  \ ! +2
,2.3,1.4,,,,,,                  \ ! +2
,2.5,1.4,,,,,,                  \ ! +2
,2.3,,,,,,,                     \ ! +1
,3.1,2.3,2.5,,,,,               \ ! +1
,3.1,,,,,,,                     \ ! +1
,1.4,,,,,,,                     \ ! +1
,2.3,2.5,2.2,,,,,               \ ! +1
,,,,,,,,                        \ ! +0
,,,,,,,,                        \ ! +0
,3.1,2.3,2.5,2.2,,,,            \ ! +0
,3.1,2.5,,,,,,                  \ ! +0
,2.3,1.4,,,,,,                  \ ! +0
,1.6,1.4,,,,,,                  \ ! +0
,2.2,2.3,2.5,,,,,]                ! -1

DO s=1,#symm

   ! HARTREE-FOCK
   {RHF;maxit,300;
     start,atden
     wf,nelectrons(s),symm(s),ss(s) 
     occ,Ag(s),B3u(s),B2u(s),B1g(s),B1u(s),B2g(s),B3g(s),Au(s)    
     open,sing_occ(8*s-8+1),sing_occ(8*s-8+2),sing_occ(8*s-8+3),sing_occ(8*s-8+4),sing_occ(8*s-8+5),sing_occ(8*s-8+6),sing_occ(8*s-8+7),sing_occ(8*s-8+8);
   }
   hf_energy(s)=energy
   hf_gap(s)=hf_energy(s)-si_hf_gs
   hf_val(s)=hf_energy(s)-si_hf_cs

   IF (cceval(s).EQ.1) THEN

      {uccsd(t)
       core,2,1,1,0,1,0,0,0
      }
      uc_energy(s)=energy
      uc_gap(s)=uc_energy(s)-si_uc_gs
      uc_val(s)=uc_energy(s)-si_hf_cs

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
    wf,28,4,2   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 si2_hf_bd=energy-2.0*si_hf_gs

 {uccsd(t)}
 si2_uc_bd=energy-2.0*si_uc_gs
 
 {uccsd(t);
    maxit,500;
    core
 }
 si2_cc_bd=energy-2.0*si_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     SiO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    wf,22,1,0   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 sio_hf_bd=energy-si_hf_gs-o_hf_gs

 pop;

 {uccsd(t)}
 sio_uc_bd=energy-si_uc_gs-o_uc_gs

 {uccsd(t);
    maxit,500;
    core
 }
 sio_cc_bd=energy-si_cc_gs-o_cc_gs

END DO