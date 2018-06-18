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
include,/remote/mbennet/projects/ECP/3p/Cl/analysis/bases/cl_aug-cc-pwCV5Z.basis

include,/remote/mbennet/projects/ECP/3p/Cl/analysis/bases/o_aug-cc-pCV5Z.basis
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LONE ATOM CORE ENERGIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! Cl !!!!!!!!!!!!!!!!!

geometry={Cl}

{RHF
  wf,10,1,0
  occ,2,1,1,0,1,0,0,0
}
cl_hf_cs=energy

{uccsd(t);
   maxit,500;
   core
}
cl_cc_cs=energy

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

!!!!!!!!!!!! Cl !!!!!!!!!!!!!!!!!

geometry={Cl}

{RHF
  wf,17,2,1
  occ,3,2,2,0,2,0,0,0
  open,2.2
}
cl_hf_gs=energy

PUT,XML,wf_ae.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
cl_cc_gs=energy

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
   core
}
o_cc_gs=energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

geometry={Cl}

!           +6 +5 +5 +4 +4 +3 +3 +3 +2 +2 +2 +1 +1 +1 +1 +0 +0 +0 -1           
nelectrons=[11,12,12,13,13,14,14,14,15,15,15,16,16,16,16,17,17,17,18]
symm      =[ 1, 1, 5, 2, 4, 4, 1, 8, 8, 5, 2, 2, 5, 7, 1, 1, 2, 8, 1]
ss        =[ 1, 0, 2, 1, 3, 2, 0, 4, 3, 5, 1, 6, 4, 2, 0, 3, 5, 7, 0]
eval      =[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
cceval    =[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

!   +6+5+5+4+4+3+3+3+2+2+2+1+1+1+1+0+0+0-1           
Ag =[3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3] ! x^2, y^2, z^2
B3u=[1,1,1,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2] ! x
B2u=[1,1,1,1,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2] ! y
B1g=[0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,1,1,1,0] ! xy
B1u=[1,1,2,1,1,1,1,2,2,2,1,2,2,2,2,2,2,2,2] ! z
B2g=[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0] ! xz
B3g=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0] ! yz
Au =[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] ! xyz

sing_occ=[3.1,,,,,,,            \ ! +6 *
,,,,,,,,                        \ ! +5 *
,3.1,2.5,,,,,,                  \ ! +5 *
,2.2,,,,,,,                     \ ! +4 *
,3.1,2.2,2.3,,,,,               \ ! +4 *
,2.2,2.3,,,,,,                  \ ! +3 *
,,,,,,,,                        \ ! +3 *
,3.1,2.2,2.3,2.5,,,,            \ ! +3 *
,2.2,2.3,2.5,,,,,               \ ! +2 *
,3.1,1.4,2.2,2.3,2.5,,,         \ ! +2 * 
,2.2,,,,,,,                     \ ! +2 *
,3.1,1.4,1.6,2.2,2.3,2.5,,      \ ! +1 *
,1.4,2.2,2.3,2.5,,,,            \ ! +1 * 
,2.3,2.5,,,,,,                  \ ! +1 *
,,,,,,,,                        \ ! +1 *
,1.4,2.2,2.3,,,,,               \ ! +0 * 
,1.4,1.6,2.2,2.3,2.5,,,         \ ! +0 
,3.1,1.4,1.6,2.2,2.3,2.5,1.7,   \ ! +0 
,,,,,,,,]                         ! -1 

DO s=1,#symm

   ! HARTREE-FOCK
   {RHF;maxit,300;
     start,atden
     wf,nelectrons(s),symm(s),ss(s) 
     occ,Ag(s),B3u(s),B2u(s),B1g(s),B1u(s),B2g(s),B3g(s),Au(s)    
     open,sing_occ(8*s-8+1),sing_occ(8*s-8+2),sing_occ(8*s-8+3),sing_occ(8*s-8+4),sing_occ(8*s-8+5),sing_occ(8*s-8+6),sing_occ(8*s-8+7),sing_occ(8*s-8+8);
   }
   hf_energy(s)=energy
   hf_gap(s)=hf_energy(s)-cl_hf_gs
   hf_val(s)=hf_energy(s)-cl_hf_cs

   IF (cceval(s).EQ.1) THEN

      {uccsd(t)
       core
      }
      cc_energy(s)=energy
      cc_gap(s)=cc_energy(s)-cl_cc_gs
      cc_val(s)=cc_energy(s)-cl_cc_cs

   END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Cl2     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO s=0,6

 y = 1.6+0.1*s
 
 geometry={
   2
   Cl Dimer
   Cl          0.0000000000        0.0000000000       -y*0.5
   Cl          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,34,1,0   ! nelec,symm,spin ( B2g * B3g = B1g )
 }
 Cl2_hf_bd=energy-2.0*cl_hf_gs

 {uccsd(t);
    maxit,500;
    core
 }
 Cl2_cc_bd=energy-2.0*cl_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Cl O     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO s=0,6

 z = 1.17+0.1*s
 
 geometry={
   2
   Cl Oxide
   Cl          0.0000000000        0.0000000000       -z*0.5
   O           0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,25,2,1   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 clo_hf_bd=energy-cl_hf_gs-o_hf_gs

 pop;

 {uccsd(t);
    maxit,500;
    core   
 }
 clo_cc_bd=energy-cl_cc_gs-o_cc_gs

END DO
