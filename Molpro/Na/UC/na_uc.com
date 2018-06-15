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
include,/remote/mbennet/projects/ECP/3p/Na/analysis/bases/na_aug-cc-pCV5Z-DK.basis

include,/remote/mbennet/projects/ECP/3p/Na/analysis/bases/o_aug-cc-pCV5Z.basis
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LONE ATOM CORE ENERGIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! Na !!!!!!!!!!!!!!!!!

geometry={Na}

{RHF
  wf,10,1,0
  occ,2,1,1,0,1,0,0,0
}
na_hf_cs=energy

na_cc_cs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

geometry={O}

{RHF
  wf,2,1,0
  occ,1,0,0,0,0,0,0,0
}
o_hf_cs=energy

o_cc_cs=energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LONE ATOM GROUNDSTATES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! Na !!!!!!!!!!!!!!!!!

geometry={Na}

{RHF
  wf,11,1,1
  occ,3,1,1,0,1,0,0,0
  open,3.1
}
na_hf_gs=energy

PUT,XML,wf_ae.xml,append;keepsph;

{uccsd(t);
   maxit,500;
}
na_cc_gs=energy

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

geometry={Na}

!           +0 +0 -1 -1           
nelectrons=[11,11,12,12]
symm      =[ 3, 4, 1, 3]
ss        =[ 1, 1, 0, 2]
eval      =[ 1, 1, 1, 1]
cceval    =[ 1, 0, 1, 1]

!   +0+0-1-1           
Ag =[2,2,3,3] ! x^2, y^2, z^2
B3u=[1,1,1,1] ! x
B2u=[2,1,1,2] ! y
B1g=[0,1,0,0] ! xy
B1u=[1,1,1,1] ! z
B2g=[0,0,0,0] ! xz
B3g=[0,0,0,0] ! yz
Au =[0,0,0,0] ! xyz

sing_occ=[2.3,,,,,,,     ,1.4,,,,,,,        ,,,,,,,,          ,3.1,2.3,,,,,,]              


DO s=1,#symm

   ! HARTREE-FOCK
   {RHF;maxit,300;
     start,atden
     wf,nelectrons(s),symm(s),ss(s) 
     occ,Ag(s),B3u(s),B2u(s),B1g(s),B1u(s),B2g(s),B3g(s),Au(s)    
     open,sing_occ(8*s-8+1),sing_occ(8*s-8+2),sing_occ(8*s-8+3),sing_occ(8*s-8+4),sing_occ(8*s-8+5),sing_occ(8*s-8+6),sing_occ(8*s-8+7),sing_occ(8*s-8+8);
   }
   hf_energy(s)=energy
   hf_gap(s)=hf_energy(s)-na_hf_gs
   hf_val(s)=hf_energy(s)-na_hf_cs

   IF (cceval(s).EQ.1) THEN

      {uccsd(t)
      }
      cc_energy(s)=energy
      cc_gap(s)=cc_energy(s)-na_cc_gs
      cc_val(s)=cc_energy(s)-na_cc_cs

   END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Na2     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO s=0,6

 y = 2.2+0.25*s
 
 geometry={
   2
   Na Dimer
   Na          0.0000000000        0.0000000000       -y*0.5
   Na          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,22,1,0   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 na2_hf_bd=energy-2.0*na_hf_gs

 {uccsd(t);
    maxit,500;
 }
 na2_cc_bd=energy-2.0*na_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     NaO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO s=0,6

 z = 1.475+0.25*s
 
 geometry={
   2
   Na Oxide
   Na          0.0000000000        0.0000000000       -z*0.5
    O          0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,19,1,1   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 nao_hf_bd=energy-na_hf_gs-o_hf_gs

 pop;

 {uccsd(t);
    maxit,500;
 }
 nao_cc_bd=energy-na_cc_gs-o_cc_gs

END DO
