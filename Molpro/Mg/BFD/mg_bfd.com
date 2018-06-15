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

!!!!!!!!!!!! Mg !!!!!!!!!!!!!!!!!

basis={
ecp,Mg,10,2,0;
3;
1,4.48537898,2.0;
3,2.24268949,8.97075796;
2,1.59710474,-7.72153408;
1;
2,1.57188656,15.07848048;
1;
2,1.42757357,12.37888383;

include,/remote/mbennet/projects/ECP/3p/Mg/analysis/bases/mg_aug-cc-pCV5Z-DK.basis
}

geometry={Mg}

{RHF
  wf,2,1,0
  occ,1,0,0,0,0,0,0,0
}
mg_hf_gs=energy

PUT,XML,wf_bfd.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
mg_cc_gs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

basis={
ecp,O,2,1,0;
3;
1,9.29793903,6.0;
3,8.86492204,55.787634;
2,8.62925665,-38.81978498;
1;
2,8.71924452,38.41914135;

include,/remote/mbennet/projects/ECP/3p/Mg/analysis/bases/o_aug-cc-pCV5Z.basis
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
ecp,Mg,10,2,0;
3;
1,4.48537898,2.0;
3,2.24268949,8.97075796;
2,1.59710474,-7.72153408;
1;
2,1.57188656,15.07848048;
1;
2,1.42757357,12.37888383;

include,/remote/mbennet/projects/ECP/3p/Mg/analysis/bases/mg_aug-cc-pCV5Z-DK.basis
}

geometry={Mg}

!           +1 +0 +0 -1           
nelectrons=[ 1, 2, 2, 3]
symm      =[ 1, 3, 4, 3]
ss        =[ 1, 2, 2, 1]
cceval    =[ 1, 1, 1, 1]
 
!   +1+0+0-1
Ag =[1,1,1,1]
B3u=[0,0,0,0]
B2u=[0,1,0,1]
B1g=[0,0,1,0]
B1u=[0,0,0,0]
B2g=[0,0,0,0]
B3g=[0,0,0,0]
Au =[0,0,0,0]

sing_occ=[1.1,,,,,,,            \ ! +1
,1.1,1.3,,,,,,                  \ ! +0
,1.1,1.4,,,,,,                  \ ! +0
,1.3,,,,,,,]                      ! -1


DO s=1,#symm

   ! HARTREE-FOCK
   {RHF;maxit,300;
     start,atden
     wf,nelectrons(s),symm(s),ss(s) 
     occ,Ag(s),B3u(s),B2u(s),B1g(s),B1u(s),B2g(s),B3g(s),Au(s)    
     open,sing_occ(8*s-8+1),sing_occ(8*s-8+2),sing_occ(8*s-8+3),sing_occ(8*s-8+4),sing_occ(8*s-8+5),sing_occ(8*s-8+6),sing_occ(8*s-8+7),sing_occ(8*s-8+8);
   }
   hf_energy(s)=energy
   hf_gap(s)=hf_energy(s)-mg_hf_gs
   hf_val(s)=hf_energy(s)-mg_hf_cs

   IF (cceval(s).EQ.1) THEN

      {uccsd(t)
       core
      }
      cc_energy(s)=energy
      cc_gap(s)=cc_energy(s)-mg_cc_gs
      cc_val(s)=cc_energy(s)-mg_cc_cs

   END IF

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ATOMIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Mg2     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO s=0,6

 y = 3.0+0.2*s
 
 geometry={
   2
   Mg Dimer
   Mg          0.0000000000        0.0000000000       -y*0.5
   Mg          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,4,1,0   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 mg2_hf_bd=energy-2.0*mg_hf_gs

 {uccsd(t);
    maxit,500;
    core
 }
 mg2_cc_bd=energy-2.0*mg_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     MgO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,Mg,10,2,0;
3;
1,4.48537898,2.0;
3,2.24268949,8.97075796;
2,1.59710474,-7.72153408;
1;
2,1.57188656,15.07848048;
1;
2,1.42757357,12.37888383;
include,/remote/mbennet/projects/ECP/3p/Mg/analysis/bases/mg_aug-cc-pCV5Z-DK.basis

ecp,O,2,1,0;
3;
1,9.29793903,6.0;
3,8.86492204,55.787634;
2,8.62925665,-38.81978498;
1;
2,8.71924452,38.41914135;
include,/remote/mbennet/projects/ECP/3p/Mg/analysis/bases/o_aug-cc-pCV5Z.basis
}

DO s=0,6

 z = 1.4+0.1*s
 
 geometry={
   2
   Mg Oxide
   Mg          0.0000000000        0.0000000000       -z*0.5
    O          0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,8,1,0   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 mgo_hf_bd=energy-mg_hf_gs-o_hf_gs

 pop;

 {uccsd(t);
    maxit,500;
    core
 }
 mgo_cc_bd=energy-mg_cc_gs-o_cc_gs

END DO
