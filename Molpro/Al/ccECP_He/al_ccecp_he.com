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

!!!!!!!!!!!! Al !!!!!!!!!!!!!!!!!

basis={
ecp,Al,10,2,0;
5;
1,0.0,-11.0;
1,0.0,3.0;
1,10.997885,11.0;
3,11.22045,120.976735;
2,11.08325,-80.393266;
2;
2,80.97577199999999,25.006358;
2,24.398239999999998,112.360799;
1;
2,1.0,0.0;

include,/remote/mbennet/projects/ECP/3p/Al/analysis_new/bases/al_aug-cc-pwCV5Z.basis
}

geometry={Al}

{RHF
  wf,8,1,0
  occ,1,1,1,0,1,0,0,0
}
al_hf_cs=energy

{uccsd(t);
   maxit,500;
   core
}
al_cc_cs=energy

{RHF
  wf,11,3,1
  occ,2,1,2,0,1,0,0,0
  open,2.3
}
al_hf_gs=energy

PUT,XML,wf_ae.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
al_cc_gs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

basis={
ecp,O,2,1,0;
3;
1,12.30997,6.0;
3,14.769620000000002,73.85984;
2,13.714189999999999,-47.876;
1;
2,13.655120000000002,85.86406;

include,/remote/mbennet/projects/ECP/3p/Al/analysis_new/bases/o_aug-cc-pCV5Z.basis
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
ecp,Al,10,2,0;
5;
1,0.0,-11.0;
1,0.0,3.0;
1,10.997885,11.0;
3,11.22045,120.976735;
2,11.08325,-80.393266;
2;
2,80.97577199999999,25.006358;
2,24.398239999999998,112.360799;
1;
2,1.0,0.0;

include,/remote/mbennet/projects/ECP/3p/Al/analysis_new/bases/al_aug-cc-pwCV5Z.basis
}

geometry={Al}

!           +2 +1 +0 +0 -1           
nelectrons=[ 9,10,11,11,12]
symm      =[ 1, 1, 7, 4, 7]
ss        =[ 1, 0, 3, 1, 2]
cceval    =[ 1, 1, 1, 1, 1]
 
!   +2+1+0+0-1
Ag =[2,2,2,2,2]
B3u=[1,1,1,1,1]
B2u=[1,1,2,1,2]
B1g=[0,0,0,1,0]
B1u=[1,1,2,1,2]
B2g=[0,0,0,0,0]
B3g=[0,0,0,0,0]
Au =[0,0,0,0,0]

sing_occ=[2.1,,,,,,,            \ ! +2
,,,,,,,,                        \ ! +1
,2.1,2.3,2.5,,,,,               \ ! +0
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
       core
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

 dzeff=(11.0*11.0-3.0*3.0)/(y/TOANG)
 
 geometry={
   2
   Al Dimer
   Al          0.0000000000        0.0000000000       -y*0.5
   Al          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,22,2,2   ! nelec,symm,spin ( B2g * B3g = B1g )
    occ,4,2,1,0,3,1,1,0
 }
 al2_hf_bd=energy-2.0*al_hf_gs+dzeff

 {uccsd(t);
    maxit,500;
    core
 }
 al2_cc_bd=energy-2.0*al_cc_gs+dzeff

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     MgO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,Al,10,2,0;
5;
1,0.0,-11.0;
1,0.0,3.0;
1,10.997885,11.0;
3,11.22045,120.976735;
2,11.08325,-80.393266;
2;
2,80.97577199999999,25.006358;
2,24.398239999999998,112.360799;
1;
2,1.0,0.0;

include,/remote/mbennet/projects/ECP/3p/Al/analysis_new/bases/al_aug-cc-pwCV5Z.basis

ecp,O,2,1,0;
3;
1,12.30997,6.0;
3,14.769620000000002,73.85984;
2,13.714189999999999,-47.876;
1;
2,13.655120000000002,85.86406;

include,/remote/mbennet/projects/ECP/3p/Al/analysis_new/bases/o_aug-cc-pCV5Z.basis
}

DO s=0,6

 z = 1.22+0.1*s

 dzeff=(11.0-3.0)*6.0/(z/TOANG) 

 geometry={
   2
   Al Oxide
   Al          0.0000000000        0.0000000000       -z*0.5
    O          0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF; !shift,-1.0,-0.5
    MAXIT,200
    wf,16,1,0   ! nelec,symm,spin ( B2g * B3g = B1g )
    sym,1,1,1,1,1,1
    occ,4,2,2,0
    save,2104.2
 }

 ! HARTREE-FOCK
 {RHF; shift,-1.0,-0.5
    start,2104.2
    wf,17,1,1   ! nelec,symm,spin ( B2g * B3g = B1g ) 
    sym,1,1,1,1,1,1
    occ,5,2,2,0
    open,5.1
 }
 alo_hf_bd=energy-al_hf_gs-o_hf_gs+dzeff

 pop;

 {uccsd(t);
    maxit,500;
    core   
 }
 alo_cc_bd=energy-al_cc_gs-o_cc_gs+dzeff

END DO
