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
3;
1,3.07301275,3.0;
3,9.97055633,9.21903825;
2,2.6413466,-9.65888637;
1;
2,1.87284747,17.16680112;
1;
2,1.79397208,14.22120694;

include,al_aug-cc-pwCV5Z.basis
}

geometry={Al}

{RHF
  wf,3,3,1
  occ,1,0,1,0,0,0,0,0
  open,1.3
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
ecp,Al,10,2,0;
3;
1,3.07301275,3.0;
3,9.97055633,9.21903825;
2,2.6413466,-9.65888637;
1;
2,1.87284747,17.16680112;
1;
2,1.79397208,14.22120694;

include,al_aug-cc-pwCV5Z.basis
}

geometry={Al}

!           +2 +1 +0 +0 -1           
nelectrons=[ 1, 2, 3, 3, 4]
symm      =[ 1, 1, 7, 4, 7]
ss        =[ 1, 0, 3, 1, 2]
cceval    =[ 1, 1, 1, 1, 1]
 
!   +2+1+0+0-1
Ag =[1,1,1,1,1]
B3u=[0,0,0,0,0]
B2u=[0,0,1,0,1]
B1g=[0,0,0,1,0]
B1u=[0,0,1,0,1]
B2g=[0,0,0,0,0]
B3g=[0,0,0,0,0]
Au =[0,0,0,0,0]

sing_occ=[1.1,,,,,,,            \ ! +2
,,,,,,,,                        \ ! +1
,1.1,1.3,1.5,,,,,               \ ! +0
,1.4,,,,,,,                     \ ! +0
,1.3,1.5,,,,,,]                   ! -1


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
 
 geometry={
   2
   Al Dimer
   Al          0.0000000000        0.0000000000       -y*0.5
   Al          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,6,2,2   ! nelec,symm,spin ( B2g * B3g = B1g )
    occ,2,1,0,0,1,0,0,0
 }
 al2_hf_bd=energy-2.0*al_hf_gs

 {uccsd(t);
    maxit,500;
    core
 }
 al2_cc_bd=energy-2.0*al_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     MgO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,Al,10,2,0;
3;
1,3.07301275,3.0;
3,9.97055633,9.21903825;
2,2.6413466,-9.65888637;
1;
2,1.87284747,17.16680112;
1;
2,1.79397208,14.22120694;

include,al_aug-cc-pwCV5Z.basis

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

 z = 1.22+0.1*s
 
 geometry={
   2
   Al Oxide
   Al          0.0000000000        0.0000000000       -z*0.5
    O          0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,9,1,1   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 alo_hf_bd=energy-al_hf_gs-o_hf_gs

 pop;

 {uccsd(t);
    maxit,500;
    core   
 }
 alo_cc_bd=energy-al_cc_gs-o_cc_gs

END DO
