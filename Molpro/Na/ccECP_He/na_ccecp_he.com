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
ecp,Na,10,2,0;
5;
1,0.0,1.0;
1,0.0,-9.0;
1,11.999835000000001,9.0;
3,11.949854,107.998515;
2,12.040046,-69.25461899999999;
1;
2,29.947121999999997,221.734689;
1;
2,1.0,0.0;

include,na_aug-cc-pCV5Z-DK.basis
}

geometry={Na}

{RHF
  wf,9,1,1
  occ,2,1,1,0,1,0,0,0
  open,2.1
}
na_hf_gs=energy

PUT,XML,wf_001_sc.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
na_cc_gs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

basis={
ecp,O,2,1,0;
3;
1,12.30997,6.0;
3,14.769620000000002,73.85984;
2,13.714189999999999,-47.876;
1;
2,13.655120000000002,85.86406;

include,o_aug-cc-pCV5Z.basis
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
ecp,Na,10,2,0;
5;
1,0.0,1.0;
1,0.0,-9.0;
1,11.999835000000001,9.0;
3,11.949854,107.998515;
2,12.040046,-69.25461899999999;
1;
2,29.947121999999997,221.734689;
1;
2,1.0,0.0;

include,na_aug-cc-pCV5Z-DK.basis
}

geometry={Na}

!           +0 +0 -1 -1           
nelectrons=[ 9, 9,10,10]
symm      =[ 3, 4, 1, 3]
ss        =[ 1, 1, 0, 2]
eval      =[ 1, 1, 1, 1]
cceval    =[ 1, 0, 1, 1]

!   +0+0-1-1           
ag =[1,1,2,2] ! x^2, y^2, z^2
b3u=[1,1,1,1] ! x
b2u=[2,1,1,2] ! y
b1g=[0,1,0,0] ! xy
b1u=[1,1,1,1] ! z
b2g=[0,0,0,0] ! xz
b3g=[0,0,0,0] ! yz
au =[0,0,0,0] ! xyz

sing_occ=[2.3,,,,,,,     ,1.4,,,,,,,        ,,,,,,,,          ,2.1,2.3,,,,,,]              


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
       core
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

 dzeff=(9.0*9.0-1.0*1.0)/y
 
 geometry={
   2
   Na Dimer
   Na          0.0000000000        0.0000000000       -y*0.5
   Na          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,18,1,0   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 na2_hf_bd=energy-2.0*na_hf_gs+dzeff

 {uccsd(t);
    maxit,500;
    core
 }
 na2_cc_bd=energy-2.0*na_cc_gs+dzeff

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     NaO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,Na,10,2,0;
5;
1,0.0,1.0;
1,0.0,-9.0;
1,11.999835000000001,9.0;
3,11.949854,107.998515;
2,12.040046,-69.25461899999999;
1;
2,29.947121999999997,221.734689;
1;
2,1.0,0.0;
include,na_aug-cc-pCV5Z-DK.basis

ecp,O,2,1,0;
3;
1,12.30997,6.0;
3,14.769620000000002,73.85984;
2,13.714189999999999,-47.876;
1;
2,13.655120000000002,85.86406;
include,o_aug-cc-pCV5Z.basis
}

DO s=0,6

 z = 1.475+0.25*s

 dzeff=(6.0*9.0-1.0*6.0)/z 

 geometry={
   2
   Na Oxide
   Na          0.0000000000        0.0000000000       -z*0.5
    O          0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,15,1,1   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 nao_hf_bd=energy-na_hf_gs-o_hf_gs+dzeff

 pop;

 {uccsd(t);
    maxit,500;
    core
 }
 nao_cc_bd=energy-na_cc_gs-o_cc_gs+dzeff

END DO
