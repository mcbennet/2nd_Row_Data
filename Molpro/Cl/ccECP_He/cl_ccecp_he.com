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

!!!!!!!!!!!! Cl !!!!!!!!!!!!!!!!!

basis={
ecp,Cl,10,2,0;
5;
1,0.0,-15.0;
1,0.0,7.0;
1,20.01071281,15.0;
3,21.44804874,300.16069215;
2,21.17991361,-151.85882234;
2;
2,44.02974125,202.03781785;
2,65.80681305,25.00902536;
1;
2,1.0,0.0;

include,cl_aug-cc-pwCV5Z.basis
}

geometry={Cl}

{RHF
  wf,8,1,0
  occ,1,1,1,0,1,0,0,0
}
cl_hf_cs=energy

{uccsd(t);
   maxit,500;
   core
}
cl_cc_cs=energy

{RHF
  wf,15,2,1
  occ,2,2,2,0,2,0,0,0
  open,2.2
}
cl_hf_gs=energy

PUT,XML,wf_003cb_mp_sc.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
cl_cc_gs=energy

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
ecp,Cl,10,2,0;
5;
1,0.0,-15.0;
1,0.0,7.0;
1,20.01071281,15.0;
3,21.44804874,300.16069215;
2,21.17991361,-151.85882234;
2;
2,44.02974125,202.03781785;
2,65.80681305,25.00902536;
1;
2,1.0,0.0;

include,cl_aug-cc-pwCV5Z.basis
}

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

sing_occ=[2.1,,,,,,,            \ ! +6
,,,,,,,,                        \ ! +5
,2.1,2.5,,,,,,                  \ ! +5
,2.2,,,,,,,                     \ ! +4
,2.1,2.2,2.3,,,,,               \ ! +4
,2.2,2.3,,,,,,                  \ ! +3
,,,,,,,,                        \ ! +3
,2.1,2.2,2.3,2.5,,,,            \ ! +3
,2.2,2.3,2.5,,,,,               \ ! +2
,2.1,1.4,2.2,2.3,2.5,,,         \ ! +2 *
,2.2,,,,,,,                     \ ! +2
,2.1,1.4,1.6,2.2,2.3,2.5,,      \ ! +1 *
,1.4,2.2,2.3,2.5,,,,            \ ! +1 *
,2.3,2.5,,,,,,                  \ ! +1
,,,,,,,,                        \ ! +1
,1.4,2.2,2.3,,,,,               \ ! +0 *
,1.4,1.6,2.2,2.3,2.5,,,         \ ! +0 *
,2.1,1.4,1.6,2.2,2.3,2.5,1.7,   \ ! +0 *
,,,,,,,,]                         ! -1 *

DO s=1,#symm

   ! HARTREE-FOCK
   {RHF;maxit,300;
     start,atden
     wf,nelectrons(s)-2,symm(s),ss(s) 
     occ,Ag(s)-1,B3u(s),B2u(s),B1g(s),B1u(s),B2g(s),B3g(s),Au(s)    
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
    wf,30,1,0   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 cl2_hf_bd=energy-2.0*cl_hf_gs

 {uccsd(t);
    maxit,500;
    core
 }
 cl2_cc_bd=energy-2.0*cl_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Cl O     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,Cl,10,2,0;
5;
1,0.0,-15.0;
1,0.0,7.0;
1,20.01071281,15.0;
3,21.44804874,300.16069215;
2,21.17991361,-151.85882234;
2;
2,44.02974125,202.03781785;
2,65.80681305,25.00902536;
1;
2,1.0,0.0;

include,cl_aug-cc-pwCV5Z.basis

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

 z = 1.17+0.1*s
 
 geometry={
   2
   Cl Oxide
   Cl          0.0000000000        0.0000000000       -z*0.5
   O           0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,21,2,1   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 clo_hf_bd=energy-cl_hf_gs-o_hf_gs

 pop;

 {uccsd(t);
    maxit,500;
    core
 }
 clo_cc_bd=energy-cl_cc_gs-o_cc_gs

END DO
