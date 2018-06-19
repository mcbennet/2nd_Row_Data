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

!!!!!!!!!!!! Si !!!!!!!!!!!!!!!!!

basis={
ecp,Si,10,2,0;
5;
1,0.0,4.0;
1,0.0,-12.0;
1,16.146033,12.0;
3,17.077556,193.752396;
2,16.425284,-104.82393499999999;
1;
2,22.549338,96.908435;
1;
2,1.0,0.0;

include,si_aug-cc-pCV5Z_2s2p2d.basis
}

geometry={Si}

{RHF
  wf,8,1,0
  occ,1,1,1,0,1,0,0,0
}
si_hf_cs=energy

{uccsd(t);
   maxit,500;
   core
}
si_cc_cs=energy

{RHF
  wf,12,7,2
  occ,2,1,2,0,2,0,0,0
  open,2.3,2.5
}
si_hf_gs=energy

PUT,XML,wf_opt_lm_small-core_initscale.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
si_cc_gs=energy

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
ecp,Si,10,2,0;
5;
1,0.0,4.0;
1,0.0,-12.0;
1,16.146033,12.0;
3,17.077556,193.752396;
2,16.425284,-104.82393499999999;
1;
2,22.549338,96.908435;
1;
2,1.0,0.0;

include,si_aug-cc-pCV5Z_2s2p2d.basis
}

geometry={Si}

!           +3 +3 +3 +2 +2 +2 +2 +2 +2 +2 +2 +1 +1 +1 +1 +1 +0 +0 +0 +0 +0 +0 -1           
nelectrons=[ 9, 9, 9,10,10,10,10,10,10,10,10,11,11,11,11,11,12,12,12,12,12,12,13]
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

sing_occ=[2.1,,,,,,,            \ ! +3
,2.3,,,,,,,                     \ ! +3
,1.4,,,,,,,                     \ ! +3
,,,,,,,,                        \ ! +2
,2.1,2.3,,,,,,                  \ ! +2
,,,,,,,,                        \ ! +2
,2.3,2.5,,,,,,                  \ ! +2
,,,,,,,,                        \ ! +2
,2.1,1.4,,,,,,                  \ ! +2
,2.3,1.4,,,,,,                  \ ! +2
,2.5,1.4,,,,,,                  \ ! +2
,2.3,,,,,,,                     \ ! +1
,2.1,2.3,2.5,,,,,               \ ! +1
,2.1,,,,,,,                     \ ! +1
,1.4,,,,,,,                     \ ! +1
,2.3,2.5,2.2,,,,,               \ ! +1
,,,,,,,,                        \ ! +0
,,,,,,,,                        \ ! +0
,2.1,2.3,2.5,2.2,,,,            \ ! +0
,2.1,2.5,,,,,,                  \ ! +0
,2.3,1.4,,,,,,                  \ ! +0
,1.6,1.4,,,,,,                  \ ! +0
,2.2,2.3,2.5,,,,,]                ! -1

DO s=1,#symm

   ! HARTREE-FOCK
   {RHF;maxit,300;
     start,atden
     wf,nelectrons(s),symm(s),ss(s) 
     occ,Ag(s)-1,B3u(s),B2u(s),B1g(s),B1u(s),B2g(s),B3g(s),Au(s)    
     open,sing_occ(8*s-8+1),sing_occ(8*s-8+2),sing_occ(8*s-8+3),sing_occ(8*s-8+4),sing_occ(8*s-8+5),sing_occ(8*s-8+6),sing_occ(8*s-8+7),sing_occ(8*s-8+8);
   }
   hf_energy(s)=energy
   hf_gap(s)=hf_energy(s)-si_hf_gs
   hf_val(s)=hf_energy(s)-si_hf_cs

   IF (cceval(s).EQ.1) THEN

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

 dzeff=(12.0*12.0-4.0*4.0)/(y*1.88973)
 
 geometry={
   2
  Si Dimer
  Si          0.0000000000        0.0000000000       -y*0.5
  Si          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,24,4,2   ! nelec,symm,spin ( B2g * B3g = B1g ) 
    sym,1,1,1,1,1
    occ,4,2,2,0,3,1,1,0
    open,2.2,2.3
 }
 si2_hf_bd=energy-2.0*si_hf_gs+dzeff

 {uccsd(t);
    maxit,500;
    core
 }
 si2_cc_bd=energy-2.0*si_cc_gs+dzeff

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     SiO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,Si,10,2,0;
5;
1,0.0,4.0;
1,0.0,-12.0;
1,16.146033,12.0;
3,17.077556,193.752396;
2,16.425284,-104.82393499999999;
1;
2,22.549338,96.908435;
1;
2,1.0,0.0;

include,si_aug-cc-pCV5Z_2s2p2d.basis

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

 z = 1.15+0.1*s
 
 dzeff=(12.0*6.0-4.0*6.0)/(z*1.88973)

 geometry={
   2
  Si Oxide
  Si          0.0000000000        0.0000000000       -z*0.5
   O          0.0000000000        0.0000000000        z*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,18,1,0   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 sio_hf_bd=energy-si_hf_gs-o_hf_gs+dzeff

 pop;

 {uccsd(t);
    maxit,500;
    core
 }
 sio_cc_bd=energy-si_cc_gs-o_cc_gs+dzeff

END DO
