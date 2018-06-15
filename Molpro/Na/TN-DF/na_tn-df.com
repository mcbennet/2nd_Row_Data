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
6;
1,0.77882055,1.0;
2,0.54760901,2.12893654;
2,1.18084101,-4.60331698;
2,0.66481312,-5.56778623;
2,1.21372677,6.56778623;
3,0.89497017,0.77882055;
6;
2,1.03099673,2078.97839708;
2,0.75861177,-2077.51151713;
2,0.56809326,15.0371604;
2,0.7403944,1084.57266409;
2,0.83152073,1590.08067606;
2,1.0046283,-2688.69023632;
6;
2,0.99343134,182.53163618;
2,0.65380201,-182.29631553;
2,0.53973673,-34.22694241;
2,0.61550632,183.42506319;
2,0.75946119,51.56329252;
2,0.97984031,-199.76011839;

include,/remote/mbennet/projects/ECP/3p/Na/analysis/bases/na_aug-cc-pCV5Z-DK.basis
}

geometry={Na}

{RHF
  wf,1,1,1
  occ,1,0,0,0,0,0,0,0
  open,1.1
}
na_hf_gs=energy

PUT,XML,wf_tn-df.xml,append;keepsph;

{uccsd(t);
   maxit,500;
   core
}
na_cc_gs=energy

!!!!!!!!!!!! O !!!!!!!!!!!!!!!!!

basis={
ecp,O,2,2,0;
6;
1,8.86932353,6.0;
2,6.05326172,-28.04199287;
2,5.51480979,11.15704031;
2,10.77878397,180.8243251;
2,10.23693413,-179.8243251;
3,7.90462675,53.21594115;
6;
2,7.28893859,-9212.20980516;
2,6.0597119,9226.8656795;
2,10.83143357,58203.26727502;
2,5.75281092,-5120.48607364;
2,10.51155711,-93321.50266843;
2,9.72227746,40239.72318888;
6;
2,7.43321349,10001.55649464;
2,5.85047476,-10012.86801601;
2,5.79011164,8554.95973537;
2,8.08750969,-20342.33136146;
2,8.4322992,11739.44079236;
2,4.71055456,48.9283704;

include,/remote/mbennet/projects/ECP/3p/Na/analysis/bases/o_aug-cc-pCV5Z.basis
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
6;
1,0.77882055,1.0;
2,0.54760901,2.12893654;
2,1.18084101,-4.60331698;
2,0.66481312,-5.56778623;
2,1.21372677,6.56778623;
3,0.89497017,0.77882055;
6;
2,1.03099673,2078.97839708;
2,0.75861177,-2077.51151713;
2,0.56809326,15.0371604;
2,0.7403944,1084.57266409;
2,0.83152073,1590.08067606;
2,1.0046283,-2688.69023632;
6;
2,0.99343134,182.53163618;
2,0.65380201,-182.29631553;
2,0.53973673,-34.22694241;
2,0.61550632,183.42506319;
2,0.75946119,51.56329252;
2,0.97984031,-199.76011839;

include,/remote/mbennet/projects/ECP/3p/Na/analysis/bases/na_aug-cc-pCV5Z-DK.basis
}

geometry={Na}

!           +0 +0 -1 -1           
nelectrons=[ 1, 1, 2, 2]
symm      =[ 3, 4, 1, 3]
ss        =[ 1, 1, 0, 2]
eval      =[ 1, 1, 1, 1]
cceval    =[ 1, 0, 1, 1]

!   +0+0-1-1           
ag =[0,0,1,1] ! x^2, y^2, z^2
b3u=[0,0,0,0] ! x
b2u=[1,0,0,1] ! y
b1g=[0,1,0,0] ! xy
b1u=[0,0,0,0] ! z
b2g=[0,0,0,0] ! xz
b3g=[0,0,0,0] ! yz
au =[0,0,0,0] ! xyz

sing_occ=[1.3,,,,,,,     ,1.4,,,,,,,        ,,,,,,,,          ,1.1,1.3,,,,,,]              


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
 
 geometry={
   2
   Na Dimer
   Na          0.0000000000        0.0000000000       -y*0.5
   Na          0.0000000000        0.0000000000        y*0.5
 }

 ! HARTREE-FOCK
 {RHF
    wf,2,1,0   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 na2_hf_bd=energy-2.0*na_hf_gs

 {uccsd(t);
    maxit,500;
    core
 }
 na2_cc_bd=energy-2.0*na_cc_gs

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     NaO     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis={
ecp,Na,10,2,0;
6;
1,0.77882055,1.0;
2,0.54760901,2.12893654;
2,1.18084101,-4.60331698;
2,0.66481312,-5.56778623;
2,1.21372677,6.56778623;
3,0.89497017,0.77882055;
6;
2,1.03099673,2078.97839708;
2,0.75861177,-2077.51151713;
2,0.56809326,15.0371604;
2,0.7403944,1084.57266409;
2,0.83152073,1590.08067606;
2,1.0046283,-2688.69023632;
6;
2,0.99343134,182.53163618;
2,0.65380201,-182.29631553;
2,0.53973673,-34.22694241;
2,0.61550632,183.42506319;
2,0.75946119,51.56329252;
2,0.97984031,-199.76011839;
include,/remote/mbennet/projects/ECP/3p/Na/analysis/bases/na_aug-cc-pCV5Z-DK.basis

ecp,O,2,2,0;
6;
1,8.86932353,6.0;
2,6.05326172,-28.04199287;
2,5.51480979,11.15704031;
2,10.77878397,180.8243251;
2,10.23693413,-179.8243251;
3,7.90462675,53.21594115;
6;
2,7.28893859,-9212.20980516;
2,6.0597119,9226.8656795;
2,10.83143357,58203.26727502;
2,5.75281092,-5120.48607364;
2,10.51155711,-93321.50266843;
2,9.72227746,40239.72318888;
6;
2,7.43321349,10001.55649464;
2,5.85047476,-10012.86801601;
2,5.79011164,8554.95973537;
2,8.08750969,-20342.33136146;
2,8.4322992,11739.44079236;
2,4.71055456,48.9283704;
include,/remote/mbennet/projects/ECP/3p/Na/analysis/bases/o_aug-cc-pCV5Z.basis
}

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
    wf,7,1,1   ! nelec,symm,spin ( B2g * B3g = B1g ) 
 }
 nao_hf_bd=energy-na_hf_gs-o_hf_gs

 pop;

 {uccsd(t);
    maxit,500;
    core
 }
 nao_cc_bd=energy-na_cc_gs-o_cc_gs

END DO