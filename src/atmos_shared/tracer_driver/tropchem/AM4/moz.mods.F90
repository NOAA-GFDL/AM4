





      module mo_grid_mod
!---------------------------------------------------------------------
! ... Basic grid point resolution parameters
!---------------------------------------------------------------------
      implicit none

      save

      integer, parameter :: &
                pcnst = 99 +1, & ! number of advected constituents including cloud water
                pcnstm1 = 99, & ! number of advected constituents excluding cloud water
                plev = 1, & ! number of vertical levels
                plevp = plev+1, & ! plev plus 1
                plevm = plev-1, & ! plev minus 1
                plon = 1, & ! number of longitudes
                plat = 1 ! number of latitudes

      integer, parameter :: &
                pnats = 0 ! number of non-advected trace species






      integer :: nodes ! mpi task count
      integer :: plonl ! longitude tile dimension
      integer :: pplon ! longitude tile count
      integer :: plnplv ! plonl * plev

      end module mo_grid_mod

      module chem_mods_mod
!--------------------------------------------------------------
! ... basic chemistry array parameters
!--------------------------------------------------------------

      use mo_grid_mod, only : pcnstm1

      use mpp_mod, only : mpp_error, FATAL

      implicit none

      save

      integer, parameter :: hetcnt = 0, & ! number of heterogeneous processes
                            phtcnt = 43, & ! number of photo processes
                            rxntot = 247, & ! number of total reactions
                            gascnt = 204, & ! number of gas phase reactions
                            nfs = 3, & ! number of "fixed" species
                            relcnt = 0, & ! number of relationship species
                            grpcnt = 0, & ! number of group members
                            imp_nzcnt = 998, & ! number of non-zero implicit matrix entries
                            rod_nzcnt = 0, & ! number of non-zero rodas matrix entries
                            extcnt = 0, & ! number of species with external forcing
                            clscnt1 = 8, & ! number of species in explicit class
                            clscnt2 = 0, & ! number of species in hov class
                            clscnt3 = 0, & ! number of species in ebi class
                            clscnt4 = 91, & ! number of species in implicit class
                            clscnt5 = 0, & ! number of species in rodas class
                            indexm = 1, & ! index of total atm density in invariant array
                            ncol_abs = 2, & ! number of column densities
                            indexh2o = 0, & ! index of water vapor density
                            clsze = 1 ! loop length for implicit chemistry

      integer :: ngrp = 0
      integer :: drydep_cnt = 0
      integer :: srfems_cnt = 0
      integer :: rxt_alias_cnt = 0
      integer, allocatable :: grp_mem_cnt(:)
      integer, allocatable :: rxt_alias_map(:)
      real :: adv_mass(pcnstm1)
      real :: nadv_mass(grpcnt)
      character(len=16), allocatable :: rxt_alias_lst(:)
      character(len=8), allocatable :: drydep_lst(:)
      character(len=8), allocatable :: srfems_lst(:)
      character(len=8), allocatable :: grp_lst(:)
      character(len=8) :: het_lst(max(1,hetcnt))
      character(len=8) :: extfrc_lst(max(1,extcnt))
      character(len=8) :: inv_lst(max(1,nfs))

      type solver_class
  integer :: clscnt
  integer :: lin_rxt_cnt
  integer :: nln_rxt_cnt
  integer :: indprd_cnt
  integer :: iter_max
         integer :: cls_rxt_cnt(4)
         integer, pointer :: permute(:)
         integer, pointer :: diag_map(:)
         integer, pointer :: clsmap(:)
      end type solver_class

      type(solver_class) :: explicit, implicit, rodas

      contains
      subroutine endrun(msg)

      implicit none

      character(len=128), intent(in), optional :: msg
      call mpp_error(FATAL, msg)

      end subroutine endrun

      subroutine chem_mods_init
!--------------------------------------------------------------
! ... intialize the class derived type
!--------------------------------------------------------------

      implicit none

      integer :: astat

      explicit%clscnt = 8
      explicit%indprd_cnt = 56

      implicit%clscnt = 91
      implicit%lin_rxt_cnt = 76
      implicit%nln_rxt_cnt = 168
      implicit%indprd_cnt = 3
      implicit%iter_max = 11

      rodas%clscnt = 0
      rodas%lin_rxt_cnt = 0
      rodas%nln_rxt_cnt = 0
      rodas%indprd_cnt = 0

      if( explicit%clscnt > 0 ) then
  allocate( explicit%clsmap(explicit%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate explicit%clsmap ; error = ',astat
     call endrun
  end if
         explicit%clsmap(:) = 0
      end if
      if( implicit%clscnt > 0 ) then
  allocate( implicit%permute(implicit%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate implicit%permute ; error = ',astat
     call endrun
  end if
         implicit%permute(:) = 0
  allocate( implicit%diag_map(implicit%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate implicit%diag_map ; error = ',astat
     call endrun
  end if
         implicit%diag_map(:) = 0
  allocate( implicit%clsmap(implicit%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate implicit%clsmap ; error = ',astat
     call endrun
  end if
         implicit%clsmap(:) = 0
      end if
      if( rodas%clscnt > 0 ) then
  allocate( rodas%permute(rodas%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate rodas%permute ; error = ',astat
     call endrun
  end if
         rodas%permute(:) = 0
  allocate( rodas%diag_map(rodas%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate rodas%diag_map ; error = ',astat
     call endrun
  end if
         rodas%diag_map(:) = 0
  allocate( rodas%clsmap(rodas%clscnt),stat=astat )
  if( astat /= 0 ) then
     write(*,*) 'chem_mods_inti: failed to allocate rodas%clsmap ; error = ',astat
     call endrun
  end if
         rodas%clsmap(:) = 0
      end if

      end subroutine chem_mods_init

      end module chem_mods_mod

      module M_SPC_ID_MOD

      implicit none

      integer, parameter :: id_O3 = 1
      integer, parameter :: id_O = 2
      integer, parameter :: id_O1D = 3
      integer, parameter :: id_N2O = 4
      integer, parameter :: id_N = 5
      integer, parameter :: id_NO = 6
      integer, parameter :: id_NO2 = 7
      integer, parameter :: id_NO3 = 8
      integer, parameter :: id_HNO3 = 9
      integer, parameter :: id_HO2NO2 = 10
      integer, parameter :: id_N2O5 = 11
      integer, parameter :: id_CH4 = 12
      integer, parameter :: id_CH3O2 = 13
      integer, parameter :: id_CH3OOH = 14
      integer, parameter :: id_CH2O = 15
      integer, parameter :: id_CO = 16
      integer, parameter :: id_OH = 17
      integer, parameter :: id_HO2 = 18
      integer, parameter :: id_H2O2 = 19
      integer, parameter :: id_C3H6 = 20
      integer, parameter :: id_ISOP = 21
      integer, parameter :: id_PO2 = 22
      integer, parameter :: id_CH3CHO = 23
      integer, parameter :: id_POOH = 24
      integer, parameter :: id_CH3CO3 = 25
      integer, parameter :: id_CH3COOOH = 26
      integer, parameter :: id_PAN = 27
      integer, parameter :: id_C2H6 = 28
      integer, parameter :: id_C2H4 = 29
      integer, parameter :: id_C4H10 = 30
      integer, parameter :: id_MPAN = 31
      integer, parameter :: id_ISOPO2 = 32
      integer, parameter :: id_MVK = 33
      integer, parameter :: id_MACR = 34
      integer, parameter :: id_MACRO2 = 35
      integer, parameter :: id_MACROOH = 36
      integer, parameter :: id_C2H5O2 = 37
      integer, parameter :: id_C2H5OOH = 38
      integer, parameter :: id_C10H16 = 39
      integer, parameter :: id_C3H8 = 40
      integer, parameter :: id_C3H7O2 = 41
      integer, parameter :: id_C3H7OOH = 42
      integer, parameter :: id_CH3COCH3 = 43
      integer, parameter :: id_CH3OH = 44
      integer, parameter :: id_C2H5OH = 45
      integer, parameter :: id_GLYALD = 46
      integer, parameter :: id_HYAC = 47
      integer, parameter :: id_EO2 = 48
      integer, parameter :: id_EO = 49
      integer, parameter :: id_ISOPOOH = 50
      integer, parameter :: id_H2 = 51
      integer, parameter :: id_SO2 = 52
      integer, parameter :: id_SO4 = 53
      integer, parameter :: id_DMS = 54
      integer, parameter :: id_NH3 = 55
      integer, parameter :: id_NH4NO3 = 56
      integer, parameter :: id_NH4 = 57
      integer, parameter :: id_HCl = 58
      integer, parameter :: id_HOCl = 59
      integer, parameter :: id_ClONO2 = 60
      integer, parameter :: id_Cl = 61
      integer, parameter :: id_ClO = 62
      integer, parameter :: id_Cl2O2 = 63
      integer, parameter :: id_Cl2 = 64
      integer, parameter :: id_HOBr = 65
      integer, parameter :: id_HBr = 66
      integer, parameter :: id_BrONO2 = 67
      integer, parameter :: id_Br = 68
      integer, parameter :: id_BrO = 69
      integer, parameter :: id_BrCl = 70
      integer, parameter :: id_LCH4 = 71
      integer, parameter :: id_H = 72
      integer, parameter :: id_H2O = 73
      integer, parameter :: id_ROH = 74
      integer, parameter :: id_RCHO = 75
      integer, parameter :: id_ISOPNB = 76
      integer, parameter :: id_ISOPNBO2 = 77
      integer, parameter :: id_MACRN = 78
      integer, parameter :: id_MVKN = 79
      integer, parameter :: id_R4N2 = 80
      integer, parameter :: id_MEK = 81
      integer, parameter :: id_R4N1 = 82
      integer, parameter :: id_IEPOX = 83
      integer, parameter :: id_IEPOXOO = 84
      integer, parameter :: id_GLYX = 85
      integer, parameter :: id_MGLY = 86
      integer, parameter :: id_MVKO2 = 87
      integer, parameter :: id_MVKOOH = 88
      integer, parameter :: id_MACRNO2 = 89
      integer, parameter :: id_MAO3 = 90
      integer, parameter :: id_MAOP = 91
      integer, parameter :: id_MAOPO2 = 92
      integer, parameter :: id_ATO2 = 93
      integer, parameter :: id_ATOOH = 94
      integer, parameter :: id_INO2 = 95
      integer, parameter :: id_INPN = 96
      integer, parameter :: id_ISNOOA = 97
      integer, parameter :: id_ISN1 = 98
      integer, parameter :: id_O3S = 99


      end module M_SPC_ID_MOD

      module M_RXT_ID_MOD

      implicit none

      integer, parameter :: rid_jo2 = 1
      integer, parameter :: rid_jo1d = 2
      integer, parameter :: rid_jo3p = 3
      integer, parameter :: rid_jn2o = 4
      integer, parameter :: rid_jno = 5
      integer, parameter :: rid_jno2 = 6
      integer, parameter :: rid_jn2o5 = 7
      integer, parameter :: rid_jhno3 = 8
      integer, parameter :: rid_jno3 = 9
      integer, parameter :: rid_jho2no2 = 10
      integer, parameter :: rid_jch3ooh = 11
      integer, parameter :: rid_jch2o_a = 12
      integer, parameter :: rid_jch2o_b = 13
      integer, parameter :: rid_jh2o = 14
      integer, parameter :: rid_jh2o2 = 15
      integer, parameter :: rid_jch3cho = 16
      integer, parameter :: rid_jpooh = 17
      integer, parameter :: rid_jch3co3h = 18
      integer, parameter :: rid_jpan = 19
      integer, parameter :: rid_jmpan = 20
      integer, parameter :: rid_jmacr_a = 21
      integer, parameter :: rid_jmacr_b = 22
      integer, parameter :: rid_jmvk = 23
      integer, parameter :: rid_jc2h5ooh = 24
      integer, parameter :: rid_jc3h7ooh = 25
      integer, parameter :: rid_jacet = 26
      integer, parameter :: rid_jmgly = 27
      integer, parameter :: rid_jglyoxal = 28
      integer, parameter :: rid_jisopooh = 29
      integer, parameter :: rid_jhyac = 30
      integer, parameter :: rid_jglyald = 31
      integer, parameter :: rid_jisopnb = 32
      integer, parameter :: rid_jmacrn = 33
      integer, parameter :: rid_jmvkn = 34
      integer, parameter :: rid_jr4n2 = 35
      integer, parameter :: rid_jclono2 = 36
      integer, parameter :: rid_jhocl = 37
      integer, parameter :: rid_jcl2o2 = 38
      integer, parameter :: rid_jbrono2 = 39
      integer, parameter :: rid_jhobr = 40
      integer, parameter :: rid_jbrcl = 41
      integer, parameter :: rid_jbro = 42
      integer, parameter :: rid_jcl2 = 43
      integer, parameter :: rid_uo_o2 = 44
      integer, parameter :: rid_uco_oha = 48
      integer, parameter :: rid_uco_ohb = 49
      integer, parameter :: rid_ol_oh = 53
      integer, parameter :: rid_ol_ho2 = 54
      integer, parameter :: rid_uho2_ho2 = 55
      integer, parameter :: rid_o1d_n2 = 60
      integer, parameter :: rid_o1d_o2 = 61
      integer, parameter :: rid_ol_o1d = 62
      integer, parameter :: rid_op_ho2 = 65
      integer, parameter :: rid_uno2_no3 = 70
      integer, parameter :: rid_un2o5 = 71
      integer, parameter :: rid_uoh_no2 = 73
      integer, parameter :: rid_uoh_hno3 = 74
      integer, parameter :: rid_uho2_no2 = 76
      integer, parameter :: rid_uhno4 = 78
      integer, parameter :: rid_op_mo2 = 81
      integer, parameter :: rid_uoh_c2h4 = 88
      integer, parameter :: rid_op_eo2 = 89
      integer, parameter :: rid_ol_c2h4 = 92
      integer, parameter :: rid_uoh_c3h6 = 93
      integer, parameter :: rid_ol_c3h6 = 94
      integer, parameter :: rid_op_po2 = 96
      integer, parameter :: rid_op_ch3co3 = 101
      integer, parameter :: rid_upan_f = 102
      integer, parameter :: rid_upan_b = 106
      integer, parameter :: rid_op_c2h5o2 = 109
      integer, parameter :: rid_op_c3h7o2 = 116
      integer, parameter :: rid_uoh_acet = 120
      integer, parameter :: rid_ol_isop = 125
      integer, parameter :: rid_op_isopo2 = 126
      integer, parameter :: rid_op_isopnbo2 = 131
      integer, parameter :: rid_ol_isopnb = 133
      integer, parameter :: rid_op_iepoxo2 = 138
      integer, parameter :: rid_ol_mvk = 140
      integer, parameter :: rid_op_mvko2 = 141
      integer, parameter :: rid_ol_macr = 148
      integer, parameter :: rid_op_macro2 = 150
      integer, parameter :: rid_op_macrno2 = 156
      integer, parameter :: rid_op_mao3 = 158
      integer, parameter :: rid_umpan_f = 162
      integer, parameter :: rid_umpan_b = 163
      integer, parameter :: rid_op_maopo2 = 168
      integer, parameter :: rid_op_ato2 = 170
      integer, parameter :: rid_op_ino2 = 182
      integer, parameter :: rid_ol_c10h16 = 195
      integer, parameter :: rid_n2o5h = 197
      integer, parameter :: rid_no3h = 198
      integer, parameter :: rid_ho2h = 199
      integer, parameter :: rid_no2h = 200
      integer, parameter :: rid_uoh_dms = 203
      integer, parameter :: rid_nh3h = 205
      integer, parameter :: rid_strat13 = 207
      integer, parameter :: rid_strat14 = 208
      integer, parameter :: rid_strat20 = 209
      integer, parameter :: rid_strat21 = 210
      integer, parameter :: rid_strat22 = 211
      integer, parameter :: rid_strat23 = 212
      integer, parameter :: rid_strat24 = 213
      integer, parameter :: rid_strat25 = 214
      integer, parameter :: rid_strat26 = 215
      integer, parameter :: rid_strat27 = 216
      integer, parameter :: rid_strat28 = 217
      integer, parameter :: rid_strat29 = 218
      integer, parameter :: rid_strat33 = 219
      integer, parameter :: rid_strat35 = 220
      integer, parameter :: rid_strat37 = 221
      integer, parameter :: rid_strat38 = 222
      integer, parameter :: rid_strat39 = 223
      integer, parameter :: rid_strat40 = 224
      integer, parameter :: rid_strat41 = 225
      integer, parameter :: rid_strat42 = 226
      integer, parameter :: rid_strat43 = 227
      integer, parameter :: rid_strat44 = 228
      integer, parameter :: rid_strat45 = 229
      integer, parameter :: rid_strat46 = 230
      integer, parameter :: rid_strat47 = 231
      integer, parameter :: rid_strat48 = 232
      integer, parameter :: rid_strat69 = 233
      integer, parameter :: rid_strat58 = 234
      integer, parameter :: rid_strat59 = 235
      integer, parameter :: rid_strat64 = 236
      integer, parameter :: rid_strat71 = 237
      integer, parameter :: rid_strat72 = 238
      integer, parameter :: rid_strat73 = 239
      integer, parameter :: rid_strat74 = 240
      integer, parameter :: rid_strat75 = 241
      integer, parameter :: rid_strat76 = 242
      integer, parameter :: rid_strat77 = 243
      integer, parameter :: rid_strat78 = 244
      integer, parameter :: rid_strat79 = 245
      integer, parameter :: rid_strat80 = 246

      integer, parameter :: rid_r0045 = 45
      integer, parameter :: rid_r0046 = 46
      integer, parameter :: rid_r0047 = 47
      integer, parameter :: rid_r0050 = 50
      integer, parameter :: rid_r0051 = 51
      integer, parameter :: rid_r0052 = 52
      integer, parameter :: rid_r0056 = 56
      integer, parameter :: rid_r0057 = 57
      integer, parameter :: rid_r0058 = 58
      integer, parameter :: rid_r0059 = 59
      integer, parameter :: rid_r0063 = 63
      integer, parameter :: rid_r0064 = 64
      integer, parameter :: rid_r0066 = 66
      integer, parameter :: rid_r0067 = 67
      integer, parameter :: rid_r0068 = 68
      integer, parameter :: rid_r0069 = 69
      integer, parameter :: rid_r0072 = 72
      integer, parameter :: rid_r0075 = 75
      integer, parameter :: rid_r0077 = 77
      integer, parameter :: rid_r0079 = 79
      integer, parameter :: rid_r0080 = 80
      integer, parameter :: rid_r0082 = 82
      integer, parameter :: rid_r0083 = 83
      integer, parameter :: rid_r0084 = 84
      integer, parameter :: rid_r0085 = 85
      integer, parameter :: rid_r0086 = 86
      integer, parameter :: rid_r0087 = 87
      integer, parameter :: rid_r0090 = 90
      integer, parameter :: rid_r0091 = 91
      integer, parameter :: rid_r0095 = 95
      integer, parameter :: rid_r0097 = 97
      integer, parameter :: rid_r0098 = 98
      integer, parameter :: rid_r0099 = 99
      integer, parameter :: rid_r0100 = 100
      integer, parameter :: rid_r0103 = 103
      integer, parameter :: rid_r0104 = 104
      integer, parameter :: rid_r0105 = 105
      integer, parameter :: rid_r0107 = 107
      integer, parameter :: rid_r0108 = 108
      integer, parameter :: rid_r0110 = 110
      integer, parameter :: rid_r0111 = 111
      integer, parameter :: rid_r0112 = 112
      integer, parameter :: rid_r0113 = 113
      integer, parameter :: rid_r0114 = 114
      integer, parameter :: rid_r0115 = 115
      integer, parameter :: rid_r0117 = 117
      integer, parameter :: rid_r0118 = 118
      integer, parameter :: rid_r0119 = 119
      integer, parameter :: rid_r0121 = 121
      integer, parameter :: rid_r0122 = 122
      integer, parameter :: rid_r0123 = 123
      integer, parameter :: rid_r0124 = 124
      integer, parameter :: rid_r0127 = 127
      integer, parameter :: rid_r0128 = 128
      integer, parameter :: rid_r0129 = 129
      integer, parameter :: rid_r0130 = 130
      integer, parameter :: rid_r0132 = 132
      integer, parameter :: rid_r0134 = 134
      integer, parameter :: rid_r0135 = 135
      integer, parameter :: rid_r0136 = 136
      integer, parameter :: rid_r0137 = 137
      integer, parameter :: rid_r0139 = 139
      integer, parameter :: rid_r0142 = 142
      integer, parameter :: rid_r0143 = 143
      integer, parameter :: rid_r0144 = 144
      integer, parameter :: rid_r0145 = 145
      integer, parameter :: rid_r0146 = 146
      integer, parameter :: rid_r0147 = 147
      integer, parameter :: rid_r0149 = 149
      integer, parameter :: rid_r0151 = 151
      integer, parameter :: rid_r0152 = 152
      integer, parameter :: rid_r0153 = 153
      integer, parameter :: rid_r0154 = 154
      integer, parameter :: rid_r0155 = 155
      integer, parameter :: rid_r0157 = 157
      integer, parameter :: rid_r0159 = 159
      integer, parameter :: rid_r0160 = 160
      integer, parameter :: rid_r0161 = 161
      integer, parameter :: rid_r0164 = 164
      integer, parameter :: rid_r0165 = 165
      integer, parameter :: rid_r0166 = 166
      integer, parameter :: rid_r0167 = 167
      integer, parameter :: rid_r0169 = 169
      integer, parameter :: rid_r0171 = 171
      integer, parameter :: rid_r0172 = 172
      integer, parameter :: rid_r0173 = 173
      integer, parameter :: rid_r0174 = 174
      integer, parameter :: rid_r0175 = 175
      integer, parameter :: rid_r0176 = 176
      integer, parameter :: rid_r0177 = 177
      integer, parameter :: rid_r0178 = 178
      integer, parameter :: rid_r0179 = 179
      integer, parameter :: rid_r0180 = 180
      integer, parameter :: rid_r0181 = 181
      integer, parameter :: rid_r0183 = 183
      integer, parameter :: rid_r0184 = 184
      integer, parameter :: rid_r0185 = 185
      integer, parameter :: rid_r0186 = 186
      integer, parameter :: rid_r0187 = 187
      integer, parameter :: rid_r0188 = 188
      integer, parameter :: rid_r0189 = 189
      integer, parameter :: rid_r0190 = 190
      integer, parameter :: rid_r0191 = 191
      integer, parameter :: rid_r0192 = 192
      integer, parameter :: rid_r0193 = 193
      integer, parameter :: rid_r0194 = 194
      integer, parameter :: rid_r0196 = 196
      integer, parameter :: rid_r0201 = 201
      integer, parameter :: rid_r0202 = 202
      integer, parameter :: rid_r0204 = 204
      integer, parameter :: rid_r0206 = 206
      integer, parameter :: rid_r0247 = 247

      end module M_RXT_ID_MOD

      module M_HET_ID_MOD

      implicit none


      end module M_HET_ID_MOD
