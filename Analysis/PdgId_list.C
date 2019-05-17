
#define ntupleClass_MC_cxx
#define NPARTICLES 560

void ntupleClass_MC::Fill_particleName(TString pId[NPARTICLES]){
    // Given a vector of strings the function initializes it with the names of particles
    
    //Particle name list
    //QUARKS
    pId[0] = "d";
    pId[1] = "u";
    pId[2] = "s";
    pId[3] = "c";
    pId[4] = "b";
    pId[5] = "t";
    //LEPTONS
    pId[6] = "e^{-}";
    pId[7] = "e^{+}";
    pId[8] = "#nu_{e}";
    pId[9] = "#bar{#nu}_{e}";
    pId[10] = "#mu^{-}";
    pId[11] = "#mu^{+}";
    pId[12] = "#nu_{#mu}";
    pId[13] = "#bar{#nu}_{#mu}";
    pId[14] = "#tau^{-}";
    pId[15] = "#tau^{+}";
    pId[16] = "#nu_{#tau}";
    pId[17] = "#bar{#nu}_{#tau}";
    //GAUGE BOSONS
    pId[18] = "g";
    pId[19] = "#gamma";
    pId[20] = "Z^{0}";
    pId[21] = "W^{+}";
    pId[22] = "W^{-}";
    //DIQUARKS
    pId[23] = "(dd)_{1}";
    pId[24] = "(ud)_{0}";
    pId[25] = "(ud)_{1}";
    pId[26] = "(uu)_{1}";
    pId[27] = "(sd)_{0}";
    pId[28] = "(sd)_{1}";
    pId[29] = "(su)_{0}";
    pId[30] = "(su)_{1}";
    pId[31] = "(ss)_{1}";
    pId[32] = "(cd)_{0}";
    pId[33] = "(cd)_{1}";
    pId[34] = "(cu)_{0}";
    pId[35] = "(cu)_{1}";
    pId[36] = "(cs)_{0}";
    pId[37] = "(cs)_{1}";
    pId[38] = "(cc)_{1}";
    pId[39] = "(bd)_{0}";
    pId[40] = "(bd)_{1}";
    pId[41] = "(bu)_{0}";
    pId[42] = "(bu)_{1}";
    pId[43] = "(bs)_{0}";
    pId[44] = "(bs)_{1}";
    pId[45] = "(bc)_{0}";
    pId[46] = "(bc)_{1}";
    pId[47] = "(bb)_{1}";
    //LIGHT I=1 MESONS
    pId[48] = "#Pi^{0}";
    pId[49] = "#bar{#Pi}^{0}";
    pId[50] = "#Pi^{+}";
    pId[51] = "#Pi^{-}";
    pId[52] = "a_{0}(980)^{0}";
    pId[53] = "#bar{a}_{0}(980)^{0}";
    pId[54] = "a_{0}(980)^{+}";
    pId[55] = "a_{0}(980)^{-}";
    pId[56] = "#Pi(1300)^{0}";
    pId[57] = "#bar{#Pi}(1300)^{0}";
    pId[58] = "#Pi(1300)^{+}";
    pId[59] = "#Pi(1300)^{-}";
    pId[60] = "a_{0}(1450)^{0}";
    pId[61] = "#bar{a}_{0}(1450)^{0}";
    pId[62] = "a_{0}(1450)^{+}";
    pId[63] = "a_{0}(1450)^{-}";
    pId[64] = "#Pi(1800)^{0}";
    pId[65] = "#bar{#Pi}(1800)^{0}";
    pId[66] = "#Pi(1800)^{+}";
    pId[67] = "#Pi(1800)^{-}";
    pId[68] = "#rho(770)^{0}";
    pId[69] = "#bar{#rho}(770)^{0}";
    pId[70] = "#rho(770)^{+}";
    pId[71] = "#rho(770)^{-}";
    pId[72] = "b_{1}(1235)^{0}";
    pId[73] = "#bar{b}_{1}(1235)^{0}";
    pId[74] = "b_{1}(1235)^{+}";
    pId[75] = "b_{1}(1235)^{-}";
    pId[76] = "a_{1}(1260)^{0}";
    pId[77] = "#bar{a}_{1}(1260)^{0}";
    pId[78] = "a_{1}(1260)^{+}";
    pId[79] = "a_{1}(1260)^{-}";
    pId[80] = "#Pi_{1}(1400)^{0}";
    pId[81] = "#bar{#Pi}_{1}(1400)^{0}";
    pId[82] = "#Pi_{1}(1400)^{+}";
    pId[83] = "#Pi_{1}(1400)^{-}";
    pId[84] = "#rho(1450)^{0}";
    pId[85] = "#bar{#rho}(1450)^{0}";
    pId[86] = "#rho(1450)^{+}";
    pId[87] = "#rho(1450)^{-}";
    pId[88] = "#Pi_{1}(1600)^{0}";
    pId[89] = "#bar{#Pi}_{1}(1600)^{0}";
    pId[90] = "#Pi_{1}(1600)^{+}";
    pId[91] = "#Pi_{1}(1600)^{-}";
    pId[92] = "a_{1}(1640)^{0}";
    pId[93] = "#bar{a}_{1}(1640)^{0}";
    pId[94] = "a_{1}(1640)^{+}";
    pId[95] = "a_{1}(1640)^{-}";
    pId[96] = "#rho(1700)^{0}";
    pId[97] = "#bar{#rho}(1700)^{0}";
    pId[98] = "#rho(1700)^{+}";
    pId[99] = "#rho(1700)^{-}";
    pId[100] = "#rho(1900)^{0}";
    pId[101] = "#bar{#rho}(1900)^{0}";
    pId[102] = "#rho(1900)^{+}";
    pId[103] = "#rho(1900)^{-}";
    pId[104] = "#rho(2150)^{0}";
    pId[105] = "#bar{#rho}(2150)^{0}";
    pId[106] = "#rho(2150)^{+}";
    pId[107] = "#rho(2150)^{-}";
    pId[108] = "a_{2}(1320)^{0}";
    pId[109] = "#bar{a}_{2}(1320)^{0}";
    pId[110] = "a_{2}(1320)^{+}";
    pId[111] = "a_{2}(1320)^{-}";
    pId[112] = "#Pi_{2}(1670)^{0}";
    pId[113] = "#bar{#Pi}_{2}(1670)^{0}";
    pId[114] = "#Pi_{2}(1670)^{+}";
    pId[115] = "#Pi_{2}(1670)^{-}";
    pId[116] = "a_{2}(1700)^{0}";
    pId[117] = "#bar{a}_{2}(1700)^{0}";
    pId[118] = "a_{2}(1700)^{+}";
    pId[119] = "a_{2}(1700)^{-}";
    pId[120] = "#Pi_{2}(2100)^{0}";
    pId[121] = "#bar{#Pi}_{2}(2100)^{0}";
    pId[122] = "#Pi_{2}(2100)^{+}";
    pId[123] = "#Pi_{2}(2100)^{-}";
    pId[124] = "#rho_{3}(1690)^{0}";
    pId[125] = "#bar{#rho}_{3}(1690)^{0}";
    pId[126] = "#rho_{3}(1690)^{+}";
    pId[127] = "#rho_{3}(1690)^{-}";
    pId[128] = "#rho_{3}(1990)^{0}";
    pId[129] = "#bar{#rho}_{3}(1990)^{0}";
    pId[130] = "#rho_{3}(1990)^{+}";
    pId[131] = "#rho_{3}(1990)^{-}";
    pId[132] = "#rho_{3}(2250)^{0}";
    pId[133] = "#bar{#rho}_{3}(12250)^{0}";
    pId[134] = "#rho_{3}(2250)^{+}";
    pId[135] = "#rho_{3}(2250)^{-}";
    pId[136] = "a_{4}(2040)^{0}";
    pId[137] = "#bar{a}_{4}(2040)^{0}";
    pId[138] = "a_{4}(2040)^{+}";
    pId[139] = "a_{4}(2040)^{-}";
    // LIGHT I=0 MESONS
    // (u\bar{u}, d\bar{d}, and s\bar{s} Admixtures)
    pId[140] = "#eta";
    pId[141] = "#eta'(958)";
    pId[142] = "f_{0}(600)";
    pId[143] = "f_{0}(980)";
    pId[144] = "#eta(1295)";
    pId[145] = "f_{0}(1370)";
    pId[146] = "#eta(1440)";
    pId[147] = "f_{0}(1500)";
    pId[148] = "f_{0}(1710)";
    pId[149] = "#eta(1760)";
    pId[150] = "f_{0}(2020)";
    pId[151] = "f_{0}(2100)";
    pId[152] = "f_{0}(2200)";
    pId[153] = "#eta(2225)";
    pId[154] = "f_{0}(2330)";
    pId[155] = "#omega(782)";
    pId[156] = "#Phi(1020)";
    pId[157] = "h_{1}(1170)";
    pId[158] = "f_{1}(1285)";
    pId[159] = "h_{1}(1380)";
    pId[160] = "f_{1}(1420)";
    pId[161] = "#omega(1420)";
    pId[162] = "f_{1}(1510)";
    pId[163] = "h_{1}(1595)";
    pId[164] = "#omega(1650)";
    pId[165] = "#Phi(1680)";
    pId[166] = "f_{2}(1270)";
    pId[167] = "f_{2}(1430)";
    pId[168] = "f'_{2}(1525)";
    pId[169] = "f_{2}(1565)";
    pId[170] = "f_{2}(1640)";
    pId[171] = "#eta_{2}(1645)";
    pId[172] = "f_{2}(1810)";
    pId[173] = "#eta_{2}(1870)";
    pId[174] = "f_{2}(1910)";
    pId[175] = "f_{2}(1950)";
    pId[176] = "f_{2}(2010)";
    pId[177] = "f_{2}(2150)";
    pId[178] = "f_{2}(2300)";
    pId[179] = "f_{2}(2340)";
    pId[180] = "#omega_{3}(1670)";
    pId[181] = "#Phi_{3}(1850)";
    pId[182] = "f_{4}(2050)";
    pId[183] = "f_{j}(2220)";
    pId[184] = "f_{4}(2300)";
    // STRANGE MESONS
    pId[185] = "K_{L}^{0}";
    pId[186] = "#bar{K}_{L}^{0}";
    pId[187] = "K_{S}^{0}";
    pId[188] = "#bar{K}_{S}^{0}";
    pId[189] = "K^{0}";
    pId[190] = "#bar{K}^{0}";
    pId[191] = "K^{+}";
    pId[192] = "K^{-}";
    pId[193] = "K_{0}^{*}(1430)^{0}";
    pId[194] = "#bar{K}_{0}^{*}(1430)^{0}";
    pId[195] = "K_{0}^{*}(1430)^{+}";
    pId[196] = "K_{0}^{*}(1430)^{-}";
    pId[197] = "K(1460)^{0}";
    pId[198] = "#bar{K}(1460)^{0}";
    pId[199] = "K(1460)^{+}";
    pId[200] = "K(1460)^{-}";
    pId[201] = "K(1830)^{0}";
    pId[202] = "#bar{K}(1830)^{0}";
    pId[203] = "K(1830)^{+}";
    pId[204] = "K(1830)^{-}";
    pId[205] = "K_{0}^{*}(1950)^{0}";
    pId[206] = "#bar{K}_{0}^{*}(1950)^{0}";
    pId[207] = "K_{0}^{*}(1950)^{+}";
    pId[208] = "K_{0}^{*}(1950)^{-}";
    pId[209] = "K^{*}(892)^{0}";
    pId[210] = "#bar{K}^{*}(892)^{0}";
    pId[211] = "K^{*}(892)^{+}";
    pId[212] = "K^{*}(892)^{-}";
    pId[213] = "K_{1}(1270)^{0}";
    pId[214] = "#bar{K}_{1}(1270)^{0}";
    pId[215] = "K_{1}(1270)^{+}";
    pId[216] = "K_{1}(1270)^{-}";
    pId[217] = "K_{1}(1400)^{0}";
    pId[218] = "#bar{K}_{1}(1400)^{0}";
    pId[219] = "K_{1}(1400)^{+}";
    pId[220] = "K_{1}(1400)^{-}";
    pId[221] = "K^{*}(1410)^{0}";
    pId[222] = "#bar{K}^{*}(1410)^{0}";
    pId[223] = "K^{*}(1410)^{+}";
    pId[224] = "K^{*}(1410)^{-}";
    pId[225] = "K_{1}(1650)^{0}";
    pId[226] = "#bar{K}_{1}(1650)^{0}";
    pId[227] = "K_{1}(1650)^{+}";
    pId[228] = "K_{1}(1650)^{-}";
    pId[229] = "K^{*}(1680)^{0}";
    pId[230] = "#bar{K}^{*}(1680)^{0}";
    pId[231] = "K^{*}(1680)^{+}";
    pId[232] = "K^{*}(1680)^{-}";
    pId[233] = "K_{2}^{*}(1430)^{0}";
    pId[234] = "#bar{K}_{2}^{*}(1430)^{0}";
    pId[235] = "K_{2}^{*}(1430)^{+}";
    pId[236] = "K_{2}^{*}(1430)^{-}";
    pId[237] = "K_{2}(1580)^{0}";
    pId[238] = "#bar{K}_{2}(1580)^{0}";
    pId[239] = "K_{2}(1580)^{+}";
    pId[240] = "K_{2}(1580)^{-}";
    pId[241] = "K_{2}(1770)^{0}";
    pId[242] = "#bar{K}_{2}(1770)^{0}";
    pId[243] = "K_{2}(1770)^{+}";
    pId[244] = "K_{2}(1770)^{-}";
    pId[245] = "K_{2}(1820)^{0}";
    pId[246] = "#bar{K}_{2}(1820)^{0}";
    pId[247] = "K_{2}(1820)^{+}";
    pId[248] = "K_{2}(1820)^{-}";
    pId[249] = "K_{2}^{∗}(1980)^{0}";
    pId[250] = "#bar{K}_{2}^{∗}(1980)^{0}";
    pId[251] = "K_{2}^{∗}(1980)^{+}";
    pId[252] = "K_{2}^{∗}(1980)^{-}";
    pId[253] = "K_{2}(2250)^{0}";
    pId[254] = "#bar{K}_{2}(2250)^{0}";
    pId[255] = "K_{2}(2250)^{+}";
    pId[256] = "K_{2}(2250)^{-}";
    pId[257] = "K_{3}^{∗}(1780)^{0}";
    pId[258] = "#bar{K}_{3}^{∗}(1780)^{0}";
    pId[259] = "K_{3}^{∗}(1780)^{+}";
    pId[260] = "K_{3}^{∗}(1780)^{-}";
    pId[261] = "K_{3}(2320)^{0}";
    pId[262] = "#bar{K}_{3}(2320)^{0}";
    pId[263] = "K_{3}(2320)^{+}";
    pId[264] = "K_{3}(2320)^{-}";
    pId[265] = "K_{4}^{∗}(2045)^{0}";
    pId[266] = "#bar{K}_{4}^{∗}(2045)^{0}";
    pId[267] = "K_{4}^{∗}(2045)^{+}";
    pId[268] = "K_{4}^{∗}(2045)^{-}";
    pId[269] = "K_{4}(2500)^{0}";
    pId[270] = "#bar{K}_{4}(2500)^{0}";
    pId[271] = "K_{4}(2500)^{+}";
    pId[272] = "K_{4}(2500)^{-}";
    // CHARMED MESONS
    pId[273] = "D^{+}";
    pId[274] = "D^{-}";
    pId[275] = "D^{0}";
    pId[276] = "#bar{D}^{0}";
    pId[277] = "D_{0}^{*+}";
    pId[278] = "D_{0}^{*-}";
    pId[279] = "D_{0}^{*0}";
    pId[280] = "#bar{D}_{0}^{*0}";
    pId[281] = "D^{*}(2010)^{+}";
    pId[282] = "D^{*}(2010)^{-}";
    pId[283] = "D^{*}(2007)^{0}";
    pId[284] = "#bar{D}^{*}(2007)^{0}";
    pId[285] = "D_{1}(2420)^{+}";
    pId[286] = "D_{1}(2420)^{-}";
    pId[287] = "D_{1}(2420)^{0}";
    pId[288] = "#bar{D}_{1}(2420)^{0}";
    pId[289] = "D_{1}(H)^{+}";
    pId[290] = "D_{1}(H)^{-}";
    pId[291] = "D_{1}(H)^{0}";
    pId[292] = "#bar{D}_{1}(H)^{0}";
    pId[293] = "D^{*}_{2}(2460)^{+}";
    pId[294] = "D^{*}_{2}(2460)^{-}";
    pId[295] = "D^{*}_{2}(2460)^{0}";
    pId[296] = "#bar{D}^{*}_{2}(2460)^{0}";
    pId[297] = "D_{s}^{+}";
    pId[298] = "D_{s}^{-}";
    pId[299] = "D_{s0}^{*+}";
    pId[300] = "D_{s0}^{*-}";
    pId[301] = "D_{s}^{*+}";
    pId[302] = "D_{s}^{*-}";
    pId[303] = "D_{s1}(2536)^{+}";
    pId[304] = "D_{s1}(2536)^{-}";
    pId[305] = "D_{s1}(H)^{+}";
    pId[306] = "D_{s1}(H)^{-}";
    pId[307] = "D_{s2}^{*+}";
    pId[308] = "D_{s2}^{*-}";
    // BOTTOM MESONS
    pId[309] = "B^{0}";
    pId[310] = "#bar{B}^{0}";
    pId[311] = "B^{+}";
    pId[312] = "B^{-}";
    pId[313] = "B_{0}^{*0}";
    pId[314] = "#bar{B}_{0}^{*0}";
    pId[315] = "B_{0}^{*+}";
    pId[316] = "B_{0}^{*-}";
    pId[317] = "B^{*0}";
    pId[318] = "#bar{B}^{*0}";
    pId[319] = "B^{*+}";
    pId[320] = "B^{*-}";
    pId[321] = "B_{1}(L)^{0}";
    pId[322] = "#bar{B}_{1}(L)^{0}";
    pId[323] = "B_{1}(L)^{+}";
    pId[324] = "B_{1}(L)^{-}";
    pId[325] = "B_{1}(H)^{0}";
    pId[326] = "#bar{B}_{1}(H)^{0}";
    pId[327] = "B_{1}(H)^{+}";
    pId[328] = "B_{1}(H)^{-}";
    pId[329] = "B_{2}^{*0}";
    pId[330] = "#bar{B}_{2}^{*0}}";
    pId[331] = "B_{2}^{*+}";
    pId[332] = "B_{2}^{*-}";
    pId[333] = "B_{s}^{0}";
    pId[334] = "#bar{B}_{s}^{0}";
    pId[335] = "B_{s0}^{*0}";
    pId[336] = "#bar{B}_{s0}^{*0}";
    pId[337] = "B_{s}^{*0}";
    pId[338] = "#bar{B}_{s}^{*0}";
    pId[339] = "B_{s1}(L)^{0}";
    pId[340] = "#bar{B}_{s1}(L)^{0}";
    pId[341] = "B_{s1}(H)^{0}";
    pId[342] = "#bar{B}_{s1}(H)^{0}";
    pId[343] = "B_{s2}^{*0}";
    pId[344] = "#bar{B}_{s2}^{*0}";
    pId[345] = "B_{c}^{+}";
    pId[346] = "B_{c}^{-}";
    pId[347] = "B_{c0}^{*+}";
    pId[348] = "B_{c0}^{*-}";
    pId[349] = "B_{c}^{*+}";
    pId[350] = "B_{c}^{*-}";
    pId[351] = "B_{c1}(L)^{+}";
    pId[352] = "B_{c1}(L)^{-}";
    pId[353] = "B_{c1}(H)^{+}";
    pId[354] = "B_{c1}(H)^{-}";
    pId[355] = "B_{c2}^{*+}";
    pId[356] = "B_{c2}^{*-}";
    // c\bar{c} MESONS
    pId[357] = "#eta_{c}(1S)";
    pId[358] = "#chi_{c0}(1P)";
    pId[359] = "#eta_{c}(2S)";
    pId[360] = "J/#psi(1S)";
    pId[361] = "h_{c}(1P)";
    pId[362] = "#chi_{c1}(1P)";
    pId[363] = "#psi(2S)";
    pId[364] = "#psi(3770)";
    pId[365] = "#psi(4040)";
    pId[366] = "#psi(4160)";
    pId[367] = "#psi(4415)";
    pId[368] = "#chi_{c2}(1P)";
    pId[369] = "#psi(3836)";
    // b\bar{b} MESONS
    pId[370] = "#eta_{b}(1S)";
    pId[371] = "#chi_{b0}(1P)";
    pId[372] = "#eta_{b}(2S)";
    pId[373] = "#chi_{b0}(2P)";
    pId[374] = "#eta_{b}(3S)";
    pId[375] = "#chi_{b0}(3P)";
    pId[376] = "#Upsilon(1S)";
    pId[377] = "h_{b}(1P)";
    pId[378] = "#chi_{b1}(1P)";
    pId[379] = "#Upsilon_{1}(1D)";
    pId[380] = "#Upsilon(2S)";
    pId[381] = "h_{b}(2P)";
    pId[382] = "#chi_{b1}(2P)";
    pId[383] = "#Upsilon_{1}(2D)";
    pId[384] = "#Upsilon(3S)";
    pId[385] = "h_{b}(3P)";
    pId[386] = "#chi_{b1}(3P)";
    pId[387] = "#Upsilon(4S)";
    pId[388] = "#Upsilon(10860)";
    pId[389] = "#Upsilon(11020)";
    pId[390] = "#chi_{b2}(1P)";
    pId[391] = "#eta_{b2}(1D)";
    pId[392] = "#Upsilon_{2}(1D)";
    pId[393] = "#chi_{b2}(2P)";
    pId[394] = "#eta_{b2}(2D)";
    pId[395] = "#Upsilon_{2}(2D)";
    pId[396] = "#chi_{b2}(3P)";
    pId[397] = "#Upsilon_{3}(1D)";
    pId[398] = "#Upsilon_{3}(2D)";
    // LIGHT BARYONS
    pId[399] = "p";
    pId[400] = "#bar{p}";
    pId[401] = "n";
    pId[402] = "#bar{n}";
    pId[403] = "#Delta^{++}";
    pId[404] = "#bar{#Delta}^{--}";
    pId[405] = "#Delta^{+}";
    pId[406] = "bar{#Delta}^{-}";
    pId[407] = "#Delta^{0}";
    pId[408] = "#bar{#Delta}^{0}";
    pId[409] = "#Delta^{-}";
    pId[410] = "#bar{#Delta}^{+}";
    // STRANGE BARYONS
    pId[411] = "#Lambda";
    pId[412] = "#bar{#Lambda}";
    pId[413] = "#Sigma^{+}";
    pId[414] = "#bar{#Sigma}^{-}";
    pId[415] = "#Sigma^{0}";
    pId[416] = "#bar{#Sigma}^{0}";
    pId[417] = "#Sigma^{-}";
    pId[418] = "#bar{#Sigma}^{+}";
    pId[419] = "#Sigma^{*+}";
    pId[420] = "#bar{#Sigma}^{*-}";
    pId[421] = "#Sigma^{*0}";
    pId[422] = "#bar{#Sigma}^{*0}";
    pId[423] = "#Sigma^{*-}";
    pId[424] = "#bar{#Sigma}^{*+}";
    pId[425] = "#Xi^{0}";
    pId[426] = "#bar{#Xi}^{0}";
    pId[427] = "#Xi^{-}";
    pId[428] = "#bar{#Xi}^{+}";
    pId[429] = "#Xi^{*0}";
    pId[430] = "#bar{#Xi}^{*0}";
    pId[431] = "#Xi^{*-}";
    pId[432] = "#bar{#Xi}^{*+}";
    pId[433] = "#Omega^{-}";
    pId[434] = "#bar{#Omega}^{-}";
    // CHARMED BARYONS
    pId[435] = "#Lambda_{c}^{+}";
    pId[436] = "#Lambda_{c}^{-}";
    pId[437] = "#Sigma_{c}^{++}";
    pId[438] = "#bar{#Sigma}_{c}^{--}";
    pId[439] = "#Sigma_{c}^{+}";
    pId[440] = "#Sigma_{c}^{-}";
    pId[441] = "#Sigma_{c}^{0}";
    pId[442] = "#bar{#Sigma}_{c}^{0}";
    pId[443] = "#Sigma_{c}^{*++}";
    pId[444] = "#bar{#Sigma}_{c}^{*--}";
    pId[445] = "#Sigma_{c}^{*+}";
    pId[446] = "#Sigma_{c}^{*-}";
    pId[447] = "#Sigma_{c}^{*0}";
    pId[448] = "#bar{#Sigma}_{c}^{*0}";
    pId[449] = "#Xi_{c}^{+}";
    pId[450] = "#Xi_{c}^{-}";
    pId[451] = "#Xi_{c}^{0}";
    pId[452] = "#bar{#Xi}_{c}^{+}";
    pId[453] = "#Xi'_{c}^{+}";
    pId[454] = "#Xi'_{c}^{-}";
    pId[455] = "#Xi'_{c}^{0}";
    pId[456] = "#bar{#Xi}'_{c}^{0}";
    pId[457] = "#Xi_{c}^{*+}";
    pId[458] = "#Xi_{c}^{*-}";
    pId[459] = "#Xi_{c}^{*0}";
    pId[460] = "#bar{#Xi}_{c}^{*0}";
    pId[461] = "#Omega_{c}^{0}";
    pId[462] = "#bar{#Omega}_{c}^{0}";
    pId[463] = "#Omega_{c}^{*0}";
    pId[464] = "#bar{#Omega}_{c}^{*0}";
    pId[465] = "#Xi_{cc}^{+}";
    pId[466] = "#Xi_{cc}^{-}";
    pId[467] = "#Xi_{cc}^{++}";
    pId[468] = "#Xi_{cc}^{--}";
    pId[469] = "#Xi_{cc}^{*+}";
    pId[470] = "#Xi_{cc}^{*-}";
    pId[471] = "#Xi_{cc}^{*++}";
    pId[472] = "#Xi_{cc}^{*--}";
    pId[473] = "#Omega_{cc}^{+}";
    pId[474] = "#Omega_{cc}^{-}";
    pId[475] = "#Omega_{cc}^{*+}";
    pId[476] = "Omega_{cc}^{*-}";
    pId[477] = "#Omega_{ccc}^{++}";
    pId[478] = "#Omega_{ccc}^{--}";
    // BOTTOM BARYONS
    pId[479] = "#Lambda_{b}^{0}";
    pId[480] = "#bar{#Lambda}_{b}^{0}";
    pId[481] = "#Sigma_{b}^{-}";
    pId[482] = "#bar{#Sigma}_{b}^{+}";
    pId[483] = "#Sigma_{b}^{0}";
    pId[484] = "#bar{#Sigma}_{b}^{0}";
    pId[485] = "#Sigma_{b}^{+}";
    pId[486] = "#bar{#Sigma}_{b}^{-}";
    pId[487] = "#Sigma_{b}^{*-}";
    pId[488] = "#bar{#Sigma}_{b}^{*+}";
    pId[489] = "#Sigma_{b}^{*0}";
    pId[490] = "#bar{#Sigma}_{b}^{*0}";
    pId[491] = "#Sigma_{b}^{*+}";
    pId[492] = "#bar{#Sigma}_{b}^{*-}";
    pId[493] = "#Xi_{b}^{-}";
    pId[494] = "#Xi_{b}^{+}";
    pId[495] = "#Xi_{b}^{0}";
    pId[496] = "#bar{#Xi}_{b}^{0}";
    pId[497] = "#Xi'_{b}^{-}";
    pId[498] = "#Xi'_{b}^{+}";
    pId[499] = "#Xi'_{b}^{0}";
    pId[500] = "#bar{#Xi}'_{b}^{0}";
    pId[501] = "#Xi_{b}^{*-}";
    pId[502] = "#Xi_{b}^{*+}";
    pId[503] = "#Xi_{b}^{*0}";
    pId[504] = "#bar{#Xi}_{b}^{*0}";
    pId[505] = "#Omega_{b}^{-}";
    pId[506] = "#Omega_{b}^{+}";
    pId[507] = "#Omega_{b}^{*-}";
    pId[508] = "#Omega_{b}^{*+}";
    pId[509] = "#Xi_{bc}^{0}";
    pId[510] = "#bar{#Xi}_{bc}^{0}";
    pId[511] = "#Xi_{bc}^{+}";
    pId[512] = "#Xi_{bc}^{-}";
    pId[513] = "#Xi'_{bc}^{0}";
    pId[514] = "#bar{#Xi}'_{bc}^{0}";
    pId[515] = "#Xi'_{bc}^{+}";
    pId[516] = "#Xi'_{bc}^{-}";
    pId[517] = "#Xi_{bc}^{*0}";
    pId[518] = "#bar{#Xi}_{bc}^{*0}";
    pId[519] = "#Xi_{bc}^{*+}";
    pId[520] = "#Xi_{bc}^{*-}";
    pId[521] = "#Omega_{bc}^{0}";
    pId[522] = "#bar{#Omega}_{bc}^{0}";
    pId[523] = "#Omega'_{bc}^{0}";
    pId[524] = "#bar{#Omega}'_{bc}^{0}";
    pId[525] = "#Omega_{bc}^{*0}";
    pId[526] = "#bar{#Omega}_{bc}^{*0}";
    pId[527] = "#Omega_{bcc}^{+}";
    pId[528] = "#Omega_{bcc}^{-}";
    pId[529] = "#Omega_{bcc}^{*+}";
    pId[530] = "#Omega_{bcc}^{*-}";
    pId[531] = "#Xi_{bb}^{-}";
    pId[532] = "#Xi_{bb}^{+}";
    pId[533] = "#Xi_{bb}^{0}";
    pId[534] = "#bar{#Xi}_{bb}^{0}";
    pId[535] = "#Xi_{bb}^{*-}";
    pId[536] = "#Xi_{bb}^{*+}";
    pId[537] = "#Xi_{bb}^{*0}";
    pId[538] = "#bar{#Xi}_{bb}^{*0}";
    pId[539] = "#Omega_{bb}^{-}";
    pId[540] = "#Omega_{bb}^{+}";
    pId[541] = "#Omega_{bb}^{*-}";
    pId[542] = "#Omega_{bb}^{*+}";
    pId[543] = "#Omega_{bbc}^{0}";
    pId[544] = "#bar{#Omega}_{bbc}^{0}";
    pId[545] = "#Omega_{bbc}^{*0}";
    pId[546] = "#bar{#Omega}_{bbc}^{*0}";
    pId[547] = "#Omega_{bbb}^{-}";
    pId[548] = "#Omega_{bbb}^{+}";
    pId[549] = "#Lambda_{c}(2593)^{+}";
    pId[550] = "#Lambda_{c}(2593)^{-}";
    pId[551] = "#eta(1440)";
    pId[552] = "#Lambda_{c}(2625)^{+}";
    pId[553] = "#Lambda_{c}(2625)^{-}";
    pId[554] = "";
    pId[555] = "";
    pId[556] = "";
    pId[557] = "";
    pId[558] = "";
    pId[559] = "NN";
}

void ntupleClass_MC::Fill_particleId(Int_t pdgId, Int_t IdSummary[NPARTICLES]){
    // Given the particle Id, it increases the relative counter
    switch (pdgId) {
            // QUARKS [6]
        case 1: // is a quark d
            IdSummary[0]++;
            break;
        case 2: // is a quark u
            IdSummary[1]++;
            break;
        case 3: // is a quark s
            IdSummary[2]++;
            break;
        case 4: // is a quark c
            IdSummary[3]++;
            break;
        case 5: // is a quark b
            IdSummary[4]++;
            break;
        case 6: // is a quark t
            IdSummary[5]++;
            break;
            // LEPTONS [12]
        case 11: // is a e-
            IdSummary[6]++;
            break;
        case -11: // is a e+
            IdSummary[7]++;
            break;
        case 12: // is a nu_e
            IdSummary[8]++;
            break;
        case -12: // is an anti nu_e
            IdSummary[9]++;
            break;
        case 13: // is a mu-
            IdSummary[10]++;
            break;
        case -13: // is a mu+
            IdSummary[11]++;
            break;
        case 14: // is a nu_mu
            IdSummary[12]++;
            break;
        case -14: // is an anti nu_mu
            IdSummary[13]++;
            break;
        case 15: // is a tau-
            IdSummary[14]++;
            break;
        case -15: // is a tau+
            IdSummary[15]++;
            break;
        case 16: // is a nu_tau
            IdSummary[16]++;
            break;
        case -16: // is an anti nu_tau
            IdSummary[17]++;
            break;
            // GAUGE & HIGGS BOSON(S) [5]
        case 21: // is a gluon
            IdSummary[18]++;
            break;
        case 22: // is a photon
            IdSummary[19]++;
            break;
        case 23: // is a Z0
            IdSummary[20]++;
            break;
        case 24: // is a W+
            IdSummary[21]++;
            break;
        case -24: // is a W-
            IdSummary[22]++;
            break;
            // DIQUARKS [25]
        case 1103: // is a (dd)1
            IdSummary[23]++;
            break;
        case 2101: // is a (ud)0
            IdSummary[24]++;
            break;
        case 2103: // is a (ud)1
            IdSummary[25]++;
            break;
        case 2203: // is a (uu)1
            IdSummary[26]++;
            break;
        case 3101: // is a (sd)0
            IdSummary[27]++;
            break;
        case 3103: // is a (sd)1
            IdSummary[28]++;
            break;
        case 3201: // is a (su)0
            IdSummary[29]++;
            break;
        case 3203: // is a (su)1
            IdSummary[30]++;
            break;
        case 3303: // is a (ss)1
            IdSummary[31]++;
            break;
        case 4101: // is a (cd)0
            IdSummary[32]++;
            break;
        case 4103: // is a (cd)1
            IdSummary[33]++;
            break;
        case 4201: // is a (cu)0
            IdSummary[34]++;
            break;
        case 4203: // is a (cu)1
            IdSummary[35]++;
            break;
        case 4301: // is a (cs)0
            IdSummary[36]++;
            break;
        case 4303: // is a (cs)1
            IdSummary[37]++;
            break;
        case 4403: // is a (cc)1
            IdSummary[38]++;
            break;
        case 5101: // is a (bd)0
            IdSummary[39]++;
            break;
        case 5103: // is a (bd)1
            IdSummary[40]++;
            break;
        case 5201: // is a (bu)0
            IdSummary[41]++;
            break;
        case 5203: // is a (bu)1
            IdSummary[42]++;
            break;
        case 5301: // is a (bs)0
            IdSummary[43]++;
            break;
        case 5303: // is a (bs)1
            IdSummary[44]++;
            break;
        case 5401: // is a (bc)0
            IdSummary[45]++;
            break;
        case 5403: // is a (bc)1
            IdSummary[46]++;
            break;
        case 5503: // is a (bb)1
            IdSummary[47]++;
            break;
            // LIGHT I=1 MESONS [92]
        case 111: // is a Pi0
            IdSummary[48]++;
            break;
        case -111: // is an anti Pi0
            IdSummary[49]++;
            break;
        case 211: // is a Pi+
            IdSummary[50]++;
            break;
        case -211: // is a Pi-
            IdSummary[51]++;
            break;
        case 9000111: // is a a0 (980)0
            IdSummary[52]++;
            break;
        case -9000111: // is an anti a0 (980)0
            IdSummary[53]++;
            break;
        case 9000211: // is a a0 (980)+
            IdSummary[54]++;
            break;
        case -9000211: // is a a0 (980)-
            IdSummary[55]++;
            break;
        case 100111: // is a π(1300)0
            IdSummary[56]++;
            break;
        case -100111: // is an anti π(1300)0
            IdSummary[57]++;
            break;
        case 100211: // is a π(1300)+
            IdSummary[58]++;
            break;
        case -100211: // is a π(1300)-
            IdSummary[59]++;
            break;
        case 10111: // is a a0 (1450)0
            IdSummary[60]++;
            break;
        case -10111: // is an anti a0 (1450)0
            IdSummary[61]++;
            break;
        case 10211: // is a a0 (1450)+
            IdSummary[62]++;
            break;
        case -10211: // is a a0 (1450)-
            IdSummary[63]++;
            break;
        case 200111: // is a π(1800)0
            IdSummary[64]++;
            break;
        case -200111: // is an anti π(1800)0
            IdSummary[65]++;
            break;
        case 200211: // is a π(1800)+
            IdSummary[66]++;
            break;
        case -200211: // is a π(1800)-
            IdSummary[67]++;
            break;
        case 113: // is a ρ(770)0
            IdSummary[68]++;
            break;
        case -113: // is an anti ρ(770)0
            IdSummary[69]++;
            break;
        case 213: // is a ρ(770)+
            IdSummary[70]++;
            break;
        case -213: // is a ρ(770)-
            IdSummary[71]++;
            break;
        case 10113: // is a b1 (1235)0
            IdSummary[72]++;
            break;
        case -10113: // is an anti b1 (1235)0
            IdSummary[73]++;
            break;
        case 10213: // is a b1 (1235)+
            IdSummary[74]++;
            break;
        case -10213: // is a b1 (1235)-
            IdSummary[75]++;
            break;
        case 20113: // is a a1 (1260)0
            IdSummary[76]++;
            break;
        case -20113: // is an anti a1 (1260)0
            IdSummary[77]++;
            break;
        case 20213: // is a a1 (1260)+
            IdSummary[78]++;
            break;
        case -20213: // is a a1 (1260)-
            IdSummary[79]++;
            break;
        case 9000113: // is a π1 (1400)0
            IdSummary[80]++;
            break;
        case -9000113: // is an anti π1 (1400)0
            IdSummary[81]++;
            break;
        case 9000213: // is a π1 (1400)+
            IdSummary[82]++;
            break;
        case -9000213: // is a π1 (1400)-
            IdSummary[83]++;
            break;
        case 100113: // is a ρ(1450)0
            IdSummary[84]++;
            break;
        case -100113: // is an anti ρ(1450)0
            IdSummary[85]++;
            break;
        case 100213: // is a ρ(1450)+
            IdSummary[86]++;
            break;
        case -100213: // is a ρ(1450)-
            IdSummary[87]++;
            break;
        case 9010113: // is a π1 (1600)0
            IdSummary[88]++;
            break;
        case -9010113: // is an anti π1 (1600)0
            IdSummary[89]++;
            break;
        case 9010213: // is a π1 (1600)+
            IdSummary[90]++;
            break;
        case -9010213: // is a π1 (1600)-
            IdSummary[91]++;
            break;
        case 9020113: // is a a1 (1640)0
            IdSummary[92]++;
            break;
        case -9020113: // is an anti a1 (1640)0
            IdSummary[93]++;
            break;
        case 9020213: // is a a1 (1640)+
            IdSummary[94]++;
            break;
        case -9020213: // is a a1 (1640)-
            IdSummary[95]++;
            break;
        case 30113: // is a ρ(1700)0
            IdSummary[96]++;
            break;
        case -30113: // is an anti ρ(1700)0
            IdSummary[97]++;
            break;
        case 30213: // is a ρ(1700)+
            IdSummary[98]++;
            break;
        case -30213: // is a ρ(1700)-
            IdSummary[99]++;
            break;
        case 9030113: // is a ρ(1900)0
            IdSummary[100]++;
            break;
        case -9030113: // is an anti ρ(1900)0
            IdSummary[101]++;
            break;
        case 9030213: // is a ρ(1900)+
            IdSummary[102]++;
            break;
        case -9030213: // is a ρ(1900)-
            IdSummary[103]++;
            break;
        case 9040113: // is a ρ(2150)0
            IdSummary[104]++;
            break;
        case -9040113: // is an anti ρ(2150)0
            IdSummary[105]++;
            break;
        case 9040213: // is a ρ(2150)+
            IdSummary[106]++;
            break;
        case -9040213: // is a ρ(2150)-
            IdSummary[107]++;
            break;
        case 115: // is a a2 (1320)0
            IdSummary[108]++;
            break;
        case -115: // is an anti a2 (1320)0
            IdSummary[109]++;
            break;
        case 215: // is a a2 (1320)+
            IdSummary[110]++;
            break;
        case -215: // is a a2 (1320)-
            IdSummary[111]++;
            break;
        case 10115: // is a π2 (1670)0
            IdSummary[112]++;
            break;
        case -10115: // is an anti π2 (1670)0
            IdSummary[113]++;
            break;
        case 10215: // is a π2 (1670)+
            IdSummary[114]++;
            break;
        case -10215: // is a π2 (1670)-
            IdSummary[115]++;
            break;
        case 100115: // is a a2(1700)0
            IdSummary[116]++;
            break;
        case -100115: // is an anti a2(1700)0
            IdSummary[117]++;
            break;
        case 100215: // is a a2(1700)+
            IdSummary[118]++;
            break;
        case -100215: // is a a2(1700)-
            IdSummary[119]++;
            break;
        case 9000115: // is a π2 (2100)0
            IdSummary[120]++;
            break;
        case -9000115: // is an anti π2 (2100)0
            IdSummary[121]++;
            break;
        case 9000215: // is a π2 (2100)+
            IdSummary[122]++;
            break;
        case -9000215: // is a π2 (2100)-
            IdSummary[123]++;
            break;
        case 117: // is a ρ3 (1690)0
            IdSummary[124]++;
            break;
        case -117: // is an anti ρ3 (1690)0
            IdSummary[125]++;
            break;
        case 217: // is a ρ3 (1690)+
            IdSummary[126]++;
            break;
        case -217: // is a ρ3 (1690)-
            IdSummary[127]++;
            break;
        case 9000117: // is a ρ3(1990)0
            IdSummary[128]++;
            break;
        case -9000117: // is an anti ρ3(1990)0
            IdSummary[129]++;
            break;
        case 9000217: // is a ρ3(1990)+
            IdSummary[130]++;
            break;
        case -9000217: // is a ρ3(1990)-
            IdSummary[131]++;
            break;
        case 9010117: // is a ρ3 (2250)0
            IdSummary[132]++;
            break;
        case -9010117: // is an anti ρ3 (2250)0
            IdSummary[133]++;
            break;
        case 9010217: // is a ρ3 (2250)+
            IdSummary[134]++;
            break;
        case -9010217: // is a ρ3 (2250)-
            IdSummary[135]++;
            break;
        case 119: // is a a4 (2040)0
            IdSummary[136]++;
            break;
        case -119: // is an anti a4 (2040)0
            IdSummary[137]++;
            break;
        case 219: // is a a4 (2040)+
            IdSummary[138]++;
            break;
        case -219: // is a a4 (2040)-
            IdSummary[139]++;
            break;
            // LIGHT I=0 MESONS [45]
            // (u\bar{u}, d\bar{d}, and s\bar{s} Admixtures)
        case 221: // is a η
            IdSummary[140]++;
            break;
        case 331: // is a η′(958)
            IdSummary[141]++;
            break;
        case 9000221: // is a f0 (600)
            IdSummary[142]++;
            break;
        case 9010221: // is a f0 (980)
            IdSummary[143]++;
            break;
        case 100221: // is a η(1295)
            IdSummary[144]++;
            break;
        case 10221: // is a f0 (1370)
            IdSummary[145]++;
            break;
        case 100331: // is a η(1440)
            IdSummary[146]++;
            break;
        case 9020221: // is a f0 (1500)
            IdSummary[147]++;
            break;
        case 10331: // is a f0 (1710)
            IdSummary[148]++;
            break;
        case 200221: // is a η(1760)
            IdSummary[149]++;
            break;
        case 9030221: // is a f0 (2020)
            IdSummary[150]++;
            break;
        case 9040221: // is a f0 (2100)
            IdSummary[151]++;
            break;
        case 9050221: // is a f0 (2200)
            IdSummary[152]++;
            break;
        case 9060221: // is a η(2225)
            IdSummary[153]++;
            break;
        case 9070221: // is a f0 (2330)
            IdSummary[154]++;
            break;
        case 223: // is a ω(782)
            IdSummary[155]++;
            break;
        case 333: // is a φ(1020)
            IdSummary[156]++;
            break;
        case 10223: // is a h1 (1170)
            IdSummary[157]++;
            break;
        case 20223: // is a f1 (1285)
            IdSummary[158]++;
            break;
        case 10333: // is a h1 (1380)
            IdSummary[159]++;
            break;
        case 20333: // is a f1 (1420)
            IdSummary[160]++;
            break;
        case 100223: // is a ω(1420)
            IdSummary[161]++;
            break;
        case 9000223: // is a f1 (1510)
            IdSummary[162]++;
            break;
        case 9010223: // is a h1 (1595)
            IdSummary[163]++;
            break;
        case 30223: // is a ω(1650)
            IdSummary[164]++;
            break;
        case 100333: // is a φ(1680)
            IdSummary[165]++;
            break;
        case 225: // is a f2 (1270)
            IdSummary[166]++;
            break;
        case 9000225: // is a f2 (1430)
            IdSummary[167]++;
            break;
        case 335: // is a f2′ (1525)
            IdSummary[168]++;
            break;
        case 9010225: // is a f2 (1565)
            IdSummary[169]++;
            break;
        case 9020225: // is a f2 (1640)
            IdSummary[170]++;
            break;
        case 10225: // is a η2 (1645)
            IdSummary[171]++;
            break;
        case 9030225: // is a f2 (1810)
            IdSummary[172]++;
            break;
        case 10335: // is a η2 (1870)
            IdSummary[173]++;
            break;
        case 9040225: // is a f2 (1910)
            IdSummary[174]++;
            break;
        case 100225: // is a f2 (1950)
            IdSummary[175]++;
            break;
        case 100335: // is a f2 (2010)
            IdSummary[176]++;
            break;
        case 9050225: // is a f2 (2150)
            IdSummary[177]++;
            break;
        case 9060225: // is a f2 (2300)
            IdSummary[178]++;
            break;
        case 9070225: // is a f2 (2340)
            IdSummary[179]++;
            break;
        case 227: // is a ω3 (1670)
            IdSummary[180]++;
            break;
        case 337: // is a φ3 (1850)
            IdSummary[181]++;
            break;
        case 229: // is a f4 (2050)
            IdSummary[182]++;
            break;
        case 9000339: // is a fJ (2220)
            IdSummary[183]++;
            break;
        case 9000229: // is a f4 (2300)
            IdSummary[184]++;
            break;
            // STRANGE MESONS [88]
        case 130: // is a KL0
            IdSummary[185]++;
            break;
        case -130: // is an anti KL0
            IdSummary[186]++;
            break;
        case 310: // is a KS0
            IdSummary[187]++;
            break;
        case -310: // is an anti KS0
            IdSummary[188]++;
            break;
        case 311: // is a K0
            IdSummary[189]++;
            break;
        case -311: // is an anti K0
            IdSummary[190]++;
            break;
        case 321: // is a K+
            IdSummary[191]++;
            break;
        case -321: // is a K-
            IdSummary[192]++;
            break;
        case 10311: // is a K0∗ (1430)0
            IdSummary[193]++;
            break;
        case -10311: // is an anti K0∗ (1430)0
            IdSummary[194]++;
            break;
        case 10321: // is a K0∗ (1430)+
            IdSummary[195]++;
            break;
        case -10321: // is a K0∗ (1430)-
            IdSummary[196]++;
            break;
        case 100311: // is a K (1460)0
            IdSummary[197]++;
            break;
        case -100311: // is an anti K (1460)0
            IdSummary[198]++;
            break;
        case 100321: // is a K (1460)+
            IdSummary[199]++;
            break;
        case -100321: // is a K (1460)-
            IdSummary[200]++;
            break;
        case 200311: // is a K (1830)0
            IdSummary[201]++;
            break;
        case -200311: // is an anti K (1830)0
            IdSummary[202]++;
            break;
        case 200321: // is a K (1830)+
            IdSummary[203]++;
            break;
        case -200321: // is a K (1830)-
            IdSummary[204]++;
            break;
        case 9000311: // is a K0∗ (1950)0
            IdSummary[205]++;
            break;
        case -9000311: // is an anti K0∗ (1950)0
            IdSummary[206]++;
            break;
        case 9000321: // is a K0∗ (1950)+
            IdSummary[207]++;
            break;
        case -9000321: // is a K0∗ (1950)-
            IdSummary[208]++;
            break;
        case 313: // is a K ∗ (892)0
            IdSummary[209]++;
            break;
        case -313: // is an anti K ∗ (892)0
            IdSummary[210]++;
            break;
        case 323: // is a K ∗ (892)+
            IdSummary[211]++;
            break;
        case -323: // is a K ∗ (892)-
            IdSummary[212]++;
            break;
        case 10313: // is a K1 (1270)0
            IdSummary[213]++;
            break;
        case -10313: // is an anti K1 (1270)0
            IdSummary[214]++;
            break;
        case 10323: // is a K1 (1270)+
            IdSummary[215]++;
            break;
        case -10323: // is a K1 (1270)-
            IdSummary[216]++;
            break;
        case 20313: // is a K1 (1400)0
            IdSummary[217]++;
            break;
        case -20313: // is an anti K1 (1400)0
            IdSummary[218]++;
            break;
        case 20323: // is a K1 (1400)+
            IdSummary[219]++;
            break;
        case -20323: // is a K1 (1400)-
            IdSummary[220]++;
            break;
        case 100313: // is a K ∗ (1410)0
            IdSummary[221]++;
            break;
        case -100313: // is an anti K ∗ (1410)0
            IdSummary[222]++;
            break;
        case 100323: // is a K ∗ (1410)+
            IdSummary[223]++;
            break;
        case -100323: // is a K ∗ (1410)-
            IdSummary[224]++;
            break;
        case 9000313: // is a K1 (1650)0
            IdSummary[225]++;
            break;
        case -9000313: // is an anti K1 (1650)0
            IdSummary[226]++;
            break;
        case 9000323: // is a K1 (1650)+
            IdSummary[227]++;
            break;
        case -9000323: // is a K1 (1650)-
            IdSummary[228]++;
            break;
        case 30313: // is a K ∗ (1680)0
            IdSummary[229]++;
            break;
        case -30313: // is an anti K ∗ (1680)0
            IdSummary[230]++;
            break;
        case 30323: // is a K ∗ (1680)+
            IdSummary[231]++;
            break;
        case -30323: // is a K ∗ (1680)-
            IdSummary[232]++;
            break;
        case 315: // is a K2∗ (1430)0
            IdSummary[233]++;
            break;
        case -315: // is an anti K2∗ (1430)0
            IdSummary[234]++;
            break;
        case 325: // is a K2∗ (1430)+
            IdSummary[235]++;
            break;
        case -325: // is a K2∗ (1430)-
            IdSummary[236]++;
            break;
        case 9000315: // is a K2 (1580)0
            IdSummary[237]++;
            break;
        case -9000315: // is an anti K2 (1580)0
            IdSummary[238]++;
            break;
        case 9000325: // is a K2 (1580)+
            IdSummary[239]++;
            break;
        case -9000325: // is a K2 (1580)-
            IdSummary[240]++;
            break;
        case 10315: // is a K2 (1770)0
            IdSummary[241]++;
            break;
        case -10315: // is an anti K2 (1770)0
            IdSummary[242]++;
            break;
        case 10325: // is a K2 (1770)+
            IdSummary[243]++;
            break;
        case -10325: // is a K2 (1770)-
            IdSummary[244]++;
            break;
        case 20315: // is a K2 (1820)0
            IdSummary[245]++;
            break;
        case -20315: // is an anti K2 (1820)0
            IdSummary[246]++;
            break;
        case 20325: // is a K2 (1820)+
            IdSummary[247]++;
            break;
        case -20325: // is a K2 (1820)-
            IdSummary[248]++;
            break;
        case 100315: // is a K2∗ (1980)0
            IdSummary[249]++;
            break;
        case -100315: // is an anti K2∗ (1980)0
            IdSummary[250]++;
            break;
        case 100325: // is a K2∗ (1980)+
            IdSummary[251]++;
            break;
        case -100325: // is a K2∗ (1980)-
            IdSummary[252]++;
            break;
        case 9010315: // is a K2 (2250)0
            IdSummary[253]++;
            break;
        case -9010315: // is an anti K2 (2250)0
            IdSummary[254]++;
            break;
        case 9010325: // is a K2 (2250)+
            IdSummary[255]++;
            break;
        case -9010325: // is a K2 (2250)-
            IdSummary[256]++;
            break;
        case 317: // is a K3∗ (1780)0
            IdSummary[257]++;
            break;
        case -317: // is an anti K3∗ (1780)0
            IdSummary[258]++;
            break;
        case 327: // is a K3∗ (1780)+
            IdSummary[259]++;
            break;
        case -327: // is a K3∗ (1780)-
            IdSummary[260]++;
            break;
        case 9010317: // is a K3 (2320)0
            IdSummary[261]++;
            break;
        case -9010317: // is an anti K3 (2320)0
            IdSummary[262]++;
            break;
        case 9010327: // is a K3 (2320)+
            IdSummary[263]++;
            break;
        case -9010327: // is a K3 (2320)-
            IdSummary[264]++;
            break;
        case 319: // is a K4∗ (2045)0
            IdSummary[265]++;
            break;
        case -319: // is an anti K4∗ (2045)0
            IdSummary[266]++;
            break;
        case 329: // is a K4∗ (2045)+
            IdSummary[267]++;
            break;
        case -329: // is a K4∗ (2045)-
            IdSummary[268]++;
            break;
        case 9000319: // is a K4 (2500)0
            IdSummary[269]++;
            break;
        case -9000319: // is an anti K4 (2500)0
            IdSummary[270]++;
            break;
        case 9000329: // is a K4 (2500)+
            IdSummary[271]++;
            break;
        case -9000329: // is a K4 (2500)-
            IdSummary[272]++;
            break;
            // CHARMED MESONS [36]
        case 411: // is a D+
            IdSummary[273]++;
            break;
        case -411: // is a D-
            IdSummary[274]++;
            break;
        case 421: // is a D0
            IdSummary[275]++;
            break;
        case -421: // is an anti D0
            IdSummary[276]++;
            break;
        case 10411: // is a D0∗+
            IdSummary[277]++;
            break;
        case -10411: // is a D0∗-
            IdSummary[278]++;
            break;
        case 10421: // is a D0∗0
            IdSummary[279]++;
            break;
        case -10421: // is an anti D0∗0
            IdSummary[280]++;
            break;
        case 413: // is a D∗(2010)+
            IdSummary[281]++;
            break;
        case -413: // is a D∗(2010)-
            IdSummary[282]++;
            break;
        case 423: // is a D∗(2007)0
            IdSummary[283]++;
            break;
        case -423: // is an anti D∗(2007)0
            IdSummary[284]++;
            break;
        case 10413: // is a D1 (2420)+
            IdSummary[285]++;
            break;
        case -10413: // is a D1 (2420)-
            IdSummary[286]++;
            break;
        case 10423: // is a D1(2420)0
            IdSummary[287]++;
            break;
        case -10423: // is an anti D1(2420)0
            IdSummary[288]++;
            break;
        case 20413: // is a D1(H)+
            IdSummary[289]++;
            break;
        case -20413: // is a D1(H)-
            IdSummary[290]++;
            break;
        case 20423: // is a D1(H)0
            IdSummary[291]++;
            break;
        case -20423: // is an anti D1(H)0
            IdSummary[292]++;
            break;
        case 415: // is a D2∗(2460)+
            IdSummary[293]++;
            break;
        case -415: // is a D2∗(2460)-
            IdSummary[294]++;
            break;
        case 425: // is a D2∗(2460)0
            IdSummary[295]++;
            break;
        case -425: // is an anti D2∗(2460)0
            IdSummary[296]++;
            break;
        case 431: // is a Ds+
            IdSummary[297]++;
            break;
        case -431: // is a Ds-
            IdSummary[298]++;
            break;
        case 10431: // is a Ds0∗+
            IdSummary[299]++;
            break;
        case -10431: // is a Ds0∗-
            IdSummary[300]++;
            break;
        case 433: // is a Ds∗+
            IdSummary[301]++;
            break;
        case -433: // is a Ds∗-
            IdSummary[302]++;
            break;
        case 10433: // is a Ds1 (2536)+
            IdSummary[303]++;
            break;
        case -10433: // is a Ds1 (2536)-
            IdSummary[304]++;
            break;
        case 20433: // is a Ds1(H)+
            IdSummary[305]++;
            break;
        case -20433: // is a Ds1(H)-
            IdSummary[306]++;
            break;
        case 435: // is a Ds2*+
            IdSummary[307]++;
            break;
        case -435: // is a Ds2*-
            IdSummary[308]++;
            break;
            // BOTTOM MESONS [48]
        case 511: // is a B0
            IdSummary[309]++;
            break;
        case -511: // is an anti B0
            IdSummary[310]++;
            break;
        case 521: // is a B+
            IdSummary[311]++;
            break;
        case -521: // is a B-
            IdSummary[312]++;
            break;
        case 10511: // is a B0∗0
            IdSummary[313]++;
            break;
        case -10511: // is an anti B0∗0
            IdSummary[314]++;
            break;
        case 10521: // is a B0∗+
            IdSummary[315]++;
            break;
        case -10521: // is a B0∗-
            IdSummary[316]++;
            break;
        case 513: // is a B∗0
            IdSummary[317]++;
            break;
        case -513: // is an anti B∗0
            IdSummary[318]++;
            break;
        case 523: // is a B∗+
            IdSummary[319]++;
            break;
        case -523: // is a B∗-
            IdSummary[320]++;
            break;
        case 10513: // is a B1 (L)0
            IdSummary[321]++;
            break;
        case -10513: // is an anti B1 (L)0
            IdSummary[322]++;
            break;
        case 10523: // is a B1(L)+
            IdSummary[323]++;
            break;
        case -10523: // is a B1(L)-
            IdSummary[324]++;
            break;
        case 20513: // is a B1(H)0
            IdSummary[325]++;
            break;
        case -20513: // is an anti B1(H)0
            IdSummary[326]++;
            break;
        case 20523: // is a B1(H)+
            IdSummary[327]++;
            break;
        case -20523: // is a B1(H)-
            IdSummary[328]++;
            break;
        case 515: // is a B2∗0
            IdSummary[329]++;
            break;
        case -515: // is an anti B2∗0
            IdSummary[330]++;
            break;
        case 525: // is a B2∗+
            IdSummary[331]++;
            break;
        case -525: // is a B2∗-
            IdSummary[332]++;
            break;
        case 531: // is a Bs0
            IdSummary[333]++;
            break;
        case -531: // is an anti Bs0
            IdSummary[334]++;
            break;
        case 10531: // is a Bs0*0
            IdSummary[335]++;
            break;
        case -10531: // is an anti Bs0*0
            IdSummary[336]++;
            break;
        case 533: // is a Bs∗0
            IdSummary[337]++;
            break;
        case -533: // is an anti Bs∗0
            IdSummary[338]++;
            break;
        case 10533: // is a Bs1(L)0
            IdSummary[339]++;
            break;
        case -10533: // is an anti Bs1(L)0
            IdSummary[340]++;
            break;
        case 20533: // is a Bs1(H)0
            IdSummary[341]++;
            break;
        case -20533: // is an anti Bs1(H)0
            IdSummary[342]++;
            break;
        case 535: // is a Bs2∗0
            IdSummary[343]++;
            break;
        case -535: // is an anti Bs2∗0
            IdSummary[344]++;
            break;
        case 541: // is a Bc+
            IdSummary[345]++;
            break;
        case -541: // is a Bc-
            IdSummary[346]++;
            break;
        case 10541: // is a Bc0∗+
            IdSummary[347]++;
            break;
        case -10541: // is a Bc0∗-
            IdSummary[348]++;
            break;
        case 543: // is a Bc∗+
            IdSummary[349]++;
            break;
        case -543: // is a Bc∗-
            IdSummary[350]++;
            break;
        case 10543: // is a Bc1(L)+
            IdSummary[351]++;
            break;
        case -10543: // is a Bc1(L)-
            IdSummary[352]++;
            break;
        case 20543: // is a Bc1(H)+
            IdSummary[353]++;
            break;
        case -20543: // is a Bc1(H)-
            IdSummary[354]++;
            break;
        case 545: // is a Bc2∗+
            IdSummary[355]++;
            break;
        case -545: // is a Bc2∗-
            IdSummary[356]++;
            break;
            // c\bar{c} MESONS [13]
        case 441: // is a ηc(1S)
            IdSummary[357]++;
            break;
        case 10441: // is a χc0(1P)
            IdSummary[358]++;
            break;
        case 100441: // is a ηc (2S)
            IdSummary[359]++;
            break;
        case 443: // is a J/ψ(1S)
            IdSummary[360]++;
            break;
        case 10443: // is a hc(1P)
            IdSummary[361]++;
            break;
        case 20443: // is a χc1(1P)
            IdSummary[362]++;
            break;
        case 100443: // is a ψ(2S)
            IdSummary[363]++;
            break;
        case 30443: // is a ψ(3770)
            IdSummary[364]++;
            break;
        case 9000443: // is a ψ(4040)
            IdSummary[365]++;
            break;
        case 9010443: // is a ψ(4160)
            IdSummary[366]++;
            break;
        case 9020443: // is a ψ(4415)
            IdSummary[367]++;
            break;
        case 445: // is a χc2(1P)
            IdSummary[368]++;
            break;
        case 9000445: // is a ψ(3836)
            IdSummary[369]++;
            break;
            //b\bar{b} MESONS [29]
        case 551: // is a ηb (1S)
            IdSummary[370]++;
            break;
        case 10551: // is a χb0(1P)
            IdSummary[371]++;
            break;
        case 100551: // is a ηb (2S)
            IdSummary[372]++;
            break;
        case 110551: // is a χb0(2P)
            IdSummary[373]++;
            break;
        case 200551: // is a ηb (3S)
            IdSummary[374]++;
            break;
        case 210551: // is a χb0(3P)
            IdSummary[375]++;
            break;
        case 553: // is a Υ(1S)
            IdSummary[376]++;
            break;
        case 10553: // is a hb(1P)
            IdSummary[377]++;
            break;
        case 20553: // is a χb1(1P)
            IdSummary[378]++;
            break;
        case 30553: // is a Υ1(1D)
            IdSummary[379]++;
            break;
        case 100553: // is a Υ(2S)
            IdSummary[380]++;
            break;
        case 110553: // is a hb(2P)
            IdSummary[381]++;
            break;
        case 120553: // is a χb1(2P)
            IdSummary[382]++;
            break;
        case 130553: // is a Υ1(2D)
            IdSummary[383]++;
            break;
        case 200553: // is a Υ(3S)
            IdSummary[384]++;
            break;
        case 210553: // is a hb(3P)
            IdSummary[385]++;
            break;
        case 220553: // is a χb1(3P)
            IdSummary[386]++;
            break;
        case 300553: // is a Υ(4S)
            IdSummary[387]++;
            break;
        case 9000553: // is a Υ (10860)
            IdSummary[388]++;
            break;
        case 9010553: // is a Υ (11020)
            IdSummary[389]++;
            break;
        case 555: // is a χb2 (1P )
            IdSummary[390]++;
            break;
        case 10555: // is a ηb2(1D)
            IdSummary[391]++;
            break;
        case 20555: // is a Υ2(1D)
            IdSummary[392]++;
            break;
        case 100555: // is a χb2(2P)
            IdSummary[393]++;
            break;
        case 110555: // is a ηb2(2D)
            IdSummary[394]++;
            break;
        case 120555: // is a Υ2(2D)
            IdSummary[395]++;
            break;
        case 200555: // is a χb2(3P)
            IdSummary[396]++;
            break;
        case 557: // is a Υ3(1D)
            IdSummary[397]++;
            break;
        case 100557: // is a Υ3(2D)
            IdSummary[398]++;
            break;
            // LIGHT BARYONS [12]
        case 2212: // is a proton
            IdSummary[399]++;
            break;
        case -2212: // is an anti proton
            IdSummary[400]++;
            break;
        case 2112: // is a neutron
            IdSummary[401]++;
            break;
        case -2112: // is an anti neutron
            IdSummary[402]++;
            break;
        case 2224: // is a Delta++
            IdSummary[403]++;
            break;
        case -2224: // is an anti Delta++
            IdSummary[404]++;
            break;
        case 2214: // is a Delta+
            IdSummary[405]++;
            break;
        case -2214: // is an anti Delta+
            IdSummary[406]++;
            break;
        case 2114: // is a Delta0
            IdSummary[407]++;
            break;
        case -2114: // is an anti Delta0
            IdSummary[408]++;
            break;
        case 1114: // is a Delta -
            IdSummary[409]++;
            break;
        case -1114: // is an AntiDelta -
            IdSummary[410]++;
            break;
            //STRANGE BARYONS [24]
        case 3122: // is a Λ
            IdSummary[411]++;
            break;
        case -3122: // is an anti Λ
            IdSummary[412]++;
            break;
        case 3222: // is a Σ+
            IdSummary[413]++;
            break;
        case -3222: // is an anti Σ+
            IdSummary[414]++;
            break;
        case 3212: // is a Σ0
            IdSummary[415]++;
            break;
        case -3212: // is an anti Σ0
            IdSummary[416]++;
            break;
        case 3112: // is a Σ-
            IdSummary[417]++;
            break;
        case -3112: // is am anti Σ-
            IdSummary[418]++;
            break;
        case 3224: // is a Σ*+
            IdSummary[419]++;
            break;
        case -3224: // is an anti Σ*+
            IdSummary[420]++;
            break;
        case 3214: // is a Σ*0
            IdSummary[421]++;
            break;
        case -3214: // is an anti Σ*0
            IdSummary[422]++;
            break;
        case 3114: // is a Σ*-
            IdSummary[423]++;
            break;
        case -3114: // is an anti Σ*-
            IdSummary[424]++;
            break;
        case 3322: // is a Ξ0
            IdSummary[425]++;
            break;
        case -3322: // is an anti Ξ0
            IdSummary[426]++;
            break;
        case 3312: // is a Ξ-
            IdSummary[427]++;
            break;
        case -3312: // is an anti Ξ-
            IdSummary[428]++;
            break;
        case 3324: // is a Ξ*0
            IdSummary[429]++;
            break;
        case -3324: // is an anti Ξ*0
            IdSummary[430]++;
            break;
        case 3314: // is a Ξ*-
            IdSummary[431]++;
            break;
        case -3314: // is an anti Ξ*-
            IdSummary[432]++;
            break;
        case 3334: // is a Ω−
            IdSummary[433]++;
            break;
        case -3334: // is an anti Ω−
            IdSummary[434]++;
            break;
            // CHARMED BARYONS [44]
        case 4122: // is a Lambda+c
            IdSummary[435]++;
            break;
        case -4122: // is a Lambda-c
            IdSummary[436]++;
            break;
        case 4222: // is a Σc++
            IdSummary[437]++;
            break;
        case -4222: //is an anti Σc++
            IdSummary[438]++;
            break;
        case 4212: // is a Σc+
            IdSummary[439]++;
            break;
        case -4212: // is a Σc-
            IdSummary[440]++;
            break;
        case 4112: // is a Σc0
            IdSummary[441]++;
            break;
        case -4112: // is an anti Σc0
            IdSummary[442]++;
            break;
        case 4224: // is a Σc∗++
            IdSummary[443]++;
            break;
        case -4224: // is an anti Σc∗++
            IdSummary[444]++;
            break;
        case 4214: // is a Σc∗+
            IdSummary[445]++;
            break;
        case -4214: // is a Σc∗-
            IdSummary[446]++;
            break;
        case 4114: // is a Σc∗0
            IdSummary[447]++;
            break;
        case -4114: // is an anti Σc∗0
            IdSummary[448]++;
            break;
        case 4232: // is Ξc+
            IdSummary[449]++;
            break;
        case -4232: // is a Ξc-
            IdSummary[450]++;
            break;
        case 4132: // is a Ξc0
            IdSummary[451]++;
            break;
        case -4132: // is an anti Ξc0
            IdSummary[452]++;
            break;
        case 4322: // is a Ξ′c+
            IdSummary[453]++;
            break;
        case -4322: // is a Ξ′c-
            IdSummary[454]++;
            break;
        case 4312: // is a Ξ′c0
            IdSummary[455]++;
            break;
        case -4312: // is an anti Ξ′c0
            IdSummary[456]++;
            break;
        case 4324: // is a Ξc∗+
            IdSummary[457]++;
            break;
        case -4324: // is a Ξc∗-
            IdSummary[458]++;
            break;
        case 4314: // is a Ξc∗0
            IdSummary[459]++;
            break;
        case -4314: // is an anti Ξc∗0
            IdSummary[460]++;
            break;
        case 4332: // is a Ωc0
            IdSummary[461]++;
            break;
        case -4332: // is an anti Ωc0
            IdSummary[462]++;
            break;
        case 4334: // is a Ω*c0
            IdSummary[463]++;
            break;
        case -4334: // is an anti Ω*c0
            IdSummary[464]++;
            break;
        case 4412: // is a Ξcc+
            IdSummary[465]++;
            break;
        case -4412: // is a Ξcc-
            IdSummary[466]++;
            break;
        case 4422: // is a Ξcc++
            IdSummary[467]++;
            break;
        case -4422: // is a Ξcc--
            IdSummary[468]++;
            break;
        case 4414: // is a Ξcc∗+
            IdSummary[469]++;
            break;
        case -4414: // is a Ξcc∗-
            IdSummary[470]++;
            break;
        case 4424: // is a Ξcc∗++
            IdSummary[471]++;
            break;
        case -4424: // is a Ξcc∗--
            IdSummary[472]++;
            break;
        case 4432: // is a Ωcc+
            IdSummary[473]++;
            break;
        case -4432: // is a Ωcc-
            IdSummary[474]++;
            break;
        case 4434: // is a Ωcc∗+
            IdSummary[475]++;
            break;
        case -4434: // is a Ωcc∗-
            IdSummary[476]++;
            break;
        case 4444: // is a Ωccc++
            IdSummary[477]++;
            break;
        case -4444: // is a Ωccc--
            IdSummary[478]++;
            break;
            //BOTTOM BARYONS [70]
        case 5122: // is a Λ0b
            IdSummary[479]++;
            break;
        case -5122: // is an anti Λ0b
            IdSummary[480]++;
            break;
        case 5112: // is a Σb-
            IdSummary[481]++;
            break;
        case -5112: // is an anti Σb-
            IdSummary[482]++;
            break;
        case 5212: // is a Σb0
            IdSummary[483]++;
            break;
        case -5212: // is an anti Σb0
            IdSummary[484]++;
            break;
        case 5222: // is a Σb+
            IdSummary[485]++;
            break;
        case -5222: // is an anti Σb+
            IdSummary[486]++;
            break;
        case 5114: // is a Σb*-
            IdSummary[487]++;
            break;
        case -5114: // is an anti Σb*-
            IdSummary[488]++;
            break;
        case 5214: // is a Σb*0
            IdSummary[489]++;
            break;
        case -5214: // is an anti Σb*0
            IdSummary[490]++;
            break;
        case 5224: // is a Σb∗+
            IdSummary[491]++;
            break;
        case -5224: // is an anti Σb∗+
            IdSummary[492]++;
            break;
        case 5132: // is a Ξb−
            IdSummary[493]++;
            break;
        case -5132: // is a Ξb+
            IdSummary[494]++;
            break;
        case 5232: // is a Ξb0
            IdSummary[495]++;
            break;
        case -5232: // is an anti Ξb0
            IdSummary[496]++;
            break;
        case 5312: // is a Ξ′b−
            IdSummary[497]++;
            break;
        case -5312: // is a Ξ′b+
            IdSummary[498]++;
            break;
        case 5322: // is a Ξb′0
            IdSummary[499]++;
            break;
        case -5322: // is an anti Ξb′0
            IdSummary[500]++;
            break;
        case 5314: // is a Ξb∗−
            IdSummary[501]++;
            break;
        case -5314: // is a Ξb∗+
            IdSummary[502]++;
            break;
        case 5324: // is a Ξb∗0
            IdSummary[503]++;
            break;
        case -5324: // is an anti Ξb∗0
            IdSummary[504]++;
            break;
        case 5332: // is a Ωb−
            IdSummary[505]++;
            break;
        case -5332: // is a Ωb+
            IdSummary[506]++;
            break;
        case 5334: // is a Ωb*−
            IdSummary[507]++;
            break;
        case -5334: // is a Ωb*+
            IdSummary[508]++;
            break;
        case 5142: // is a Ξbc0
            IdSummary[509]++;
            break;
        case -5142: // is an anti Ξbc0
            IdSummary[510]++;
            break;
        case 5242: // is a Ξbc+
            IdSummary[511]++;
            break;
        case -5242: // is a Ξbc-
            IdSummary[512]++;
            break;
        case 5412: // is a Ξ′0bc
            IdSummary[513]++;
            break;
        case -5412: // is an anti Ξ′0bc
            IdSummary[514]++;
            break;
        case 5422: // is a Ξ′bc+
            IdSummary[515]++;
            break;
        case -5422: // is a Ξ′bc-
            IdSummary[516]++;
            break;
        case 5414: // is a Ξbc∗0
            IdSummary[517]++;
            break;
        case -5414: // is an anti Ξbc∗0
            IdSummary[518]++;
            break;
        case 5424: // is a Ξbc∗+
            IdSummary[519]++;
            break;
        case -5424: // is a Ξbc∗-
            IdSummary[520]++;
            break;
        case 5342: // is a Ωbc0
            IdSummary[521]++;
            break;
        case -5342: // is an anti Ωbc0
            IdSummary[522]++;
            break;
        case 5432: // is a Ωbc′0
            IdSummary[523]++;
            break;
        case -5432: // is an anti Ωbc′0
            IdSummary[524]++;
            break;
        case 5434: // is a Ωbc∗0
            IdSummary[525]++;
            break;
        case -5434: // is an anti Ωbc∗0
            IdSummary[526]++;
            break;
        case 5442: // is a Ωbcc+
            IdSummary[527]++;
            break;
        case -5442: // is a Ωbcc-
            IdSummary[528]++;
            break;
        case 5444: // is a Ωbcc∗+
            IdSummary[529]++;
            break;
        case -5444: // is a Ωbcc∗-
            IdSummary[530]++;
            break;
        case 5512: // is a Ξbb−
            IdSummary[531]++;
            break;
        case -5512: // is a Ξbb+
            IdSummary[532]++;
            break;
        case 5522: // is a Ξbb0
            IdSummary[533]++;
            break;
        case -5522: // is an anti Ξbb0
            IdSummary[534]++;
            break;
        case 5514: // is a Ξbb∗−
            IdSummary[535]++;
            break;
        case -5514: // is a Ξbb∗+
            IdSummary[536]++;
            break;
        case 5524: // is a Ξbb∗0
            IdSummary[537]++;
            break;
        case -5524: // is an anti Ξbb∗0
            IdSummary[538]++;
            break;
        case 5532: // is a Ωbb−
            IdSummary[539]++;
            break;
        case -5532: // is a Ωbb+
            IdSummary[540]++;
            break;
        case 5534: // is a Ωbb∗−
            IdSummary[541]++;
            break;
        case -5534: // is a Ωbb∗+
            IdSummary[542]++;
            break;
        case 5542: // is a Ωbbc0
            IdSummary[543]++;
            break;
        case -5542: // is an anti Ωbbc0
            IdSummary[544]++;
            break;
        case 5544: // is a Ωbbc∗0
            IdSummary[545]++;
            break;
        case -5544: // is an anti Ωbbc∗0
            IdSummary[546]++;
            break;
        case 5554: // is a Ωbbb−
            IdSummary[547]++;
            break;
        case -5554: // is a Ωbbb+
            IdSummary[548]++;
            break;
        case 14122: // is a Lambda_c(2593)+
            IdSummary[549]++;
            break;
        case -14122: // is a Lambda_c(2593)-
            IdSummary[550]++;
            break;
        case 30221: // is a eta(1440)
            IdSummary[551]++;
            break;
        case 4124: // is a Lambda_c(2625)+
            IdSummary[552]++;
            break;
        case -4124: // is a Lambda_c(2625)-
             IdSummary[553]++;
             break;
            /*
             case ??: // is a
             IdSummary[554]++;
             break;
             case ??: // is a
             IdSummary[555]++;
             break;
             case 0: // is a UNKNOWN
             IdSummary[556]++;
             break;
             case 0: // is a UNKNOWN
             IdSummary[557]++;
             break;
             case 0: // is a UNKNOWN
             IdSummary[558]++;
             break;
             */
            
        case 0: // is a UNKNOWN
            IdSummary[559]++;
            break;
            
        default: cout << "Unknown particle! PdgId = " << pdgId << endl;
            break;
    }
}

void ntupleClass_MC::Fill_particleId_2D(Int_t pdgId, Int_t pdgIdMother, Int_t IdSummary[NPARTICLES][NPARTICLES]){
    // Given the particle Id, it checks at which particle corresponds and increases the relative counter 2D
    switch (pdgId) {
            // QUARKS [6]
        case 1: // is a quark d
            Fill_particleId(pdgIdMother, IdSummary[0]);
            break;
        case 2: // is a quark u
            Fill_particleId(pdgIdMother, IdSummary[1]);
            break;
        case 3: // is a quark s
            Fill_particleId(pdgIdMother, IdSummary[2]);
            break;
        case 4: // is a quark c
            Fill_particleId(pdgIdMother, IdSummary[3]);
            break;
        case 5: // is a quark b
            Fill_particleId(pdgIdMother, IdSummary[4]);
            break;
        case 6: // is a quark t
            Fill_particleId(pdgIdMother, IdSummary[5]);
            break;
            // LEPTONS [12]
        case 11: // is a e-
            Fill_particleId(pdgIdMother, IdSummary[6]);
            break;
        case -11: // is a e+
            Fill_particleId(pdgIdMother, IdSummary[7]);
            break;
        case 12: // is a nu_e
            Fill_particleId(pdgIdMother, IdSummary[8]);
            break;
        case -12: // is an anti nu_e
            Fill_particleId(pdgIdMother, IdSummary[9]);
            break;
        case 13: // is a mu-
            Fill_particleId(pdgIdMother, IdSummary[10]);
            break;
        case -13: // is a mu+
            Fill_particleId(pdgIdMother, IdSummary[11]);
            break;
        case 14: // is a nu_mu
            Fill_particleId(pdgIdMother, IdSummary[12]);
            break;
        case -14: // is an anti nu_mu
            Fill_particleId(pdgIdMother, IdSummary[13]);
            break;
        case 15: // is a tau-
            Fill_particleId(pdgIdMother, IdSummary[14]);
            break;
        case -15: // is a tau+
            Fill_particleId(pdgIdMother, IdSummary[15]);
            break;
        case 16: // is a nu_tau
            Fill_particleId(pdgIdMother, IdSummary[16]);
            break;
        case -16: // is an anti nu_tau
            Fill_particleId(pdgIdMother, IdSummary[17]);
            break;
            // GAUGE & HIGGS BOSON(S) [5]
        case 21: // is a gluon
            Fill_particleId(pdgIdMother, IdSummary[18]);
            break;
        case 22: // is a photon
            Fill_particleId(pdgIdMother, IdSummary[19]);
            break;
        case 23: // is a Z0
            Fill_particleId(pdgIdMother, IdSummary[20]);
            break;
        case 24: // is a W+
            Fill_particleId(pdgIdMother, IdSummary[21]);
            break;
        case -24: // is a W-
            Fill_particleId(pdgIdMother, IdSummary[22]);
            break;
            // DIQUARKS [25]
        case 1103: // is a (dd)1
            Fill_particleId(pdgIdMother, IdSummary[23]);
            break;
        case 2101: // is a (ud)0
            Fill_particleId(pdgIdMother, IdSummary[24]);
            break;
        case 2103: // is a (ud)1
            Fill_particleId(pdgIdMother, IdSummary[25]);
            break;
        case 2203: // is a (uu)1
            Fill_particleId(pdgIdMother, IdSummary[26]);
            break;
        case 3101: // is a (sd)0
            Fill_particleId(pdgIdMother, IdSummary[27]);
            break;
        case 3103: // is a (sd)1
            Fill_particleId(pdgIdMother, IdSummary[28]);
            break;
        case 3201: // is a (su)0
            Fill_particleId(pdgIdMother, IdSummary[29]);
            break;
        case 3203: // is a (su)1
            Fill_particleId(pdgIdMother, IdSummary[30]);
            break;
        case 3303: // is a (ss)1
            Fill_particleId(pdgIdMother, IdSummary[31]);
            break;
        case 4101: // is a (cd)0
            Fill_particleId(pdgIdMother, IdSummary[32]);
            break;
        case 4103: // is a (cd)1
            Fill_particleId(pdgIdMother, IdSummary[33]);
            break;
        case 4201: // is a (cu)0
            Fill_particleId(pdgIdMother, IdSummary[34]);
            break;
        case 4203: // is a (cu)1
            Fill_particleId(pdgIdMother, IdSummary[35]);
            break;
        case 4301: // is a (cs)0
            Fill_particleId(pdgIdMother, IdSummary[36]);
            break;
        case 4303: // is a (cs)1
            Fill_particleId(pdgIdMother, IdSummary[37]);
            break;
        case 4403: // is a (cc)1
            Fill_particleId(pdgIdMother, IdSummary[38]);
            break;
        case 5101: // is a (bd)0
            Fill_particleId(pdgIdMother, IdSummary[39]);
            break;
        case 5103: // is a (bd)1
            Fill_particleId(pdgIdMother, IdSummary[40]);
            break;
        case 5201: // is a (bu)0
            Fill_particleId(pdgIdMother, IdSummary[41]);
            break;
        case 5203: // is a (bu)1
            Fill_particleId(pdgIdMother, IdSummary[42]);
            break;
        case 5301: // is a (bs)0
            Fill_particleId(pdgIdMother, IdSummary[43]);
            break;
        case 5303: // is a (bs)1
            Fill_particleId(pdgIdMother, IdSummary[44]);
            break;
        case 5401: // is a (bc)0
            Fill_particleId(pdgIdMother, IdSummary[45]);
            break;
        case 5403: // is a (bc)1
            Fill_particleId(pdgIdMother, IdSummary[46]);
            break;
        case 5503: // is a (bb)1
            Fill_particleId(pdgIdMother, IdSummary[47]);
            break;
            // LIGHT I=1 MESONS [92]
        case 111: // is a Pi0
            Fill_particleId(pdgIdMother, IdSummary[48]);
            break;
        case -111: // is an anti Pi0
            Fill_particleId(pdgIdMother, IdSummary[49]);
            break;
        case 211: // is a Pi+
            Fill_particleId(pdgIdMother, IdSummary[50]);
            break;
        case -211: // is a Pi-
            Fill_particleId(pdgIdMother, IdSummary[51]);
            break;
        case 9000111: // is a a0 (980)0
            Fill_particleId(pdgIdMother, IdSummary[52]);
            break;
        case -9000111: // is an anti a0 (980)0
            Fill_particleId(pdgIdMother, IdSummary[53]);
            break;
        case 9000211: // is a a0 (980)+
            Fill_particleId(pdgIdMother, IdSummary[54]);
            break;
        case -9000211: // is a a0 (980)-
            Fill_particleId(pdgIdMother, IdSummary[55]);
            break;
        case 100111: // is a π(1300)0
            Fill_particleId(pdgIdMother, IdSummary[56]);
            break;
        case -100111: // is an anti π(1300)0
            Fill_particleId(pdgIdMother, IdSummary[57]);
            break;
        case 100211: // is a π(1300)+
            Fill_particleId(pdgIdMother, IdSummary[58]);
            break;
        case -100211: // is a π(1300)-
            Fill_particleId(pdgIdMother, IdSummary[59]);
            break;
        case 10111: // is a a0 (1450)0
            Fill_particleId(pdgIdMother, IdSummary[60]);
            break;
        case -10111: // is an anti a0 (1450)0
            Fill_particleId(pdgIdMother, IdSummary[61]);
            break;
        case 10211: // is a a0 (1450)+
            Fill_particleId(pdgIdMother, IdSummary[62]);
            break;
        case -10211: // is a a0 (1450)-
            Fill_particleId(pdgIdMother, IdSummary[63]);
            break;
        case 200111: // is a π(1800)0
            Fill_particleId(pdgIdMother, IdSummary[64]);
            break;
        case -200111: // is an anti π(1800)0
            Fill_particleId(pdgIdMother, IdSummary[65]);
            break;
        case 200211: // is a π(1800)+
            Fill_particleId(pdgIdMother, IdSummary[66]);
            break;
        case -200211: // is a π(1800)-
            Fill_particleId(pdgIdMother, IdSummary[67]);
            break;
        case 113: // is a ρ(770)0
            Fill_particleId(pdgIdMother, IdSummary[68]);
            break;
        case -113: // is an anti ρ(770)0
            Fill_particleId(pdgIdMother, IdSummary[69]);
            break;
        case 213: // is a ρ(770)+
            Fill_particleId(pdgIdMother, IdSummary[70]);
            break;
        case -213: // is a ρ(770)-
            Fill_particleId(pdgIdMother, IdSummary[71]);
            break;
        case 10113: // is a b1 (1235)0
            Fill_particleId(pdgIdMother, IdSummary[72]);
            break;
        case -10113: // is an anti b1 (1235)0
            Fill_particleId(pdgIdMother, IdSummary[73]);
            break;
        case 10213: // is a b1 (1235)+
            Fill_particleId(pdgIdMother, IdSummary[74]);
            break;
        case -10213: // is a b1 (1235)-
            Fill_particleId(pdgIdMother, IdSummary[75]);
            break;
        case 20113: // is a a1 (1260)0
            Fill_particleId(pdgIdMother, IdSummary[76]);
            break;
        case -20113: // is an anti a1 (1260)0
            Fill_particleId(pdgIdMother, IdSummary[77]);
            break;
        case 20213: // is a a1 (1260)+
            Fill_particleId(pdgIdMother, IdSummary[78]);
            break;
        case -20213: // is a a1 (1260)-
            Fill_particleId(pdgIdMother, IdSummary[79]);
            break;
        case 9000113: // is a π1 (1400)0
            Fill_particleId(pdgIdMother, IdSummary[80]);
            break;
        case -9000113: // is an anti π1 (1400)0
            Fill_particleId(pdgIdMother, IdSummary[81]);
            break;
        case 9000213: // is a π1 (1400)+
            Fill_particleId(pdgIdMother, IdSummary[82]);
            break;
        case -9000213: // is a π1 (1400)-
            Fill_particleId(pdgIdMother, IdSummary[83]);
            break;
        case 100113: // is a ρ(1450)0
            Fill_particleId(pdgIdMother, IdSummary[84]);
            break;
        case -100113: // is an anti ρ(1450)0
            Fill_particleId(pdgIdMother, IdSummary[85]);
            break;
        case 100213: // is a ρ(1450)+
            Fill_particleId(pdgIdMother, IdSummary[86]);
            break;
        case -100213: // is a ρ(1450)-
            Fill_particleId(pdgIdMother, IdSummary[87]);
            break;
        case 9010113: // is a π1 (1600)0
            Fill_particleId(pdgIdMother, IdSummary[88]);
            break;
        case -9010113: // is an anti π1 (1600)0
            Fill_particleId(pdgIdMother, IdSummary[89]);
            break;
        case 9010213: // is a π1 (1600)+
            Fill_particleId(pdgIdMother, IdSummary[90]);
            break;
        case -9010213: // is a π1 (1600)-
            Fill_particleId(pdgIdMother, IdSummary[91]);
            break;
        case 9020113: // is a a1 (1640)0
            Fill_particleId(pdgIdMother, IdSummary[92]);
            break;
        case -9020113: // is an anti a1 (1640)0
            Fill_particleId(pdgIdMother, IdSummary[93]);
            break;
        case 9020213: // is a a1 (1640)+
            Fill_particleId(pdgIdMother, IdSummary[94]);
            break;
        case -9020213: // is a a1 (1640)-
            Fill_particleId(pdgIdMother, IdSummary[95]);
            break;
        case 30113: // is a ρ(1700)0
            Fill_particleId(pdgIdMother, IdSummary[96]);
            break;
        case -30113: // is an anti ρ(1700)0
            Fill_particleId(pdgIdMother, IdSummary[97]);
            break;
        case 30213: // is a ρ(1700)+
            Fill_particleId(pdgIdMother, IdSummary[98]);
            break;
        case -30213: // is a ρ(1700)-
            Fill_particleId(pdgIdMother, IdSummary[99]);
            break;
        case 9030113: // is a ρ(1900)0
            Fill_particleId(pdgIdMother, IdSummary[100]);
            break;
        case -9030113: // is an anti ρ(1900)0
            Fill_particleId(pdgIdMother, IdSummary[101]);
            break;
        case 9030213: // is a ρ(1900)+
            Fill_particleId(pdgIdMother, IdSummary[102]);
            break;
        case -9030213: // is a ρ(1900)-
            Fill_particleId(pdgIdMother, IdSummary[103]);
            break;
        case 9040113: // is a ρ(2150)0
            Fill_particleId(pdgIdMother, IdSummary[104]);
            break;
        case -9040113: // is an anti ρ(2150)0
            Fill_particleId(pdgIdMother, IdSummary[105]);
            break;
        case 9040213: // is a ρ(2150)+
            Fill_particleId(pdgIdMother, IdSummary[106]);
            break;
        case -9040213: // is a ρ(2150)-
            Fill_particleId(pdgIdMother, IdSummary[107]);
            break;
        case 115: // is a a2 (1320)0
            Fill_particleId(pdgIdMother, IdSummary[108]);
            break;
        case -115: // is an anti a2 (1320)0
            Fill_particleId(pdgIdMother, IdSummary[109]);
            break;
        case 215: // is a a2 (1320)+
            Fill_particleId(pdgIdMother, IdSummary[110]);
            break;
        case -215: // is a a2 (1320)-
            Fill_particleId(pdgIdMother, IdSummary[111]);
            break;
        case 10115: // is a π2 (1670)0
            Fill_particleId(pdgIdMother, IdSummary[112]);
            break;
        case -10115: // is an anti π2 (1670)0
            Fill_particleId(pdgIdMother, IdSummary[113]);
            break;
        case 10215: // is a π2 (1670)+
            Fill_particleId(pdgIdMother, IdSummary[114]);
            break;
        case -10215: // is a π2 (1670)-
            Fill_particleId(pdgIdMother, IdSummary[115]);
            break;
        case 100115: // is a a2(1700)0
            Fill_particleId(pdgIdMother, IdSummary[116]);
            break;
        case -100115: // is an anti a2(1700)0
            Fill_particleId(pdgIdMother, IdSummary[117]);
            break;
        case 100215: // is a a2(1700)+
            Fill_particleId(pdgIdMother, IdSummary[118]);
            break;
        case -100215: // is a a2(1700)-
            Fill_particleId(pdgIdMother, IdSummary[119]);
            break;
        case 9000115: // is a π2 (2100)0
            Fill_particleId(pdgIdMother, IdSummary[120]);
            break;
        case -9000115: // is an anti π2 (2100)0
            Fill_particleId(pdgIdMother, IdSummary[121]);
            break;
        case 9000215: // is a π2 (2100)+
            Fill_particleId(pdgIdMother, IdSummary[122]);
            break;
        case -9000215: // is a π2 (2100)-
            Fill_particleId(pdgIdMother, IdSummary[123]);
            break;
        case 117: // is a ρ3 (1690)0
            Fill_particleId(pdgIdMother, IdSummary[124]);
            break;
        case -117: // is an anti ρ3 (1690)0
            Fill_particleId(pdgIdMother, IdSummary[125]);
            break;
        case 217: // is a ρ3 (1690)+
            Fill_particleId(pdgIdMother, IdSummary[126]);
            break;
        case -217: // is a ρ3 (1690)-
            Fill_particleId(pdgIdMother, IdSummary[127]);
            break;
        case 9000117: // is a ρ3(1990)0
            Fill_particleId(pdgIdMother, IdSummary[128]);
            break;
        case -9000117: // is an anti ρ3(1990)0
            Fill_particleId(pdgIdMother, IdSummary[129]);
            break;
        case 9000217: // is a ρ3(1990)+
            Fill_particleId(pdgIdMother, IdSummary[130]);
            break;
        case -9000217: // is a ρ3(1990)-
            Fill_particleId(pdgIdMother, IdSummary[131]);
            break;
        case 9010117: // is a ρ3 (2250)0
            Fill_particleId(pdgIdMother, IdSummary[132]);
            break;
        case -9010117: // is an anti ρ3 (2250)0
            Fill_particleId(pdgIdMother, IdSummary[133]);
            break;
        case 9010217: // is a ρ3 (2250)+
            Fill_particleId(pdgIdMother, IdSummary[134]);
            break;
        case -9010217: // is a ρ3 (2250)-
            Fill_particleId(pdgIdMother, IdSummary[135]);
            break;
        case 119: // is a a4 (2040)0
            Fill_particleId(pdgIdMother, IdSummary[136]);
            break;
        case -119: // is an anti a4 (2040)0
            Fill_particleId(pdgIdMother, IdSummary[137]);
            break;
        case 219: // is a a4 (2040)+
            Fill_particleId(pdgIdMother, IdSummary[138]);
            break;
        case -219: // is a a4 (2040)-
            Fill_particleId(pdgIdMother, IdSummary[139]);
            break;
            // LIGHT I=0 MESONS [45]
            // (u\bar{u}, d\bar{d}, and s\bar{s} Admixtures)
        case 221: // is a η
            Fill_particleId(pdgIdMother, IdSummary[140]);
            break;
        case 331: // is a η′(958)
            Fill_particleId(pdgIdMother, IdSummary[141]);
            break;
        case 9000221: // is a f0 (600)
            Fill_particleId(pdgIdMother, IdSummary[142]);
            break;
        case 9010221: // is a f0 (980)
            Fill_particleId(pdgIdMother, IdSummary[143]);
            break;
        case 100221: // is a η(1295)
            Fill_particleId(pdgIdMother, IdSummary[144]);
            break;
        case 10221: // is a f0 (1370)
            Fill_particleId(pdgIdMother, IdSummary[145]);
            break;
        case 100331: // is a η(1440)
            Fill_particleId(pdgIdMother, IdSummary[146]);
            break;
        case 9020221: // is a f0 (1500)
            Fill_particleId(pdgIdMother, IdSummary[147]);
            break;
        case 10331: // is a f0 (1710)
            Fill_particleId(pdgIdMother, IdSummary[148]);
            break;
        case 200221: // is a η(1760)
            Fill_particleId(pdgIdMother, IdSummary[149]);
            break;
        case 9030221: // is a f0 (2020)
            Fill_particleId(pdgIdMother, IdSummary[150]);
            break;
        case 9040221: // is a f0 (2100)
            Fill_particleId(pdgIdMother, IdSummary[151]);
            break;
        case 9050221: // is a f0 (2200)
            Fill_particleId(pdgIdMother, IdSummary[152]);
            break;
        case 9060221: // is a η(2225)
            Fill_particleId(pdgIdMother, IdSummary[153]);
            break;
        case 9070221: // is a f0 (2330)
            Fill_particleId(pdgIdMother, IdSummary[154]);
            break;
        case 223: // is a ω(782)
            Fill_particleId(pdgIdMother, IdSummary[155]);
            break;
        case 333: // is a φ(1020)
            Fill_particleId(pdgIdMother, IdSummary[156]);
            break;
        case 10223: // is a h1 (1170)
            Fill_particleId(pdgIdMother, IdSummary[157]);
            break;
        case 20223: // is a f1 (1285)
            Fill_particleId(pdgIdMother, IdSummary[158]);
            break;
        case 10333: // is a h1 (1380)
            Fill_particleId(pdgIdMother, IdSummary[159]);
            break;
        case 20333: // is a f1 (1420)
            Fill_particleId(pdgIdMother, IdSummary[160]);
            break;
        case 100223: // is a ω(1420)
            Fill_particleId(pdgIdMother, IdSummary[161]);
            break;
        case 9000223: // is a f1 (1510)
            Fill_particleId(pdgIdMother, IdSummary[162]);
            break;
        case 9010223: // is a h1 (1595)
            Fill_particleId(pdgIdMother, IdSummary[163]);
            break;
        case 30223: // is a ω(1650)
            Fill_particleId(pdgIdMother, IdSummary[164]);
            break;
        case 100333: // is a φ(1680)
            Fill_particleId(pdgIdMother, IdSummary[165]);
            break;
        case 225: // is a f2 (1270)
            Fill_particleId(pdgIdMother, IdSummary[166]);
            break;
        case 9000225: // is a f2 (1430)
            Fill_particleId(pdgIdMother, IdSummary[167]);
            break;
        case 335: // is a f2′ (1525)
            Fill_particleId(pdgIdMother, IdSummary[168]);
            break;
        case 9010225: // is a f2 (1565)
            Fill_particleId(pdgIdMother, IdSummary[169]);
            break;
        case 9020225: // is a f2 (1640)
            Fill_particleId(pdgIdMother, IdSummary[170]);
            break;
        case 10225: // is a η2 (1645)
            Fill_particleId(pdgIdMother, IdSummary[171]);
            break;
        case 9030225: // is a f2 (1810)
            Fill_particleId(pdgIdMother, IdSummary[172]);
            break;
        case 10335: // is a η2 (1870)
            Fill_particleId(pdgIdMother, IdSummary[173]);
            break;
        case 9040225: // is a f2 (1910)
            Fill_particleId(pdgIdMother, IdSummary[174]);
            break;
        case 100225: // is a f2 (1950)
            Fill_particleId(pdgIdMother, IdSummary[175]);
            break;
        case 100335: // is a f2 (2010)
            Fill_particleId(pdgIdMother, IdSummary[176]);
            break;
        case 9050225: // is a f2 (2150)
            Fill_particleId(pdgIdMother, IdSummary[177]);
            break;
        case 9060225: // is a f2 (2300)
            Fill_particleId(pdgIdMother, IdSummary[178]);
            break;
        case 9070225: // is a f2 (2340)
            Fill_particleId(pdgIdMother, IdSummary[179]);
            break;
        case 227: // is a ω3 (1670)
            Fill_particleId(pdgIdMother, IdSummary[180]);
            break;
        case 337: // is a φ3 (1850)
            Fill_particleId(pdgIdMother, IdSummary[181]);
            break;
        case 229: // is a f4 (2050)
            Fill_particleId(pdgIdMother, IdSummary[182]);
            break;
        case 9000339: // is a fJ (2220)
            Fill_particleId(pdgIdMother, IdSummary[183]);
            break;
        case 9000229: // is a f4 (2300)
            Fill_particleId(pdgIdMother, IdSummary[184]);
            break;
            // STRANGE MESONS [88]
        case 130: // is a KL0
            Fill_particleId(pdgIdMother, IdSummary[185]);
            break;
        case -130: // is an anti KL0
            Fill_particleId(pdgIdMother, IdSummary[186]);
            break;
        case 310: // is a KS0
            Fill_particleId(pdgIdMother, IdSummary[187]);
            break;
        case -310: // is an anti KS0
            Fill_particleId(pdgIdMother, IdSummary[188]);
            break;
        case 311: // is a K0
            Fill_particleId(pdgIdMother, IdSummary[189]);
            break;
        case -311: // is an anti K0
            Fill_particleId(pdgIdMother, IdSummary[190]);
            break;
        case 321: // is a K+
            Fill_particleId(pdgIdMother, IdSummary[191]);
            break;
        case -321: // is a K-
            Fill_particleId(pdgIdMother, IdSummary[192]);
            break;
        case 10311: // is a K0∗ (1430)0
            Fill_particleId(pdgIdMother, IdSummary[193]);
            break;
        case -10311: // is an anti K0∗ (1430)0
            Fill_particleId(pdgIdMother, IdSummary[194]);
            break;
        case 10321: // is a K0∗ (1430)+
            Fill_particleId(pdgIdMother, IdSummary[195]);
            break;
        case -10321: // is a K0∗ (1430)-
            Fill_particleId(pdgIdMother, IdSummary[196]);
            break;
        case 100311: // is a K (1460)0
            Fill_particleId(pdgIdMother, IdSummary[197]);
            break;
        case -100311: // is an anti K (1460)0
            Fill_particleId(pdgIdMother, IdSummary[198]);
            break;
        case 100321: // is a K (1460)+
            Fill_particleId(pdgIdMother, IdSummary[199]);
            break;
        case -100321: // is a K (1460)-
            Fill_particleId(pdgIdMother, IdSummary[200]);
            break;
        case 200311: // is a K (1830)0
            Fill_particleId(pdgIdMother, IdSummary[201]);
            break;
        case -200311: // is an anti K (1830)0
            Fill_particleId(pdgIdMother, IdSummary[202]);
            break;
        case 200321: // is a K (1830)+
            Fill_particleId(pdgIdMother, IdSummary[203]);
            break;
        case -200321: // is a K (1830)-
            Fill_particleId(pdgIdMother, IdSummary[204]);
            break;
        case 9000311: // is a K0∗ (1950)0
            Fill_particleId(pdgIdMother, IdSummary[205]);
            break;
        case -9000311: // is an anti K0∗ (1950)0
            Fill_particleId(pdgIdMother, IdSummary[206]);
            break;
        case 9000321: // is a K0∗ (1950)+
            Fill_particleId(pdgIdMother, IdSummary[207]);
            break;
        case -9000321: // is a K0∗ (1950)-
            Fill_particleId(pdgIdMother, IdSummary[208]);
            break;
        case 313: // is a K ∗ (892)0
            Fill_particleId(pdgIdMother, IdSummary[209]);
            break;
        case -313: // is an anti K ∗ (892)0
            Fill_particleId(pdgIdMother, IdSummary[210]);
            break;
        case 323: // is a K ∗ (892)+
            Fill_particleId(pdgIdMother, IdSummary[211]);
            break;
        case -323: // is a K ∗ (892)-
            Fill_particleId(pdgIdMother, IdSummary[212]);
            break;
        case 10313: // is a K1 (1270)0
            Fill_particleId(pdgIdMother, IdSummary[213]);
            break;
        case -10313: // is an anti K1 (1270)0
            Fill_particleId(pdgIdMother, IdSummary[214]);
            break;
        case 10323: // is a K1 (1270)+
            Fill_particleId(pdgIdMother, IdSummary[215]);
            break;
        case -10323: // is a K1 (1270)-
            Fill_particleId(pdgIdMother, IdSummary[216]);
            break;
        case 20313: // is a K1 (1400)0
            Fill_particleId(pdgIdMother, IdSummary[217]);
            break;
        case -20313: // is an anti K1 (1400)0
            Fill_particleId(pdgIdMother, IdSummary[218]);
            break;
        case 20323: // is a K1 (1400)+
            Fill_particleId(pdgIdMother, IdSummary[219]);
            break;
        case -20323: // is a K1 (1400)-
            Fill_particleId(pdgIdMother, IdSummary[220]);
            break;
        case 100313: // is a K ∗ (1410)0
            Fill_particleId(pdgIdMother, IdSummary[221]);
            break;
        case -100313: // is an anti K ∗ (1410)0
            Fill_particleId(pdgIdMother, IdSummary[222]);
            break;
        case 100323: // is a K ∗ (1410)+
            Fill_particleId(pdgIdMother, IdSummary[223]);
            break;
        case -100323: // is a K ∗ (1410)-
            Fill_particleId(pdgIdMother, IdSummary[224]);
            break;
        case 9000313: // is a K1 (1650)0
            Fill_particleId(pdgIdMother, IdSummary[225]);
            break;
        case -9000313: // is an anti K1 (1650)0
            Fill_particleId(pdgIdMother, IdSummary[226]);
            break;
        case 9000323: // is a K1 (1650)+
            Fill_particleId(pdgIdMother, IdSummary[227]);
            break;
        case -9000323: // is a K1 (1650)-
            Fill_particleId(pdgIdMother, IdSummary[228]);
            break;
        case 30313: // is a K ∗ (1680)0
            Fill_particleId(pdgIdMother, IdSummary[229]);
            break;
        case -30313: // is an anti K ∗ (1680)0
            Fill_particleId(pdgIdMother, IdSummary[230]);
            break;
        case 30323: // is a K ∗ (1680)+
            Fill_particleId(pdgIdMother, IdSummary[231]);
            break;
        case -30323: // is a K ∗ (1680)-
            Fill_particleId(pdgIdMother, IdSummary[232]);
            break;
        case 315: // is a K2∗ (1430)0
            Fill_particleId(pdgIdMother, IdSummary[233]);
            break;
        case -315: // is an anti K2∗ (1430)0
            Fill_particleId(pdgIdMother, IdSummary[234]);
            break;
        case 325: // is a K2∗ (1430)+
            Fill_particleId(pdgIdMother, IdSummary[235]);
            break;
        case -325: // is a K2∗ (1430)-
            Fill_particleId(pdgIdMother, IdSummary[236]);
            break;
        case 9000315: // is a K2 (1580)0
            Fill_particleId(pdgIdMother, IdSummary[237]);
            break;
        case -9000315: // is an anti K2 (1580)0
            Fill_particleId(pdgIdMother, IdSummary[238]);
            break;
        case 9000325: // is a K2 (1580)+
            Fill_particleId(pdgIdMother, IdSummary[239]);
            break;
        case -9000325: // is a K2 (1580)-
            Fill_particleId(pdgIdMother, IdSummary[240]);
            break;
        case 10315: // is a K2 (1770)0
            Fill_particleId(pdgIdMother, IdSummary[241]);
            break;
        case -10315: // is an anti K2 (1770)0
            Fill_particleId(pdgIdMother, IdSummary[242]);
            break;
        case 10325: // is a K2 (1770)+
            Fill_particleId(pdgIdMother, IdSummary[243]);
            break;
        case -10325: // is a K2 (1770)-
            Fill_particleId(pdgIdMother, IdSummary[244]);
            break;
        case 20315: // is a K2 (1820)0
            Fill_particleId(pdgIdMother, IdSummary[245]);
            break;
        case -20315: // is an anti K2 (1820)0
            Fill_particleId(pdgIdMother, IdSummary[246]);
            break;
        case 20325: // is a K2 (1820)+
            Fill_particleId(pdgIdMother, IdSummary[247]);
            break;
        case -20325: // is a K2 (1820)-
            Fill_particleId(pdgIdMother, IdSummary[248]);
            break;
        case 100315: // is a K2∗ (1980)0
            Fill_particleId(pdgIdMother, IdSummary[249]);
            break;
        case -100315: // is an anti K2∗ (1980)0
            Fill_particleId(pdgIdMother, IdSummary[250]);
            break;
        case 100325: // is a K2∗ (1980)+
            Fill_particleId(pdgIdMother, IdSummary[251]);
            break;
        case -100325: // is a K2∗ (1980)-
            Fill_particleId(pdgIdMother, IdSummary[252]);
            break;
        case 9010315: // is a K2 (2250)0
            Fill_particleId(pdgIdMother, IdSummary[253]);
            break;
        case -9010315: // is an anti K2 (2250)0
            Fill_particleId(pdgIdMother, IdSummary[254]);
            break;
        case 9010325: // is a K2 (2250)+
            Fill_particleId(pdgIdMother, IdSummary[255]);
            break;
        case -9010325: // is a K2 (2250)-
            Fill_particleId(pdgIdMother, IdSummary[256]);
            break;
        case 317: // is a K3∗ (1780)0
            Fill_particleId(pdgIdMother, IdSummary[257]);
            break;
        case -317: // is an anti K3∗ (1780)0
            Fill_particleId(pdgIdMother, IdSummary[258]);
            break;
        case 327: // is a K3∗ (1780)+
            Fill_particleId(pdgIdMother, IdSummary[259]);
            break;
        case -327: // is a K3∗ (1780)-
            Fill_particleId(pdgIdMother, IdSummary[260]);
            break;
        case 9010317: // is a K3 (2320)0
            Fill_particleId(pdgIdMother, IdSummary[261]);
            break;
        case -9010317: // is an anti K3 (2320)0
            Fill_particleId(pdgIdMother, IdSummary[262]);
            break;
        case 9010327: // is a K3 (2320)+
            Fill_particleId(pdgIdMother, IdSummary[263]);
            break;
        case -9010327: // is a K3 (2320)-
            Fill_particleId(pdgIdMother, IdSummary[264]);
            break;
        case 319: // is a K4∗ (2045)0
            Fill_particleId(pdgIdMother, IdSummary[265]);
            break;
        case -319: // is an anti K4∗ (2045)0
            Fill_particleId(pdgIdMother, IdSummary[266]);
            break;
        case 329: // is a K4∗ (2045)+
            Fill_particleId(pdgIdMother, IdSummary[267]);
            break;
        case -329: // is a K4∗ (2045)-
            Fill_particleId(pdgIdMother, IdSummary[268]);
            break;
        case 9000319: // is a K4 (2500)0
            Fill_particleId(pdgIdMother, IdSummary[269]);
            break;
        case -9000319: // is an anti K4 (2500)0
            Fill_particleId(pdgIdMother, IdSummary[270]);
            break;
        case 9000329: // is a K4 (2500)+
            Fill_particleId(pdgIdMother, IdSummary[271]);
            break;
        case -9000329: // is a K4 (2500)-
            Fill_particleId(pdgIdMother, IdSummary[272]);
            break;
            // CHARMED MESONS [36]
        case 411: // is a D+
            Fill_particleId(pdgIdMother, IdSummary[273]);
            break;
        case -411: // is a D-
            Fill_particleId(pdgIdMother, IdSummary[274]);
            break;
        case 421: // is a D0
            Fill_particleId(pdgIdMother, IdSummary[275]);
            break;
        case -421: // is an anti D0
            Fill_particleId(pdgIdMother, IdSummary[276]);
            break;
        case 10411: // is a D0∗+
            Fill_particleId(pdgIdMother, IdSummary[277]);
            break;
        case -10411: // is a D0∗-
            Fill_particleId(pdgIdMother, IdSummary[278]);
            break;
        case 10421: // is a D0∗0
            Fill_particleId(pdgIdMother, IdSummary[279]);
            break;
        case -10421: // is an anti D0∗0
            Fill_particleId(pdgIdMother, IdSummary[280]);
            break;
        case 413: // is a D∗(2010)+
            Fill_particleId(pdgIdMother, IdSummary[281]);
            break;
        case -413: // is a D∗(2010)-
            Fill_particleId(pdgIdMother, IdSummary[282]);
            break;
        case 423: // is a D∗(2007)0
            Fill_particleId(pdgIdMother, IdSummary[283]);
            break;
        case -423: // is an anti D∗(2007)0
            Fill_particleId(pdgIdMother, IdSummary[284]);
            break;
        case 10413: // is a D1 (2420)+
            Fill_particleId(pdgIdMother, IdSummary[285]);
            break;
        case -10413: // is a D1 (2420)-
            Fill_particleId(pdgIdMother, IdSummary[286]);
            break;
        case 10423: // is a D1(2420)0
            Fill_particleId(pdgIdMother, IdSummary[287]);
            break;
        case -10423: // is an anti D1(2420)0
            Fill_particleId(pdgIdMother, IdSummary[288]);
            break;
        case 20413: // is a D1(H)+
            Fill_particleId(pdgIdMother, IdSummary[289]);
            break;
        case -20413: // is a D1(H)-
            Fill_particleId(pdgIdMother, IdSummary[290]);
            break;
        case 20423: // is a D1(H)0
            Fill_particleId(pdgIdMother, IdSummary[291]);
            break;
        case -20423: // is an anti D1(H)0
            Fill_particleId(pdgIdMother, IdSummary[292]);
            break;
        case 415: // is a D2∗(2460)+
            Fill_particleId(pdgIdMother, IdSummary[293]);
            break;
        case -415: // is a D2∗(2460)-
            Fill_particleId(pdgIdMother, IdSummary[294]);
            break;
        case 425: // is a D2∗(2460)0
            Fill_particleId(pdgIdMother, IdSummary[295]);
            break;
        case -425: // is an anti D2∗(2460)0
            Fill_particleId(pdgIdMother, IdSummary[296]);
            break;
        case 431: // is a Ds+
            Fill_particleId(pdgIdMother, IdSummary[297]);
            break;
        case -431: // is a Ds-
            Fill_particleId(pdgIdMother, IdSummary[298]);
            break;
        case 10431: // is a Ds0∗+
            Fill_particleId(pdgIdMother, IdSummary[299]);
            break;
        case -10431: // is a Ds0∗-
            Fill_particleId(pdgIdMother, IdSummary[300]);
            break;
        case 433: // is a Ds∗+
            Fill_particleId(pdgIdMother, IdSummary[301]);
            break;
        case -433: // is a Ds∗-
            Fill_particleId(pdgIdMother, IdSummary[302]);
            break;
        case 10433: // is a Ds1 (2536)+
            Fill_particleId(pdgIdMother, IdSummary[303]);
            break;
        case -10433: // is a Ds1 (2536)-
            Fill_particleId(pdgIdMother, IdSummary[304]);
            break;
        case 20433: // is a Ds1(H)+
            Fill_particleId(pdgIdMother, IdSummary[305]);
            break;
        case -20433: // is a Ds1(H)-
            Fill_particleId(pdgIdMother, IdSummary[306]);
            break;
        case 435: // is a Ds2*+
            Fill_particleId(pdgIdMother, IdSummary[307]);
            break;
        case -435: // is a Ds2*-
            Fill_particleId(pdgIdMother, IdSummary[308]);
            break;
            // BOTTOM MESONS [48]
        case 511: // is a B0
            Fill_particleId(pdgIdMother, IdSummary[309]);
            break;
        case -511: // is an anti B0
            Fill_particleId(pdgIdMother, IdSummary[310]);
            break;
        case 521: // is a B+
            Fill_particleId(pdgIdMother, IdSummary[311]);
            break;
        case -521: // is a B-
            Fill_particleId(pdgIdMother, IdSummary[312]);
            break;
        case 10511: // is a B0∗0
            Fill_particleId(pdgIdMother, IdSummary[313]);
            break;
        case -10511: // is an anti B0∗0
            Fill_particleId(pdgIdMother, IdSummary[314]);
            break;
        case 10521: // is a B0∗+
            Fill_particleId(pdgIdMother, IdSummary[315]);
            break;
        case -10521: // is a B0∗-
            Fill_particleId(pdgIdMother, IdSummary[316]);
            break;
        case 513: // is a B∗0
            Fill_particleId(pdgIdMother, IdSummary[317]);
            break;
        case -513: // is an anti B∗0
            Fill_particleId(pdgIdMother, IdSummary[318]);
            break;
        case 523: // is a B∗+
            Fill_particleId(pdgIdMother, IdSummary[319]);
            break;
        case -523: // is a B∗-
            Fill_particleId(pdgIdMother, IdSummary[320]);
            break;
        case 10513: // is a B1 (L)0
            Fill_particleId(pdgIdMother, IdSummary[321]);
            break;
        case -10513: // is an anti B1 (L)0
            Fill_particleId(pdgIdMother, IdSummary[322]);
            break;
        case 10523: // is a B1(L)+
            Fill_particleId(pdgIdMother, IdSummary[323]);
            break;
        case -10523: // is a B1(L)-
            Fill_particleId(pdgIdMother, IdSummary[324]);
            break;
        case 20513: // is a B1(H)0
            Fill_particleId(pdgIdMother, IdSummary[325]);
            break;
        case -20513: // is an anti B1(H)0
            Fill_particleId(pdgIdMother, IdSummary[326]);
            break;
        case 20523: // is a B1(H)+
            Fill_particleId(pdgIdMother, IdSummary[327]);
            break;
        case -20523: // is a B1(H)-
            Fill_particleId(pdgIdMother, IdSummary[328]);
            break;
        case 515: // is a B2∗0
            Fill_particleId(pdgIdMother, IdSummary[329]);
            break;
        case -515: // is an anti B2∗0
            Fill_particleId(pdgIdMother, IdSummary[330]);
            break;
        case 525: // is a B2∗+
            Fill_particleId(pdgIdMother, IdSummary[331]);
            break;
        case -525: // is a B2∗-
            Fill_particleId(pdgIdMother, IdSummary[332]);
            break;
        case 531: // is a Bs0
            Fill_particleId(pdgIdMother, IdSummary[333]);
            break;
        case -531: // is an anti Bs0
            Fill_particleId(pdgIdMother, IdSummary[334]);
            break;
        case 10531: // is a Bs0*0
            Fill_particleId(pdgIdMother, IdSummary[335]);
            break;
        case -10531: // is an anti Bs0*0
            Fill_particleId(pdgIdMother, IdSummary[336]);
            break;
        case 533: // is a Bs∗0
            Fill_particleId(pdgIdMother, IdSummary[337]);
            break;
        case -533: // is an anti Bs∗0
            Fill_particleId(pdgIdMother, IdSummary[338]);
            break;
        case 10533: // is a Bs1(L)0
            Fill_particleId(pdgIdMother, IdSummary[339]);
            break;
        case -10533: // is an anti Bs1(L)0
            Fill_particleId(pdgIdMother, IdSummary[340]);
            break;
        case 20533: // is a Bs1(H)0
            Fill_particleId(pdgIdMother, IdSummary[341]);
            break;
        case -20533: // is an anti Bs1(H)0
            Fill_particleId(pdgIdMother, IdSummary[342]);
            break;
        case 535: // is a Bs2∗0
            Fill_particleId(pdgIdMother, IdSummary[343]);
            break;
        case -535: // is an anti Bs2∗0
            Fill_particleId(pdgIdMother, IdSummary[344]);
            break;
        case 541: // is a Bc+
            Fill_particleId(pdgIdMother, IdSummary[345]);
            break;
        case -541: // is a Bc-
            Fill_particleId(pdgIdMother, IdSummary[346]);
            break;
        case 10541: // is a Bc0∗+
            Fill_particleId(pdgIdMother, IdSummary[347]);
            break;
        case -10541: // is a Bc0∗-
            Fill_particleId(pdgIdMother, IdSummary[348]);
            break;
        case 543: // is a Bc∗+
            Fill_particleId(pdgIdMother, IdSummary[349]);
            break;
        case -543: // is a Bc∗-
            Fill_particleId(pdgIdMother, IdSummary[350]);
            break;
        case 10543: // is a Bc1(L)+
            Fill_particleId(pdgIdMother, IdSummary[351]);
            break;
        case -10543: // is a Bc1(L)-
            Fill_particleId(pdgIdMother, IdSummary[352]);
            break;
        case 20543: // is a Bc1(H)+
            Fill_particleId(pdgIdMother, IdSummary[353]);
            break;
        case -20543: // is a Bc1(H)-
            Fill_particleId(pdgIdMother, IdSummary[354]);
            break;
        case 545: // is a Bc2∗+
            Fill_particleId(pdgIdMother, IdSummary[355]);
            break;
        case -545: // is a Bc2∗-
            Fill_particleId(pdgIdMother, IdSummary[356]);
            break;
            // c\bar{c} MESONS [13]
        case 441: // is a ηc(1S)
            Fill_particleId(pdgIdMother, IdSummary[357]);
            break;
        case 10441: // is a χc0(1P)
            Fill_particleId(pdgIdMother, IdSummary[358]);
            break;
        case 100441: // is a ηc (2S)
            Fill_particleId(pdgIdMother, IdSummary[359]);
            break;
        case 443: // is a J/ψ(1S)
            Fill_particleId(pdgIdMother, IdSummary[360]);
            break;
        case 10443: // is a hc(1P)
            Fill_particleId(pdgIdMother, IdSummary[361]);
            break;
        case 20443: // is a χc1(1P)
            Fill_particleId(pdgIdMother, IdSummary[362]);
            break;
        case 100443: // is a ψ(2S)
            Fill_particleId(pdgIdMother, IdSummary[363]);
            break;
        case 30443: // is a ψ(3770)
            Fill_particleId(pdgIdMother, IdSummary[364]);
            break;
        case 9000443: // is a ψ(4040)
            Fill_particleId(pdgIdMother, IdSummary[365]);
            break;
        case 9010443: // is a ψ(4160)
            Fill_particleId(pdgIdMother, IdSummary[366]);
            break;
        case 9020443: // is a ψ(4415)
            Fill_particleId(pdgIdMother, IdSummary[367]);
            break;
        case 445: // is a χc2(1P)
            Fill_particleId(pdgIdMother, IdSummary[368]);
            break;
        case 9000445: // is a ψ(3836)
            Fill_particleId(pdgIdMother, IdSummary[369]);
            break;
            //b\bar{b} MESONS [29]
        case 551: // is a ηb (1S)
            Fill_particleId(pdgIdMother, IdSummary[370]);
            break;
        case 10551: // is a χb0(1P)
            Fill_particleId(pdgIdMother, IdSummary[371]);
            break;
        case 100551: // is a ηb (2S)
            Fill_particleId(pdgIdMother, IdSummary[372]);
            break;
        case 110551: // is a χb0(2P)
            Fill_particleId(pdgIdMother, IdSummary[373]);
            break;
        case 200551: // is a ηb (3S)
            Fill_particleId(pdgIdMother, IdSummary[374]);
            break;
        case 210551: // is a χb0(3P)
            Fill_particleId(pdgIdMother, IdSummary[375]);
            break;
        case 553: // is a Υ(1S)
            Fill_particleId(pdgIdMother, IdSummary[376]);
            break;
        case 10553: // is a hb(1P)
            Fill_particleId(pdgIdMother, IdSummary[377]);
            break;
        case 20553: // is a χb1(1P)
            Fill_particleId(pdgIdMother, IdSummary[378]);
            break;
        case 30553: // is a Υ1(1D)
            Fill_particleId(pdgIdMother, IdSummary[379]);
            break;
        case 100553: // is a Υ(2S)
            Fill_particleId(pdgIdMother, IdSummary[380]);
            break;
        case 110553: // is a hb(2P)
            Fill_particleId(pdgIdMother, IdSummary[381]);
            break;
        case 120553: // is a χb1(2P)
            Fill_particleId(pdgIdMother, IdSummary[382]);
            break;
        case 130553: // is a Υ1(2D)
            Fill_particleId(pdgIdMother, IdSummary[383]);
            break;
        case 200553: // is a Υ(3S)
            Fill_particleId(pdgIdMother, IdSummary[384]);
            break;
        case 210553: // is a hb(3P)
            Fill_particleId(pdgIdMother, IdSummary[385]);
            break;
        case 220553: // is a χb1(3P)
            Fill_particleId(pdgIdMother, IdSummary[386]);
            break;
        case 300553: // is a Υ(4S)
            Fill_particleId(pdgIdMother, IdSummary[387]);
            break;
        case 9000553: // is a Υ (10860)
            Fill_particleId(pdgIdMother, IdSummary[388]);
            break;
        case 9010553: // is a Υ (11020)
            Fill_particleId(pdgIdMother, IdSummary[389]);
            break;
        case 555: // is a χb2 (1P )
            Fill_particleId(pdgIdMother, IdSummary[390]);
            break;
        case 10555: // is a ηb2(1D)
            Fill_particleId(pdgIdMother, IdSummary[391]);
            break;
        case 20555: // is a Υ2(1D)
            Fill_particleId(pdgIdMother, IdSummary[392]);
            break;
        case 100555: // is a χb2(2P)
            Fill_particleId(pdgIdMother, IdSummary[393]);
            break;
        case 110555: // is a ηb2(2D)
            Fill_particleId(pdgIdMother, IdSummary[394]);
            break;
        case 120555: // is a Υ2(2D)
            Fill_particleId(pdgIdMother, IdSummary[395]);
            break;
        case 200555: // is a χb2(3P)
            Fill_particleId(pdgIdMother, IdSummary[396]);
            break;
        case 557: // is a Υ3(1D)
            Fill_particleId(pdgIdMother, IdSummary[397]);
            break;
        case 100557: // is a Υ3(2D)
            Fill_particleId(pdgIdMother, IdSummary[398]);
            break;
            // LIGHT BARYONS [12]
        case 2212: // is a proton
            Fill_particleId(pdgIdMother, IdSummary[399]);
            break;
        case -2212: // is an anti proton
            Fill_particleId(pdgIdMother, IdSummary[400]);
            break;
        case 2112: // is a neutron
            Fill_particleId(pdgIdMother, IdSummary[401]);
            break;
        case -2112: // is an anti neutron
            Fill_particleId(pdgIdMother, IdSummary[402]);
            break;
        case 2224: // is a Delta++
            Fill_particleId(pdgIdMother, IdSummary[403]);
            break;
        case -2224: // is an anti Delta++
            Fill_particleId(pdgIdMother, IdSummary[404]);
            break;
        case 2214: // is a Delta+
            Fill_particleId(pdgIdMother, IdSummary[405]);
            break;
        case -2214: // is an anti Delta+
            Fill_particleId(pdgIdMother, IdSummary[406]);
            break;
        case 2114: // is a Delta0
            Fill_particleId(pdgIdMother, IdSummary[407]);
            break;
        case -2114: // is an anti Delta0
            Fill_particleId(pdgIdMother, IdSummary[408]);
            break;
        case 1114: // is a Delta -
            Fill_particleId(pdgIdMother, IdSummary[409]);
            break;
        case -1114: // is an AntiDelta -
            Fill_particleId(pdgIdMother, IdSummary[410]);
            break;
            //STRANGE BARYONS [24]
        case 3122: // is a Λ
            Fill_particleId(pdgIdMother, IdSummary[411]);
            break;
        case -3122: // is an anti Λ
            Fill_particleId(pdgIdMother, IdSummary[412]);
            break;
        case 3222: // is a Σ+
            Fill_particleId(pdgIdMother, IdSummary[413]);
            break;
        case -3222: // is an anti Σ+
            Fill_particleId(pdgIdMother, IdSummary[414]);
            break;
        case 3212: // is a Σ0
            Fill_particleId(pdgIdMother, IdSummary[415]);
            break;
        case -3212: // is an anti Σ0
            Fill_particleId(pdgIdMother, IdSummary[416]);
            break;
        case 3112: // is a Σ-
            Fill_particleId(pdgIdMother, IdSummary[417]);
            break;
        case -3112: // is am anti Σ-
            Fill_particleId(pdgIdMother, IdSummary[418]);
            break;
        case 3224: // is a Σ*+
            Fill_particleId(pdgIdMother, IdSummary[419]);
            break;
        case -3224: // is an anti Σ*+
            Fill_particleId(pdgIdMother, IdSummary[420]);
            break;
        case 3214: // is a Σ*0
            Fill_particleId(pdgIdMother, IdSummary[421]);
            break;
        case -3214: // is an anti Σ*0
            Fill_particleId(pdgIdMother, IdSummary[422]);
            break;
        case 3114: // is a Σ*-
            Fill_particleId(pdgIdMother, IdSummary[423]);
            break;
        case -3114: // is an anti Σ*-
            Fill_particleId(pdgIdMother, IdSummary[424]);
            break;
        case 3322: // is a Ξ0
            Fill_particleId(pdgIdMother, IdSummary[425]);
            break;
        case -3322: // is an anti Ξ0
            Fill_particleId(pdgIdMother, IdSummary[426]);
            break;
        case 3312: // is a Ξ-
            Fill_particleId(pdgIdMother, IdSummary[427]);
            break;
        case -3312: // is an anti Ξ-
            Fill_particleId(pdgIdMother, IdSummary[428]);
            break;
        case 3324: // is a Ξ*0
            Fill_particleId(pdgIdMother, IdSummary[429]);
            break;
        case -3324: // is an anti Ξ*0
            Fill_particleId(pdgIdMother, IdSummary[430]);
            break;
        case 3314: // is a Ξ*-
            Fill_particleId(pdgIdMother, IdSummary[431]);
            break;
        case -3314: // is an anti Ξ*-
            Fill_particleId(pdgIdMother, IdSummary[432]);
            break;
        case 3334: // is a Ω−
            Fill_particleId(pdgIdMother, IdSummary[433]);
            break;
        case -3334: // is an anti Ω−
            Fill_particleId(pdgIdMother, IdSummary[434]);
            break;
            // CHARMED BARYONS [44]
        case 4122: // is a Lambda+c
            Fill_particleId(pdgIdMother, IdSummary[435]);
            break;
        case -4122: // is a Lambda-c
            Fill_particleId(pdgIdMother, IdSummary[436]);
            break;
        case 4222: // is a Σc++
            Fill_particleId(pdgIdMother, IdSummary[437]);
            break;
        case -4222: //is an anti Σc++
            Fill_particleId(pdgIdMother, IdSummary[438]);
            break;
        case 4212: // is a Σc+
            Fill_particleId(pdgIdMother, IdSummary[439]);
            break;
        case -4212: // is a Σc-
            Fill_particleId(pdgIdMother, IdSummary[440]);
            break;
        case 4112: // is a Σc0
            Fill_particleId(pdgIdMother, IdSummary[441]);
            break;
        case -4112: // is an anti Σc0
            Fill_particleId(pdgIdMother, IdSummary[442]);
            break;
        case 4224: // is a Σc∗++
            Fill_particleId(pdgIdMother, IdSummary[443]);
            break;
        case -4224: // is an anti Σc∗++
            Fill_particleId(pdgIdMother, IdSummary[444]);
            break;
        case 4214: // is a Σc∗+
            Fill_particleId(pdgIdMother, IdSummary[445]);
            break;
        case -4214: // is a Σc∗-
            Fill_particleId(pdgIdMother, IdSummary[446]);
            break;
        case 4114: // is a Σc∗0
            Fill_particleId(pdgIdMother, IdSummary[447]);
            break;
        case -4114: // is an anti Σc∗0
            Fill_particleId(pdgIdMother, IdSummary[448]);
            break;
        case 4232: // is Ξc+
            Fill_particleId(pdgIdMother, IdSummary[449]);
            break;
        case -4232: // is a Ξc-
            Fill_particleId(pdgIdMother, IdSummary[450]);
            break;
        case 4132: // is a Ξc0
            Fill_particleId(pdgIdMother, IdSummary[451]);
            break;
        case -4132: // is an anti Ξc0
            Fill_particleId(pdgIdMother, IdSummary[452]);
            break;
        case 4322: // is a Ξ′c+
            Fill_particleId(pdgIdMother, IdSummary[453]);
            break;
        case -4322: // is a Ξ′c-
            Fill_particleId(pdgIdMother, IdSummary[454]);
            break;
        case 4312: // is a Ξ′c0
            Fill_particleId(pdgIdMother, IdSummary[455]);
            break;
        case -4312: // is an anti Ξ′c0
            Fill_particleId(pdgIdMother, IdSummary[456]);
            break;
        case 4324: // is a Ξc∗+
            Fill_particleId(pdgIdMother, IdSummary[457]);
            break;
        case -4324: // is a Ξc∗-
            Fill_particleId(pdgIdMother, IdSummary[458]);
            break;
        case 4314: // is a Ξc∗0
            Fill_particleId(pdgIdMother, IdSummary[459]);
            break;
        case -4314: // is an anti Ξc∗0
            Fill_particleId(pdgIdMother, IdSummary[460]);
            break;
        case 4332: // is a Ωc0
            Fill_particleId(pdgIdMother, IdSummary[461]);
            break;
        case -4332: // is an anti Ωc0
            Fill_particleId(pdgIdMother, IdSummary[462]);
            break;
        case 4334: // is a Ω*c0
            Fill_particleId(pdgIdMother, IdSummary[463]);
            break;
        case -4334: // is an anti Ω*c0
            Fill_particleId(pdgIdMother, IdSummary[464]);
            break;
        case 4412: // is a Ξcc+
            Fill_particleId(pdgIdMother, IdSummary[465]);
            break;
        case -4412: // is a Ξcc-
            Fill_particleId(pdgIdMother, IdSummary[466]);
            break;
        case 4422: // is a Ξcc++
            Fill_particleId(pdgIdMother, IdSummary[467]);
            break;
        case -4422: // is a Ξcc--
            Fill_particleId(pdgIdMother, IdSummary[468]);
            break;
        case 4414: // is a Ξcc∗+
            Fill_particleId(pdgIdMother, IdSummary[469]);
            break;
        case -4414: // is a Ξcc∗-
            Fill_particleId(pdgIdMother, IdSummary[470]);
            break;
        case 4424: // is a Ξcc∗++
            Fill_particleId(pdgIdMother, IdSummary[471]);
            break;
        case -4424: // is a Ξcc∗--
            Fill_particleId(pdgIdMother, IdSummary[472]);
            break;
        case 4432: // is a Ωcc+
            Fill_particleId(pdgIdMother, IdSummary[473]);
            break;
        case -4432: // is a Ωcc-
            Fill_particleId(pdgIdMother, IdSummary[474]);
            break;
        case 4434: // is a Ωcc∗+
            Fill_particleId(pdgIdMother, IdSummary[475]);
            break;
        case -4434: // is a Ωcc∗-
            Fill_particleId(pdgIdMother, IdSummary[476]);
            break;
        case 4444: // is a Ωccc++
            Fill_particleId(pdgIdMother, IdSummary[477]);
            break;
        case -4444: // is a Ωccc--
            Fill_particleId(pdgIdMother, IdSummary[478]);
            break;
            //BOTTOM BARYONS [70]
        case 5122: // is a Λ0b
            Fill_particleId(pdgIdMother, IdSummary[479]);
            break;
        case -5122: // is an anti Λ0b
            Fill_particleId(pdgIdMother, IdSummary[480]);
            break;
        case 5112: // is a Σb-
            Fill_particleId(pdgIdMother, IdSummary[481]);
            break;
        case -5112: // is an anti Σb-
            Fill_particleId(pdgIdMother, IdSummary[482]);
            break;
        case 5212: // is a Σb0
            Fill_particleId(pdgIdMother, IdSummary[483]);
            break;
        case -5212: // is an anti Σb0
            Fill_particleId(pdgIdMother, IdSummary[484]);
            break;
        case 5222: // is a Σb+
            Fill_particleId(pdgIdMother, IdSummary[485]);
            break;
        case -5222: // is an anti Σb+
            Fill_particleId(pdgIdMother, IdSummary[486]);
            break;
        case 5114: // is a Σb*-
            Fill_particleId(pdgIdMother, IdSummary[487]);
            break;
        case -5114: // is an anti Σb*-
            Fill_particleId(pdgIdMother, IdSummary[488]);
            break;
        case 5214: // is a Σb*0
            Fill_particleId(pdgIdMother, IdSummary[489]);
            break;
        case -5214: // is an anti Σb*0
            Fill_particleId(pdgIdMother, IdSummary[490]);
            break;
        case 5224: // is a Σb∗+
            Fill_particleId(pdgIdMother, IdSummary[491]);
            break;
        case -5224: // is an anti Σb∗+
            Fill_particleId(pdgIdMother, IdSummary[492]);
            break;
        case 5132: // is a Ξb−
            Fill_particleId(pdgIdMother, IdSummary[493]);
            break;
        case -5132: // is a Ξb+
            Fill_particleId(pdgIdMother, IdSummary[494]);
            break;
        case 5232: // is a Ξb0
            Fill_particleId(pdgIdMother, IdSummary[495]);
            break;
        case -5232: // is an anti Ξb0
            Fill_particleId(pdgIdMother, IdSummary[496]);
            break;
        case 5312: // is a Ξ′b−
            Fill_particleId(pdgIdMother, IdSummary[497]);
            break;
        case -5312: // is a Ξ′b+
            Fill_particleId(pdgIdMother, IdSummary[498]);
            break;
        case 5322: // is a Ξb′0
            Fill_particleId(pdgIdMother, IdSummary[499]);
            break;
        case -5322: // is an anti Ξb′0
            Fill_particleId(pdgIdMother, IdSummary[500]);
            break;
        case 5314: // is a Ξb∗−
            Fill_particleId(pdgIdMother, IdSummary[501]);
            break;
        case -5314: // is a Ξb∗+
            Fill_particleId(pdgIdMother, IdSummary[502]);
            break;
        case 5324: // is a Ξb∗0
            Fill_particleId(pdgIdMother, IdSummary[503]);
            break;
        case -5324: // is an anti Ξb∗0
            Fill_particleId(pdgIdMother, IdSummary[504]);
            break;
        case 5332: // is a Ωb−
            Fill_particleId(pdgIdMother, IdSummary[505]);
            break;
        case -5332: // is a Ωb+
            Fill_particleId(pdgIdMother, IdSummary[506]);
            break;
        case 5334: // is a Ωb*−
            Fill_particleId(pdgIdMother, IdSummary[507]);
            break;
        case -5334: // is a Ωb*+
            Fill_particleId(pdgIdMother, IdSummary[508]);
            break;
        case 5142: // is a Ξbc0
            Fill_particleId(pdgIdMother, IdSummary[509]);
            break;
        case -5142: // is an anti Ξbc0
            Fill_particleId(pdgIdMother, IdSummary[510]);
            break;
        case 5242: // is a Ξbc+
            Fill_particleId(pdgIdMother, IdSummary[511]);
            break;
        case -5242: // is a Ξbc-
            Fill_particleId(pdgIdMother, IdSummary[512]);
            break;
        case 5412: // is a Ξ′0bc
            Fill_particleId(pdgIdMother, IdSummary[513]);
            break;
        case -5412: // is an anti Ξ′0bc
            Fill_particleId(pdgIdMother, IdSummary[514]);
            break;
        case 5422: // is a Ξ′bc+
            Fill_particleId(pdgIdMother, IdSummary[515]);
            break;
        case -5422: // is a Ξ′bc-
            Fill_particleId(pdgIdMother, IdSummary[516]);
            break;
        case 5414: // is a Ξbc∗0
            Fill_particleId(pdgIdMother, IdSummary[517]);
            break;
        case -5414: // is an anti Ξbc∗0
            Fill_particleId(pdgIdMother, IdSummary[518]);
            break;
        case 5424: // is a Ξbc∗+
            Fill_particleId(pdgIdMother, IdSummary[519]);
            break;
        case -5424: // is a Ξbc∗-
            Fill_particleId(pdgIdMother, IdSummary[520]);
            break;
        case 5342: // is a Ωbc0
            Fill_particleId(pdgIdMother, IdSummary[521]);
            break;
        case -5342: // is an anti Ωbc0
            Fill_particleId(pdgIdMother, IdSummary[522]);
            break;
        case 5432: // is a Ωbc′0
            Fill_particleId(pdgIdMother, IdSummary[523]);
            break;
        case -5432: // is an anti Ωbc′0
            Fill_particleId(pdgIdMother, IdSummary[524]);
            break;
        case 5434: // is a Ωbc∗0
            Fill_particleId(pdgIdMother, IdSummary[525]);
            break;
        case -5434: // is an anti Ωbc∗0
            Fill_particleId(pdgIdMother, IdSummary[526]);
            break;
        case 5442: // is a Ωbcc+
            Fill_particleId(pdgIdMother, IdSummary[527]);
            break;
        case -5442: // is a Ωbcc-
            Fill_particleId(pdgIdMother, IdSummary[528]);
            break;
        case 5444: // is a Ωbcc∗+
            Fill_particleId(pdgIdMother, IdSummary[529]);
            break;
        case -5444: // is a Ωbcc∗-
            Fill_particleId(pdgIdMother, IdSummary[530]);
            break;
        case 5512: // is a Ξbb−
            Fill_particleId(pdgIdMother, IdSummary[531]);
            break;
        case -5512: // is a Ξbb+
            Fill_particleId(pdgIdMother, IdSummary[532]);
            break;
        case 5522: // is a Ξbb0
            Fill_particleId(pdgIdMother, IdSummary[533]);
            break;
        case -5522: // is an anti Ξbb0
            Fill_particleId(pdgIdMother, IdSummary[534]);
            break;
        case 5514: // is a Ξbb∗−
            Fill_particleId(pdgIdMother, IdSummary[535]);
            break;
        case -5514: // is a Ξbb∗+
            Fill_particleId(pdgIdMother, IdSummary[536]);
            break;
        case 5524: // is a Ξbb∗0
            Fill_particleId(pdgIdMother, IdSummary[537]);
            break;
        case -5524: // is an anti Ξbb∗0
            Fill_particleId(pdgIdMother, IdSummary[538]);
            break;
        case 5532: // is a Ωbb−
            Fill_particleId(pdgIdMother, IdSummary[539]);
            break;
        case -5532: // is a Ωbb+
            Fill_particleId(pdgIdMother, IdSummary[540]);
            break;
        case 5534: // is a Ωbb∗−
            Fill_particleId(pdgIdMother, IdSummary[541]);
            break;
        case -5534: // is a Ωbb∗+
            Fill_particleId(pdgIdMother, IdSummary[542]);
            break;
        case 5542: // is a Ωbbc0
            Fill_particleId(pdgIdMother, IdSummary[543]);
            break;
        case -5542: // is an anti Ωbbc0
            Fill_particleId(pdgIdMother, IdSummary[544]);
            break;
        case 5544: // is a Ωbbc∗0
            Fill_particleId(pdgIdMother, IdSummary[545]);
            break;
        case -5544: // is an anti Ωbbc∗0
            Fill_particleId(pdgIdMother, IdSummary[546]);
            break;
        case 5554: // is a Ωbbb−
            Fill_particleId(pdgIdMother, IdSummary[547]);
            break;
        case -5554: // is a Ωbbb+
            Fill_particleId(pdgIdMother, IdSummary[548]);
            break;
        case 14122: // is a Lambda_c(2593)+
            Fill_particleId(pdgIdMother, IdSummary[549]);
            break;
        case -14122: // is a Lambda_c(2593)-
            Fill_particleId(pdgIdMother, IdSummary[550]);
            break;
        case 30221: // is a eta(1440)
            Fill_particleId(pdgIdMother, IdSummary[551]);
            break;
        case 4124: // is a Lambda_c(2625)+
            Fill_particleId(pdgIdMother, IdSummary[552]);
            break;
        case -4124: // is a Lambda_c(2625)-
             Fill_particleId(pdgIdMother, IdSummary[553]);
             break;
            /*
             case ??: // is a
             Fill_particleId(pdgIdMother, IdSummary[554]);
             break;
             case ??: // is a
             Fill_particleId(pdgIdMother, IdSummary[555]);
             break;
             case 0: // is a UNKNOWN
             Fill_particleId(pdgIdMother, IdSummary[556]);
             break;
             case 0: // is a UNKNOWN
             Fill_particleId(pdgIdMother, IdSummary[557]);
             break;
             case 0: // is a UNKNOWN
             Fill_particleId(pdgIdMother, IdSummary[558]);
             break;
             */
            
        case 0: // is a UNKNOWN
            Fill_particleId(pdgIdMother, IdSummary[559]);
            break;
            
        default: cout << "Unknown particle! PdgId = " << pdgId << endl;
            break;
    }
}
