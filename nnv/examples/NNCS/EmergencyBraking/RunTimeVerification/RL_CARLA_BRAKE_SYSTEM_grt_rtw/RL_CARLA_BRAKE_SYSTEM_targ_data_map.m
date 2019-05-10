  function targMap = targDataMap(),

  ;%***********************
  ;% Create Parameter Map *
  ;%***********************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 1;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc paramMap
    ;%
    paramMap.nSections           = nTotSects;
    paramMap.sectIdxOffset       = sectIdxOffset;
      paramMap.sections(nTotSects) = dumSection; %prealloc
    paramMap.nTotData            = -1;
    
    ;%
    ;% Auto data (RL_CARLA_BRAKE_SYSTEM_P)
    ;%
      section.nData     = 189;
      section.data(189)  = dumData; %prealloc
      
	  ;% RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_A
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_B
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 3;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_C
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 5;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.DiscreteStateSpace_InitialCondi
	  section.data(4).logicalSrcIdx = 4;
	  section.data(4).dtTransOffset = 8;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Gain_Gain
	  section.data(5).logicalSrcIdx = 5;
	  section.data(5).dtTransOffset = 11;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Constant_Value
	  section.data(6).logicalSrcIdx = 6;
	  section.data(6).dtTransOffset = 12;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Gain1_Gain
	  section.data(7).logicalSrcIdx = 7;
	  section.data(7).dtTransOffset = 13;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Gain2_Gain
	  section.data(8).logicalSrcIdx = 8;
	  section.data(8).dtTransOffset = 14;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW111_Value
	  section.data(9).logicalSrcIdx = 9;
	  section.data(9).dtTransOffset = 15;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1110_Value
	  section.data(10).logicalSrcIdx = 10;
	  section.data(10).dtTransOffset = 18;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1111_Value
	  section.data(11).logicalSrcIdx = 11;
	  section.data(11).dtTransOffset = 21;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1112_Value
	  section.data(12).logicalSrcIdx = 12;
	  section.data(12).dtTransOffset = 24;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1113_Value
	  section.data(13).logicalSrcIdx = 13;
	  section.data(13).dtTransOffset = 27;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1114_Value
	  section.data(14).logicalSrcIdx = 14;
	  section.data(14).dtTransOffset = 30;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1115_Value
	  section.data(15).logicalSrcIdx = 15;
	  section.data(15).dtTransOffset = 33;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1116_Value
	  section.data(16).logicalSrcIdx = 16;
	  section.data(16).dtTransOffset = 36;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1117_Value
	  section.data(17).logicalSrcIdx = 17;
	  section.data(17).dtTransOffset = 39;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1118_Value
	  section.data(18).logicalSrcIdx = 18;
	  section.data(18).dtTransOffset = 42;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1119_Value
	  section.data(19).logicalSrcIdx = 19;
	  section.data(19).dtTransOffset = 45;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW112_Value
	  section.data(20).logicalSrcIdx = 20;
	  section.data(20).dtTransOffset = 48;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1120_Value
	  section.data(21).logicalSrcIdx = 21;
	  section.data(21).dtTransOffset = 51;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1121_Value
	  section.data(22).logicalSrcIdx = 22;
	  section.data(22).dtTransOffset = 54;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1122_Value
	  section.data(23).logicalSrcIdx = 23;
	  section.data(23).dtTransOffset = 57;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1123_Value
	  section.data(24).logicalSrcIdx = 24;
	  section.data(24).dtTransOffset = 60;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1124_Value
	  section.data(25).logicalSrcIdx = 25;
	  section.data(25).dtTransOffset = 63;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1125_Value
	  section.data(26).logicalSrcIdx = 26;
	  section.data(26).dtTransOffset = 66;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1126_Value
	  section.data(27).logicalSrcIdx = 27;
	  section.data(27).dtTransOffset = 69;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1127_Value
	  section.data(28).logicalSrcIdx = 28;
	  section.data(28).dtTransOffset = 72;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1128_Value
	  section.data(29).logicalSrcIdx = 29;
	  section.data(29).dtTransOffset = 75;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1129_Value
	  section.data(30).logicalSrcIdx = 30;
	  section.data(30).dtTransOffset = 78;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW113_Value
	  section.data(31).logicalSrcIdx = 31;
	  section.data(31).dtTransOffset = 81;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1130_Value
	  section.data(32).logicalSrcIdx = 32;
	  section.data(32).dtTransOffset = 84;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1131_Value
	  section.data(33).logicalSrcIdx = 33;
	  section.data(33).dtTransOffset = 87;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1132_Value
	  section.data(34).logicalSrcIdx = 34;
	  section.data(34).dtTransOffset = 90;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1133_Value
	  section.data(35).logicalSrcIdx = 35;
	  section.data(35).dtTransOffset = 93;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1134_Value
	  section.data(36).logicalSrcIdx = 36;
	  section.data(36).dtTransOffset = 96;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1135_Value
	  section.data(37).logicalSrcIdx = 37;
	  section.data(37).dtTransOffset = 99;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1136_Value
	  section.data(38).logicalSrcIdx = 38;
	  section.data(38).dtTransOffset = 102;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1137_Value
	  section.data(39).logicalSrcIdx = 39;
	  section.data(39).dtTransOffset = 105;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1138_Value
	  section.data(40).logicalSrcIdx = 40;
	  section.data(40).dtTransOffset = 108;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1139_Value
	  section.data(41).logicalSrcIdx = 41;
	  section.data(41).dtTransOffset = 111;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW114_Value
	  section.data(42).logicalSrcIdx = 42;
	  section.data(42).dtTransOffset = 114;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1140_Value
	  section.data(43).logicalSrcIdx = 43;
	  section.data(43).dtTransOffset = 117;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1141_Value
	  section.data(44).logicalSrcIdx = 44;
	  section.data(44).dtTransOffset = 120;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1142_Value
	  section.data(45).logicalSrcIdx = 45;
	  section.data(45).dtTransOffset = 123;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1143_Value
	  section.data(46).logicalSrcIdx = 46;
	  section.data(46).dtTransOffset = 126;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1144_Value
	  section.data(47).logicalSrcIdx = 47;
	  section.data(47).dtTransOffset = 129;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1145_Value
	  section.data(48).logicalSrcIdx = 48;
	  section.data(48).dtTransOffset = 132;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1146_Value
	  section.data(49).logicalSrcIdx = 49;
	  section.data(49).dtTransOffset = 135;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1147_Value
	  section.data(50).logicalSrcIdx = 50;
	  section.data(50).dtTransOffset = 138;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1148_Value
	  section.data(51).logicalSrcIdx = 51;
	  section.data(51).dtTransOffset = 141;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1149_Value
	  section.data(52).logicalSrcIdx = 52;
	  section.data(52).dtTransOffset = 144;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW115_Value
	  section.data(53).logicalSrcIdx = 53;
	  section.data(53).dtTransOffset = 147;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1150_Value
	  section.data(54).logicalSrcIdx = 54;
	  section.data(54).dtTransOffset = 150;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW116_Value
	  section.data(55).logicalSrcIdx = 55;
	  section.data(55).dtTransOffset = 153;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW117_Value
	  section.data(56).logicalSrcIdx = 56;
	  section.data(56).dtTransOffset = 156;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW118_Value
	  section.data(57).logicalSrcIdx = 57;
	  section.data(57).dtTransOffset = 159;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW119_Value
	  section.data(58).logicalSrcIdx = 58;
	  section.data(58).dtTransOffset = 162;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.b1_Value
	  section.data(59).logicalSrcIdx = 59;
	  section.data(59).dtTransOffset = 165;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat
	  section.data(60).logicalSrcIdx = 60;
	  section.data(60).dtTransOffset = 215;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat
	  section.data(61).logicalSrcIdx = 61;
	  section.data(61).dtTransOffset = 216;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW211_Value
	  section.data(62).logicalSrcIdx = 62;
	  section.data(62).dtTransOffset = 217;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2110_Value
	  section.data(63).logicalSrcIdx = 63;
	  section.data(63).dtTransOffset = 267;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2111_Value
	  section.data(64).logicalSrcIdx = 64;
	  section.data(64).dtTransOffset = 317;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2112_Value
	  section.data(65).logicalSrcIdx = 65;
	  section.data(65).dtTransOffset = 367;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2113_Value
	  section.data(66).logicalSrcIdx = 66;
	  section.data(66).dtTransOffset = 417;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2114_Value
	  section.data(67).logicalSrcIdx = 67;
	  section.data(67).dtTransOffset = 467;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2115_Value
	  section.data(68).logicalSrcIdx = 68;
	  section.data(68).dtTransOffset = 517;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2116_Value
	  section.data(69).logicalSrcIdx = 69;
	  section.data(69).dtTransOffset = 567;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2117_Value
	  section.data(70).logicalSrcIdx = 70;
	  section.data(70).dtTransOffset = 617;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2118_Value
	  section.data(71).logicalSrcIdx = 71;
	  section.data(71).dtTransOffset = 667;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2119_Value
	  section.data(72).logicalSrcIdx = 72;
	  section.data(72).dtTransOffset = 717;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW212_Value
	  section.data(73).logicalSrcIdx = 73;
	  section.data(73).dtTransOffset = 767;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2120_Value
	  section.data(74).logicalSrcIdx = 74;
	  section.data(74).dtTransOffset = 817;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2121_Value
	  section.data(75).logicalSrcIdx = 75;
	  section.data(75).dtTransOffset = 867;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2122_Value
	  section.data(76).logicalSrcIdx = 76;
	  section.data(76).dtTransOffset = 917;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2123_Value
	  section.data(77).logicalSrcIdx = 77;
	  section.data(77).dtTransOffset = 967;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2124_Value
	  section.data(78).logicalSrcIdx = 78;
	  section.data(78).dtTransOffset = 1017;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2125_Value
	  section.data(79).logicalSrcIdx = 79;
	  section.data(79).dtTransOffset = 1067;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2126_Value
	  section.data(80).logicalSrcIdx = 80;
	  section.data(80).dtTransOffset = 1117;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2127_Value
	  section.data(81).logicalSrcIdx = 81;
	  section.data(81).dtTransOffset = 1167;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2128_Value
	  section.data(82).logicalSrcIdx = 82;
	  section.data(82).dtTransOffset = 1217;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2129_Value
	  section.data(83).logicalSrcIdx = 83;
	  section.data(83).dtTransOffset = 1267;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW213_Value
	  section.data(84).logicalSrcIdx = 84;
	  section.data(84).dtTransOffset = 1317;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2130_Value
	  section.data(85).logicalSrcIdx = 85;
	  section.data(85).dtTransOffset = 1367;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW214_Value
	  section.data(86).logicalSrcIdx = 86;
	  section.data(86).dtTransOffset = 1417;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW215_Value
	  section.data(87).logicalSrcIdx = 87;
	  section.data(87).dtTransOffset = 1467;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW216_Value
	  section.data(88).logicalSrcIdx = 88;
	  section.data(88).dtTransOffset = 1517;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW217_Value
	  section.data(89).logicalSrcIdx = 89;
	  section.data(89).dtTransOffset = 1567;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW218_Value
	  section.data(90).logicalSrcIdx = 90;
	  section.data(90).dtTransOffset = 1617;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW219_Value
	  section.data(91).logicalSrcIdx = 91;
	  section.data(91).dtTransOffset = 1667;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.b2_Value
	  section.data(92).logicalSrcIdx = 92;
	  section.data(92).dtTransOffset = 1717;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_e
	  section.data(93).logicalSrcIdx = 93;
	  section.data(93).dtTransOffset = 1747;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_p
	  section.data(94).logicalSrcIdx = 94;
	  section.data(94).dtTransOffset = 1748;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW321_Value
	  section.data(95).logicalSrcIdx = 95;
	  section.data(95).dtTransOffset = 1749;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.b3_Value
	  section.data(96).logicalSrcIdx = 96;
	  section.data(96).dtTransOffset = 1779;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_j
	  section.data(97).logicalSrcIdx = 97;
	  section.data(97).dtTransOffset = 1780;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_n
	  section.data(98).logicalSrcIdx = 98;
	  section.data(98).dtTransOffset = 1781;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Gain1_Gain_m
	  section.data(99).logicalSrcIdx = 99;
	  section.data(99).dtTransOffset = 1782;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW321_Value_m
	  section.data(100).logicalSrcIdx = 100;
	  section.data(100).dtTransOffset = 1783;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW211_Value_a
	  section.data(101).logicalSrcIdx = 101;
	  section.data(101).dtTransOffset = 1813;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW111_Value_e
	  section.data(102).logicalSrcIdx = 102;
	  section.data(102).dtTransOffset = 1863;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW112_Value_i
	  section.data(103).logicalSrcIdx = 103;
	  section.data(103).dtTransOffset = 1865;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW113_Value_p
	  section.data(104).logicalSrcIdx = 104;
	  section.data(104).dtTransOffset = 1867;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW114_Value_p
	  section.data(105).logicalSrcIdx = 105;
	  section.data(105).dtTransOffset = 1869;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW115_Value_m
	  section.data(106).logicalSrcIdx = 106;
	  section.data(106).dtTransOffset = 1871;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW116_Value_k
	  section.data(107).logicalSrcIdx = 107;
	  section.data(107).dtTransOffset = 1873;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW117_Value_p
	  section.data(108).logicalSrcIdx = 108;
	  section.data(108).dtTransOffset = 1875;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW118_Value_f
	  section.data(109).logicalSrcIdx = 109;
	  section.data(109).dtTransOffset = 1877;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW119_Value_n
	  section.data(110).logicalSrcIdx = 110;
	  section.data(110).dtTransOffset = 1879;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1110_Value_d
	  section.data(111).logicalSrcIdx = 111;
	  section.data(111).dtTransOffset = 1881;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1111_Value_g
	  section.data(112).logicalSrcIdx = 112;
	  section.data(112).dtTransOffset = 1883;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1112_Value_b
	  section.data(113).logicalSrcIdx = 113;
	  section.data(113).dtTransOffset = 1885;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1113_Value_l
	  section.data(114).logicalSrcIdx = 114;
	  section.data(114).dtTransOffset = 1887;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1114_Value_b
	  section.data(115).logicalSrcIdx = 115;
	  section.data(115).dtTransOffset = 1889;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1115_Value_c
	  section.data(116).logicalSrcIdx = 116;
	  section.data(116).dtTransOffset = 1891;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1116_Value_p
	  section.data(117).logicalSrcIdx = 117;
	  section.data(117).dtTransOffset = 1893;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1117_Value_i
	  section.data(118).logicalSrcIdx = 118;
	  section.data(118).dtTransOffset = 1895;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1118_Value_c
	  section.data(119).logicalSrcIdx = 119;
	  section.data(119).dtTransOffset = 1897;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1119_Value_b
	  section.data(120).logicalSrcIdx = 120;
	  section.data(120).dtTransOffset = 1899;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1120_Value_m
	  section.data(121).logicalSrcIdx = 121;
	  section.data(121).dtTransOffset = 1901;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1121_Value_i
	  section.data(122).logicalSrcIdx = 122;
	  section.data(122).dtTransOffset = 1903;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1122_Value_k
	  section.data(123).logicalSrcIdx = 123;
	  section.data(123).dtTransOffset = 1905;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1123_Value_b
	  section.data(124).logicalSrcIdx = 124;
	  section.data(124).dtTransOffset = 1907;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1124_Value_i
	  section.data(125).logicalSrcIdx = 125;
	  section.data(125).dtTransOffset = 1909;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1125_Value_l
	  section.data(126).logicalSrcIdx = 126;
	  section.data(126).dtTransOffset = 1911;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1126_Value_h
	  section.data(127).logicalSrcIdx = 127;
	  section.data(127).dtTransOffset = 1913;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1127_Value_c
	  section.data(128).logicalSrcIdx = 128;
	  section.data(128).dtTransOffset = 1915;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1128_Value_a
	  section.data(129).logicalSrcIdx = 129;
	  section.data(129).dtTransOffset = 1917;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1129_Value_p
	  section.data(130).logicalSrcIdx = 130;
	  section.data(130).dtTransOffset = 1919;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1130_Value_h
	  section.data(131).logicalSrcIdx = 131;
	  section.data(131).dtTransOffset = 1921;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1131_Value_l
	  section.data(132).logicalSrcIdx = 132;
	  section.data(132).dtTransOffset = 1923;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1132_Value_i
	  section.data(133).logicalSrcIdx = 133;
	  section.data(133).dtTransOffset = 1925;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1133_Value_g
	  section.data(134).logicalSrcIdx = 134;
	  section.data(134).dtTransOffset = 1927;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1134_Value_d
	  section.data(135).logicalSrcIdx = 135;
	  section.data(135).dtTransOffset = 1929;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1135_Value_n
	  section.data(136).logicalSrcIdx = 136;
	  section.data(136).dtTransOffset = 1931;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1136_Value_k
	  section.data(137).logicalSrcIdx = 137;
	  section.data(137).dtTransOffset = 1933;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1137_Value_f
	  section.data(138).logicalSrcIdx = 138;
	  section.data(138).dtTransOffset = 1935;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1138_Value_o
	  section.data(139).logicalSrcIdx = 139;
	  section.data(139).dtTransOffset = 1937;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1139_Value_o
	  section.data(140).logicalSrcIdx = 140;
	  section.data(140).dtTransOffset = 1939;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1140_Value_l
	  section.data(141).logicalSrcIdx = 141;
	  section.data(141).dtTransOffset = 1941;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1141_Value_c
	  section.data(142).logicalSrcIdx = 142;
	  section.data(142).dtTransOffset = 1943;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1142_Value_l
	  section.data(143).logicalSrcIdx = 143;
	  section.data(143).dtTransOffset = 1945;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1143_Value_i
	  section.data(144).logicalSrcIdx = 144;
	  section.data(144).dtTransOffset = 1947;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1144_Value_d
	  section.data(145).logicalSrcIdx = 145;
	  section.data(145).dtTransOffset = 1949;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1145_Value_c
	  section.data(146).logicalSrcIdx = 146;
	  section.data(146).dtTransOffset = 1951;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1146_Value_p
	  section.data(147).logicalSrcIdx = 147;
	  section.data(147).dtTransOffset = 1953;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1147_Value_d
	  section.data(148).logicalSrcIdx = 148;
	  section.data(148).dtTransOffset = 1955;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1148_Value_l
	  section.data(149).logicalSrcIdx = 149;
	  section.data(149).dtTransOffset = 1957;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1149_Value_d
	  section.data(150).logicalSrcIdx = 150;
	  section.data(150).dtTransOffset = 1959;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW1150_Value_f
	  section.data(151).logicalSrcIdx = 151;
	  section.data(151).dtTransOffset = 1961;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.b1_Value_i
	  section.data(152).logicalSrcIdx = 152;
	  section.data(152).dtTransOffset = 1963;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_f
	  section.data(153).logicalSrcIdx = 153;
	  section.data(153).dtTransOffset = 2013;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_b
	  section.data(154).logicalSrcIdx = 154;
	  section.data(154).dtTransOffset = 2014;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW212_Value_j
	  section.data(155).logicalSrcIdx = 155;
	  section.data(155).dtTransOffset = 2015;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW213_Value_b
	  section.data(156).logicalSrcIdx = 156;
	  section.data(156).dtTransOffset = 2065;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW214_Value_a
	  section.data(157).logicalSrcIdx = 157;
	  section.data(157).dtTransOffset = 2115;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW215_Value_d
	  section.data(158).logicalSrcIdx = 158;
	  section.data(158).dtTransOffset = 2165;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW216_Value_p
	  section.data(159).logicalSrcIdx = 159;
	  section.data(159).dtTransOffset = 2215;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW217_Value_d
	  section.data(160).logicalSrcIdx = 160;
	  section.data(160).dtTransOffset = 2265;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW218_Value_j
	  section.data(161).logicalSrcIdx = 161;
	  section.data(161).dtTransOffset = 2315;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW219_Value_f
	  section.data(162).logicalSrcIdx = 162;
	  section.data(162).dtTransOffset = 2365;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2110_Value_m
	  section.data(163).logicalSrcIdx = 163;
	  section.data(163).dtTransOffset = 2415;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2111_Value_a
	  section.data(164).logicalSrcIdx = 164;
	  section.data(164).dtTransOffset = 2465;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2112_Value_a
	  section.data(165).logicalSrcIdx = 165;
	  section.data(165).dtTransOffset = 2515;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2113_Value_c
	  section.data(166).logicalSrcIdx = 166;
	  section.data(166).dtTransOffset = 2565;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2114_Value_e
	  section.data(167).logicalSrcIdx = 167;
	  section.data(167).dtTransOffset = 2615;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2115_Value_j
	  section.data(168).logicalSrcIdx = 168;
	  section.data(168).dtTransOffset = 2665;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2116_Value_j
	  section.data(169).logicalSrcIdx = 169;
	  section.data(169).dtTransOffset = 2715;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2117_Value_e
	  section.data(170).logicalSrcIdx = 170;
	  section.data(170).dtTransOffset = 2765;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2118_Value_n
	  section.data(171).logicalSrcIdx = 171;
	  section.data(171).dtTransOffset = 2815;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2119_Value_h
	  section.data(172).logicalSrcIdx = 172;
	  section.data(172).dtTransOffset = 2865;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2120_Value_j
	  section.data(173).logicalSrcIdx = 173;
	  section.data(173).dtTransOffset = 2915;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2121_Value_f
	  section.data(174).logicalSrcIdx = 174;
	  section.data(174).dtTransOffset = 2965;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2122_Value_k
	  section.data(175).logicalSrcIdx = 175;
	  section.data(175).dtTransOffset = 3015;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2123_Value_j
	  section.data(176).logicalSrcIdx = 176;
	  section.data(176).dtTransOffset = 3065;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2124_Value_a
	  section.data(177).logicalSrcIdx = 177;
	  section.data(177).dtTransOffset = 3115;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2125_Value_i
	  section.data(178).logicalSrcIdx = 178;
	  section.data(178).dtTransOffset = 3165;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2126_Value_g
	  section.data(179).logicalSrcIdx = 179;
	  section.data(179).dtTransOffset = 3215;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2127_Value_p
	  section.data(180).logicalSrcIdx = 180;
	  section.data(180).dtTransOffset = 3265;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2128_Value_j
	  section.data(181).logicalSrcIdx = 181;
	  section.data(181).dtTransOffset = 3315;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2129_Value_h
	  section.data(182).logicalSrcIdx = 182;
	  section.data(182).dtTransOffset = 3365;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.IW2130_Value_f
	  section.data(183).logicalSrcIdx = 183;
	  section.data(183).dtTransOffset = 3415;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.b2_Value_h
	  section.data(184).logicalSrcIdx = 184;
	  section.data(184).dtTransOffset = 3465;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Saturation_UpperSat_c
	  section.data(185).logicalSrcIdx = 185;
	  section.data(185).dtTransOffset = 3495;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Saturation_LowerSat_bl
	  section.data(186).logicalSrcIdx = 186;
	  section.data(186).dtTransOffset = 3496;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.b3_Value_b
	  section.data(187).logicalSrcIdx = 187;
	  section.data(187).dtTransOffset = 3497;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Gain_Gain_h
	  section.data(188).logicalSrcIdx = 188;
	  section.data(188).dtTransOffset = 3498;
	
	  ;% RL_CARLA_BRAKE_SYSTEM_P.Gain2_Gain_o
	  section.data(189).logicalSrcIdx = 189;
	  section.data(189).dtTransOffset = 3499;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(1) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (parameter)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    paramMap.nTotData = nTotData;
    


  ;%**************************
  ;% Create Block Output Map *
  ;%**************************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 1;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc sigMap
    ;%
    sigMap.nSections           = nTotSects;
    sigMap.sectIdxOffset       = sectIdxOffset;
      sigMap.sections(nTotSects) = dumSection; %prealloc
    sigMap.nTotData            = -1;
    
    ;%
    ;% Auto data (RL_CARLA_BRAKE_SYSTEM_B)
    ;%
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% RL_CARLA_BRAKE_SYSTEM_B.DiscreteStateSpace
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(1) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (signal)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    sigMap.nTotData = nTotData;
    


  ;%*******************
  ;% Create DWork Map *
  ;%*******************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 2;
    sectIdxOffset = 1;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc dworkMap
    ;%
    dworkMap.nSections           = nTotSects;
    dworkMap.sectIdxOffset       = sectIdxOffset;
      dworkMap.sections(nTotSects) = dumSection; %prealloc
    dworkMap.nTotData            = -1;
    
    ;%
    ;% Auto data (RL_CARLA_BRAKE_SYSTEM_DW)
    ;%
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% RL_CARLA_BRAKE_SYSTEM_DW.DiscreteStateSpace_DSTATE
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(1) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% RL_CARLA_BRAKE_SYSTEM_DW.Scope_PWORK.LoggedData
	  section.data(1).logicalSrcIdx = 1;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(2) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (dwork)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    dworkMap.nTotData = nTotData;
    


  ;%
  ;% Add individual maps to base struct.
  ;%

  targMap.paramMap  = paramMap;    
  targMap.signalMap = sigMap;
  targMap.dworkMap  = dworkMap;
  
  ;%
  ;% Add checksums to base struct.
  ;%


  targMap.checksum0 = 112044904;
  targMap.checksum1 = 1786873236;
  targMap.checksum2 = 1306172202;
  targMap.checksum3 = 2742050055;

