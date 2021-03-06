c
c  These contain bitmaps for the ascii character set (some
c  characters, such as "g", have a bug). These bitmaps where pinched
c  from Tex pxl files for amtt5 (?).
c
      integer maxvoff,maxv
      parameter(maxv=23,maxvoff=5)
      integer voff(33:127),hoff(33:127),pattn(16,33:127),i
      data hoff( 33),voff( 33),(pattn(i, 33),i=1,16)/ -5,  0,
     *      7,    7,    7,    0,    0,    6,    6,    6,
     *      6,    6,    7,    7,    7,    7,    7,    7/
      data hoff( 34),voff( 34),(pattn(i, 34),i=1,16)/ -2, -8,
     *    231,  231,  231,  231,  231,  231,  231,  231,
     *      0,    0,    0,    0,    0,    0,    0,    0/
      data hoff( 35),voff( 35),(pattn(i, 35),i=1,16)/ -1,  0,
     *    102,  110,  108,  108,  236, 1023, 1023,  204,
     *    204,  204, 1023, 1023,  216,  216,  472,  408/
      data hoff( 36),voff( 36),(pattn(i, 36),i=1,16)/ -1, -2,
     *    951,  823,  823,  823,  944,  496,  252,   62,
     *    951,  947,  947,  951,  510,  252,   48,   48/
      data hoff( 37),voff( 37),(pattn(i, 37),i=1,16)/ -1, -2,
     *    972,  988,  984,  408,   56,   48,   48,  112,
     *    102,  111,  239,  207,  207,  463,  399,  390/
      data hoff( 38),voff( 38),(pattn(i, 38),i=1,16)/ -1,  0,
     *    478, 1023,  883,  883,  243,  251,  479,  414,
     *    924, 2012, 2044,   60,   60,   60,   60,   24/
      data hoff( 39),voff( 39),(pattn(i, 39),i=1,16)/ -5,-10,
     *      3,    7,    6,    7,    7,    7,    0,    0,
     *      0,    0,    0,    0,    0,    0,    0,    0/
      data hoff( 40),voff( 40),(pattn(i, 40),i=1,16)/ -4, -2,
     *      6,    6,    7,    3,    3,    3,    3,    3,
     *      3,    7,    6,    6,   14,   28,   56,   48/
      data hoff( 41),voff( 41),(pattn(i, 41),i=1,16)/ -2, -2,
     *     24,   24,   56,   48,   48,   48,   48,   48,
     *     48,   56,   24,   24,   28,   14,    7,    3/
      data hoff( 42),voff( 42),(pattn(i, 42),i=1,16)/ -1, -3,
     *     48,   48,  819,  951,  252,   48,  252,  951,
     *    819,   48,   48,    0,    0,    0,    0,    0/
      data hoff( 43),voff( 43),(pattn(i, 43),i=1,16)/ -1, -3,
     *     48,   48,   48,   48, 1023, 1023,   48,   48,
     *     48,   48,    0,    0,    0,    0,    0,    0/
      data hoff( 44),voff( 44),(pattn(i, 44),i=1,16)/ -5,  3,
     *      3,    7,    6,    7,    7,    7,    0,    0,
     *      0,    0,    0,    0,    0,    0,    0,    0/
      data hoff( 45),voff( 45),(pattn(i, 45),i=1,16)/ -2, -7,
     *    255,  255,    0,    0,    0,    0,    0,    0,
     *      0,    0,    0,    0,    0,    0,    0,    0/
      data hoff( 46),voff( 46),(pattn(i, 46),i=1,16)/ -5,  0,
     *      7,    7,    7,    0,    0,    0,    0,    0,
     *      0,    0,    0,    0,    0,    0,    0,    0/
      data hoff( 47),voff( 47),(pattn(i, 47),i=1,16)/ -1, -2,
     *     14,   12,   28,   24,   56,   48,   48,  112,
     *     96,  224,  192,  448,  384,  896,  768,  768/
      data hoff( 48),voff( 48),(pattn(i, 48),i=1,16)/ -1,  0,
     *    120,  252,  462,  390,  903,  771,  771,  771,
     *    771,  771,  771,  903,  390,  462,  252,  120/
      data hoff( 49),voff( 49),(pattn(i, 49),i=1,16)/ -2,  0,
     *    255,  255,   24,   24,   24,   24,   24,   24,
     *     24,   24,   24,   31,   31,   28,   24,   24/
      data hoff( 50),voff( 50),(pattn(i, 50),i=1,16)/ -1,  0,
     *   1023, 1023,  782,  796,   56,  112,  224,  448,
     *    384,  896,  775,  775,  775,  967,  510,  124/
      data hoff( 51),voff( 51),(pattn(i, 51),i=1,16)/ -1,  0,
     *    124,  510,  455,  903,  775,  775,  768,  384,
     *    120,  248,  192,  398,  398,  462,  254,  124/
      data hoff( 52),voff( 52),(pattn(i, 52),i=1,16)/ -1,  0,
     *   1008, 1008,  192,  192, 1023, 1023,  194,  198,
     *    196,  204,  200,  200,  208,  208,  240,  224/
      data hoff( 53),voff( 53),(pattn(i, 53),i=1,16)/ -1,  0,
     *    124,  510,  455,  903,  775,  775,  768,  902,
     *    510,  254,    6,    6,    6,    6,  510,  510/
      data hoff( 54),voff( 54),(pattn(i, 54),i=1,16)/ -1,  0,
     *    248,  508,  396,  910,  774,  775,  775,  911,
     *    399,  511,  251,    7,  462,  476,  508,  240/
      data hoff( 55),voff( 55),(pattn(i, 55),i=1,16)/ -1,  0,
     *     24,   24,   24,   24,   24,   56,   48,   48,
     *    112,   96,   96,  224,  451,  387, 1023, 1023/
      data hoff( 56),voff( 56),(pattn(i, 56),i=1,16)/ -1,  0,
     *    252,  510,  903,  771,  771,  771,  771,  390,
     *    120,  252,  390,  771,  771,  903,  510,  252/
      data hoff( 57),voff( 57),(pattn(i, 57),i=1,16)/ -1,  0,
     *     60,  254,  494,  398,  896,  892, 1022,  966,
     *    903,  771,  899,  899,  455,  486,  254,  124/
      data hoff( 58),voff( 58),(pattn(i, 58),i=1,16)/ -5,  0,
     *      7,    7,    7,    0,    0,    0,    0,    0,
     *      7,    7,    7,    0,    0,    0,    0,    0/
      data hoff( 59),voff( 59),(pattn(i, 59),i=1,16)/ -5,  3,
     *      3,    7,    6,    7,    7,    7,    0,    0,
     *      0,    0,    0,    7,    7,    7,    0,    0/
      data hoff( 60),voff( 60),(pattn(i, 60),i=1,16)/ -1, -2,
     *    768,  960,  496,  120,   30,   15,   15,   62,
     *    120,  480,  960,  768,    0,    0,    0,    0/
      data hoff( 61),voff( 61),(pattn(i, 61),i=1,16)/ -1, -5,
     *   1023, 1023,    0,    0, 1023, 1023,    0,    0,
     *      0,    0,    0,    0,    0,    0,    0,    0/
      data hoff( 62),voff( 62),(pattn(i, 62),i=1,16)/ -1, -2,
     *      3,   15,   62,  120,  480,  960,  960,  496,
     *    120,   30,   15,    3,    0,    0,    0,    0/
      data hoff( 63),voff( 63),(pattn(i, 63),i=1,16)/ -1,  0,
     *    112,  112,  112,    0,    0,   48,   48,   48,
     *    112,   96,  480,  967,  775,  903,  510,  252/
      data hoff( 64),voff( 64),(pattn(i, 64),i=1,16)/ -1,  0,
     *    496, 1020,  796,  238,  503,  955,  795,  795,
     *    795,  795,  955, 1015, 1006,  412,  508,  240/
      data hoff( 65),voff( 65),(pattn(i, 65),i=1,16)/  0, -1,
     *   3855, 1542,  780, 1020, 1020,  780,  392,  408,
     *    408,  408,  208,  208,  240,  240,  240,   96/
      data hoff( 66),voff( 66),(pattn(i, 66),i=1,16)/ -1,  0,
     *    255,  511,  902,  774,  774,  774,  774,  902,
     *    510,  510,  902,  774,  774,  902,  511,  255/
      data hoff( 67),voff( 67),(pattn(i, 67),i=1,16)/ -1,  0,
     *    248,  508,  974,  902,  775,  771,    3,    3,
     *      3,    3,  771,  775,  902,  974, 1020,  888/
      data hoff( 68),voff( 68),(pattn(i, 68),i=1,16)/ -1,  0,
     *    127,  255,  454,  390,  902,  774,  774,  774,
     *    774,  774,  774,  902,  390,  454,  255,  127/
      data hoff( 69),voff( 69),(pattn(i, 69),i=1,16)/ -1,  0,
     *   1023, 1023,  774,  774,  774,  102,  102,  126,
     *    126,  102,  102,    6,  774,  774, 1023, 1023/
      data hoff( 70),voff( 70),(pattn(i, 70),i=1,16)/  0,  0,
     *     63,   63,   12,   12,   12,  204,  204,  252,
     *    252,  204,  204,   12, 1548, 1548, 2047, 2047/
      data hoff( 71),voff( 71),(pattn(i, 71),i=1,16)/ -1,  0,
     *    504,  508,  462,  390,  391,  995,  995,    3,
     *      3,    3,  387,  391,  454,  462,  508,  504/
      data hoff( 72),voff( 72),(pattn(i, 72),i=1,16)/ -1,  0,
     *    975,  975,  390,  390,  390,  390,  390,  510,
     *    510,  390,  390,  390,  390,  390,  975,  975/
      data hoff( 73),voff( 73),(pattn(i, 73),i=1,16)/ -2,  0,
     *    255,  255,   24,   24,   24,   24,   24,   24,
     *     24,   24,   24,   24,   24,   24,  255,  255/
      data hoff( 74),voff( 74),(pattn(i, 74),i=1,16)/ -1,  0,
     *    124,  255,  455,  391,  391,  384,  384,  384,
     *    384,  384,  384,  384,  384,  384,  992,  992/
      data hoff( 75),voff( 75),(pattn(i, 75),i=1,16)/ -1,  0,
     *   1999, 1999,  390,  390,  198,  198,  110,  110,
     *     62,   62,   54,  102,  102,  198, 1007, 1007/
      data hoff( 76),voff( 76),(pattn(i, 76),i=1,16)/  0,  0,
     *   2047, 2047, 1548, 1548, 1548,   12,   12,   12,
     *     12,   12,   12,   12,   12,   12,   63,   63/
      data hoff( 77),voff( 77),(pattn(i, 77),i=1,16)/  0,  0,
     *   3855, 3855, 1542, 1766, 1766, 2022, 2038, 1974,
     *   1974, 1974, 1942, 1950, 1950, 1950, 3855, 3855/
      data hoff( 78),voff( 78),(pattn(i, 78),i=1,16)/ -1,  0,
     *    463,  463,  486,  486,  486,  422,  422,  422,
     *    406,  406,  406,  414,  414,  414,  975,  975/
      data hoff( 79),voff( 79),(pattn(i, 79),i=1,16)/ -1,  0,
     *    252,  510,  903,  771,  771,  771,  771,  771,
     *    771,  771,  771,  771,  771,  903,  510,  252/
      data hoff( 80),voff( 80),(pattn(i, 80),i=1,16)/ -1,  0,
     *     15,   15,    6,    6,    6,    6,    6,  254,
     *    510,  902,  774,  774,  774,  902,  511,  255/
      data hoff( 81),voff( 81),(pattn(i, 81),i=1,16)/ -1,  0,
     *    508,  510,  999,  883,  819,  819,  771,  771,
     *    771,  771,  771,  771,  771,  903,  510,  252/
      data hoff( 82),voff( 82),(pattn(i, 82),i=1,16)/ -1,  0,
     *    911, 1999, 1734, 1734,  198,  198,  230,  126,
     *    254,  454,  390,  390,  390,  454,  255,  127/
      data hoff( 83),voff( 83),(pattn(i, 83),i=1,16)/ -1,  0,
     *    255,  511,  903,  771,  771,  771,  896,  496,
     *    252,   30,    7,  771,  771,  903, 1022, 1020/
      data hoff( 84),voff( 84),(pattn(i, 84),i=1,16)/ -1,  0,
     *    252,  252,   48,   48,   48,   48,   48,   48,
     *     48,   48,   48,   48,  819,  819, 1023, 1023/
      data hoff( 85),voff( 85),(pattn(i, 85),i=1,16)/  0,  0,
     *    240, 1020,  924, 1806, 1542, 1542, 1542, 1542,
     *   1542, 1542, 1542, 1542, 1542, 1542, 3855, 3855/
      data hoff( 86),voff( 86),(pattn(i, 86),i=1,16)/  0,  0,
     *    240,  240,  240,  176,  176,  408,  408,  408,
     *    280,  780,  780,  780,  780, 1542, 3855, 3855/
      data hoff( 87),voff( 87),(pattn(i, 87),i=1,16)/ -1,  0,
     *    510,  510,  510,  510,  478,  338,  338,  338,
     *    338,  338,  851,  851,  883,  883,  803,  771/
      data hoff( 88),voff( 88),(pattn(i, 88),i=1,16)/  0,  0,
     *   3855, 3855,  780,  408,  408,  240,  240,   96,
     *     96,  240,  240,  408,  408,  780, 3855, 3855/
      data hoff( 89),voff( 89),(pattn(i, 89),i=1,16)/ -1,  0,
     *    120,  120,   48,   48,   48,   48,   48,   48,
     *    120,   88,   88,  204,  204,  204,  975,  975/
      data hoff( 90),voff( 90),(pattn(i, 90),i=1,16)/ -1,  0,
     *   1023, 1023,  774,  782,  780,   28,   24,   56,
     *    112,   96,  224,  192,  454,  390, 1022, 1022/
      data hoff( 91),voff( 91),(pattn(i, 91),i=1,16)/ -5, -2,
     *      3,    3,    3,    3,    3,    3,    3,    3,
     *      3,    3,    3,    3,    3,    3,   63,   63/
      data hoff( 92),voff( 92),(pattn(i, 92),i=1,16)/ -1, -2,
     *    448,  192,  224,   96,  112,   48,   48,   56,
     *     24,   28,   12,   14,    6,    7,    3,    3/
      data hoff( 93),voff( 93),(pattn(i, 93),i=1,16)/ -1, -2,
     *     48,   48,   48,   48,   48,   48,   48,   48,
     *     48,   48,   48,   48,   48,   48,   63,   63/
      data hoff( 94),voff( 94),(pattn(i, 94),i=1,16)/ -2,-12,
     *    231,  255,   60,   24,    0,    0,    0,    0,
     *      0,    0,    0,    0,    0,    0,    0,    0/
      data hoff( 95),voff( 95),(pattn(i, 95),i=1,16)/ -1,  2,
     *   1023, 1023,    0,    0,    0,    0,    0,    0,
     *      0,    0,    0,    0,    0,    0,    0,    0/
      data hoff( 96),voff( 96),(pattn(i, 96),i=1,16)/ -5,-10,
     *      7,    7,    7,    3,    7,    6,    0,    0,
     *      0,    0,    0,    0,    0,    0,    0,    0/
      data hoff( 97),voff( 97),(pattn(i, 97),i=1,16)/ -1,  0,
     *   1980, 2047,  231,  195,  207,  254,  248,  192,
     *    238,  126,   62,    0,    0,    0,    0,    0/
      data hoff( 98),voff( 98),(pattn(i, 98),i=1,16)/  0,  0,
     *    508, 1020,  796, 1820, 1548, 1548, 1548, 1804,
     *    796, 1020,  508,   12,   12,   12,   15,   15/
      data hoff( 99),voff( 99),(pattn(i, 99),i=1,16)/ -2,  0,
     *    248,  508,  398,    7,    3,    3,    3,    7,
     *    462,  508,  504,    0,    0,    0,    0,    0/
      data hoff(100),voff(100),(pattn(i,100),i=1,16)/ -1,  0,
     *   2044, 2046,  454,  455,  387,  387,  387,  391,
     *    454,  510,  508,  384,  384,  384,  480,  480/
      data hoff(101),voff(101),(pattn(i,101),i=1,16)/ -2,  0,
     *    248,  510,  398,    7,    3,  511,  511,  391,
     *    454,  254,  124,    0,    0,    0,    0,    0/
      data hoff(102),voff(102),(pattn(i,102),i=1,16)/  0,  0,
     *    510,  510,   48,   48,   48,   48,   48,   48,
     *     48, 1023, 1023,   48,   48,  944, 1008,  992/
      data hoff(103),voff(103),(pattn(i,103),i=1,16)/ -1,  5,
     *    252,  510,  903,  771,  903,  510,  254,    6,
     *     62,  126,  238,  198,  198, 2030, 2044, 2040/
      data hoff(104),voff(104),(pattn(i,104),i=1,16)/  0,  0,
     *   3999, 3999,  780,  780,  780,  780,  780,  780,
     *    924, 1020,  508,   12,   12,   12,   15,   15/
      data hoff(105),voff(105),(pattn(i,105),i=1,16)/ -1,  0,
     *    511,  511,   48,   48,   48,   48,   48,   48,
     *     48,   62,   62,    0,    0,   56,   56,   56/
      data hoff(106),voff(106),(pattn(i,106),i=1,16)/ -2,  0,
     *     48,   48,   48,   48,   48,   48,   48,   48,
     *     48,   62,   62,    0,    0,   56,   56,   56/
      data hoff(107),voff(107),(pattn(i,107),i=1,16)/  0,  0,
     *   3903, 3903,  780,  396,  396,  252,  108,  204,
     *    396, 4044, 4044,   12,   12,   12,   15,   15/
      data hoff(108),voff(108),(pattn(i,108),i=1,16)/ -1,  0,
     *   1023, 1023,   48,   48,   48,   48,   48,   48,
     *     48,   48,   48,   48,   48,   48,   63,   63/
      data hoff(109),voff(109),(pattn(i,109),i=1,16)/  0,  0,
     *   3687, 3687, 1638, 1638, 1638, 1638, 1638, 1638,
     *   1774, 2047,  959,    0,    0,    0,    0,    0/
      data hoff(110),voff(110),(pattn(i,110),i=1,16)/  0,  0,
     *   3999, 3999,  780,  780,  780,  780,  780,  780,
     *    924, 1023,  511,    0,    0,    0,    0,    0/
      data hoff(111),voff(111),(pattn(i,111),i=1,16)/ -2,  0,
     *    124,  254,  198,  455,  387,  387,  387,  455,
     *    198,  254,  124,    0,    0,    0,    0,    0/
      data hoff(112),voff(112),(pattn(i,112),i=1,16)/  0,  5,
     *     63,   63,   12,   12,   12,  508, 1020,  796,
     *   1820, 1548, 1548, 1548, 1804,  796, 1023,  511/
      data hoff(113),voff(113),(pattn(i,113),i=1,16)/ -1,  5,
     *   4032, 4032,  768,  768,  768,  888, 1022,  974,
     *    903,  771,  771,  771,  903,  974, 1020,  888/
      data hoff(114),voff(114),(pattn(i,114),i=1,16)/ -1,  0,
     *    127,  127,   12,   12,   12,   12,   12,   28,
     *   1852, 2047, 2031,    0,    0,    0,    0,    0/
      data hoff(115),voff(115),(pattn(i,115),i=1,16)/ -2,  0,
     *    127,  255,  455,  387,  499,  254,   31,  195,
     *    199,  255,  252,    0,    0,    0,    0,    0/
      data hoff(116),voff(116),(pattn(i,116),i=1,16)/ -1,  0,
     *    248,  508,  396,  396,  396,   12,   12,   12,
     *     12,  511,  511,   12,   12,   12,   12,    0/
      data hoff(117),voff(117),(pattn(i,117),i=1,16)/  0,  0,
     *   4088, 4092,  908,  780,  780,  780,  780,  780,
     *    780,  975,  975,    0,    0,    0,    0,    0/
      data hoff(118),voff(118),(pattn(i,118),i=1,16)/ -1,  0,
     *    112,  112,  248,  248,  152,  408,  396,  396,
     *    268, 2015, 2015,    0,    0,    0,    0,    0/
      data hoff(119),voff(119),(pattn(i,119),i=1,16)/  0,  0,
     *   1020, 1020, 1020,  956,  676,  676,  676, 1702,
     *   1766, 3823, 3663,    0,    0,    0,    0,    0/
      data hoff(120),voff(120),(pattn(i,120),i=1,16)/ -1,  0,
     *   2015, 2015,  396,  216,  112,  112,  112,  216,
     *    396, 2015, 2015,    0,    0,    0,    0,    0/
      data hoff(121),voff(121),(pattn(i,121),i=1,16)/ -1,  5,
     *     14,   31,   55,   55,   96,   96,   96,  240,
     *    240,  152,  408,  408,  396,  268, 2015, 2015/
      data hoff(122),voff(122),(pattn(i,122),i=1,16)/ -1,  0,
     *   1023, 1023,  782,  796,   56,   48,  112,  227,
     *    451, 1023, 1023,    0,    0,    0,    0,    0/
      data hoff(123),voff(123),(pattn(i,123),i=1,16)/ -1, -2,
     *     48,   48,   48,   48,   56,   31,   31,   56,
     *     48,   48,   48,   48,   48,  112,  992,  960/
      data hoff(124),voff(124),(pattn(i,124),i=1,16)/ -5, -2,
     *      3,    3,    3,    3,    3,    3,    3,    3,
     *      3,    3,    3,    3,    3,    3,    3,    3/
      data hoff(125),voff(125),(pattn(i,125),i=1,16)/ -1, -2,
     *     48,   48,   48,   48,  112,  992,  992,  112,
     *     48,   48,   48,   48,   48,   56,   31,   15/
      data hoff(126),voff(126),(pattn(i,126),i=1,16)/ -2,-13,
     *    123,  255,  206,    0,    0,    0,    0,    0,
     *      0,    0,    0,    0,    0,    0,    0,    0/
      data hoff(127),voff(127),(pattn(i,127),i=1,16)/ -2,-13,
     *    231,  231,  231,    0,    0,    0,    0,    0,
     *      0,    0,    0,    0,    0,    0,    0,    0/
