|Algorithm #|4th Ed. Page # (PDF)|          Book Name           |                    Function Name                      | Test Program | Example Page # (If Applicable) | Verified Working |
| :---: | :---: | :---: | :---: | :---: | :---: |:---: |
|  1   |  pg 63 (90)    |  Find C2C3                                   | findc2c3() (used within lambertu())  | ex7_5.py & testlam.py (dungeon) |           -                |  Y/N  |
|  2   |  pg 65 (92)    | KepEqtnE                                     | newtonm()                            | ex2_1.py                        | pg 66 (93)                 |   Y   |
|  3   |  pg 69 (96)    | KepEqtnP                                     | newtonm()                            | ex2_2.py                        | pg 69 (96)                 |   Y   |
|  4   |  pg 71 (98)    | KepEqtnH                                     | newtonm()                            | ex2_3.py                        | pg 71 (98)                 |   Y   |
|  5   |  pg 77 (104)   | v to Anomaly                                 | newtonnu()                           | ex2_2.py                        |           -                |   Y   |
|  6   |  pg 77 (104)   | Anomaly to v                                 | newtone()                            | ex2_2.py                        |           -                |   Y   |
|  7   |  pg 81 (108)   | KeplerCOE                                    | N/A                                  |                -                |           -                |   -   |
|  8   |  pg 93 (120)   | Kepler                                       | kepler()                             | ex2_4.py                        | pg 94 (121)                |   Y   |
|  9   |  pg 113 (140)  | RV2COE                                       | rv2coe(), rv2coeh(), rv2coeS()       | ex2_5.py                        | pg 114 (141)               |   Y   |
|  10  |  pg 118 (145)  | COE2RV                                       | coe2rv(), coe2rvh(), coe2rvS()       | ex2_6.py                        | pg 119 (146)               |   Y   |
|  11  |  pg 126 (153)  | FindTOF                                      | findtof()                            |                -                |           -                |   N   |
|  12  |  pg 172 (199)  | ECEF To LatLon                               | ecef2ll()                            | ex3_3.py                        | pg 173 (200)               |   Y   |
|  13  |  pg 173 (200)  | ECEF To LatLon (B)                           | ecef2llb()                           | ex3_3.py                        | pg 173 (200)               |   Y   |
|  14  |  pg 183 (210)  | Julian Date                                  | jday()                               | ex3_4.py                        | pg 184 (211)               |   Y   |
|  15  |  pg 188 (215)  | LSTime                                       | lstime()                             | ex3_5.py & ex3_6.py             | pg 188-189 (215-216)       |   Y   |
|  16  |  pg 195 (222)  | ConvTime (Calculate Dynamic Time)            | convtime()                           | ex3_7.py                        | pg 195 (222)               |   Y   |
|  17  |  pg 197 (224)  | DMStoRAD                                     | dms2rad()                            | ex3_8.py                        | pg 198 (225)               |   Y   |
|  18  |  pg 197 (224)  | RADtoDMS                                     | rad2dms()                            | ex3_8.py                        | pg 198 (225)               |   Y   |
|  19  |  pg 198 (225)  | HMStoRAD                                     | hms2rad()                            | ex3_9.py                        | pg 199 (226)               |   Y   |
|  20  |  pg 198 (225)  | RADtoHMS                                     | rad2hms()                            | ex3_9.py                        | pg 199 (226)               |   Y   |
|  21  |  pg 199 (226)  | TimeToHMS                                    | sec2hms()                            | ex3_10.py                       | pg 199 (226)               |   Y   |
|  22  |  pg 202 (229)  | JDtoGregorianDate                            | invjday()                            | ex3_13.py                       | pg 203 (230)               |   Y   |
|  23  |  pg 219 (246)  | IAU-2000, CIO-Based                          | ecef2eci() (IAU-2006)                | ex3_14.py                       | pg 220 (247)               |   -   |
|  24  |  pg 228 (255)  | IAU-76/FK5 Reduction                         | eci2ecef()                           | ex3_15.py                       | pg 230 (257)               |   -   |
|  25  |  pg 259 (286)  | Geocentric RaDec                             | rv2radec()                           | ex4_1.py                        | pg 273 (300)               |   Y   |
|  26  |  pg 260 (287)  | Topocentric                                  | rv2tradec()                          | ex4_1.py                        | pg 273 (300)               |   Y   |
|  27  |  pg 265 (292)  | RAZEL                                        | rv2razel()                           | ex4_1.py                        | pg 273 (300)               |   N   |
|  28  |  pg 268 (295)  | AzELToRaDec                                  | azel2radec()                         |                -                |           -                |   -   |
|  29  |  pg 279 (306)  | Sun                                          | sun(), sunalmanac()                  | ex5_1.py                        | pg 280 (307)               |   Y   |
|  30  |  pg 283 (310)  | SunriseSet                                   | sunriset()                           | ex5_2.py                        | pg 284 (311)               |   Y   |
|  31  |  pg 288 (315)  | Moon                                         | moon()                               | ex5_3.py                        | pg 288 (315)               |   Y   |
|  32  |  pg 290 (317)  | MoonRiseSet                                  | moonrise(), moonrise2(), moonrise3() | ex5_4.py                        | pg 292 (319)               |   Y   |
|  33  |  pg 296 (323)  | PlanetRV                                     | N/A                                  | ex5_5.py                        | pg 297 (324)               |   -   |
|  34  |  pg 301 (328)  | Shadow                                       | shadow()                             | shadowtest.py                   |           -                |   N   |
|  35  |  pg 308 (335)  | Sight                                        | sight(), light()                     | ex5_6.py                        | pg 309 (336)               |   Y   |
|  36  |  pg 325 (352)  | Hohmann Transfer                             | hohmann()                            | ex6_1.py                        | pg 326 (353)               |   Y   |
|  37  |  pg 326 (353)  | Bi-elliptic Transfer                         | biellip()                            | ex6_2.py                        | pg 327 (354)               |   Y   |
|  38  |  pg 333 (360)  | One-Tangent Burn                             | onetang()                            | ex6_3.py                        | pg 334 (360)               |   Y   |
|  39  |  pg 344 (371)  | Inclination Only                             | ionlychg()                           | ex6_4.py                        | pg 344 (371)               |   Y   |
|  40  |  pg 347 (374)  | Change in Ascending Node-Circular            | nodeonly()                           | ex6_5.py                        | pg 347 (374)               |   Y   |
|  41  |  pg 348 (375)  | Combined Changes (i & Omega)-Circular        | iandnode()                           | ex6_6.py                        | pg 349 (376)               |   Y   |
|  42  |  pg 353 (380)  | Minimum Combined Plane Change                | combined()                           | ex6_7.py                        | pg 356 (383)               |   Y   |
|  43  |  pg 356 (383)  | Fixed-delta v Maneuvers                      | combined()                           | ex6_7.py                        | pg 356 (383)               |   Y   |
|  44  |  pg 362 (389)  | Circular Coplanar Phasing (Same Orbits)      | rendz()                              | ex6_8.py                        | pg 363 (390)               |   Y   |
|  45  |  pg 363 (390)  | Circular Coplanar Phasing (Different Orbits) | rendz()                              | ex6_9.py                        | pg 364 (391)               |   Y   |
|  46  |  pg 368 (395)  | Noncoplanar Phasing                          | noncoplr()                           | ex6_10.py                       | pg 369 (396)               |   Y   |
|  47  |  pg 387 (414)  | Low Thrust Transfer                          | lowthrust()                          | ex6_12.py & ex6_13.py           | pg 381 (408) & 388 (415)   |   -   |
|  48  |  pg 396 (423)  | Hill's Equation                              | hillsr(), hillsv()                   | ex6_14.py and ex6_15.py         | pg 397 (424) & 410 (437)   |   Y   |
|  49  |  pg 414 (441)  | HillEQCM to ECI                              | EQCM_to_ECI_RTN_sal()  (dungeon)     | EQCM_to_ECI_RTN_sal.py          |           -                |   N   |
|  50  |  pg 415 (442)  | ECI to HillEQCM                              | ECI_to_EQCM_RTN_sal()  (dungeon)     | ECI_to_EQCM_RTN_sal.py          |           -                |   N   |
|  51  |  pg 430 (457)  | Site-Track                                   | site()                               | ex7_1.py                        | pg 431 (458)               |   -   |
|  52  |  pg 442 (469)  | Angles-Only Gauss                            | anglesg()                            | ex7_2.py                        | pg 447 (474)               |   N   |
|  53  |  pg 460 (487)  | Angles-Double-R                              | anglesdr()                           | ex7_2.py                        | pg 447 (474)               |   N   |
|  54  |  pg 460 (487)  | Gibbs                                        | gibbs(), gibbsh()                    | ex7_34.py                       | pg 461 (488)               |   Y   |
|  55  |  pg 466 (493)  | Herrick-Gibbs                                | hgibbs()                             | ex7_34.py                       | pg 467 (494)               |   Y   |
|  56  |  pg 475 (502)  | Lambert's Problem - Minimum Energy           | lambertmin()                         | ex7_5.py & testlam.py (dungeon) | pg 497 (524)               |  Y/N  |
|  57  |  pg 478 (505)  | Lambert - Gauss's Solution                   | N/A                                  |                -                |           -                |   -   |
|  58  |  pg 492 (519)  | Lambert - Universal Variables                | lambertu()                           | ex7_5.py & testlam.py (dungeon) | pg 497 (524)               |  Y/N  |
|  59  |  pg 494 (521)  | Lambert - Battin Method                      | lambertb()                           | ex7_5.py & testlam.py (dungeon) | pg 497 (524)               |  Y/N  |
|  60  |  pg 503 (530)  | Hit Earth                                    | checkhitearth()                      |                -                |           -                |   N   |
|  61  |  pg 503 (530)  | Target                                       | target()                             |                -                |           -                |   N   |
|  62  |  pg 525 (552)  | ENCKE                                        | N/A                                  |                -                |           -                |   -   |
|  63  |  pg 558 (585)  | Ap2Kp                                        | ap2kp()                              |                -                |           -                |   N   |
|  63  |  pg 558 (585)  | Ap2Kp                                        | kp2ap()                              |                -                |           -                |   N   |
|  64  |  pg 591 (618)  | Numerical Integration                        | N/A                                  |                -                |           -                |   -   |
|  65  |  pg 691 (718)  | PKepler                                      | pkepler()                            | ex10_4.py & ex10_5.py           | pg 769 (796) & 775 (802)   |  Y/N  |
|  66  |  pg 766 (793)  | Nominal State                                | N/A                                  | ex10_4.py?                      | pg 769 (796)               |   -   |
|  67  |  pg 768 (795)  | Differential Correction                      | N/A                                  | ex10_4.py                       | pg 769 (796)               |   -   |
|  68  |  pg 785 (812)  | Kalman Filter - Linear System                | N/A                                  | ex10_6.py                       | pg 786 (813)               |   -   |
|  69  |  pg 791 (818)  | Linearized                                   | N/A                                  | ex10_6.py                       | pg 786 (813)               |   -   |
|  70  |  pg 793 (820)  | Extended Kalman Filter                       | N/A                                  | ex10_6.py                       | pg 786 (813)               |   -   |
|  71  |  pg 873 (900)  | Repeat Ground Track                          | repeatgt.py                          | ex11_3.py                       |       -                    |   N   |
|  72  |  pg 879 (906)  | Main Repeat Groundtrack                      | N/A                                  | ex11_4.py                       | pg 893 (920)               |   -   |
|  73  |  pg 885 (912)  | Minimum Altitude Variation                   | N/A                                  | ex11_5.py                       | pg 895 (922)               |   -   |
|  74  |  pg 911 (938)  | Predict                                      | predict()                            | ex11_6.py                       | pg 912 (939)               |   N   |
|  75  |  pg 916 (943)  | Rise/Set                                     | N/A                                  | ex11_7.py (DNE)                 | pg 918 (945)               |   -   |
|  76  |  pg 935 (962)  | Close Approach                               | N/A                                  |                -                |           -                |   -   |
|  77  |  pg 953 (980)  | Patched Conic                                | N/A                                  |                -                |           -                |   -   |
|  78  |  pg 962 (989)  | Algorithm B-plane I                          | N/A                                  |                -                |           -                |   -   |
|  79  |  pg 963 (990)  | Algorithm B-plane II                         | N/A                                  |                -                |           -                |   -   |

