from Pygrad.Pygrad cimport Pygrad
import numpy as np
#This file takes the place of the CASCDAT subroutine.  It defines several 
#arrays related to cascade data.
cpdef loadData(Pygrad object):
    object.cascdata['NSD'] = [2,2,2,4,2,2,4,4,6,2,2,4,4,6,2,2,4]
    object.cascdata['SCRPT'] = [' K  ',' L1 ',' L2 ',' L3 ',' M1',' M2 ',' M3 ',' M4 ', ' M5 ',' N1 ',' N2 ',' N3 ',' N4 ',' N5 ',' O1 ',' O2 ',' O3 ']
    object.cascdata['SCRPT1'] = [' 1s ',' 2s ',' 2p1/2',' 2p3/2',' 3s ',' 3p1/2',' 3p3/2',' 3d3/2',' 3d5/2',' 4s ',' 4p1/2',' 4p3/2',' 4d3/2',' 4d5/2',' 5s ',' 5p1/2',' 5p3/2']
    object.cascdata['AA'] = [0.0,0.0,0.25,0.25,0.0,0.25,0.25,0.50,0.50,0.0,0.25,0.25,0.50,0.50,0.0,0.25,0.25]
    object.cascdata['BB'] = [1.5,1.5,1.25,1.25,1.5,1.25,1.25,0.75,0.75,1.5,1.25,1.25,0.75,0.75,1.5,1.25,1.25]
    object.cascdata['ELEV'] = np.array([[13.598, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 24.587, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 54.7, 5.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 111.5, 9.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 188.0, 12.6, 4.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 284.2, 18.0, 6.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 401.6, 24.4, 14.534, 14.524, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 532.0, 28.5, 13.618, 13.618, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 685.4, 34.0, 16.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 870.2, 48.475, 21.661, 21.565, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1070.8, 63.5, 30.65, 30.81, 5.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1303.0, 88.7, 49.78, 49.5, 7.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1559.6, 117.8, 72.95, 72.55, 10.6, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1839.0, 149.7, 99.82, 99.42, 13.5, 8.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2145.5, 189.0, 136.0, 135.0, 16.1, 10.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2472.0, 230.9, 163.6, 162.5, 20.2, 10.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2822.4, 270.0, 202.0, 200.0, 24.5, 12.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3205.9, 326.3, 250.6, 248.4, 29.239, 15.937, 15.76, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3608.4, 378.6, 297.3, 294.6, 34.8, 18.3, 18.2, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4038.5, 438.4, 349.7, 346.2, 44.3, 25.4, 25.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4492.0, 498.0, 403.6, 398.7, 51.1, 28.3, 28.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4966.0, 560.9, 460.2, 453.8, 58.7, 32.6, 32.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5465.0, 626.7, 519.8, 512.1, 66.3, 37.2, 37.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5989.0, 696.0, 583.8, 574.1], [74.1, 42.2, 42.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6539.0, 769.1, 649.9, 638.7, 82.3, 47.2, 47.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7112.0, 844.6, 719.9, 706.8, 91.3, 52.7, 52.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7709.0, 925.1, 793.2, 778.1, 101.0, 59.8, 58.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8333.0, 1008.6, 870.0, 852.7, 110.8, 68.0, 66.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 8979.0, 1096.7, 952.3, 932.7, 122.5, 77.3, 75.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9659.0, 1196.2, 1044.9, 1021.8, 139.8, 91.4, 88.6, 10.2, 10.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10367.0, 1299.0, 1143.2, 1116.4, 159.5, 103.5, 100.0, 18.7, 18.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11103.0, 1414.6, 1248.1, 1217.0, 180.1, 124.9, 120.8, 29.8, 29.2, 14.3, 7.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11867.0, 1527.0, 1359.1, 1323.6, 204.7, 146.2, 141.2, 41.7, 41.7], [17.0, 9.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12658.0, 1652.0, 1474.3, 1433.9, 229.6, 166.5, 160.7, 55.5, 54.6, 20.1, 9.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13474.0, 1782.0, 1596.0, 1550.0, 257.0, 189.0, 182.0, 70.0, 69.0, 23.8, 11.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14327.26, 1921.0, 1730.9, 1678.4, 292.8, 222.2, 214.4, 95.0, 93.8, 27.5, 14.666, 13.9996, 0.0, 0.0, 0.0, 0.0, 0.0, 15200.0, 2065.0, 1864.0, 1804.0, 326.7, 248.7, 239.1, 113.0, 112.0, 30.5, 16.3, 15.3, 0.0, 0.0, 0.0, 0.0, 0.0, 16105.0, 2216.0, 2007.0], [1940.0, 358.7, 280.3, 270.0, 136.0, 134.2, 38.9, 21.3, 20.1, 0.0, 0.0, 0.0, 0.0, 0.0, 17038.0, 2373.0, 2156.0, 2080.0, 392.0, 310.6, 298.8, 157.7, 155.8, 43.8, 24.4, 23.1, 0.0, 0.0, 0.0, 0.0, 0.0, 17998.0, 2532.0, 2307.0, 2223.0, 430.3, 343.5, 329.8, 181.1, 178.8, 50.6, 28.5, 27.1, 0.0, 0.0, 0.0, 0.0, 0.0, 4118986.0, 2698.0, 2465.0, 2371.0, 466.6, 376.1, 360.6, 205.0, 202.3, 56.4, 32.6, 30.8, 0.0, 0.0, 0.0, 0.0, 0.0, 20000.0, 2866.0, 2625.0, 2520.0, 506.3, 411.6, 394.0, 231.1, 227.9, 63.2, 37.6, 35.5, 0.0, 0.0], [0.0, 0.0, 0.0, 21044.0, 3043.0, 2793.0, 2677.0, 544.0, 447.6, 417.7, 257.6, 253.9, 69.5, 42.3, 39.9, 0.0, 0.0, 0.0, 0.0, 0.0, 22117.0, 3224.0, 2967.0, 2838.0, 586.1, 483.5, 461.4, 284.2, 280.0, 75.0, 46.3, 43.2, 0.0, 0.0, 0.0, 0.0, 0.0, 23220.0, 3412.0, 3146.0, 3004.0, 628.1, 521.3, 496.5, 311.9, 307.2, 81.4, 50.5, 47.3, 0.0, 0.0, 0.0, 0.0, 0.0, 24350.0, 3604.0, 3330.0, 3173.0, 671.6, 559.9, 532.3, 340.5, 335.2, 87.1, 55.7, 50.9, 0.0, 0.0, 0.0, 0.0, 0.0, 25514.0, 3806.0, 3524.0, 3351.0, 719.0, 603.8, 573.0, 374.0], [368.3, 97.0, 63.7, 58.3, 0.0, 0.0, 0.0, 0.0, 0.0, 26711.0, 4018.0, 3727.0, 3538.0, 772.0, 652.6, 618.4, 411.9, 405.2, 109.8, 63.9, 63.8, 11.7, 10.7, 0.0, 0.0, 0.0, 27940.0, 4238.0, 3938.0, 3730.0, 827.2, 703.2, 665.3, 451.4, 443.9, 122.9, 73.6, 73.5, 17.7, 16.9, 0.0, 0.0, 0.0, 29200.0, 4465.0, 4156.0, 3929.0, 884.7, 756.5, 714.6, 493.2, 484.9, 137.1, 83.6, 83.5, 24.9, 23.9, 0.0, 0.0, 0.0, 30491.0, 4698.0, 4380.0, 4132.0, 946.0, 812.7, 766.4, 537.5, 528.2, 153.2, 95.6, 95.5, 33.3, 32.1, 0.0, 0.0, 0.0, 31814.0, 4939.0], [4612.0, 4341.0, 1006.0, 870.8, 820.0, 583.4, 573.0, 169.4, 103.3, 103.2, 41.9, 40.4, 0.0, 0.0, 0.0, 33169.0, 5188.0, 4852.0, 4557.0, 1072.0, 931.0, 875.0, 630.8, 619.3, 186.0, 123.0, 122.9, 50.6, 48.9, 0.0, 0.0, 0.0, 34561.0, 5453.0, 5107.0, 4786.0, 1148.7, 1002.1, 940.6, 689.0, 676.4, 213.2, 146.7, 145.5, 69.5, 67.5, 23.3, 13.43, 12.129843, 35985.0, 5714.0, 5359.0, 5012.0, 1211.0, 1071.0, 1003.0, 740.5, 726.6, 232.3, 172.4, 161.3, 79.8, 77.5, 23.7, 14.2, 12.6, 37441.0, 5989.0, 5624.0, 5247.0, 1293.0, 1137.0, 1063.0, 795.7, 780.5, 253.5, 192.0, 178.6, 92.6], [89.9, 30.3, 17.0, 14.8, 38925.0, 6266.0, 5891.0, 5483.0, 1362.0, 1209.0, 1128.0, 853.0, 836.0, 274.7, 205.8, 196.0, 105.3, 102.5, 34.3, 19.3, 16.8, 40443.0, 6549.0, 6164.0, 5723.0, 1436.0, 1274.0, 1187.0, 902.4, 883.8, 291.0, 223.2, 206.5, 109.0, 107.0, 37.2, 19.8, 17.0, 41991.0, 6835.0, 6440.0, 5964.0, 1511.0, 1337.0, 1242.0, 948.3, 928.8, 304.5, 236.3, 217.6, 115.1, 115.0, 37.4, 21.0, 20.9, 43569.0, 7126.0, 6722.0, 6208.0, 1575.0, 1403.0, 1297.0, 1003.3, 980.4, 319.2, 243.3, 224.6, 120.5, 120.4, 37.5, 21.1, 21.0, 45184.0, 7428.0, 7013.0, 6459.0, 1650.0, 1471.0, 1357.0], [1052.0, 1027.0, 332.0, 251.0, 231.0, 124.0, 123.0, 37.6, 21.4, 21.3, 46834.0, 7737.0, 7312.0, 6716.0, 1723.0, 1541.0, 1420.0, 1110.9, 1083.4, 347.2, 265.6, 247.4, 128.0, 127.0, 37.7, 21.4, 21.3, 48519.0, 8052.0, 7617.0, 6977.0, 1800.0, 1614.0, 1481.0, 1158.6, 1127.5, 360.0, 284.0, 257.0, 132.0, 127.7, 37.8, 22.0, 21.9, 50239.0, 8376.0, 7930.0, 7243.0, 1881.0, 1688.0, 1544.0, 1221.9, 1189.6, 378.6, 286.0, 271.0, 143.0, 142.6, 36.0, 28.0, 22.0, 51996.0, 8708.0, 8252.0, 7514.0, 1968.0, 1768.0, 1611.0, 1276.9, 1241.1, 396.0, 322.4, 284.1, 150.5, 150.4, 45.6, 28.7, 22.6, 53789.0], [9046.0, 8581.0, 7790.0, 2047.0, 1842.0, 1676.0, 1333.0, 1292.6, 414.2, 333.5, 293.2, 153.6, 153.5, 48.9, 29.5, 23.3, 55618.0, 9394.0, 8918.0, 8071.0, 2128.0, 1923.0, 1741.0, 1392.0, 1351.0, 432.4, 343.5, 308.2, 160.1, 160.0, 49.3, 30.8, 24.1, 57486.0, 9751.0, 9264.0, 8358.0, 2207.0, 2006.0, 1812.0, 1453.0, 1409.0, 449.8, 366.2, 320.2, 167.6, 167.5, 50.6, 31.4, 24.7, 59390.0, 10116.0, 9617.0, 8648.0, 2307.0, 2090.0, 1885.0, 1515.0, 1468.0, 470.9, 385.9, 332.6, 175.5, 175.4, 54.7, 31.8, 25.0, 61332.0, 10486.0, 9978.0, 8944.0, 2398.0, 2173.0, 1950.0, 1576.0, 1528.0, 480.5, 388.7, 339.7], [191.2, 182.4, 55.0, 32.5, 25.8, 63314.0, 10870.0, 10349.0, 9244.0, 2491.0, 2264.0, 2024.0, 1639.0, 1589.0, 506.8, 412.4, 359.2, 206.1, 196.3, 57.3, 33.6, 26.7, 65351.0, 11271.0, 10739.0, 9561.0, 2601.0, 2365.0, 2108.0, 1716.0, 1662.0, 538.0, 438.2, 380.7, 220.0, 211.5, 64.2, 38.0, 29.9, 67416.0, 11682.0, 11136.0, 9881.0, 2708.0, 2469.0, 2194.0, 1793.0, 1735.0, 563.4, 463.4, 400.9, 237.9, 226.4, 69.7, 42.2, 32.7, 69525.0, 12100.0, 11544.0, 10207.0, 2820.0, 2575.0, 2281.0, 1872.0, 1809.0, 594.1, 490.4, 423.6, 255.9, 243.5, 75.6, 45.3, 36.8, 71676.0, 12527.0, 11959.0, 10535.0, 2932.0, 2682.0], [2367.0, 1949.0, 1883.0, 625.4, 518.7, 446.8, 273.9, 260.5, 83.0, 45.6, 38.0, 73871.0, 12968.0, 12385.0, 10871.0, 3049.0, 2792.0, 2457.0, 2031.0, 1960.0, 658.2, 549.1, 470.7, 293.1, 278.5, 84.0, 58.0, 45.0, 76111.0, 13419.0, 12824.0, 11215.0, 3174.0, 2909.0, 2551.0, 2116.0, 2040.0, 691.1, 577.8, 495.8, 311.9, 296.3, 95.2, 63.0, 49.0, 78395.0, 13880.0, 13273.0, 11564.0, 3296.0, 3027.0, 2645.0, 2202.0, 2122.0, 725.4, 609.1, 519.4, 331.6, 314.6, 101.7, 65.3, 52.0, 80725.0, 14353.0, 13734.0, 11919.0, 3425.0, 3148.0, 2743.0, 2291.0, 2206.0, 762.1, 642.7, 546.3, 353.2, 335.1, 107.2, 74.2, 57.2]])



