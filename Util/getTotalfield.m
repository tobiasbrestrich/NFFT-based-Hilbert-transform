function [a0, b0, c0, X0, Y0, Z0] = getTotalfield(year, level)
%GETTOTALFIELD Returns the geomagnetic total field around Kirchberg, Saxony.
%Calculated according to the International Geomagnetic Reference Field.
%   year       date (in decimals : yyy.yy) of the calculated field
%   level       
%RETURN:
%   a0          transform parameter for the x-component  
%   b0          transform parameter for the y-component    
%   c0          transform parameter for the z-component  
%   X0          x-component of the total field  
%   Y0          y-component of the total field    
%   Z0          z-component of the total field    


assert(year > 1964 & year < 2018, 'Year must be 1964 < year < 2018.');
% data = textread('igrfwmmData.csv','', ...
%     'commentstyle', 'shell', 'delimiter',',');

if nargin < 2
    level = 0;
end
if level == 0
    data = [
        1965.00000	-1.60732	65.90598	19412.8	47553.0	19405.2	-544.5	43410.0
        1966.00000	-1.56275	65.89852	19429.9	47581.0	19422.7	-529.9	43433.0
        1967.00000	-1.51827	65.89105	19447.0	47609.0	19440.2	-515.3	43456.1
        1968.00000	-1.47385	65.88358	19464.1	47637.0	19457.7	-500.6	43479.1
        1969.00000	-1.42952	65.87610	19481.2	47665.0	19475.2	-486.0	43502.1
        1970.00000	-1.38527	65.86862	19498.4	47693.1	19492.7	-471.4	43525.2
        1971.00000	-1.34304	65.86180	19518.6	47729.8	19513.2	-457.5	43556.3
        1972.00000	-1.30090	65.85497	19538.8	47766.5	19533.7	-443.6	43587.5
        1973.00000	-1.25884	65.84814	19559.0	47803.2	19554.3	-429.7	43618.7
        1974.00000	-1.21687	65.84131	19579.2	47839.9	19574.8	-415.8	43649.8
        1975.00000	-1.17499	65.83448	19599.4	47876.6	19595.3	-401.9	43681.0
        1976.00000	-1.05209	65.83963	19607.9	47906.9	19604.6	-360.0	43710.4
        1977.00000	-0.92930	65.84468	19616.5	47937.3	19613.9	-318.2	43739.8
        1978.00000	-0.80662	65.84963	19625.1	47967.7	19623.2	-276.3	43769.3
        1979.00000	-0.68404	65.85447	19633.9	47998.1	19632.5	-234.4	43798.7
        1980.00000	-0.56158	65.85920	19642.7	48028.6	19641.8	-192.5	43828.1
        1981.00000	-0.45337	65.86900	19642.9	48047.3	19642.3	-155.4	43848.6
        1982.00000	-0.34517	65.87870	19643.1	48066.0	19642.8	-118.3	43869.0
        1983.00000	-0.23698	65.88833	19643.4	48084.8	19643.3	-81.2	43889.5
        1984.00000	-0.12878	65.89787	19643.8	48103.6	19643.8	-44.2	43909.9
        1985.00000	-0.02059	65.90732	19644.3	48122.5	19644.3	-7.1	43930.4
        1986.00000	0.05690	65.92819	19637.6	48145.5	19637.6	19.5	43958.5
        1987.00000	0.13444	65.94901	19631.0	48168.5	19631.0	46.1	43986.6
        1988.00000	0.21203	65.96976	19624.5	48191.5	19624.4	72.6	44014.8
        1989.00000	0.28968	65.99046	19618.0	48214.6	19617.7	99.2	44042.9
        1990.00000	0.36737	66.01109	19611.5	48237.7	19611.1	125.7	44071.1
        1991.00000	0.45564	66.02544	19608.7	48257.9	19608.1	155.9	44094.5
        1992.00000	0.54394	66.03972	19605.9	48278.2	19605.0	186.1	44117.9
        1993.00000	0.63226	66.05393	19603.2	48298.5	19602.0	216.3	44141.3
        1994.00000	0.72061	66.06809	19600.5	48318.8	19599.0	246.5	44164.7
        1995.00000	0.80898	66.08218	19597.9	48339.1	19596.0	276.7	44188.1
        1996.00000	0.91241	66.09311	19603.6	48374.0	19601.1	312.2	44223.7
        1997.00000	1.01578	66.10394	19609.4	48408.9	19606.3	347.6	44259.4
        1998.00000	1.11908	66.11470	19615.2	48443.8	19611.5	383.1	44295.0
        1999.00000	1.22233	66.12537	19621.1	48478.8	19616.7	418.6	44330.6
        2000.00000	1.32551	66.13595	19627.1	48513.8	19621.9	454.0	44366.2
        2001.00000	1.42864	66.14563	19633.5	48548.0	19627.4	489.5	44400.8
        2002.00000	1.53170	66.15522	19639.9	48582.2	19632.8	525.0	44435.4
        2003.00000	1.63470	66.16473	19646.3	48616.4	19638.3	560.4	44470.0
        2004.00000	1.73762	66.17416	19652.8	48650.7	19643.8	595.9	44504.6
        2005.00000	1.84048	66.18351	19659.4	48685.0	19649.3	631.4	44539.2
        2006.00000	1.95677	66.18810	19668.0	48715.1	19656.5	671.6	44568.2
        2007.00000	2.07297	66.19259	19676.6	48745.1	19663.8	711.7	44597.3
        2008.00000	2.18906	66.19700	19685.4	48775.3	19671.0	751.9	44626.4
        2009.00000	2.30505	66.20131	19694.2	48805.4	19678.2	792.1	44655.4
        2010.00000	2.42093	66.20553	19703.1	48835.6	19685.5	832.3	44684.5
        2011.00000	2.56036	66.21250	19709.7	48865.6	19690.0	880.5	44714.3
        2012.00000	2.69970	66.21934	19716.5	48895.6	19694.6	928.7	44744.2
        2013.00000	2.83895	66.22605	19723.4	48925.7	19699.2	976.9	44774.0
        2014.00000	2.97809	66.23262	19730.4	48955.8	19703.7	1025.1	44803.8
        2015.00000	3.11714	66.23906	19737.5	48986.0	19708.3	1073.3	44833.7
        2016.00000	3.25917	66.24182	19748.5	49018.5	19716.5	1122.7	44864.4
        2017.00000	3.40104	66.24445	19759.5	49051.2	19724.7	1172.2	44895.2
        ];
else
    data = [
        1965.00000 -1.60868 65.90531 19409.0 47542.4 19401.3 -544.9 43400.1
1966.00000 -1.56411 65.89785 19426.1 47570.4 19418.8 -530.2 43423.1
1967.00000 -1.51962 65.89038 19443.2 47598.4 19436.3 -515.6 43446.2
1968.00000 -1.47521 65.88290 19460.3 47626.4 19453.8 -501.0 43469.2
1969.00000 -1.43088 65.87543 19477.4 47654.4 19471.3 -486.4 43492.2
1970.00000 -1.38662 65.86794 19494.5 47682.4 19488.8 -471.7 43515.2
1971.00000 -1.34439 65.86111 19514.7 47719.1 19509.3 -457.9 43546.4
1972.00000 -1.30224 65.85428 19534.9 47755.8 19529.9 -444.0 43577.5
1973.00000 -1.26017 65.84745 19555.1 47792.4 19550.4 -430.1 43608.7
1974.00000 -1.21820 65.84061 19575.3 47829.1 19570.9 -416.2 43639.8
1975.00000 -1.17631 65.83378 19595.6 47865.8 19591.4 -402.3 43670.9
1976.00000 -1.05341 65.83892 19604.0 47896.1 19600.7 -360.4 43700.4
1977.00000 -0.93063 65.84397 19612.6 47926.5 19610.0 -318.5 43729.8
1978.00000 -0.80795 65.84891 19621.3 47956.9 19619.3 -276.7 43759.2
1979.00000 -0.68538 65.85374 19630.0 47987.3 19628.6 -234.8 43788.6
1980.00000 -0.56291 65.85848 19638.9 48017.7 19637.9 -192.9 43818.0
1981.00000 -0.45471 65.86826 19639.0 48036.4 19638.4 -155.9 43838.4
1982.00000 -0.34651 65.87797 19639.3 48055.2 19638.9 -118.8 43858.9
1983.00000 -0.23831 65.88759 19639.6 48074.0 19639.4 -81.7 43879.3
1984.00000 -0.13011 65.89713 19639.9 48092.8 19639.9 -44.6 43899.7
1985.00000 -0.02192 65.90658 19640.4 48111.6 19640.4 -7.5 43920.2
1986.00000 0.05557 65.92745 19633.8 48134.6 19633.7 19.0 43948.3
1987.00000 0.13312 65.94825 19627.2 48157.6 19627.1 45.6 43976.4
1988.00000 0.21073 65.96900 19620.6 48180.6 19620.5 72.2 44004.5
1989.00000 0.28838 65.98969 19614.1 48203.6 19613.9 98.7 44032.7
1990.00000 0.36608 66.01032 19607.6 48226.7 19607.2 125.3 44060.8
1991.00000 0.45436 66.02466 19604.8 48246.9 19604.2 155.5 44084.2
1992.00000 0.54266 66.03893 19602.1 48267.2 19601.2 185.7 44107.6
1993.00000 0.63098 66.05315 19599.4 48287.5 19598.2 215.8 44131.0
1994.00000 0.71933 66.06730 19596.7 48307.8 19595.1 246.0 44154.4
1995.00000 0.80771 66.08139 19594.1 48328.1 19592.1 276.2 44177.8
1996.00000 0.91115 66.09231 19599.8 48362.9 19597.3 311.7 44213.4
1997.00000 1.01453 66.10315 19605.5 48397.8 19602.5 347.1 44249.0
1998.00000 1.11784 66.11390 19611.4 48432.7 19607.6 382.6 44284.6
1999.00000 1.22110 66.12457 19617.3 48467.7 19612.8 418.1 44320.2
2000.00000 1.32429 66.13515 19623.2 48502.7 19618.0 453.5 44355.8
2001.00000 1.42743 66.14483 19629.6 48536.9 19623.5 489.0 44390.4
2002.00000 1.53050 66.15442 19636.0 48571.1 19629.0 524.5 44425.0
2003.00000 1.63350 66.16393 19642.4 48605.3 19634.4 559.9 44459.5
2004.00000 1.73644 66.17336 19649.0 48639.6 19639.9 595.4 44494.1
2005.00000 1.83930 66.18271 19655.5 48673.8 19645.4 630.9 44528.7
2006.00000 1.95561 66.18730 19664.1 48703.9 19652.6 671.0 44557.7
2007.00000 2.07181 66.19179 19672.7 48733.9 19659.9 711.2 44586.8
2008.00000 2.18791 66.19619 19681.5 48764.0 19667.1 751.4 44615.8
2009.00000 2.30390 66.20051 19690.3 48794.2 19674.4 791.5 44644.9
2010.00000 2.41979 66.20472 19699.2 48824.3 19681.6 831.7 44673.9
2011.00000 2.55923 66.21170 19705.8 48854.3 19686.2 879.9 44703.8
2012.00000 2.69857 66.21853 19712.6 48884.3 19690.7 928.1 44733.6
2013.00000 2.83782 66.22524 19719.5 48914.4 19695.3 976.3 44763.4
2014.00000 2.97697 66.23181 19726.5 48944.5 19699.8 1024.5 44793.2
2015.00000 3.11602 66.23825 19733.6 48974.7 19704.4 1072.7 44823.0
2016.00000 3.25804 66.24101 19744.5 49007.2 19712.6 1122.1 44853.8
2017.00000 3.39991 66.24364 19755.6 49039.9 19720.8 1171.6 44884.5 
        ];
end

YY = data(:, 1);
DE = data(:, 2);
IN = data(:, 3);
X = data(:, 6);
Y = data(:, 7);
Z = data(:, 8);

X0 = spline(YY, X, year);
Y0 = spline(YY, Y, year);
Z0 = spline(YY, Z, year);
inkl = pi / 180 * spline(YY, IN, year);
dekl = pi / 180 * spline(YY, DE, year);


% X0 = 18400;
% Y0 = 0; %#ok
% Z0 = 43500;

% inkl = atan(Z0 / X0);

a0 = cos(inkl) * cos(dekl);
b0 = cos(inkl) * sin(dekl);
c0 = sin(inkl);

if nargout < 1
    fprintf('X, Y, Z: %5.0f %5.0f %5.0f\n', X0, Y0, Z0);
end
