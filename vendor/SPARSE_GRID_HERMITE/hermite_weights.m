function weight = hermite_weights ( order )

%*****************************************************************************80
%
%% HERMITE_WEIGHTS returns weights for certain Gauss-Hermite quadrature rules.
%
%  Discussion:
%
%    The allowed orders are 1, 3, 7, 15, 31, 63 or 127.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 October 2007
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Milton Abramowitz, Irene Stegun,
%    Handbook of Mathematical Functions,
%    National Bureau of Standards, 1964,
%    ISBN: 0-486-61272-4,
%    LC: QA47.A34.
%
%    Arthur Stroud, Don Secrest,
%    Gaussian Quadrature Formulas,
%    Prentice Hall, 1966,
%    LC: QA299.4G3S7.
%
%  Parameters:
%
%    Input, integer ORDER, the order of the rule.
%    ORDER must be 1, 3, 7, 15, 31, 63, or 127.
%
%    Output, real WEIGHT(ORDER), the weights.
%    The weights are positive, symmetric and should sum to sqrt(PI).
%
  weight = zeros ( 1, order );

  if ( order == 1 )

    weight(1) = 1.77245385090551602729816748334E+00;

  elseif ( order == 3 )

    weight(1) = 0.295408975150919337883027913890E+00;
    weight(2) = 0.118163590060367735153211165556E+01;
    weight(3) = 0.295408975150919337883027913890E+00;

  elseif ( order == 7 )

    weight(1) = 0.971781245099519154149424255939E-03;
    weight(2) = 0.545155828191270305921785688417E-01;
    weight(3) = 0.425607252610127800520317466666E+00;
    weight(4) = 0.810264617556807326764876563813E+00;
    weight(5) = 0.425607252610127800520317466666E+00;
    weight(6) = 0.545155828191270305921785688417E-01;
    weight(7) = 0.971781245099519154149424255939E-03;

  elseif ( order == 15 )

    weight(1) =  0.152247580425351702016062666965E-08;
    weight(2) =  0.105911554771106663577520791055E-05;
    weight(3) =  0.100004441232499868127296736177E-03;
    weight(4) =  0.277806884291277589607887049229E-02;
    weight(5) =  0.307800338725460822286814158758E-01;
    weight(6) =  0.158488915795935746883839384960E+00;
    weight(7) =  0.412028687498898627025891079568E+00;
    weight(8) =  0.564100308726417532852625797340E+00;
    weight(9) =  0.412028687498898627025891079568E+00;
    weight(10) = 0.158488915795935746883839384960E+00;
    weight(11) = 0.307800338725460822286814158758E-01;
    weight(12) = 0.277806884291277589607887049229E-02;
    weight(13) = 0.100004441232499868127296736177E-03;
    weight(14) = 0.105911554771106663577520791055E-05;
    weight(15) = 0.152247580425351702016062666965E-08;

  elseif ( order == 31 )

    weight(  1) =   0.46189683944498305857470556847735E-21;
    weight(  2) =   0.51106090079112519643027197715274E-17;
    weight(  3) =   0.58995564987355133075257722133966E-14;
    weight(  4) =   0.18603735214463569590294465062239E-11;
    weight(  5) =   0.23524920032013205739850619940094E-09;
    weight(  6) =   0.14611988344865057576066495091513E-07;
    weight(  7) =   0.50437125589241034841778074689627E-06;
    weight(  8) =   0.10498602757642934202945441341697E-04;
    weight(  9) =   0.13952090395003623854995664958146E-03;
    weight( 10) =   0.12336833073030489880608311394968E-02;
    weight( 11) =   0.74827999140119116765002499116934E-02;
    weight( 12) =   0.31847230731201222775249585776902E-01;
    weight( 13) =   0.96717948160569462991143316029341E-01;
    weight( 14) =   0.21213278866810461318136114862419E+00;
    weight( 15) =   0.33877265789305344906000174083214E+00;
    weight( 16) =   0.39577855609737786462923720809676E+00;
    weight( 17) =   0.33877265789305344906000174083214E+00;
    weight( 18) =   0.21213278866810461318136114862419E+00;
    weight( 19) =   0.96717948160569462991143316029341E-01;
    weight( 20) =   0.31847230731201222775249585776902E-01;
    weight( 21) =   0.74827999140119116765002499116934E-02;
    weight( 22) =   0.12336833073030489880608311394968E-02;
    weight( 23) =   0.13952090395003623854995664958146E-03;
    weight( 24) =   0.10498602757642934202945441341697E-04;
    weight( 25) =   0.50437125589241034841778074689627E-06;
    weight( 26) =   0.14611988344865057576066495091513E-07;
    weight( 27) =   0.23524920032013205739850619940094E-09;
    weight( 28) =   0.18603735214463569590294465062239E-11;
    weight( 29) =   0.58995564987355133075257722133966E-14;
    weight( 30) =   0.51106090079112519643027197715274E-17;
    weight( 31) =   0.46189683944498305857470556847735E-21;

  elseif ( order == 63 )

    weight(  1) =   0.37099206434787551197827130470031E-47;
    weight(  2) =   0.10400778615192299534481914814892E-41;
    weight(  3) =   0.19796804708258311251124226474396E-37;
    weight(  4) =   0.84687478191640015120141181138947E-34;
    weight(  5) =   0.13071305930779945903630127634063E-30;
    weight(  6) =   0.93437837175367456929765381518998E-28;
    weight(  7) =   0.36027426635173044862245783257252E-25;
    weight(  8) =   0.82963863115951789374753323156164E-23;
    weight(  9) =   0.12266629909105281472971700203949E-20;
    weight( 10) =   0.12288435628797061539461585325494E-18;
    weight( 11) =   0.86925536958188009075932426691516E-17;
    weight( 12) =   0.44857058689176221240330804981619E-15;
    weight( 13) =   0.17335817955735154599902643794700E-13;
    weight( 14) =   0.51265062385038307838565047455223E-12;
    weight( 15) =   0.11808921844532942490513037158404E-10;
    weight( 16) =   0.21508698297808025739828859845140E-09;
    weight( 17) =   0.31371929535285447801497640621672E-08;
    weight( 18) =   0.37041625984781705796752840204084E-07;
    weight( 19) =   0.35734732949879669663960738150956E-06;
    weight( 20) =   0.28393114498380927832990899215541E-05;
    weight( 21) =   0.18709113003730498008961134765721E-04;
    weight( 22) =   0.10284880800653635546698378640623E-03;
    weight( 23) =   0.47411702610173128107201781718693E-03;
    weight( 24) =   0.18409222622384813438539657470055E-02;
    weight( 25) =   0.60436044551187631655712178246467E-02;
    weight( 26) =   0.16829299199599730926458559757600E-01;
    weight( 27) =   0.39858264027692992170237391875317E-01;
    weight( 28) =   0.80467087993950415219587554532823E-01;
    weight( 29) =   0.13871950817615293377792092082674E+00;
    weight( 30) =   0.20448695346833761570957197160475E+00;
    weight( 31) =   0.25799889943058042204920467417642E+00;
    weight( 32) =   0.27876694884838411919175686949858E+00;
    weight( 33) =   0.25799889943058042204920467417642E+00;
    weight( 34) =   0.20448695346833761570957197160475E+00;
    weight( 35) =   0.13871950817615293377792092082674E+00   ;
    weight( 36) =   0.80467087993950415219587554532823E-01;
    weight( 37) =   0.39858264027692992170237391875317E-01;
    weight( 38) =   0.16829299199599730926458559757600E-01;
    weight( 39) =   0.60436044551187631655712178246467E-02;
    weight( 40) =   0.18409222622384813438539657470055E-02;
    weight( 41) =   0.47411702610173128107201781718693E-03;
    weight( 42) =   0.10284880800653635546698378640623E-03;
    weight( 43) =   0.18709113003730498008961134765721E-04;
    weight( 44) =   0.28393114498380927832990899215541E-05;
    weight( 45) =   0.35734732949879669663960738150956E-06;
    weight( 46) =   0.37041625984781705796752840204084E-07;
    weight( 47) =   0.31371929535285447801497640621672E-08;
    weight( 48) =   0.21508698297808025739828859845140E-09;
    weight( 49) =   0.11808921844532942490513037158404E-10;
    weight( 50) =   0.51265062385038307838565047455223E-12;
    weight( 51) =   0.17335817955735154599902643794700E-13;
    weight( 52) =   0.44857058689176221240330804981619E-15;
    weight( 53) =   0.86925536958188009075932426691516E-17;
    weight( 54) =   0.12288435628797061539461585325494E-18;
    weight( 55) =   0.12266629909105281472971700203949E-20;
    weight( 56) =   0.82963863115951789374753323156164E-23;
    weight( 57) =   0.36027426635173044862245783257252E-25;
    weight( 58) =   0.93437837175367456929765381518998E-28;
    weight( 59) =   0.13071305930779945903630127634063E-30;
    weight( 60) =   0.84687478191640015120141181138947E-34;
    weight( 61) =   0.19796804708258311251124226474396E-37;
    weight( 62) =   0.10400778615192299534481914814892E-41;
    weight( 63) =   0.37099206434787551197827130470031E-47;

  elseif ( order == 127 )

    weight(  1) =   0.12504497577050595552677230002883E-100;
    weight(  2) =   0.17272798059419131415318615789672E-93;
    weight(  3) =   0.89321681571986548608031150791499E-88;
    weight(  4) =   0.77306185240893578449625186483810E-83;
    weight(  5) =   0.20143957652648255497735460506196E-78;
    weight(  6) =   0.21503714733610239701351039429345E-74;
    weight(  7) =   0.11341924208594594813715533569504E-70;
    weight(  8) =   0.33489139011795051950683388483136E-67;
    weight(  9) =   0.60486548964016681064424451668405E-64;
    weight( 10) =   0.71375092946352177824971347343892E-61;
    weight( 11) =   0.57884563374885556636801095624030E-58;
    weight( 12) =   0.33581166223858230300409326551248E-55;
    weight( 13) =   0.14394641949253923568603163698953E-52;
    weight( 14) =   0.46821808383216117724080263903889E-50;
    weight( 15) =   0.11817054440684264071348471955361E-47;
    weight( 16) =   0.23581659156008927203181682045005E-45;
    weight( 17) =   0.37814427940797540210712758405540E-43;
    weight( 18) =   0.49411031115771638145610738414006E-41;
    weight( 19) =   0.53255303775425059266087298458297E-39;
    weight( 20) =   0.47854390680131484999315199332765E-37;
    weight( 21) =   0.36191883445952356128627543209554E-35;
    weight( 22) =   0.23232083386343554805352497446119E-33;
    weight( 23) =   0.12753331411008716683688974281454E-31;
    weight( 24) =   0.60277753850758742112436095241270E-30;
    weight( 25) =   0.24679773241777200207460855084439E-28;
    weight( 26) =   0.88019567691698482573264198727415E-27;
    weight( 27) =   0.27482489212040561315005725890593E-25;
    weight( 28) =   0.75468218903085486125222816438456E-24;
    weight( 29) =   0.18303134636280466270545996891835E-22;
    weight( 30) =   0.39355990860860813085582448449811E-21;
    weight( 31) =   0.75293161638581191068419292570042E-20;
    weight( 32) =   0.12857997786722855037584105682618E-18;
    weight( 33) =   0.19659326888445857792541925311450E-17;
    weight( 34) =   0.26986511907214101894995783364250E-16;
    weight( 35) =   0.33344414303198856330118301113874E-15;
    weight( 36) =   0.37173303125150639885726463109574E-14;
    weight( 37) =   0.37473954472839737091885387788983E-13;
    weight( 38) =   0.34230094493397259538669512076007E-12;
    weight( 39) =   0.28385303724993373166810860630552E-11;
    weight( 40) =   0.21406920290454669208938772802828E-10;
    weight( 41) =   0.14706331273431716244229273183839E-09;
    weight( 42) =   0.92173940967434659264335883218167E-09;
    weight( 43) =   0.52781663936972714041837056042506E-08;
    weight( 44) =   0.27650497044951117835905283127679E-07;
    weight( 45) =   0.13267855842539464770913063113371E-06;
    weight( 46) =   0.58380944276113062188573331195042E-06;
    weight( 47) =   0.23581561724775629112332165335800E-05;
    weight( 48) =   0.87524468034280444703919485644809E-05;
    weight( 49) =   0.29876790535909012274846532159647E-04;
    weight( 50) =   0.93874435720072545206729594267039E-04;
    weight( 51) =   0.27170762627931172053444716883938E-03;
    weight( 52) =   0.72493929742498358979684249380921E-03;
    weight( 53) =   0.17841208326763432884316727108264E-02;
    weight( 54) =   0.40524855186046131499765636276283E-02;
    weight( 55) =   0.85000263041544110385806705526917E-02;
    weight( 56) =   0.16471142241609687824005585301760E-01;
    weight( 57) =   0.29499296248213632269675010319119E-01;
    weight( 58) =   0.48847387114300011006959603975676E-01;
    weight( 59) =   0.74807989768583731416517226905270E-01;
    weight( 60) =   0.10598520508090929403834368934301E+00;
    weight( 61) =   0.13893945309051540832066283010510E+00;
    weight( 62) =   0.16856236074207929740526975049765E+00;
    weight( 63) =   0.18927849580120432177170145550076E+00;
    weight( 64) =   0.19673340688823289786163676995151E+00;
    weight( 65) =   0.18927849580120432177170145550076E+00;
    weight( 66) =   0.16856236074207929740526975049765E+00;
    weight( 67) =   0.13893945309051540832066283010510E+00;
    weight( 68) =   0.10598520508090929403834368934301E+00;
    weight( 69) =   0.74807989768583731416517226905270E-01;
    weight( 70) =   0.48847387114300011006959603975676E-01;
    weight( 71) =   0.29499296248213632269675010319119E-01;
    weight( 72) =   0.16471142241609687824005585301760E-01;
    weight( 73) =   0.85000263041544110385806705526917E-02;
    weight( 74) =   0.40524855186046131499765636276283E-02;
    weight( 75) =   0.17841208326763432884316727108264E-02;
    weight( 76) =   0.72493929742498358979684249380921E-03;
    weight( 77) =   0.27170762627931172053444716883938E-03;
    weight( 78) =   0.93874435720072545206729594267039E-04;
    weight( 79) =   0.29876790535909012274846532159647E-04;
    weight( 80) =   0.87524468034280444703919485644809E-05;
    weight( 81) =   0.23581561724775629112332165335800E-05;
    weight( 82) =   0.58380944276113062188573331195042E-06;
    weight( 83) =   0.13267855842539464770913063113371E-06;
    weight( 84) =   0.27650497044951117835905283127679E-07;
    weight( 85) =   0.52781663936972714041837056042506E-08;
    weight( 86) =   0.92173940967434659264335883218167E-09;
    weight( 87) =   0.14706331273431716244229273183839E-09;
    weight( 88) =   0.21406920290454669208938772802828E-10;
    weight( 89) =   0.28385303724993373166810860630552E-11;
    weight( 90) =   0.34230094493397259538669512076007E-12;
    weight( 91) =   0.37473954472839737091885387788983E-13;
    weight( 92) =   0.37173303125150639885726463109574E-14;
    weight( 93) =   0.33344414303198856330118301113874E-15;
    weight( 94) =   0.26986511907214101894995783364250E-16;
    weight( 95) =   0.19659326888445857792541925311450E-17;
    weight( 96) =   0.12857997786722855037584105682618E-18;
    weight( 97) =   0.75293161638581191068419292570042E-20;
    weight( 98) =   0.39355990860860813085582448449811E-21;
    weight( 99) =   0.18303134636280466270545996891835E-22;
    weight(100) =   0.75468218903085486125222816438456E-24;
    weight(101) =   0.27482489212040561315005725890593E-25;
    weight(102) =   0.88019567691698482573264198727415E-27;
    weight(103) =   0.24679773241777200207460855084439E-28;
    weight(104) =   0.60277753850758742112436095241270E-30;
    weight(105) =   0.12753331411008716683688974281454E-31;
    weight(106) =   0.23232083386343554805352497446119E-33;
    weight(107) =   0.36191883445952356128627543209554E-35;
    weight(108) =   0.47854390680131484999315199332765E-37;
    weight(109) =   0.53255303775425059266087298458297E-39;
    weight(110) =   0.49411031115771638145610738414006E-41;
    weight(111) =   0.37814427940797540210712758405540E-43;
    weight(112) =   0.23581659156008927203181682045005E-45;
    weight(113) =   0.11817054440684264071348471955361E-47;
    weight(114) =   0.46821808383216117724080263903889E-50;
    weight(115) =   0.14394641949253923568603163698953E-52;
    weight(116) =   0.33581166223858230300409326551248E-55;
    weight(117) =   0.57884563374885556636801095624030E-58;
    weight(118) =   0.71375092946352177824971347343892E-61;
    weight(119) =   0.60486548964016681064424451668405E-64;
    weight(120) =   0.33489139011795051950683388483136E-67;
    weight(121) =   0.11341924208594594813715533569504E-70;
    weight(122) =   0.21503714733610239701351039429345E-74;
    weight(123) =   0.20143957652648255497735460506196E-78;
    weight(124) =   0.77306185240893578449625186483810E-83;
    weight(125) =   0.89321681571986548608031150791499E-88;
    weight(126) =   0.17272798059419131415318615789672E-93;
    weight(127) =   0.12504497577050595552677230002883E-100;

  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'HERMITE_WEIGHTS - Fatal error!\n' );
    fprintf ( 1, '  Illegal value of ORDER = %d\n', order );
    fprintf ( 1, '  Legal values are 1, 3, 7, 15, 31, 63 or 127.\n' );
    error ( 'HERMITE_WEIGHTS - Fatal error!' );

  end

  return
end
