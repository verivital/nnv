function [rf, a11, a12, a21, a22, a31, a32, a41, a42, b11, b12, b21, b22, b31, b32, b41, b42, c11, c12, c21, c22, c31, c32, c41, c42] = get_verifyTime(VT, VT_impr, i)
    % vt: verification results
    % i : row index  
    vt = VT(i,:);
    vt_impr = VT_impr(i, :);
    rf = vt.RelaxFactor;
    a11 = vt.de_005(1);
    a21 = vt.de_005(2);
    a31 = vt.de_005(3);
    a41 = vt.de_005(4);
    b11 = vt.de_01(1);
    b21 = vt.de_01(2);
    b31 = vt.de_01(3);
    b41 = vt.de_01(4);
    c11 = vt.de_02(1);
    c21 = vt.de_02(2);
    c31 = vt.de_02(3);
    c41 = vt.de_02(4);
    
    a12 = vt_impr.de_005(1);
    a22 = vt_impr.de_005(2);
    a32 = vt_impr.de_005(3);
    a42 = vt_impr.de_005(4);
    b12 = vt_impr.de_01(1);
    b22 = vt_impr.de_01(2);
    b32 = vt_impr.de_01(3);
    b42 = vt_impr.de_01(4);
    c12 = vt_impr.de_02(1);
    c22 = vt_impr.de_02(2);
    c32 = vt_impr.de_02(3);
    c42 = vt_impr.de_02(4);
end